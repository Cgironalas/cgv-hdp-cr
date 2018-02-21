#include "cache_manager.h"
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <vol_data/volume.h>
#include <cgv/media/image/image_writer.h>
#include <chrono>
#include <thread>
#include "volume_slicer.h"

//toggles multithreading in individual block retrieval
bool in_thread = false;

//path to the slices (.bvx files from blockgen.cxx project) 	
std::string input_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices";

// output to generate .tiff blocks for validation if the write_tiff function is called
//std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";


cache_manager::cache_manager(volume_slicer &f) : vs(f) {
	cache_size_blocks = 15000; //should be calculated from the GPU capacity if possible
	thread_limit = 16;
	process_running = false;
	signal_restart = false;
}

void cache_manager::init_listener() {
	
	thread_report.lock();
	thread_amount += 1;
	std::thread t([&]() { main_loop(); });
	t.detach();
	thread_report.unlock();
}

void cache_manager::main_loop() {
	
	while (true) {
		
		blocks_in_progress_lock.lock();
		bool something_to_do = !blocks_in_progress.empty();
		blocks_in_progress_lock.unlock();

		thread_report.lock();
		if (!in_thread && thread_amount > 2) {
			std::cout << "SOMETHING IS WRONG" << std::endl;
		}

		if (thread_amount < thread_limit && something_to_do && !process_running) {
			process_running = true;
			thread_amount += 1;
			std::thread t([&]() { retrieve_blocks_in_plane(); });
			t.detach();
		}
		thread_report.unlock();
	}
}

void cache_manager::request_blocks(std::vector<ivec3> blocks_batch) {
	blocks_in_progress_lock.lock();
	for (auto it : blocks_batch) {

		blocks_cache_lock.lock();
		bool block_in_cache = block_cache_map.find(it) != block_cache_map.end();
		blocks_cache_lock.unlock();

		if (!block_in_cache) {
			blocks_in_progress.push_back(it);
		}
	}
	signal_restart = true;
	blocks_in_progress_lock.unlock();
	
}

void cache_manager::retrieve_blocks_in_plane() {
	const vec3 dimensions = vs.dimensions;
	const vec3 block_dimensions = vs.block_dimensions;
	const vec3 overlap = vs.overlap;

	blocks_in_progress_lock.lock();
	std::vector<ivec3> blocks_batch_frame;
	std::copy(blocks_in_progress.begin(), blocks_in_progress.end(), std::back_inserter(blocks_batch_frame));
	blocks_in_progress.clear();
	blocks_in_progress_lock.unlock();

	//these three values need to be recalculated only if the GUI options have changed
	ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	size_t block_size = (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size()) * (block_dimensions(2) + overlap(2));
	ivec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));

	// traverse the batch
	bool cancelled = false;
	for (int i = 0; i < blocks_batch_frame.size(); i++){

		blocks_in_progress_lock.lock();
		bool queued_blocks = !blocks_in_progress.empty();
		blocks_in_progress_lock.unlock();

		if (signal_restart && queued_blocks) {
			signal_restart = false;
			cancelled = true;
			std::cout << "cancelled remaining " << blocks_batch_frame.size() - i << " blocks." << std::endl;
			break;
		}
		
		if (in_thread) {
			thread_report.lock();

			if (thread_amount < thread_limit) {
				thread_amount += 1;
				thread_report.unlock();
				if (i % 1000 == 0) {
					std::cout << "blocks read_: " << i << " " << std::endl;
				}
				std::thread t([&]() { retrieve_block(blocks_batch_frame[i], nr_blocks, block_size, df_dim); });
				t.detach();
			}
			else {
				i--;
				thread_report.unlock();
			}
		} else {
			retrieve_block(blocks_batch_frame[i], nr_blocks, block_size, df_dim);
			if (i % 1000 == 0) {
				std::cout << "blocks read: " << i << " " << std::endl;
			}
		}
	}

	if (!cancelled) {
		std::cout << "finished batch of  " << blocks_batch_frame.size() << " blocks." << std::endl;

		// no real effect for now but this should update the frame with the new block textures loaded.
		vs.post_redraw();
	}
		
	thread_report.lock();
	thread_amount -= 1;
	process_running = false;
	thread_report.unlock();

}

bool cache_manager::retrieve_block(ivec3& block, const ivec3& nr_blocks, const size_t& block_size, const vec3& df_dim) {
	
	//std::cout << "retrieving block at: " << block(0) << ", " << block(1) << ", " << block(2) << ", " << std::endl;
	try {

		std::stringstream ss;
		ss << input_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(2) << ".bvx";

		// compute block index and pointer to block data
		unsigned bi = block(1) * nr_blocks(0) + block(0);
		char* block_ptr = new char[block_size];

		// throws exception at fread in some points 
		FILE* fp = fopen(ss.str().c_str(), "rb");
		fseek(fp, bi * block_size, SEEK_SET);

		bool result = fread(block_ptr, block_size, 1, fp) == 1;

		if (result) {
			lru_refer(block, block_ptr);
			fclose(fp);
		}
		else {
			std::cout << "failed to retrieve block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
		}

		if (in_thread) {
			thread_report.lock();
			thread_amount -= 1;
			thread_report.unlock();
		}
		
		return result; 

		/// to ensure correct block being loaded
		//return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	} catch (...) {
		std::cout << "exception thrown at block retrieval due to existing block slice dimensions" << std::endl;
	
		if (in_thread) {
			thread_report.lock();
			thread_amount -= 1;
			thread_report.unlock();
		}

		return false;
	}
}

void cache_manager::lru_refer(ivec3 block, char* block_ptr) {
	
	intersected_blocks_queue_lock.lock();
	std::list<ivec3>::iterator it = std::find(intersected_blocks_queue.begin(), intersected_blocks_queue.end(), block);
	intersected_blocks_queue_lock.unlock();

	// not in cache
	if (it == intersected_blocks_queue.end()) {
		// cache is full
		if (intersected_blocks_queue.size() == cache_size_blocks) {

			intersected_blocks_queue_lock.lock();
			ivec3 last = intersected_blocks_queue.back();
			intersected_blocks_queue.pop_back();
			intersected_blocks_queue_lock.unlock();

			blocks_cache_lock.lock();
			block_cache_map.erase(last);
			blocks_cache_lock.unlock();
		}
	} else {
		intersected_blocks_queue_lock.lock();
		intersected_blocks_queue.erase(it);
		intersected_blocks_queue_lock.unlock();

		blocks_cache_lock.lock();
		block_cache_map.erase(block);
		blocks_cache_lock.unlock();
	}

	intersected_blocks_queue_lock.lock();
	blocks_cache_lock.lock();

	intersected_blocks_queue.push_front(block);
	block_cache_map[block] = intersected_blocks_queue.begin();
	
	blocks_cache_lock.unlock();
	intersected_blocks_queue_lock.unlock();
}


//! write a block located at \c data_ptr to a tiff file
bool cache_manager::write_tiff_block(const std::string& file_name, const char* data_ptr, const vec3 df_dim, const std::string& options)
{
	// setup output format 
	cgv::data::data_format df;
	df.set_component_format(cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB));
	df.set_width(df_dim(0));
	df.set_height(df_dim(1));

	// setup image writer
	cgv::media::image::image_writer iw(file_name + ".tif");
	iw.multi_set(options);

	// write the file
	size_t slice_size = df.get_nr_bytes();
	for (unsigned i = 0; i<(unsigned)(df_dim(2)); ++i) {
		cgv::data::const_data_view dv(&df, data_ptr + slice_size*i);
		if (!iw.write_image(dv)) {
			std::cerr << "could not write slice " << i << " to tiff file " << file_name << std::endl;
			return false;
		}
	}
	return iw.close();
}
