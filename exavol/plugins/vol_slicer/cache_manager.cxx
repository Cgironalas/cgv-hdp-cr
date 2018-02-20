#include "cache_manager.h"
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <vol_data/volume.h>
#include <cgv/media/image/image_writer.h>
#include <chrono>
#include <thread>
#include "volume_slicer.h"

cache_manager::cache_manager(volume_slicer &f) : vs(f) {
	//should be calculated from the GPU capacity
	cache_size_blocks = 5000;
	thread_limit = 8;
	process_running = false;
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
		
		thread_report.lock();

		intersected_blocks_queue_lock.lock();
		if (thread_amount < thread_limit && !blocks_in_progress.empty() && !process_running) {
			thread_amount += 1;
			process_running = true;
			std::cout << "frame started " << std::endl;
			std::thread t([&]() { retrieve_blocks_in_plane(); });
			t.detach();
		}

		intersected_blocks_queue_lock.unlock();
		thread_report.unlock();

		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		
	}
}

void cache_manager::process_blocks(std::vector<ivec3> pIntersected_blocks) {
	intersected_blocks_queue_lock.lock();
	for (auto it : pIntersected_blocks) {
		blocks_cache_lock.lock();
		if (block_cache_map.find(it) == block_cache_map.end()) {
			blocks_cache_lock.unlock();
			blocks_in_progress.push_back(it);
		} else {
			blocks_cache_lock.unlock();
		}
	}
	intersected_blocks_queue_lock.unlock();
}

void cache_manager::retrieve_blocks_in_plane() {

	const vec3 dimensions = vs.dimensions;
	const vec3 block_dimensions = vs.block_dimensions;
	const vec3 overlap = vs.overlap;

	std::vector<ivec3> blocks_to_retrieve;
	blocks_in_progress_lock.lock();
	std::copy(blocks_in_progress.begin(), blocks_in_progress.end(), std::back_inserter(blocks_to_retrieve));
	blocks_in_progress.clear();
	blocks_in_progress_lock.unlock();

	//these three values need to be recalculated only if the GUI options have changed
	ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	size_t block_size = (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size()) * (block_dimensions(2) + overlap(2));
	ivec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));


	//this code loads from the slices (.bvx files from blockgen.cxx project) 
	// and uses the output to generate .tiff blocks for validation if the code is accesible

	//interaction with the GPU cache should occur here

	bool are_blocks_in_queue = true;

	for (int i = 0; i < blocks_to_retrieve.size(); i++){

		thread_report.lock();

		if (thread_amount < thread_limit) {
			thread_amount += 1;
			std::thread t([&]() { retrieve_block(blocks_to_retrieve[i], nr_blocks, block_size, df_dim); });
			//blocks_to_retrieve.erase(it);
			t.detach();
		} else {
			i--;
		}

		thread_report.unlock();
		are_blocks_in_queue = blocks_to_retrieve.empty();
	}

	process_running = false;
	std::cout << "frame loaded " << std::endl;


	thread_report.lock();
	thread_amount -= 1;
	thread_report.unlock();

	//vs.post_redraw();
}

bool cache_manager::retrieve_block(ivec3& block, const ivec3& nr_blocks, const size_t& block_size, const vec3& df_dim) {

	lru_refer(block);

	thread_report.lock();
	thread_amount -= 1;
	thread_report.unlock();
	return false;

	std::string input_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices";
	std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";

	std::cout << "retrieving block at: " << block(0) << ", " << block(1) << ", " << block(2) << ", " << std::endl;

	std::stringstream ss;
	ss << input_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(2) << ".bvx";

	// compute block index and pointer to block data
	unsigned bi = block(1) * nr_blocks(0) + block(0);
	char* block_ptr = new char[block_size];

	// throws exception at fread in some points 
	FILE* fp = fopen(ss.str().c_str(), "rb");
	fseek(fp, bi * block_size, SEEK_SET);
	try {

		bool result = fread(block_ptr, block_size, 1, fp) == 1;

		if (result) {
			fclose(fp);
		}
		else {
			std::cout << "failed to retrieve block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
		}

		thread_report.lock();
		thread_amount -= 1;
		thread_report.unlock();

		return result;
			
		//return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	}
	catch (const std::exception& e) {
		std::cout << "Exception thrown: " << e.what() << std::endl;
		thread_report.lock();
		thread_amount -= 1;
		thread_report.unlock();
		return false;
	}
}

void cache_manager::lru_refer(ivec3 block) {
	

	std::list<ivec3>::iterator it = std::find(intersected_blocks_queue.begin(), intersected_blocks_queue.end(), block);

	if (it == intersected_blocks_queue.end()) {
		if (intersected_blocks_queue.size() == cache_size_blocks) {

			ivec3 last = intersected_blocks_queue.back();

			intersected_blocks_queue_lock.lock();
			intersected_blocks_queue.pop_back();
			intersected_blocks_queue_lock.unlock();

			// this should actually erase the block from GPU memory 
			blocks_cache_lock.lock();
			block_cache_map.erase(last);
			blocks_cache_lock.unlock();
		}
	}
	else {
		intersected_blocks_queue_lock.lock();
		intersected_blocks_queue.erase(it);
		intersected_blocks_queue_lock.unlock();

		// this should actually erase the block from GPU memory
		blocks_cache_lock.lock();
		block_cache_map.erase(block);
		blocks_cache_lock.unlock();
	}

	intersected_blocks_queue_lock.lock();
	intersected_blocks_queue.push_front(block);
	intersected_blocks_queue_lock.unlock();

	// this should store the gpu position, not the queue position
	blocks_cache_lock.lock();
	block_cache_map[block] = intersected_blocks_queue.begin();
	blocks_cache_lock.unlock();
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
