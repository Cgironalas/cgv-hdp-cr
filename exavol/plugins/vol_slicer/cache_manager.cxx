#include "cache_manager.h"
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <vol_data/volume.h>
#include <cgv/media/image/image_writer.h>
#include <chrono>
#include <thread>
#include "volume_slicer.h"

//path to the slices (.bvx files from blockgen.cxx project) 	
std::string input_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices";

// output to generate .tiff blocks for validation if the write_tiff function is called
//std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";


cache_manager::cache_manager(volume_slicer &f) : vs(f) {
	cpu_cache_size_blocks = 20000; //should be calculated from the CPU capacity if possible
	gpu_cache_size_blocks = 13000; //should be calculated from the GPU capacity if possible

	// extra validation for threaded approaches
	thread_limit = 8; 

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
		bool something_to_do = !disk_queue_blocks.empty();
		blocks_in_progress_lock.unlock();

		thread_report.lock();

		if (thread_amount > 2) {
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
	bool new_blocks = false;
	for (auto it : blocks_batch) {

		cpu_cache_lock.lock();
		bool block_in_cache = cpu_block_cache_map.find(it) != cpu_block_cache_map.end();
		cpu_cache_lock.unlock();

		if (!block_in_cache) {
			disk_queue_blocks.push_back(it);
			new_blocks = true;
			
		} else {
			cache_queue_blocks.push_back(it);
		}
	}

	if (new_blocks) {
		signal_restart = true;
	}
		
	blocks_in_progress_lock.unlock();
}

void cache_manager::retrieve_blocks_in_plane() {

	const vec3 dimensions = vs.dimensions;
	const vec3 block_dimensions = vs.block_dimensions;
	const vec3 overlap = vs.overlap;

	blocks_in_progress_lock.lock();
	
	std::vector<ivec3> blocks_batch_frame;
	std::copy(disk_queue_blocks.begin(), disk_queue_blocks.end(), std::back_inserter(blocks_batch_frame));
	disk_queue_blocks.clear();
	
	std::vector<ivec3> cached_batch_frame;
	std::copy(cache_queue_blocks.begin(), cache_queue_blocks.end(), std::back_inserter(cached_batch_frame));
	cache_queue_blocks.clear();

	blocks_in_progress_lock.unlock();

	//these three values need to be recalculated only if the GUI options have changed
	ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	size_t block_size = (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size()) * (block_dimensions(2) + overlap(2));
	ivec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));

	// traverse the batch
	bool cancelled = false;

	std::cout << "\n\n total blocks to retrieve: " << blocks_batch_frame.size() + cached_batch_frame.size() << " \n" << std::endl;


	std::string sources[2] = { "disk", "cache" };
	std::vector<ivec3> frame_sources[2] = { blocks_batch_frame, cached_batch_frame };

	for (int s = 0; s < 2; s++) {
		if (!cancelled) {

			std::cout << "\nsize batch from " << sources[s] << ": " << frame_sources[s].size() << " \n" << std::endl;

			for (int i = 0; i < frame_sources[s].size(); i++) {

				blocks_in_progress_lock.lock();
				bool queued_blocks = s == 0 ? !disk_queue_blocks.empty() : !cache_queue_blocks.empty();
				blocks_in_progress_lock.unlock();

				if (signal_restart && queued_blocks) {
					signal_restart = false;
					cancelled = true;

					std::cout << "cancelled remaining " << frame_sources[s].size() - i << " blocks from " << sources[s] << std::endl;
					break;
				}

				if (s == 0) {
					lru_refer(blocks_batch_frame[i], nr_blocks, block_size, df_dim);

				} else {
					cpu_cache_lock.lock();
					char* block_ptr = cpu_block_cache_map[cached_batch_frame[i]];
					cpu_cache_lock.unlock();

					gpu_lru_refer(cached_batch_frame[i], block_ptr);
				}

				if (i % 1000 == 0) {
					std::cout << "blocks read from " << sources[s] << ":" << i << " " << std::endl;
				}
			}
		}
	}

	if (!cancelled) {
		std::cout << "finished batch of  " << blocks_batch_frame.size() << " blocks. Currently  " << cpu_block_cache_map.size() << " blocks in cpu and " << gpu_block_cache_map.size() << " in gpu " << std::endl;
		
		// no real effect for now but this should update the frame with the new block textures loaded.
		vs.post_redraw();
	}
		
	thread_report.lock();
	thread_amount -= 1;
	process_running = false;
	thread_report.unlock();

}

char* cache_manager::retrieve_block(ivec3& block, const ivec3& nr_blocks, const size_t& block_size, const vec3& df_dim) {
	
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
			fclose(fp);
			return block_ptr;
		}
		else {
			std::cout << "failed to retrieve block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
			return NULL;
		}

		/// to ensure correct block being loaded
		//return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	} catch (...) {
		std::cout << "exception thrown at block retrieval due to existing block slice dimensions" << std::endl;
		return NULL;
	}
}

void cache_manager::lru_refer(ivec3& block, const ivec3& nr_blocks, const size_t& block_size, const vec3& df_dim) {
	
	// note: because of the request_blocks filtering of blocks not in cache, we assume this function will never be called with a block already present in cpu cache
	
	//lock mutexes
	cpu_blocks_queue_lock.lock();
	cpu_cache_lock.lock();
	
	if (cpu_blocks_queue.size() == cpu_cache_size_blocks) {

		ivec3 last = cpu_blocks_queue.back();
		cpu_blocks_queue.pop_back();
		cpu_block_cache_map.erase(last);
	}

	cpu_block_cache_map[block] = retrieve_block(block, nr_blocks, block_size, df_dim);
	cpu_blocks_queue.push_front(block);

	char* block_ptr = cpu_block_cache_map[block];

	//unlock mutexes
	cpu_cache_lock.unlock();
	cpu_blocks_queue_lock.unlock();
	
	gpu_lru_refer(block, block_ptr);
}

void cache_manager::gpu_lru_refer(ivec3& block, char* block_ptr) {

	//lock mutexes
	gpu_blocks_queue_lock.lock();
	gpu_cache_lock.lock();

	std::list<ivec3>::iterator it = std::find(gpu_blocks_queue.begin(), gpu_blocks_queue.end(), block);

	// not in cache
	if (it == gpu_blocks_queue.end()) {
		// cache is full
		if (gpu_blocks_queue.size() == gpu_cache_size_blocks) {

			ivec3 last = gpu_blocks_queue.back();
			gpu_blocks_queue.pop_back();
			gpu_block_cache_map.erase(last);
		}

	} else { // in cache
		gpu_blocks_queue.erase(it);
	}

	gpu_blocks_queue.push_front(block);
	gpu_block_cache_map[block] = "";
	// call upload to gpu should go here

	//unlock mutexes
	gpu_cache_lock.unlock();
	gpu_blocks_queue_lock.unlock();
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
