#include "cache_manager.h"
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <vol_data/volume.h>
#include <cgv/media/image/image_writer.h>
#include <chrono>
#include <thread>
#include <string>
#include "volume_slicer.h"


// output to generate .tiff blocks for validation if the write_tiff function is called
//std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";

cache_manager::cache_manager(volume_slicer &f) : vs(f) {

	slices_path = "";
	slices_dimensions = ivec3(0, 0, 0);
	cpu_cache_size_blocks = 150000; //should be calculated from the CPU capacity if possible
	gpu_cache_size_blocks = 15000; //should be calculated from the GPU capacity if possible

	// extra validation for threaded approaches
	thread_limit = 8; 

	process_running = false;
	signal_restart = false;
}

void cache_manager::close_slices_files() {

	slices_files_names.clear();
	for (auto it : slices_files)
		fclose(it);
	slices_files.clear();
}

void cache_manager::set_block_folder(std::string folder_path) {
	restart_lock.lock();
	gpu_cache_lock.lock();
	cpu_cache_lock.lock();
	blocks_in_progress_lock.lock();

	slices_path = folder_path;

	close_slices_files();
	if (slices_path != "") {
	
		void* handle = cgv::utils::file::find_first(folder_path + "*.*");
		while (handle) {
			if (!cgv::utils::file::find_directory(handle)) {
				std::string file_name = folder_path + cgv::utils::file::find_name(handle);
				slices_files_names.push_back(file_name);

				FILE* fp = fopen(file_name.c_str(), "rb");
				if (fp != NULL) 
					slices_files.push_back(fp);
				else {
					slices_path = folder_path;
					break;
				}
			}
				
			handle = cgv::utils::file::find_next(handle);
		}
	}
	
	cpu_blocks_queue.clear();
	cpu_block_cache_map.clear();

	gpu_blocks_queue.clear();
	gpu_block_cache_map.clear();

	disk_queue_blocks.clear();
	cache_queue_blocks.clear();
	
	restart_lock.unlock();
	cpu_cache_lock.unlock();
	gpu_cache_lock.unlock();
	blocks_in_progress_lock.unlock();
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

		restart_lock.lock();
		bool start_retrieval = thread_amount < thread_limit && something_to_do && !process_running && slices_path != "";
		restart_lock.unlock();

		if (start_retrieval) {
			process_running = true;
			retrieve_blocks_in_plane();
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
		restart_lock.lock();
		signal_restart = true;
		restart_lock.unlock();
	}
		
	blocks_in_progress_lock.unlock();
}

void cache_manager::retrieve_blocks_in_plane() {

	vec3 dimensions = vs.slices_dimensions;
	vec3 block_dimensions = vs.block_dimensions;
	std::cout << "\nblock sizes " << block_dimensions(0) << " " << block_dimensions(1) << " " << block_dimensions(1) << std::endl;
	vec3 overlap = vs.overlap;

	blocks_in_progress_lock.lock();
	
	std::vector<ivec3> blocks_batch_frame;
	std::copy(disk_queue_blocks.begin(), disk_queue_blocks.end(), std::back_inserter(blocks_batch_frame));
	disk_queue_blocks.clear();
	
	std::vector<ivec3> cached_batch_frame;
	std::copy(cache_queue_blocks.begin(), cache_queue_blocks.end(), std::back_inserter(cached_batch_frame));
	cache_queue_blocks.clear();

	blocks_in_progress_lock.unlock();

	
	ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	size_t block_size = (long) ( (block_dimensions(1) + overlap(1))*
		((block_dimensions(0) + overlap(0))*
			 (cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size())) 
		* (block_dimensions(2) + overlap(2)) );
	vec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));

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

				restart_lock.lock();
				bool restart = (signal_restart && queued_blocks) | slices_path == "";
				restart_lock.unlock();

				if (restart) {
					restart_lock.lock();
					signal_restart = false;
					restart_lock.unlock();

					cancelled = true;
					std::cout << "cancelled remaining " << frame_sources[s].size() - i << " blocks from " << sources[s] << std::endl;
					break;
				}

				if (s == 0) {
					fifo_refer(blocks_batch_frame[i], nr_blocks, block_size, df_dim);

				} else {
					cpu_cache_lock.lock();
					char* block_ptr = cpu_block_cache_map[cached_batch_frame[i]];
					cpu_cache_lock.unlock();

					gpu_fifo_refer(cached_batch_frame[i], block_ptr);
				}

				if (i % 1000 == 0 || i % (frame_sources[s].size()-1) == 0) {
					std::cout << "blocks read from " << sources[s] << ":" << i+1 << " " << std::endl;
				}
			}
		}
	}

	if (!cancelled) {
		std::cout << "finished batch of " << blocks_batch_frame.size() + cached_batch_frame.size() << " blocks. Currently  " << cpu_block_cache_map.size() << " blocks in cpu and " << gpu_block_cache_map.size() << " in gpu " << std::endl;
		
		// no real effect for now but this should update the frame with the new block textures loaded.
		vs.post_redraw();
	}
	process_running = false;
}

char* cache_manager::retrieve_block(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {
	
	//std::cout << "retrieving block at: " << block(0) << ", " << block(1) << ", " << block(2) << ", " << std::endl;
	try {

		std::stringstream ss;
		ss << slices_path << "level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(2) << ".bvx";

		// compute block index and pointer to block data
		unsigned bi = block(1) * nr_blocks(0) + block(0);
		char* block_ptr = new char[block_size];

		// throws exception at fread in some points 
		if (block(2) <= slices_files.size()) {
			FILE* fp = slices_files.at(block(2));
			if (fp != NULL) {
				long offset = (long)(bi * block_size);

				// checks size of the block slice
				fseek(fp, 0, SEEK_END);
				off_t file_length = ftell(fp);

				if (file_length >= offset) 
					fseek(fp, offset, SEEK_SET); 
				else {
					std::cout << "failed to locate pointer for block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slice " << ss.str() << std::endl;
					return "";
				}

				bool result = fread(block_ptr, block_size, 1, fp) == 1;

				if (result) 
					return block_ptr;
				else {
					std::cout << "failed to read block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
					return "";
				}
			}
		}
		std::cout << "slice isn't open: " << ss.str() << " for block :" << block(0) << ", " << block(1) << ", " << block(2) << std::endl;
		return "";
	
		/// to ensure correct block being loaded
		//return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	} catch (...) {
		std::cout << "unhandled exception thrown at block retrieval" << std::endl;
		return "";
	}
}

void cache_manager::fifo_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {
	
	// note: because of the request_blocks filtering of blocks not in cache, we assume this function will never be called with a block already present in cpu cache
	
	//lock mutexes
	cpu_blocks_queue_lock.lock();
	cpu_cache_lock.lock();
	
	if (cpu_blocks_queue.size() == cpu_cache_size_blocks) {

		ivec3 last = cpu_blocks_queue.back();
		cpu_blocks_queue.pop_back();
		cpu_block_cache_map.erase(last);
	}

	restart_lock.lock();
	bool slices_available = slices_path != "";
	restart_lock.unlock();

	cpu_block_cache_map[block] = slices_available ? retrieve_block(block, nr_blocks, block_size, df_dim) : "";
	cpu_blocks_queue.push_front(block);

	char* block_ptr = cpu_block_cache_map[block];

	//unlock mutexes
	cpu_cache_lock.unlock();
	cpu_blocks_queue_lock.unlock();
	
	gpu_fifo_refer(block, block_ptr);
}

void cache_manager::gpu_fifo_refer(ivec3& block, char* block_ptr) {

	//lock mutexes
	gpu_blocks_queue_lock.lock();
	gpu_cache_lock.lock();

	// cache is full
	if (gpu_blocks_queue.size() == gpu_cache_size_blocks) {
		ivec3 last = gpu_blocks_queue.back();
		gpu_blocks_queue.pop_back();
		gpu_block_cache_map.erase(last);
	}

	gpu_block_cache_map[block] = "";
	gpu_blocks_queue.push_front(block);
	// call upload to gpu should go here

	//unlock mutexes
	gpu_cache_lock.unlock();
	gpu_blocks_queue_lock.unlock();
}


//! write a block located at \c data_ptr to a tiff file
bool cache_manager::write_tiff_block(std::string& file_name, char* data_ptr, vec3 df_dim, std::string& options)
{
	// setup output format 
	cgv::data::data_format df;
	df.set_component_format(cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB));
	df.set_width( (int) df_dim(0));
	df.set_height( (int) df_dim(1));

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
