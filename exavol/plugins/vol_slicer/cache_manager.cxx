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

// Constructor
//
cache_manager::cache_manager(volume_slicer &f) : vs(f) {

	slices_path = "";
	cpu_cache_size_blocks = 20000; //could be calculated from the CPU capacity if possible
	gpu_cache_size_blocks = 5000; //could be calculated from the GPU capacity if possible

	// extra validation for threaded approaches
	thread_limit = 8; 
	
	selected_cache_policy = 1;
	std::cout << "selected cache policy " << selected_cache_policy << std::endl;

	process_running = false;
	signal_restart = false;
	test_mode = 1;
	if (test_mode == 1) {
		init_test_array();
	}
}

// Sets up a thread with the main loop, separates cache manager execution from vol slicer. 
//
void cache_manager::init_listener() {

	thread_report.lock();
	thread_amount += 1;
	std::thread t([&]() { main_loop(); });
	t.detach();
	thread_report.unlock();
}

// Closes the listener when program execution is finished
//
void cache_manager::kill_listener() {
	kill_lock.lock();
	signal_kill = true;
	kill_lock.unlock();
}

/* ================================== */
// global structures (files array and path from configuration files) : 
/* ================================== */

// Sets the required variables and mutex locks, unloads and closes active files and replaces them with the new files. 
//
void cache_manager::set_block_folder(std::string folder_path) {
	restart_lock.lock();
	gpu_cache_lock.lock();
	cpu_cache_lock.lock();
	blocks_in_progress_lock.lock();

	slices_path = folder_path;

	close_slices_files();
	if (slices_path != "") {
		open_slices_files(folder_path + "_x/");
		open_slices_files(folder_path + "_y/");
		open_slices_files(folder_path + "_z/");
	}
	
	manager_cpu_blocks_queue.clear();
	cpu_blocks_queue.clear();
	cpu_block_cache_map.clear();

	manager_gpu_blocks_queue.clear();
	gpu_blocks_queue.clear();
	gpu_block_cache_map.clear();

	disk_queue_blocks.clear();
	cache_queue_blocks.clear();
	
	restart_lock.unlock();
	cpu_cache_lock.unlock();
	gpu_cache_lock.unlock();
	blocks_in_progress_lock.unlock();
}

// Opens all files from the respective block folder path and keeps them in the "slices_files" array to avoid overhead of 
// opening and closing the files everytime a block must be retrieved. 
void cache_manager::open_slices_files(std::string folder_path) {
	void* handle = cgv::utils::file::find_first(folder_path + "*.*");
	std::cout << "\n\n opening slices files at: '" << folder_path << "' \n\n" << std::endl;

	bool success = false;
	while (handle) {

		if (!cgv::utils::file::find_directory(handle)) {
			std::string file_name = folder_path + cgv::utils::file::find_name(handle);

			std::cout << file_name << std::endl;
			slices_files_names.push_back(file_name);

			FILE* fp = fopen(file_name.c_str(), "rb");
			if (fp != NULL) {
				slices_files.push_back(fp);
				success = true;
			} else {
				success = false;
				slices_path = folder_path;
				break;
			}
		}
		handle = cgv::utils::file::find_next(handle);
	}

	if (success) 
		std::cout << "\n\n successfully opened all files at: '" << folder_path << "' \n\n" << std::endl;
	else 
		std::cout << "\n\n something went wrong when opening the files at: '" << folder_path <<  "'  please check your configuration file.\n\n" << std::endl;
}

// Prevents unclosed files after the program's execution. 
//
void cache_manager::close_slices_files() {

	slices_files_names.clear();
	for (auto it : slices_files)
		fclose(it);
	slices_files.clear();
}

/* ================================== */
// main program methods : 
/* ================================== */

// The main loop of the cache manager listens to requests from the vol slicer, and reports successfully loaded block batches when they are finished. 
// If the test mode is activated, the loop also calls the test utils to force the vol slicer into processing different slice configurations
void cache_manager::main_loop() {
	
	while (true) {

		// handle termination from the vol slicer when it didnt die on close
		kill_lock.lock();
		if (signal_kill) {
			kill_lock.unlock();
			close_slices_files();
			break;
		}
		kill_lock.unlock();
		
		blocks_in_progress_lock.lock();
		bool something_to_do = !disk_queue_blocks.empty() || !cache_queue_blocks.empty();
		blocks_in_progress_lock.unlock();

		thread_report.lock();
		if (thread_amount > 1) 
			std::cout << "SOMETHING IS WRONG" << std::endl;
		thread_report.unlock();

		restart_lock.lock();
		bool start_retrieval = thread_amount < thread_limit && something_to_do && !process_running && slices_path != "";
		restart_lock.unlock();

		if (start_retrieval) {
			process_running = true;
			retrieve_blocks_in_plane();
			test_handle();
		} else {
			if (test_mode == 1 && automatic_test_ongoing == true) {
				inactivity_counter++;
				std::cout << "inactivity..." << inactivity_counter <<" restarting at: " << inactivity_limit << std::endl;
				if (inactivity_counter == inactivity_limit) {
					test_handle();
				}
			}
		}
	}
}

// This method processes the batch of blocks depending of their status: loaded in cache (cpu or gpu) or unloaded. 
// It has mutex locks on some structures for thread implementations. 
void cache_manager::request_blocks(std::vector<ivec3> blocks_batch) {
	if (selected_cache_policy == 0) {
		blocks_in_progress_lock.lock();
		bool new_blocks = false;
		for (auto it : blocks_batch) {

			cpu_cache_lock.lock();
			bool block_in_cache = cpu_block_cache_map.find(it) != cpu_block_cache_map.end();
			cpu_cache_lock.unlock();

			if (!block_in_cache) {
				disk_queue_blocks.push_back(it);
				new_blocks = true;

			}
			else {
				gpu_cache_lock.lock();
				bool block_in_gpu = gpu_block_cache_map.find(it) != gpu_block_cache_map.end();
				gpu_cache_lock.unlock();

				if (!block_in_gpu) {
					cache_queue_blocks.push_back(it);
					new_blocks = true;
				}
			}
		}

		restart_lock.lock();
		signal_restart = true;
		restart_lock.unlock();

		blocks_in_progress_lock.unlock();
	} else { // plain lru or fifo

		blocks_in_progress_lock.lock();
		
		for (auto it : blocks_batch) {
			disk_queue_blocks.push_back(it);
		}

		restart_lock.lock();
		signal_restart = true;
		restart_lock.unlock();

		blocks_in_progress_lock.unlock();
	}
}

/* ================================== */
// block management and processing : 
/* ================================== */

// This method handles the requested batch of blocks in a queue and checks for interruptions caused by user interaction and subsquent change in the required blocks. 
//
int disk_read_count = 0;
void cache_manager::retrieve_blocks_in_plane() {

	if (test_mode == 1) {
		automatic_test_ongoing = true;
	}

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
	size_t block_size = (long) ( (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*(cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size())) * (block_dimensions(2) + overlap(2)) );
	vec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));

	// cancel flag for cpu -> gpu loading: not used right now, but could be helpful if cpu -> gpu takes long
	bool cancelled = false;

	std::cout << "\n\n total blocks to retrieve: " << blocks_batch_frame.size() + cached_batch_frame.size() << " \n" << std::endl;

	std::string sources[2] = { "cpu", "gpu" };
	std::vector<ivec3> frame_sources[2] = { blocks_batch_frame, cached_batch_frame };

	// plain fifo or lru never use cached_batch_frame
	int s = frame_sources[0].empty() ? 1 : 0;
	int end_batch = frame_sources[1].empty() ? 1 : 2;

	disk_read_count = 0;

	std::clock_t start = std::clock();
	for (s; s < end_batch; s++) {
		std::cout << "\nsize batch for " << sources[s] << " cache: " << frame_sources[s].size() << " \n" << std::endl;

		for (int i = 0; i < frame_sources[s].size(); i++) {

			// premature termination ========= 
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
				std::cout << "cancelled remaining " << frame_sources[s].size() - i << " blocks for: " << sources[s] << std::endl;
				break;
			}
			// ==============================


			if (s == 0) {
				cpu_refer(blocks_batch_frame[i], nr_blocks, block_size, df_dim);
			} else {
				cpu_cache_lock.lock();
				std::shared_ptr<char> block_ptr = cpu_block_cache_map[cached_batch_frame[i]];
				cpu_cache_lock.unlock();

				gpu_refer(cached_batch_frame[i], block_ptr);
			}

			int ten_percentage = int(ceil(0.1*frame_sources[s].size()));
			ten_percentage = ten_percentage == 0 ? 1 : ten_percentage;

			if (i % ten_percentage == 0 || i % (frame_sources[s].size()-1) == 0) {
				std::cout << "blocks processed for " << sources[s] << ": " << i + 1 << ", total read from disk : " << disk_read_count << std::endl;
			}
		}
	}

	if (!cancelled) {
		std::cout << "finished batch of " << blocks_batch_frame.size() + cached_batch_frame.size() << " blocks. Currently  " << cpu_block_cache_map.size() << " blocks in cpu and " << gpu_block_cache_map.size() << " in gpu " << std::endl;
		
		// update the frame with the new block textures loaded.
		vs.post_redraw();
		
		if (test_mode == 1) {
			time_record_lock.lock();
			duration = static_cast<double>(std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
			durations.push_back(duration);
			time_record_lock.unlock();
			std::cout << " duration of refers" << duration << "\n\n" << std::endl;
		}
	}

	process_running = false;
}

// util method for list sorting
//
static bool abs_compare(float a, float b) { return (std::abs(a) < std::abs(b)); }

// This method computes the location of the block in the open files array, then retrieves it as a shared_ptr<char>
//
std::shared_ptr<char> cache_manager::retrieve_block(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {
	
	//std::cout << "retrieving block at: " << block(0) << ", " << block(1) << ", " << block(2) << ", " << std::endl;
	try {

		vec3 world_slice_normal = vs.world_from_texture_normals(vs.slice_normal_tex);
		int d = (int)std::distance(std::begin(world_slice_normal), std::max_element(std::begin(world_slice_normal), std::end(world_slice_normal), abs_compare));

		// compute block index and pointer to block data
		std::string path_to_orientation;
		unsigned bi;
		unsigned slices_offset;

		if (d == 0) {
			slices_offset = 0; //x are opened first
			path_to_orientation = "_x/";
			bi = block(1) * nr_blocks(2) + block(2);
		}

		if (d == 1) {
			slices_offset = nr_blocks(0); //ignore files from _x
			path_to_orientation = "_y/";
			bi = block(2) * nr_blocks(1) + block(0);
		}

		if (d == 2) {
			slices_offset = nr_blocks(0) + nr_blocks(1); //ignore files from _x and _y
			path_to_orientation = "_z/";
			bi = block(0) * nr_blocks(1) + block(1);
		}

		std::stringstream ss;
		ss << slices_path << path_to_orientation << "level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(d) << ".bvx";

		char* block_ptr = new char[block_size];

		// throws exception at fread in some points 
		if (block(d) <= slices_files.size()) {
			FILE* fp = slices_files.at(slices_offset + block(d));
			if (fp != NULL) {
				long offset = (long)(bi * block_size);

				// checks size of the block slice
				fseek(fp, 0, SEEK_END);
				off_t file_length = ftell(fp);

				if (file_length >= offset) 
					fseek(fp, offset, SEEK_SET); 
				else {
					std::cout << "failed to locate pointer for block at: [" << block(0) << ", " << block(1) << ", " << block(2) << "] from block slice " << ss.str() << std::endl;
					std::cout << "block size " << block_size << ", block slice size: " << file_length << ", bi: " << bi << ", nr_blocks: [" << nr_blocks(0) << ", " << nr_blocks(1) << ", " << nr_blocks(2) <<"]\n" << std::endl;
					return NULL;
				}

				bool result = fread(block_ptr, block_size, 1, fp) == 1;

				if (result) {
					disk_read_count++;
					return std::shared_ptr<char>(block_ptr);
				} else {
					std::cout << "failed to read block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
					return NULL;
				}
			}
		}
		std::cout << "slice isn't open: " << ss.str() << " for block :" << block(0) << ", " << block(1) << ", " << block(2) << std::endl;
		return NULL;
	
		/// to ensure correct block being loaded
		//return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	} catch (...) {
		std::cout << "unhandled exception thrown at block retrieval" << std::endl;
		return NULL;
	}
}

//Calls the respective method depending on the selected cache policy selected, for CPU
//
void cache_manager::cpu_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {
	
	if (selected_cache_policy == 0 || selected_cache_policy == 2) {
		cpu_fifo_refer(block, nr_blocks, block_size, df_dim);
	} else {
		cpu_lru_refer(block, nr_blocks, block_size, df_dim);
	}
}

// Calls the respective method depending on the selected cache policy selected, for GPU
//
void cache_manager::gpu_refer(ivec3& block, std::shared_ptr<char> block_ptr) {
	if (selected_cache_policy == 0 || selected_cache_policy == 2) {
		gpu_fifo_refer(block, block_ptr);
	} else {
		gpu_lru_refer(block, block_ptr);
	}
}

/* ================================== */
// fifo methods :
/* ================================== */

// Manages the data structures for an FIFO queue refer in the structure that handles CPU storage of the blocks
//
void cache_manager::cpu_fifo_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {

	// note: because of the request_blocks filtering of blocks not in cache, we can assume this function will never be called with a block already present in cpu cache

	cpu_cache_lock.lock();

	if (cpu_blocks_queue.size() == cpu_cache_size_blocks) {

		ivec3 last = cpu_blocks_queue.back();
		cpu_blocks_queue.pop_back();

		cpu_block_cache_map.erase(last);
	}

	restart_lock.lock();
	bool slices_available = slices_path != "";
	restart_lock.unlock();

	cpu_block_cache_map[block] = slices_available ? retrieve_block(block, nr_blocks, block_size, df_dim) : NULL;
	cpu_blocks_queue.push_front(block);

	std::shared_ptr<char> block_ptr = cpu_block_cache_map[block];

	cpu_cache_lock.unlock();

	gpu_refer(block, block_ptr);
}

// Manages the data structures for an FIFO queue refer in the structure that mirrors storage in the GPU
//
void cache_manager::gpu_fifo_refer(ivec3& block, std::shared_ptr<char> block_ptr) {

	gpu_cache_lock.lock();

	// cache is full
	if (gpu_blocks_queue.size() == gpu_cache_size_blocks) {
		ivec3 last = gpu_blocks_queue.back();
		gpu_blocks_queue.pop_back();
		gpu_block_cache_map.erase(last);
	}

	// call upload to gpu should go here, using block_ptr
	gpu_block_cache_map[block] = NULL;
	gpu_blocks_queue.push_front(block);

	gpu_cache_lock.unlock();
}

/* ================================== */
// lru methods :
/* ================================== */

// Manages the data structures for an LRU queue refer in the structure that handles CPU storage of the blocks
//
void cache_manager::cpu_lru_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim) {

	// note: because of the request_blocks filtering of blocks not in cache, we can assume this function will never be called with a block already present in cpu cache

	cpu_cache_lock.lock();

	//not in cache
	if (manager_cpu_blocks_queue.find(block) == manager_cpu_blocks_queue.end()) {
		// cache is full
		if (cpu_blocks_queue.size() == cpu_cache_size_blocks) {

			ivec3 last = cpu_blocks_queue.back();
			cpu_blocks_queue.pop_back();
			manager_cpu_blocks_queue.erase(last);
			cpu_block_cache_map.erase(last);
		}
	} else {
		cpu_blocks_queue.erase(manager_cpu_blocks_queue[block]);
	}

	restart_lock.lock();
	bool slices_available = slices_path != "";
	restart_lock.unlock();

	if (cpu_block_cache_map.find(block) == cpu_block_cache_map.end()) {
		cpu_block_cache_map[block] = slices_available ? retrieve_block(block, nr_blocks, block_size, df_dim) : NULL;
	}

	cpu_blocks_queue.push_front(block);
	manager_cpu_blocks_queue[block] = cpu_blocks_queue.begin();
	std::shared_ptr<char> block_ptr = cpu_block_cache_map[block];

	cpu_cache_lock.unlock();

	gpu_refer(block, block_ptr);
}

// Manages the data structures for an LRU queue refer in the structure that mirrors storage in the GPU
// 
void cache_manager::gpu_lru_refer(ivec3& block, std::shared_ptr<char> block_ptr) {

	gpu_cache_lock.lock();

	//not in cache
	if (manager_gpu_blocks_queue.find(block) == manager_gpu_blocks_queue.end()) {
		// cache is full
		if (gpu_blocks_queue.size() == gpu_cache_size_blocks) {
			ivec3 last = gpu_blocks_queue.back();
			gpu_blocks_queue.pop_back();
			manager_gpu_blocks_queue.erase(last);
			gpu_block_cache_map.erase(last);
		}
	} else {
		gpu_blocks_queue.erase(manager_gpu_blocks_queue[block]);
	}
	
	// call upload to gpu should go here, using block_ptr
	gpu_block_cache_map[block] = NULL;
	gpu_blocks_queue.push_front(block);
	manager_gpu_blocks_queue[block] = gpu_blocks_queue.begin();

	gpu_cache_lock.unlock();
}

/* ================================== */
// testing utils :
/* ================================== */

// Writes a block located at \c data_ptr to a tiff file
// Requires an output path. This is not used currently but it's useful for debugging. 
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

// Fills "tests" with different slice normal values. Each one of these is then fed to the vol slicer for block retrieval
// The tests are both random and deliberate, trying to emulate some possible user behaviors
void cache_manager::init_test_array() {
	
	double values[3];
	tests.push_back(vec3(0, 0, 1)); //default configuration
	tests.push_back(vec3(0, 1, 0));
	tests.push_back(vec3(1, 0, 0));

	// creates random positions
	for (int i = 0; i < 30; i++) {
		double x = std::rand() % 100;
		double y = std::rand() % 100;
		double z = std::rand() % 100;
		double sum = (x + y + z);
		x /= sum;
		y /= sum;
		z /= sum;
		tests.push_back(vec3(x, y, z));
	}
	
	// creates rotations in all axis
	for (int i = 0; i < 3; i++) {
		for (double r = 0; r < 30; r++) {

			values[i] = 0;
			values[(i + 1) % 3] = r / 10;
			values[(i + 2) % 3] = 1 - (r / 10);
			
			tests.push_back(vec3(values[0], values[1], values[2]));

			values[(i + 1) % 3] += 0.005;
			values[(i + 2) % 3] -= 0.005;

			tests.push_back(vec3(values[0], values[1], values[2]));

			values[(i + 1) % 3] -= 0.01;
			values[(i + 2) % 3] += 0.01;

			tests.push_back(vec3(values[0], values[1], values[2]));
		}
	}


	for (auto it : tests) {
		std::cout << "[" << it[0] << "," << it[1] << "," << it[2] << "]" << std::endl;
	}

	//block configs (all tests in the "tests" are done to each configuration)
	block_configs.push_back("D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/male/fullbody/slices/multiple/block_config64.bsdc");
	block_configs.push_back("D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/male/fullbody/slices/multiple/block_config32.bsdc");
}

// Records the test duration in the position x of "durations" array, where x can be mapped to 
// the position of the respective test in the "tests" array. 
void::cache_manager::test_handle() {

	inactivity_counter = 0;
	int currently_executed = tests_executed + (tests.size()*(2 - block_configs.size()));
	std::cout << "[" << durations.size() << "vs" << currently_executed << "]" << std::endl;

	// empty test detected by inactivity counter
	if (durations.size() == currently_executed ) {
		std::cout << "[ === inserting empty test time == ]" << std::endl;
		time_record_lock.lock();
		duration = 0;
		durations.push_back(duration);
		time_record_lock.unlock();
	}


	if (durations.size() > tests_executed) {
		next_test();
	} else {
		automatic_test_ongoing = false;
		std::cout << "[ thread desynchronization ... ]" << std::endl;
	}

	std::cout << std::fixed << " duration of refers" << duration << std::endl;
}

// Consumes the current position in the "tests" array and moves into the next configuration. 
// If the tests are over, it moves on to the next block configuration and resets the array. 
void cache_manager::next_test() {
	
	tests_executed += 1;
	if (tests_executed >= tests.size()) {

		if (block_configs.size() == 0) { //testing is finished
			automatic_test_ongoing = false; //kill automatic execution

			std::cout << " \n\nTEST results: \n " << std::endl;

			for (int j = 0; j < 3; j++) { //repeats tests for the 3 available block config 
				double avg_duration = 0;
				std::cout << " \n\n new block size \n\n" << std::endl;

				for (int i = 0; i < tests.size(); i++) {
					std::cout << " test: [" << tests[i] << "] duration : " << durations[i + (j*tests.size())] << std::endl;
					avg_duration += durations[i + (j*tests.size())];
				}

				avg_duration /= tests.size();
				std::cout << " average duration of tests " << avg_duration << std::endl;
			}
			
		} else { //move on to next block config

			vs.block_config = block_configs[block_configs.size() - 1];
			block_configs.pop_back();
			vs.on_set(&vs.block_config);
			
			tests_executed = 0;
			vs.slice_normal_tex = tests[tests_executed];

			vs.on_set(&vs.slice_normal_tex[0]);
			vs.on_set(&vs.slice_normal_tex[1]);
			vs.on_set(&vs.slice_normal_tex[2]);
		}
	} else { //change slice_normal values to those of new test

		std::cout << "automatically executing next test [" << tests[tests_executed] <<"]" << std::endl;
		vs.slice_normal_tex = tests[tests_executed];

		vs.on_set(&vs.slice_normal_tex[0]);
		vs.on_set(&vs.slice_normal_tex[1]);
		vs.on_set(&vs.slice_normal_tex[2]);
		
	}
}
