#pragma once
#include <queue>
#include <vector>
#include <iterator>
#include <string>
#include <mutex>
#include <cgv/render/drawable.h>
#include <cgv/utils/file.h>
#include <vol_data/volume.h>
#include <unordered_map>
#include <unordered_set>

class volume_slicer;

// hash function for vector based on the one from boost libraries. Used for maps with ivec3 keys
template <typename Container>
struct container_hash {
	std::size_t operator()(Container const& c) const {
		std::size_t seed = c.size();
		for (auto& i : c) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

class cache_manager {

/**@ general */
//@{
	private:
		// reference to the main program
		volume_slicer &vs;

		// used also by the vol slicer
		typedef cgv::math::fvec<float, 3> vec3;
		typedef cgv::math::fvec<int, 3> ivec3;
		typedef cgv::media::axis_aligned_box<float, 3> box3;

		// folder with the slices files to load
		std::string slices_path;

		// limit the size of the cache structures 
		int cpu_cache_size_blocks;
		int gpu_cache_size_blocks;

		// limit the amount of worker threads for block retrieval
		int thread_amount;
		int thread_limit;
		std::mutex thread_report;

		// cache policy selection (0 = custom, 1 = lru, 2 = fifo)
		int selected_cache_policy;
	
		bool signal_kill = false;
		bool signal_restart;
		bool process_running;

		std::vector<std::string> slices_files_names;
		std::vector<FILE*> slices_files;

		// receives the stream of intersected blocks
		std::vector<ivec3> disk_queue_blocks;
		std::vector<ivec3> cache_queue_blocks;
	
		// general mutexes
		std::mutex blocks_in_progress_lock;
		std::mutex slices_files_lock;
		std::mutex restart_lock;
		std::mutex kill_lock;

	public:
		///constructor receives volume_slicer reference 
		/// to use the post_redraw method when a batch of blocks is ready
		cache_manager(volume_slicer &f);

		// starts infinite loop to receive requests and handle them
		void cache_manager::set_block_folder(std::string);

		// applies fopen to all files in folder
		void cache_manager::open_slices_files(std::string folder_path);
		
		// applies fclose to all slice files open
		void cache_manager::close_slices_files();

		// starts infinite loop to receive requests and handle them
		void cache_manager::init_listener();

		// handles termination of the listener when the main execution triggers it
		void cache_manager::kill_listener();

		// endless loop to check for new blocks to load
		void cache_manager::main_loop();

		// fills the queues with unloaded blocks
		void cache_manager::request_blocks(std::vector<ivec3> blocks_batch);
	
		// processes blocks that need to be retrieved from disk to cpu or cpu to gpu
		void cache_manager::retrieve_blocks_in_plane();

		// loads a block from disk
		char* cache_manager::retrieve_block(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

		// for testing purposes
		bool cache_manager::write_tiff_block(std::string& file_name, char* data_ptr, vec3 df_dim, std::string& options);
//@}


/**@ common cache structures*/
//@{
	private:
		// cpu cache queue
		std::list<ivec3> cpu_blocks_queue;
		// cpu level cache
		std::unordered_map < ivec3, char*, container_hash<ivec3 >> cpu_block_cache_map;
		// cpu mutex
		std::mutex cpu_cache_lock;

		// cpu cache queue
		std::list<ivec3> gpu_blocks_queue;
		// gpu level cache mirror (typically smaller than cpu)
		std::unordered_map < ivec3, char*, container_hash<ivec3 >> gpu_block_cache_map;
		// gpu cache mutex
		std::mutex gpu_cache_lock;
	public:

		// updates cache and the manager queue with fifo policy
		void cache_manager::cpu_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

		// updates the gpu 
		void cache_manager::gpu_refer(ivec3& block, char* block_ptr);

//@}


/**@fifo cache policy related*/
//@{
	private:

	public:
		// updates cache and the manager queue with fifo policy
		void cache_manager::cpu_fifo_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

		// updates the gpu 
		void cache_manager::gpu_fifo_refer(ivec3& block, char* block_ptr);
//@}

/**@lru cache policy related*/
//@{

	private:
	
		//save iterators from list for faster lookup
		std::unordered_map < ivec3, std::list<ivec3>::iterator, container_hash<ivec3 >> manager_cpu_blocks_queue;
		std::unordered_map < ivec3, std::list<ivec3>::iterator, container_hash<ivec3 >> manager_gpu_blocks_queue;

	public:
		// updates an lru controlled list with the size limit
		void cache_manager::cpu_lru_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

		// updates the gpu 
		void cache_manager::gpu_lru_refer(ivec3& block, char* block_ptr);

//@}


/**@configuration testing and time recording*/
//@{

	private:
		int tests_executed = 0;
		double duration = 0.0;

		std::vector<vec3> tests;
		std::vector<std::string> block_configs;
		std::vector<double> durations;

		std::mutex time_record_lock;

	public:
		// 0 = inactive, 1 = active
		int test_mode = 1; 
		int inactivity_counter = 0;
		int inactivity_limit = 1000;
		bool automatic_test_ongoing = false;
		std::mutex test_lock;
		// moves the slice and records loading times
		void cache_manager::test_handle();
		void cache_manager::next_test();
		void cache_manager::init_test_array();
//@}

};

