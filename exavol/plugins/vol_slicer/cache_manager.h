#pragma once
#include <queue>
#include <vector>
#include <string>
#include <mutex>
#include <cgv/render/drawable.h>
#include <vol_data/volume.h>
#include <unordered_map>
#include <unordered_set>

class volume_slicer;

// hash function for vector like the one from boost
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

	// reference to the main program
	volume_slicer &vs;

	typedef cgv::math::fvec<float, 3> vec3;
	typedef cgv::math::fvec<int, 3> ivec3;
	typedef cgv::media::axis_aligned_box<float, 3> box3;

	std::string slices_path;
	ivec3 slices_dimensions;

	// limit the size of the cache structures 
	int cpu_cache_size_blocks;
	int gpu_cache_size_blocks;

	// limit the amount of worker threads for block retrieval
	int thread_amount;
	int thread_limit;

	bool signal_restart;
	bool process_running;

	// receives the stream of intersected blocks
	std::vector<ivec3> disk_queue_blocks;
	std::vector<ivec3> cache_queue_blocks;

	// manages cpu cache blocks
	std::list<ivec3> cpu_blocks_queue;
	// cpu level cache
	std::unordered_map < ivec3, char*, container_hash<ivec3 >> cpu_block_cache_map;
	// cpu mutexes
	std::mutex cpu_cache_lock;
	std::mutex cpu_blocks_queue_lock;

	// manages gpu cache blocks
	std::list<ivec3> gpu_blocks_queue;
	// gpu level cache mirror (typically smaller than cpu)
	std::unordered_map < ivec3, char*, container_hash<ivec3 >> gpu_block_cache_map;
	// gpu mutexes
	std::mutex gpu_cache_lock;
	std::mutex gpu_blocks_queue_lock;

	// general mutexes
	std::mutex blocks_in_progress_lock;
	std::mutex thread_report;
	std::mutex slices_files_lock;
	std::mutex restart_lock;

public:
	///constructor receives volume_slicer reference 
	/// to use the post_redraw method when a batch of blocks is ready
	cache_manager(volume_slicer &f);

	// starts infinite loop to receive requests and handle them
	void cache_manager::set_block_folder(std::string);

	// starts infinite loop to receive requests and handle them
	void cache_manager::init_listener();

	// endless loop to check for new blocks to load
	void cache_manager::main_loop();

	// starts queue traversal of the blocks to refresh screen
	void cache_manager::request_blocks(std::vector<ivec3> blocks_batch);
	
	// updates an lru controlled list with the size limit
	void cache_manager::fifo_refer(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

	// updates the gpu 
	void cache_manager::gpu_fifo_refer(ivec3& block, char* block_ptr);

	// spawns threads to retrieve all blocks
	void cache_manager::retrieve_blocks_in_plane();

	// loads a block for the lru_refer
	char* cache_manager::retrieve_block(ivec3& block, ivec3& nr_blocks, size_t& block_size, vec3& df_dim);

	// for testing purposes
	bool cache_manager::write_tiff_block(std::string& file_name, char* data_ptr, vec3 df_dim, std::string& options);
};

