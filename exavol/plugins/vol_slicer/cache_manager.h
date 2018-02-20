#pragma once
#include <queue>
#include <vector>
#include <string>
#include <mutex>
#include <cgv/render/drawable.h>
#include <vol_data/volume.h>
#include <unordered_map>


class volume_slicer;

///unordered map functions like the ones from boost libraries
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

class cache_manager  {
	
	volume_slicer &vs;

	typedef cgv::math::fvec<float, 3> vec3;
	typedef cgv::math::fvec<int, 3> ivec3;
	typedef cgv::media::axis_aligned_box<float, 3> box3;

	// limit the size of the cache structure. 
	int cache_size_blocks;

	// limit the amount of worker threads for block retrieval
	int thread_amount;
	int thread_limit;

	bool process_running;

	// receives the stream of intersected blocks.
	std::list<ivec3> blocks_in_progress;
	std::list<ivec3> intersected_blocks_queue;
	// manages gpu cache blocks
	std::unordered_map < ivec3, std::list<ivec3>::iterator, container_hash<ivec3 >> block_cache_map;


	// mutexes
	std::mutex blocks_cache_lock;
	std::mutex blocks_in_progress_lock;
	std::mutex intersected_blocks_queue_lock;
	std::mutex thread_report;

	

public:
	///constructor receives volume_slicer reference 
	/// to use the post_redraw method when a batch of blocks is ready
	cache_manager(volume_slicer &f);

	// starts infinite loop to receive requests and handle them
	void cache_manager::init_listener();

	void cache_manager::main_loop();

	// starts queue consumption of the blocks to refresh screen
	void cache_manager::process_blocks(std::vector<ivec3> pIntersected_blocks);

	// thread safely updates the list with the reported blocks from the UI.
	void cache_manager::update_intersected_blocks_queue();
	
	// updates an lru controlled list with the size limit
	void cache_manager::lru_refer(ivec3 block);

	// block loading functions
	bool cache_manager::retrieve_block(ivec3& block, const ivec3& nr_blocks, const size_t& block_size, const vec3& df_dim);
	void cache_manager::retrieve_blocks_in_plane();
	bool cache_manager::write_tiff_block(const std::string& file_name, const char* data_ptr, const vec3 df_dim, const std::string& options);
};

