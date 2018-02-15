#pragma once
#include <queue>
#include <vector>
#include <string>
#include <mutex>
#include <vol_data/volume.h>

class cache_manager {
	int thread_amount;
	std::queue < std::vector<int>> blocks_to_load;

public:
	/// copied from volume_slicer
	typedef cgv::math::fvec<float, 3> vec3;
	typedef cgv::math::fvec<int, 3> ivec3;
	typedef cgv::media::axis_aligned_box<float, 3> box3;
	
	std::mutex intersected_blocks_lock;

	//constructor
	cache_manager();
	void test(int);

	void cache_manager::update_intersected_blocks_to_retrieve(std::vector<ivec3> pIntersected_blocks);
	void cache_manager::update_volume_dim(std::vector<ivec3> pIntersected_blocks);

	bool cache_manager::retrieve_block(const std::string& input_path, const std::string& output_path, ivec3& block, const ivec3 nr_blocks, const size_t block_size, const vec3 df_dim);
	void cache_manager::retrieve_blocks_in_plane(const vec3 dimensions, const ivec3 block_dimensions, const ivec3 overlap);
	bool cache_manager::write_tiff_block(const std::string& file_name, const char* data_ptr, const vec3 df_dim, const std::string& options);

private: 
	std::vector<ivec3> intersected_blocks;
};

