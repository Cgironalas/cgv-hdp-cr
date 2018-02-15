#include "cache_manager.h"
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <cgv/media/image/image_writer.h>

cache_manager::cache_manager() {}

void cache_manager::test(int value){
	std::cout << "value = " << value << std::endl;
}

bool cache_manager::retrieve_block(const std::string& input_path, const std::string& output_path, ivec3& block, const ivec3 nr_blocks, const size_t block_size, const vec3 df_dim) {

	std::cout << "Trying to get block at: " << block(0) << ", " << block(1) << ", " << block(2) << ", " << std::endl;

	std::stringstream ss;
	ss << input_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(2) << ".bvx";

	// compute block index and pointer to block data
	unsigned bi = block(1) * nr_blocks(0) + block(0);
	char* block_ptr = new char[block_size];

	// throws exception at fread in some points 
	FILE* fp = fopen(ss.str().c_str(), "rb");
	fseek(fp, bi * block_size, SEEK_SET);
	try {

		if (fread(block_ptr, block_size, 1, fp) == 1) {
			fclose(fp);
			return true;
		}
		else
			return false;

		return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, df_dim, "");

	}
	catch (const std::exception& e) {
		std::cout << "Exception thrown: " << e.what() << std::endl;
		return false;
	}

}

void cache_manager::retrieve_blocks_in_plane(const vec3 dimensions, const ivec3 block_dimensions, const ivec3 overlap) {

	//intersected_blocks_lock.lock();

	if (intersected_blocks.empty()) {
		//intersected_blocks_lock.unlock();
	}
	else {
		std::vector<ivec3> blocks_to_retrieve = intersected_blocks;
		//intersected_blocks_lock.unlock();

		ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
		size_t block_size = (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size()) * (block_dimensions(2) + overlap(2));
		ivec3 df_dim(block_dimensions(0) + overlap(0), block_dimensions(1) + overlap(1), block_dimensions(2) + overlap(2));

		for (auto block : blocks_to_retrieve) {
			std::string input_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices";
			std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";

			if (!retrieve_block(input_path, output_path, block, nr_blocks, block_size, df_dim))
				std::cout << "failed to retrieve block at: " << block(0) << ", " << block(1) << ", " << block(2) << " from block slices" << std::endl;
		}
	}
}

//! write a block located at \c data_ptr to a tiff file
bool cache_manager::write_tiff_block(const std::string& file_name, const char* data_ptr, const vec3 df_dim, const std::string& options)
{
	// setup output format 
	cgv::data::data_format df;
	df.set_component_format(cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB));
	df.set_width(df_dim(0));
	df.set_height(df_dim(1));

	cgv::media::image::image_writer iw(file_name + ".tif");
	iw.multi_set(options);

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
