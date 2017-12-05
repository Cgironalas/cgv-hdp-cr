#include <iostream>
#include <cmath>
#include <cgv/utils/file.h>
#include <3rd/zlib/zlib.h>
#include <cgv/utils/dir.h>
#include <cgv/utils/scan.h>
#include <cgv/media/image/image_writer.h>
#include <cgv/media/image/image_reader.h>
#include <cgv/media/image/image_proc.h>
#include <vol_data/volume_io.h>

/// structure to hold the information needed to build a hierarchy of blocks
struct block_generation_info : public volume_info
{
	/// size of block in voxel, for example (16,16,16) or (32,32,16)
	volume::dimension_type block_dimensions;
	/// overlap between adjacent blocks in voxel, typically (1,1,1)
	volume::dimension_type overlap;
	/// for each dimension a factor for downsampling, typically (2,2,2)
	volume::dimension_type subsampling_factor;
	//! hierarchy level, when to start downsampling, (0,0,3) means for example that x and y 
	/*! are subdivided from the start and z only starting with generation of the 3rd hierarchy level*/
	volume::dimension_type subsampling_offset;

	/// return the size of a block in bytes
	size_t get_block_size() const {
		return get_block_slice_size() * (block_dimensions(2) + overlap(2));
	}
	/// return the size of a slice of a block in bytes
	size_t get_block_slice_size() const {
		return (block_dimensions(1) + overlap(1))*get_block_line_size();
	}
	/// return the size of a line of a block in bytes
	size_t get_block_line_size() const {
		return (block_dimensions(0) + overlap(0))*get_block_pixel_size();
	}
	/// return the size of a pixel in a slice of a block in bytes
	size_t get_block_pixel_size() const {
		return cgv::data::component_format(type_id, components).get_entry_size();
	}
};

/// init the block generation structure for the male visible human data set
void init_to_visible_human_male(block_generation_info& bgi)
{
	bgi.type_id = cgv::type::info::TI_UINT8;
	bgi.components = cgv::data::CF_RGB;
	bgi.dimensions.set(4096, 2700, 1878);
	bgi.extent.set(4096 * 0.114f, 2700 * 0.114f, 1878 * 1.0f);
	bgi.block_dimensions.set(16, 16, 16);
	bgi.overlap.set(1, 1, 1);
	bgi.subsampling_factor.set(2, 2, 2);
	bgi.subsampling_offset.set(0, 0, 3);
}

/// init the block generation structure for the female visible human data set
void init_to_visible_human_female(block_generation_info& bgi)
{
	bgi.type_id = cgv::type::info::TI_UINT8;
	bgi.components = cgv::data::CF_RGB;
	bgi.dimensions.set(2048, 1216, 5189);
	bgi.extent.set(2048 * 0.33f, 1216 * 0.33f, 5189 * 0.33f);
	bgi.block_dimensions.set(16, 16, 16);
	bgi.overlap.set(1, 1, 1);
	bgi.subsampling_factor.set(2, 2, 2);
	bgi.subsampling_offset.set(0, 0, 0);
}

/// currently only prints out information about each level in a block hierarchy
void build_level_infos(const block_generation_info& bgi)
{
	cgv::type::uint64_type original_size = cgv::type::uint64_type(bgi.dimensions[0]) * bgi.dimensions[1] * bgi.dimensions[2] * 3;
	std::cout << "original data: <" << bgi.dimensions << "> = " << original_size << " = " << double(original_size) / (1024 * 1024 * 1024) << std::endl;
	unsigned level = 0;
	volume::dimension_type level_dimensions = bgi.dimensions;
	volume::dimension_type nr_blocks;
	do {
		nr_blocks = volume::dimension_type(
			unsigned(ceil(float(level_dimensions[0]) / bgi.block_dimensions[0])),
			unsigned(ceil(float(level_dimensions[1]) / bgi.block_dimensions[1])),
			unsigned(ceil(float(level_dimensions[2]) / bgi.block_dimensions[2])));
		cgv::type::uint64_type block_size = cgv::type::uint64_type(bgi.block_dimensions[0] + 1) * (bgi.block_dimensions[1] + 1) * (bgi.block_dimensions[2] + 1) * 3;
		cgv::type::uint64_type nr_all_blocks = cgv::type::uint64_type(nr_blocks[0]) * nr_blocks[1] * nr_blocks[2];
		std::cout << level << ": <" << level_dimensions << ">, nr blocks: <" << nr_blocks << "> = " << nr_all_blocks * block_size << " = " << double(nr_all_blocks * block_size) / (1024 * 1024 * 1024) << std::endl;

		for (int c = 0; c < 3; ++c) {
			if (bgi.subsampling_offset[c] <= int(level))
				level_dimensions[c] = unsigned(ceil(float(level_dimensions[c]) / bgi.subsampling_factor[c]));
		}
		++level;
	} while (std::min(std::min(nr_blocks[0], nr_blocks[1]), nr_blocks[2]) >= 2);
}

/// ensure the existance of a directory
bool ensure_directory(const std::string& dir)
{
	if (cgv::utils::dir::exists(dir))
		return true;

	std::string path(cgv::utils::file::get_path(dir));
	if (path.empty())
		return false;
	if (!ensure_directory(path))
		return false;
	cgv::utils::dir::mkdir(dir);
	return true;
}

/// read an image file either from a gzipped raw file or a png file
bool read_image(const std::string& file_name, std::vector<char>& data)
{
	if (cgv::utils::to_upper(cgv::utils::file::get_extension(file_name)) == "GZ") {
		gzFile fp = gzopen(file_name.c_str(), "rb");
		int n;
		do {
			char buffer[4096];
			n = gzread(fp, buffer, 4096);
			if (n == -1)
				return false;
			size_t old_size = data.size();
			data.resize(data.size() + unsigned(n));
			std::copy(buffer, buffer + n, data.begin() + old_size);
		} while (n == 4096);
		gzclose(fp);
		return true;
	}

	else {
		cgv::data::data_format file_format;
		cgv::media::image::image_reader img_reader(file_format);

		if (!img_reader.open(file_name)) {
			return false;
		}
			
		data.resize(file_format.get_nr_bytes());
		cgv::data::data_view dv(&file_format, &data[0]);
		bool success = img_reader.read_image(dv);
		img_reader.close();
		return success;
	}
}

//! write a block located at \c data_ptr to a tiff file
/*! in options one can specify compression level through
- "" ... defaults to no compression
- "tiff_compression=none" ... same as default
- "tiff_compression=jpeg;jpeg_quality=50" ... use jpeg compression with quality value of 50 [lossless jpeg compression for quality 100]
- "tiff_compression=lzw" ... use lzw compression */
bool write_tiff_block(const std::string& file_name, const block_generation_info& bgi, const char* data_ptr, const std::string& options)
{
	// setup output format 
	cgv::data::data_format df;
	df.set_component_format(cgv::data::component_format(bgi.type_id, bgi.components));
	df.set_width(bgi.block_dimensions(0) + 1);
	df.set_height(bgi.block_dimensions(1) + 1);
	cgv::media::image::image_writer iw(file_name + ".tif");
	iw.multi_set(options);
	size_t slice_size = df.get_nr_bytes();
	for (unsigned i = 0; i<(unsigned)(bgi.block_dimensions(2) + bgi.overlap(2)); ++i) {
		cgv::data::const_data_view dv(&df, data_ptr + slice_size*i);
		if (!iw.write_image(dv)) {
			std::cerr << "could not write slice " << i << " to tiff file " << file_name << std::endl;
			return false;
		}
	}
	return iw.close();
}

/// write a single image of given component format and dimensions, detect file format from extension (bmp,png,tiff are supported)
bool write_image(const std::string& file_name, const cgv::data::component_format& cf, int w, int h, const char* data_ptr, const std::string& options = "")
{
	cgv::data::data_format df;
	df.set_component_format(cf);
	df.set_width(w);
	df.set_height(h);
	cgv::media::image::image_writer iw(file_name);
	if (!options.empty())
		iw.multi_set(options);
	cgv::data::const_data_view dv(&df, data_ptr);
	if (!iw.write_image(dv)) {
		std::cerr << "could not write image to tiff file " << file_name << std::endl;
		return false;
	}
	return iw.close();
}

/// single resolution generation of block files implemented during one sequential scan of files
bool build_blocks_from_directory(const std::string& input_path, const std::string& output_path, const block_generation_info& bgi)
{
	// ensure that output directory exists
	if (!ensure_directory(output_path)) {
		std::cout << "could not create output path <" << output_path << ">" << std::endl;
		return false;
	}

	// collect input files
	std::vector<std::string> file_names;
	void* handle = cgv::utils::file::find_first(input_path + "/*.*");
	int file_count = 0;
	while (handle) {
		//ignore files like Thumbs.db
		if (!cgv::utils::file::find_directory(handle) && file_count < bgi.dimensions(2)) {
			file_names.push_back(cgv::utils::file::find_name(handle));
			file_count += 1;
		}
		handle = cgv::utils::file::find_next(handle);
	}
	if (file_names.empty()) {
		std::cerr << "did not find files in directory <" << input_path << ">" << std::endl;
		return false;
	}
	// compute number of blocks on base levels
	volume::dimension_type level_dimensions = bgi.dimensions;
	volume::dimension_type nr_blocks(
		unsigned(ceil(float(level_dimensions(0)) / bgi.block_dimensions(0))),
		unsigned(ceil(float(level_dimensions(1)) / bgi.block_dimensions(1))),
		unsigned(ceil(float(level_dimensions(2)) / bgi.block_dimensions(2)))
	);

	// allocate data for one slice of blocks and initialize all values to 0
	size_t pixel_size = bgi.get_block_pixel_size();
	size_t block_line_size = bgi.get_block_line_size();
	size_t block_slice_size = bgi.get_block_slice_size();
	size_t block_size = bgi.get_block_size();
	size_t image_line_size = level_dimensions(0)*pixel_size;
	size_t image_size = level_dimensions(1)*image_line_size;
	std::vector<char> block_data(nr_blocks(0)*nr_blocks(1)*block_size, 0);
	size_t next_block_completion = bgi.block_dimensions(2) + bgi.overlap(2) - 1;

	// iterate over slices once
	unsigned block_z = 0;
	unsigned k = 0;
	for (unsigned z = 0; z < file_names.size(); ++z) {
		std::cout << "read slice " << z << ": <" << file_names[z] << ">" << std::endl;
		std::vector<char> data;

		if (!read_image(input_path + "/" + file_names[z], data)) {
			std::cerr << "failed to read <" << file_names[z] << ">" << std::endl;
			continue;
		}
		if (data.size() < image_size) {
			size_t old_size = data.size();
			data.resize(image_size);
			std::fill(data.begin() + old_size, data.end(), 0);
		}
		std::cout << "slice size = " << data.size() << ", block_z = " << block_z << std::endl;

		// insert data into current blocks
		for (int j = 0; j < nr_blocks(1); ++j) {
			for (int i = 0; i < nr_blocks(0); ++i) {
				// compute rectangle in input image corresponding to block 
				int x0 = i*bgi.block_dimensions(0);
				int x1 = x0 + bgi.block_dimensions(0) + bgi.overlap(0);
				if (x1 > level_dimensions(0))
					x1 = level_dimensions(0);
				size_t copy_line_size = (x1 - x0) * pixel_size;

				int y0 = j*bgi.block_dimensions(1);
				int y1 = y0 + bgi.block_dimensions(1) + bgi.overlap(1);
				if (y1 > level_dimensions(1))
					y1 = level_dimensions(1);
				const char* image_data_ptr = &data[y0 * image_line_size + x0 * pixel_size];

				// compute block index and pointer to block data
				unsigned bi = j * nr_blocks(0) + i;
				char* block_data_ptr = &block_data[bi * block_size + block_z * block_slice_size];

				// copy data from input image into current slice of block
				for (int y = y0; y < y1; ++y) {
					std::copy(image_data_ptr, image_data_ptr + copy_line_size, block_data_ptr);
					image_data_ptr += image_line_size;
					block_data_ptr += block_line_size;
				}
			}
		}
		// if blocks are full or no more input slices left, write blocks to file
		block_z += 1;
		if (block_z == bgi.block_dimensions(2) + bgi.overlap(2) || z + 1 >= file_names.size()) {
			// generate out file name
			std::stringstream ss;
			ss << output_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << k << ".bvx";
			if (!cgv::utils::file::exists(ss.str())) {
				if (!cgv::utils::file::write(ss.str(), &block_data.front(), block_data.size())) {
					std::cerr << "could not write block slice " << k << " of level 0" << std::endl;
					return false;
				}
				std::cout << "wrote block slice " << k << " of level 0" << std::endl;
			}
			else
				std::cout << "skipped block slice " << k << " of level 0" << std::endl;
		}
		// reuse overlap for next block slice
		if (block_z == bgi.block_dimensions(2) + bgi.overlap(2)) {
			// iterate all blocks
			for (int j = 0; j < nr_blocks(1); ++j) {
				for (int i = 0; i < nr_blocks(0); ++i) {
					// copy overlap to block begin and initialize remaining of block to 0
					int bi = j * nr_blocks(0) + i;
					char* block_begin_ptr = &block_data[bi * block_size];
					char* block_overlap_ptr = &block_data[bi * block_size + bgi.block_dimensions(2) * block_slice_size];
					std::copy(block_overlap_ptr, block_overlap_ptr + bgi.overlap(2) * block_slice_size, block_begin_ptr);
					std::fill(block_begin_ptr + bgi.overlap(2) * block_slice_size, block_begin_ptr + block_size, 0);
				}
			}
			block_z = bgi.overlap(2);
			k += 1;
		}
	}
	return true;
}

/*
bool create_levels(const std::string& path, unsigned nr_levels, const block_generation_info& bgi)
{
cgv::type::uint64_type block_size = cgv::type::uint64_type(bgi.block_dimensions[0] + 1) * (bgi.block_dimensions[1] + 1) * (bgi.block_dimensions[2] + 1) * 3;

// build level information
std::vector<volume::dimension_type> nr_blocks;
std::vector<volume::dimension_type> level_dimensions;
volume::dimension_type dimensions = bgi.dimensions;
level_dimensions.push_back(bgi.dimensions);
for (unsigned level = 0; level <= nr_levels; ++level) {
level_dimensions.push_back(dimensions);
nr_blocks.push_back(volume::dimension_type(
unsigned(ceil(float(level_dimensions[level](0)) / bgi.block_dimensions(0))),
unsigned(ceil(float(level_dimensions[level](1)) / bgi.block_dimensions(1))),
unsigned(ceil(float(level_dimensions[level](2)) / bgi.block_dimensions(2)))));

for (int c = 0; c < 3; ++c) {
if (bgi.subsampling_offset[c] <= int(level))
dimensions(c) = unsigned(ceil(float(dimensions(c)) / bgi.subsampling_factor[c]));
}
};

// iterate
std::stringstream ss;
ss << output_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << k << ".bvx";
if (!cgv::utils::file::exists(ss.str())) {
if (!cgv::utils::file::write(ss.str(), &block_data.front(), block_data.size())) {
std::cerr << "could not write block slice " << k << " of level 0" << std::endl;
return false;
}
std::cout << "wrote block slice " << k << " of level 0" << std::endl;
}
else
std::cout << "skipped block slice " << k << " of level 0" << std::endl;

}
*/

class slice_processor;

typedef cgv::data::ref_ptr<slice_processor, true> slice_processor_ptr;



class slice_processor : public cgv::data::ref_counted
{
public:
	typedef cgv::data::data_view slice_type;
	typedef cgv::data::ref_ptr<slice_type> slice_ptr;
	typedef cgv::data::ref_ptr<const slice_type> const_slice_ptr;
protected:
	slice_processor_ptr child_processor;
	void handle_children(int z, const_slice_ptr slice)
	{
		if (child_processor)
			child_processor->handle_slice(z, slice);
	}
public:
	slice_processor() {}
	virtual ~slice_processor() {}
	void set_child_processor(slice_processor_ptr _child_ptr) { child_processor = _child_ptr; }
	virtual void handle_slice(int z, const_slice_ptr slice) {
		handle_children(z, slice);
	}
};

class block_level_writer : public slice_processor
{
protected:
	unsigned level;
	volume::dimension_type level_dimensions;
	volume::dimension_type block_dimensions;
	volume::dimension_type block_overlap;
	volume::dimension_type nr_blocks;
	size_t pixel_size;
	size_t block_line_size;
	size_t block_slice_size;
	size_t block_size;
	size_t image_line_size;
	size_t image_size;
	std::vector<char> block_data;
	size_t next_block_completion;
	unsigned block_z, k;
	std::string output_path;
public:
	block_level_writer(const std::string& _output_path, unsigned _level, volume::dimension_type _level_dimensions, volume::dimension_type _block_dimensions, volume::dimension_type _block_overlap, size_t _pixel_size)
	{
		level = _level;
		output_path = _output_path;
		level_dimensions = _level_dimensions;
		block_dimensions = _block_dimensions;
		block_overlap = _block_overlap;
		nr_blocks = volume::dimension_type(
			unsigned(ceil(float(level_dimensions(0)) / block_dimensions(0))),
			unsigned(ceil(float(level_dimensions(1)) / block_dimensions(1))),
			unsigned(ceil(float(level_dimensions(2)) / block_dimensions(2)))
		);
		// allocate data for one slice of blocks and initialize all values to 0
		pixel_size = _pixel_size;
		block_line_size = (block_dimensions(0) + block_overlap(0)) * pixel_size;
		block_slice_size = (block_dimensions(1) + block_overlap(1)) * block_line_size;
		block_size = (block_dimensions(2) + block_overlap(2)) * block_slice_size;
		image_line_size = level_dimensions(0)*pixel_size;
		image_size = level_dimensions(1)*image_line_size;
		block_data.resize(nr_blocks(0)*nr_blocks(1)*block_size, 0);
		next_block_completion = block_dimensions(2) + block_overlap(2) - 1;
		block_z = 0;
		k = 0;
	}
	void handle_slice(int z, const_slice_ptr slice)
	{
		// insert data into current blocks
		for (int j = 0; j < nr_blocks(1); ++j) {
			for (int i = 0; i < nr_blocks(0); ++i) {
				// compute rectangle in input image corresponding to block 
				int x0 = i*block_dimensions(0);
				int x1 = x0 + block_dimensions(0) + block_overlap(0);
				if (x1 > level_dimensions(0))
					x1 = level_dimensions(0);
				size_t copy_line_size = (x1 - x0) * pixel_size;

				int y0 = j*block_dimensions(1);
				int y1 = y0 + block_dimensions(1) + block_overlap(1);
				if (y1 > level_dimensions(1))
					y1 = level_dimensions(1);
				const char* image_data_ptr = slice->get_ptr<const char>() + (y0 * image_line_size + x0 * pixel_size);

				// compute block index and pointer to block data
				unsigned bi = j * nr_blocks(0) + i;
				char* block_data_ptr = &block_data[bi * block_size + block_z * block_slice_size];

				// copy data from input image into current slice of block
				for (int y = y0; y < y1; ++y) {
					std::copy(image_data_ptr, image_data_ptr + copy_line_size, block_data_ptr);
					image_data_ptr += image_line_size;
					block_data_ptr += block_line_size;
				}
			}
		}
		// if blocks are full or no more input slices left, write blocks to file
		block_z += 1;
		if (block_z == block_dimensions(2) + block_overlap(2) || z + 1 >= level_dimensions(2)) {
			// generate out file name
			std::stringstream ss;
			ss << output_path << "/level_" << std::setw(2) << std::setfill('0') << level << "_blockslice_" << std::setw(3) << std::setfill('0') << k << ".bvx";
			if (!cgv::utils::file::exists(ss.str())) {
				if (!cgv::utils::file::write(ss.str(), &block_data.front(), block_data.size())) {
					std::cerr << "could not write block slice " << k << " of level 0" << std::endl;
					handle_children(z, slice);
					return;
				}
				std::cout << "wrote block slice " << k << " of level " << level << std::endl;
			}
			else
				std::cout << "skipped block slice " << k << " of level " << level << std::endl;
		}
		// reuse overlap for next block slice
		if (block_z == block_dimensions(2) + block_overlap(2)) {
			// iterate all blocks
			for (int j = 0; j < nr_blocks(1); ++j) {
				for (int i = 0; i < nr_blocks(0); ++i) {
					// copy overlap to block begin and initialize remaining of block to 0
					int bi = j * nr_blocks(0) + i;
					char* block_begin_ptr = &block_data[bi * block_size];
					char* block_overlap_ptr = &block_data[bi * block_size + block_dimensions(2) * block_slice_size];
					std::copy(block_overlap_ptr, block_overlap_ptr + block_overlap(2) * block_slice_size, block_begin_ptr);
					std::fill(block_begin_ptr + block_overlap(2) * block_slice_size, block_begin_ptr + block_size, 0);
				}
			}
			block_z = block_overlap(2);
			k += 1;
		}
		handle_children(z, slice);
	}
};

template <typename T> struct calc_type { typedef T type; };
template <> struct calc_type<cgv::type::uint8_type> { typedef cgv::type::int16_type type; };

template <typename T, unsigned N>
class slice_subsampler : public slice_processor
{
public:
	typedef cgv::math::fvec<T, N> pixel_type;
protected:
	int k;
	bool is_even;
	bool reduce_slices;
	const_slice_ptr prev_slice;
	cgv::data::data_format df;
	slice_ptr subsampled_slice;
public:
	slice_subsampler(bool _reduce_slices) : reduce_slices(_reduce_slices)
	{
		is_even = true;
		k = 0;
	}
	void handle_slice(int z, const_slice_ptr slice)
	{
		if (!subsampled_slice) {
			df = *slice->get_format();
			df.set_width((df.get_width() + 1) / 2);
			df.set_height((df.get_height() + 1) / 2);
			subsampled_slice = new cgv::data::data_view(&df);
		}
		if (reduce_slices) {
			if (is_even) {
				prev_slice = slice;
				is_even = false;
				return;
			}
			else {
				cgv::media::image::subsample_slice<calc_type<T>::type>(prev_slice->get_ptr<pixel_type>(), slice->get_ptr<pixel_type>(),
					subsampled_slice->get_ptr<pixel_type>(), slice->get_format()->get_width(), slice->get_format()->get_height(), N);
				is_even = true;
			}
		}
		else {
			cgv::media::image::subsample_image<calc_type<T>::type>(slice->get_ptr<pixel_type>(),
				subsampled_slice->get_ptr<pixel_type>(), slice->get_format()->get_width(), slice->get_format()->get_height(), N);
		}
		handle_children(k, subsampled_slice);
		++k;
	}
};

class slice_writer : public slice_processor
{
protected:
	std::string output_path;
	std::string prefix;
	std::string postfix;
	std::string options;
public:
	slice_writer(const std::string& _output_path, const std::string& _prefix, const std::string& _postfix, const std::string& _options = "") :
		output_path(_output_path), prefix(_prefix), postfix(_postfix), options(_options)
	{
		// ensure that output directory exists
		if (!ensure_directory(output_path)) {
			std::cerr << "could not create output path <" << output_path << ">" << std::endl;
			abort();
		}
	}
	void handle_slice(int z, const_slice_ptr slice)
	{
		// generate out file name
		std::stringstream ss;
		ss << output_path << "/" << prefix << std::setw(4) << std::setfill('0') << z << postfix;

		if (!cgv::utils::file::exists(ss.str())) {
			cgv::media::image::image_writer iw(ss.str());
			if (!options.empty())
				iw.multi_set(options);
			if (!iw.write_image(*slice))
				std::cerr << "could not write slice " << z << " to path <" << ss.str() << ">" << std::endl;
			else {
				std::cerr << "wrote slice " << z << " to path <" << ss.str() << ">" << std::endl;
			}
			iw.close();
		}
		handle_children(z, slice);
	}
};

slice_processor_ptr build_processor_tree(const block_generation_info& bgi, const std::string& block_output_path, const std::string& image_ouput_path = "")
{
	cgv::type::uint64_type original_size = cgv::type::uint64_type(bgi.dimensions[0]) * bgi.dimensions[1] * bgi.dimensions[2] * 3;
	std::cout << "original data: <" << bgi.dimensions << "> = " << original_size << " = " << double(original_size) / (1024 * 1024 * 1024) << std::endl;
	unsigned level = 0;
	volume::dimension_type level_dimensions = bgi.dimensions;
	volume::dimension_type nr_blocks;

	slice_processor_ptr main_sp_ptr;
	slice_processor_ptr last_sp_ptr;
	do {
		nr_blocks = volume::dimension_type(
			unsigned(ceil(float(level_dimensions[0]) / bgi.block_dimensions[0])),
			unsigned(ceil(float(level_dimensions[1]) / bgi.block_dimensions[1])),
			unsigned(ceil(float(level_dimensions[2]) / bgi.block_dimensions[2])));

		slice_processor_ptr blw_ptr = new block_level_writer(block_output_path, level, level_dimensions, bgi.block_dimensions, bgi.overlap, bgi.get_block_pixel_size());
		if (last_sp_ptr)
			last_sp_ptr->set_child_processor(blw_ptr);
		else
			main_sp_ptr = blw_ptr;
		slice_processor_ptr ssp_ptr = new slice_subsampler<cgv::type::uint8_type, 3>(int(level) + 1 >= bgi.subsampling_offset[2]);

		if (!image_ouput_path.empty() && level > 0) {
			std::stringstream ss;
			ss << "level_" << std::setw(2) << std::setfill('0') << level << "_";
			slice_processor_ptr swp_ptr = new slice_writer(image_ouput_path, ss.str(), ".png");
			blw_ptr->set_child_processor(swp_ptr);
			swp_ptr->set_child_processor(ssp_ptr);
		}
		else {
			blw_ptr->set_child_processor(ssp_ptr);
		}
		last_sp_ptr = ssp_ptr;

		cgv::type::uint64_type block_size = cgv::type::uint64_type(bgi.block_dimensions[0] + 1) * (bgi.block_dimensions[1] + 1) * (bgi.block_dimensions[2] + 1) * 3;
		cgv::type::uint64_type nr_all_blocks = cgv::type::uint64_type(nr_blocks[0]) * nr_blocks[1] * nr_blocks[2];
		std::cout << level << ": <" << level_dimensions << ">, nr blocks: <" << nr_blocks << "> = " << nr_all_blocks * block_size << " = " << double(nr_all_blocks * block_size) / (1024 * 1024 * 1024) << std::endl;

		for (int c = 0; c < 3; ++c) {
			if (bgi.subsampling_offset[c] <= int(level))
				level_dimensions[c] = unsigned(ceil(float(level_dimensions[c]) / bgi.subsampling_factor[c]));
		}
		++level;
	} while (std::min(std::min(nr_blocks[0], nr_blocks[1]), nr_blocks[2]) >= 2);
	return main_sp_ptr;
}

bool process_directory(const std::string& input_path, cgv::type::info::TypeId component_type, cgv::data::ComponentFormat component_format, int w, int h, slice_processor_ptr processor)
{
	// collect input files
	std::vector<std::string> file_names;
	void* handle = cgv::utils::file::find_first(input_path + "/*.*");
	while (handle) {
		if (!cgv::utils::file::find_directory(handle))
			file_names.push_back(cgv::utils::file::find_name(handle));
		handle = cgv::utils::file::find_next(handle);
	}
	if (file_names.empty()) {
		std::cerr << "did not find files in directory <" << input_path << ">" << std::endl;
		return false;
	}

	cgv::data::data_format df(w, h, component_type, component_format);
	size_t image_size = df.get_nr_entries()*df.get_entry_size();
	for (unsigned z = 0; z < file_names.size(); ++z) {
		// read and if necessary extent next slice
		std::cout << "read slice " << z << ": <" << file_names[z] << ">" << std::endl;
		std::vector<char> input_slice;
		if (!read_image(input_path + "/" + file_names[z], input_slice)) {
			std::cerr << "failed to read <" << file_names[z] << ">" << std::endl;
			continue;
		}
		if (input_slice.size() < image_size) {
			size_t old_size = input_slice.size();
			input_slice.resize(image_size);
			std::fill(input_slice.begin() + old_size, input_slice.end(), 0);
		}
		// construct data view from slice
		slice_processor::const_slice_ptr slice = slice_processor::const_slice_ptr(new slice_processor::slice_type(&df));
		std::copy(input_slice.begin(), input_slice.end(), slice->get_ptr<char>());

		// process slice
		processor->handle_slice(z, slice);
	}
	return true;
}

bool subsample_directory(const std::string& input_path, const std::string& output_path, const block_generation_info& bgi, bool subsample_slices)
{
	// ensure that output directory exists
	if (!ensure_directory(output_path)) {
		std::cout << "could not create output path <" << output_path << ">" << std::endl;
		return false;
	}

	// collect input files
	std::vector<std::string> file_names;
	void* handle = cgv::utils::file::find_first(input_path + "/*.*");
	while (handle) {
		if (!cgv::utils::file::find_directory(handle))
			file_names.push_back(cgv::utils::file::find_name(handle));
		handle = cgv::utils::file::find_next(handle);
	}
	if (file_names.empty()) {
		std::cerr << "did not find files in directory <" << input_path << ">" << std::endl;
		return false;
	}

	// compute number of blocks on base levels
	volume::dimension_type level_dimensions = bgi.dimensions;
	size_t pixel_size = bgi.get_block_pixel_size();
	size_t image_line_size = level_dimensions(0)*pixel_size;
	size_t image_size = level_dimensions(1)*image_line_size;

	volume::dimension_type next_level_dimensions = volume::dimension_type((level_dimensions(0) + 1) / 2, (level_dimensions(1) + 1) / 2, subsample_slices ? (level_dimensions(2) + 1) / 2 : level_dimensions(2));
	size_t subsampled_image_line_size = next_level_dimensions(0)*pixel_size;
	size_t subsampled_image_size = next_level_dimensions(1)*subsampled_image_line_size;

	// iterate over slices once
	bool is_even = true;
	std::vector<char> prev_slice;
	std::vector<char> subsampled_slice;
	subsampled_slice.resize(subsampled_image_size);
	cgv::math::fvec<cgv::type::uint8_type, 3>* structured_subsampled_slice_ptr = reinterpret_cast<cgv::math::fvec<cgv::type::uint8_type, 3>*>(&subsampled_slice.front());
	unsigned k = 0;
	for (unsigned z = 0; z < file_names.size(); ++z) {

		// read and if necessary extent next slice
		std::cout << "read slice " << z << ": <" << file_names[z] << ">" << std::endl;
		std::vector<char> slice;
		if (!read_image(input_path + "/" + file_names[z], slice)) {
			std::cerr << "failed to read <" << file_names[z] << ">" << std::endl;
			continue;
		}
		if (slice.size() < image_size) {
			size_t old_size = slice.size();
			slice.resize(image_size);
			std::fill(slice.begin() + old_size, slice.end(), 0);
		}

		// gather two slices if also slices should be subsampled
		cgv::math::fvec<cgv::type::uint8_type, 3>* structured_slice_ptr = reinterpret_cast<cgv::math::fvec<cgv::type::uint8_type, 3>*>(&slice.front());
		if (subsample_slices) {
			if (is_even) {
				prev_slice = slice;
				is_even = false;
				continue;
			}
			else {
				cgv::math::fvec<cgv::type::uint8_type, 3>* structured_prev_slice_ptr = reinterpret_cast<cgv::math::fvec<cgv::type::uint8_type, 3>*>(&prev_slice.front());
				cgv::media::image::subsample_slice<cgv::type::int16_type>(structured_prev_slice_ptr, structured_slice_ptr, structured_subsampled_slice_ptr, level_dimensions(0), level_dimensions(1), 3);
				is_even = true;
			}
		}
		else
			cgv::media::image::subsample_image<cgv::type::int16_type>(structured_slice_ptr, structured_subsampled_slice_ptr, level_dimensions(0), level_dimensions(1), 3);

		//		cgv::math::fvec<cgv::type::uint8_type, 3>* structured_data_ptr = reinterpret_cast<cgv::math::fvec<cgv::type::uint8_type, 3>*>(&data.front());
		//		cgv::math::fvec<cgv::type::int8_type, 3>* signed_data_ptr = reinterpret_cast<cgv::math::fvec<cgv::type::int8_type, 3>*>(&data.front());
		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(structured_data_ptr, level_dimensions(0), level_dimensions(1), 1, level_dimensions(0), 3, 255, true, 2, 1);
		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(1), level_dimensions(0), level_dimensions(0), 1, 3, 255, true, 1, 1);

		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(0) / 2, level_dimensions(1) / 2, 1, level_dimensions(0), 3, 255, true, 1, 2);
		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(1) / 2, level_dimensions(0) / 2, level_dimensions(0), 1, 3, 255, true, 1, 1);
		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(0) / 4, level_dimensions(1) / 4, 1, level_dimensions(0), 3, 255, true, 1, 2);
		//		integer_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(1) / 4, level_dimensions(0) / 4, level_dimensions(0), 1, 3, 255, true, 1, 1);

		//integer_inverse_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(signed_data_ptr, level_dimensions(1), level_dimensions(0), level_dimensions(0), 1, 3, 255, true, 1, 1);
		//integer_inverse_wavelet_transform<cgv::math::fvec<cgv::type::int16_type, 3>, cgv::math::fvec<cgv::type::int8_type, 3> >(structured_data_ptr, level_dimensions(0), level_dimensions(1), 1, level_dimensions(0), 3, 255, true, 2, 1);

		// generate out file name
		std::stringstream ss;
		ss << output_path << "/level_" << std::setw(2) << std::setfill('0') << 1 << "_" << std::setw(2) << std::setfill('0') << k << ".tif";
		if (!write_image(ss.str(), cgv::data::component_format(bgi.type_id, bgi.components), next_level_dimensions(0), next_level_dimensions(1), &subsampled_slice.front())) {
			std::cerr << "could not write transformed image " << ss.str() << std::endl;
		}
		else
			std::cout << "wrote transformed image " << ss.str() << std::endl;
		++k;
	}
	return true;
}





//Testing:
void init_to_visible_human_male_png(block_generation_info& bgi)
{
	int width = 2048;
	int height = 1216;
	int amount = 377;

	bgi.type_id = cgv::type::info::TI_UINT8;
	bgi.components = cgv::data::CF_RGB;
	bgi.dimensions.set(width, height, amount);
	bgi.extent.set(width * 0.114f, height * 0.114f, amount * 1.0f);

	bgi.block_dimensions.set(16, 16, 13);
	bgi.overlap.set(1, 1, 1);
	bgi.subsampling_factor.set(2, 2, 2);
	bgi.subsampling_offset.set(0, 0, 3);
}

#include <cgv/base/register.h>
int main(int argc, char** argv)
{
	// needed for plugin stuff
	cgv::base::enable_registration();
	// read image io support for png and tiff
	cgv::base::load_plugin("cmi_io");

	/// init block generation info structures
	block_generation_info visible_human_male;
	init_to_visible_human_male(visible_human_male);
	block_generation_info visible_human_female;
	init_to_visible_human_female(visible_human_female);

	/// test
	block_generation_info visible_human_male_png;
	init_to_visible_human_male_png(visible_human_male_png);

	block_generation_info visible_human = visible_human_male_png;

	// single resolution approach
	build_level_infos(visible_human);
	if (build_blocks_from_directory("E:/data/visual_human/male/PNG_format/head", "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/male/block16/head", visible_human))
		std::cout << "success" << std::endl;
	else
		std::cout << "failed" << std::endl;
	std::cin.get();
	return 1;

	// advanced multi resolution approach not validated tested!!!
	slice_processor_ptr sp_ptr = build_processor_tree(visible_human, "E:/visible human/male/block16", "E:/visible human/male/mipmap");
	process_directory("E:/visible human/male/70mm/fullbody", visible_human.type_id, visible_human.components, visible_human.dimensions(0), visible_human.dimensions(1), sp_ptr);

	if (subsample_directory("E:/visible human/male/70mm/fullbody", "E:/visible human/male/transformed", visible_human, true))
		std::cout << "success" << std::endl;
	else
		std::cout << "failed" << std::endl;
}