#pragma once

#include <string>
#include <vector>
#include <cgv/reflect/reflect_enum.h>

#include <cgv/base/node.h>
#include <cgv/render/drawable.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/event_handler.h>

#include <cgv/render/texture.h>
#include <cgv/render/shader_program.h>
#include <cgv_gl/gl/gl_view.h>

#include "lib_begin.h"

class CGV_API volume_slicer :
	public cgv::base::node,
	public cgv::render::drawable,
	public cgv::gui::provider,
	public cgv::gui::event_handler
{
public:
	/// type for 3d vectors
	typedef cgv::math::fvec<float, 3> vec3;
	/// type for axis aligned 3d boxes
	typedef cgv::media::axis_aligned_box<float, 3> box3;
	/// type for index tripels
	typedef cgv::math::fvec<int, 3> ivec3;
	/// standard constructor
	volume_slicer();
	/// returns "volume_slicer"
	std::string get_type_name() const;
	/// callback used to notify instance of member changes
	void on_set(void* member_ptr);
	/// self reflection is used to publish standard member variables to the set and get mechanism
	bool self_reflect(cgv::reflect::reflection_handler& rh);

	/**@name volume*/
	//@{
protected:
	/// filename of the volume dataset
	std::string file_name;
	/// for each major direction the number of voxels in the volume
	ivec3 dimensions;
	/// for each major direction the extent of the volume  
	vec3 extent;
	/// whether to use tri-linear texture interpolation
	bool interpolate;
public:
	// read volume from sliced or single file representations
	bool open_block_volume(const std::string& directory_name);
	/// read regular volume file
	bool open_volume(const std::string& file_name);
	//@}

	/**@name slicing*/
	//@{
protected:
	/// slice normal in texture coordinates
	vec3 slice_normal_tex;
	/// distance to center of volume in texture coordinates
	float slice_distance_tex;
	//@}

	/**@name rendering*/
	//@{
protected:
	/// flag that tells whether to show the bounding box of the volume
	bool show_box;
	/// material used to render the bounding box faces
	cgv::media::illum::phong_material box_mat;
private:
	/// volume slicing shader program not using blocks
	cgv::render::shader_program tex_slice_prog;
	/// used later to access block based volume data through indices
	cgv::render::shader_program indexed_tex_slice_prog;
	/// texture storing the volume data
	cgv::render::texture data_texture;
	/// index texture will be used later to store for each block the index where it is stored in the data texture
	cgv::render::texture index_texture;
public:
	/// called once per context for initialization of rendering stuff
	bool init(cgv::render::context& ctx);
	/// called once per frame before drawing
	void init_frame(cgv::render::context& ctx);
	/// called to draw the frame
	void draw(cgv::render::context& ctx);
	/// called to destruct the rendering objects
	void clear(cgv::render::context& ctx);
	//@}

	/**@name interaction*/
	//@{
private:
	/// store pointer to view which is controlled by the user
	cgv::render::gl::gl_view* view_ptr;
	/// traverse object hierarchy to extract view_ptr
	bool ensure_view_ptr();
	/// sensitivity used to map mouse movement to rotation of slice normal
	double rotate_sensitivity;
	/// sensitivity used to map mouse movement to slice translations
	double translation_sensitivity;
public:
	/// adjusts view to bounding box of all instances
	void auto_adjust_view();
	/// overload and implement this method to handle events
	bool handle(cgv::gui::event& e);
	/// overload to stream help information to the given output stream
	void stream_help(std::ostream& os);
	/// stream statistical information about volume
	void stream_stats(std::ostream& os);
	/// called to create gui elements
	void create_gui();
	//@}

};

#include <cgv/config/lib_end.h>
