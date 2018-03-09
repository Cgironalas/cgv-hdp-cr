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
#include <cgv/math/mfunc.h>
#include <cgv/media/mesh/streaming_mesh.h>

#include <vol_data/volume.h>
#include "cache_manager.h"

#include "lib_begin.h"


enum DiffuseColorMapping
{
	DCM_MATERIAL,
	DCM_COLOR,
	DCM_TEXTURE,
	DCM_NORMAL,
	DCM_TEXTURE_COORDINATES,
	DCM_POSITION
};

struct surface_render_style
{
	/// diffuse color mapping mode
	DiffuseColorMapping diffuse_mapping;
	/// whether to use illumination
	bool illumination;
	/// whether to show surface in wire frame mode
	bool show_wireframe;
	/// line width used for wireframe rendering
	float wireframe_width;
	/// culling mode
	cgv::render::CullingMode culling_mode;
	/// whether to compute per face normals in the shader pipeline
	bool face_normals;
	/// whether to negate surface normals
	bool negate_normals;
	/// whether to show surface with two sided illumination
	bool two_sided;
	/// material used to render the surface faces
	cgv::media::illum::phong_material material;
	/// set default values (no wireframe, no culling, no face normals, no negation, no two sided)
	surface_render_style();
};

class CGV_API volume_slicer :
	public cgv::base::node, 
	public cgv::render::drawable, 
	public cgv::gui::provider, 
	public cgv::gui::event_handler,
	public cgv::media::mesh::streaming_mesh_callback_handler,
	public cgv::math::v3_func<float, float>
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
public:
	/// filename of the volume dataset
	std::string file_name;
	/// filename to open volume dataset for iso surface extraction
	std::string iso_file_name;
	/// folder path with the .bvx slices
	std::string slices_path;
	// volume data structure stored for iso surface extaction
	volume V;
	/// return the value of a given voxel scaled to [0,1]
	float get_voxel_value(const ivec3& voxel) const;
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
	bool open_volume(const std::string& file_name, bool is_iso);
	/// convert world to texture coordinates
	vec3 texture_from_world_coordinates(const vec3& p_world) const;
	/// convert texture to voxel coordinates
	vec3 voxel_from_texture_coordinates(const vec3& p_texture) const;
	vec3 voxel_from_texture_coordinates_blocks(const vec3& p_texture) const;
	/// convert texture to world coordinates
	vec3 world_from_texture_coordinates(const vec3& p_texture) const;
	/// convert voxel to texture coordinates
	vec3 texture_from_voxel_coordinates(const vec3& p_voxel) const;
	vec3 texture_from_voxel_coordinates_blocks(const vec3& p_voxel) const;

	vec3 world_from_texture_normals(const vec3& n_texture) const;

	//@}

	/**@name blocks*/
	//@{
public:
	/// size of block in voxel, for example (16,16,16) or (32,32,16)
	ivec3 block_dimensions;
	/// size of block in voxel, for example (16,16,16) or (32,32,16)
	ivec3 slices_dimensions;
	/// overlap between adjacent blocks in voxel, typically (1,1,1)
	ivec3 overlap;
	/// store the intersected blocks
	std::vector<ivec3> intersected_blocks;
	/// identifies wether the batch has been updated
	bool new_batch = false;
	/// store the previous iteration of intersected blocks
	std::unordered_set<ivec3, container_hash<ivec3>> previous_intersected_blocks;
	/// logic for putting blocks in cache
	cache_manager threaded_cache_manager;
public:
	/// convert voxel to block coordinates
	vec3 block_from_voxel_coordinates(const vec3& p_voxel) const;
	/// convert voxel to block coordinates
	vec3 voxel_from_block_coordinates(const vec3& p_block) const;
	/// computes intersected blocks in the vector
	void update_intersected_blocks(cgv::render::context& ctx);
	/// compute distance of point in texture coordinates to current slice
	float compute_distance_to_slice_tex(const vec3& p) const;
	/// determine if a block is intersected by the plane
	bool is_block_intersected(const box3& B);
	/// return a block with the var_dim in the max_ind position, and the other two according to the place of the var_dim
	vec3 make_vec_d_max(int max_ind, float var_dim, float dim0, float dim1);

	//@}

	/**@name iso surface*/
	//@{
protected:
	/// store the iso surface in form of a triangle mesh with positions in world coordinates and texture coordinates of the 3D volume 
	std::vector<vec3> positions;
	std::vector<vec3> normals;
	std::vector<vec3> texture_coords;
	std::vector<cgv::type::uint32_type>  corner_indices;
	/// number of vertices per face
	unsigned face_vertex_count;
	/// downsampling factor allows to extract iso-surface on a coarser resolution
	unsigned down_sampling_factor;
	/// the iso_value to extract the iso surface - always in [0,1]
	float iso_value;
	std::map<int, std::string> *value_name_map;
	std::string last_name;
	std::map<std::string, std::vector<int> > *value_map;
	std::string value_name_map_file_name;
	bool read_value_name_map(const std::string& file_name);
private:
	/// whether to access iso volume with trilinear interpolation
	bool interpolated_evaluation;
	/// pointer to marching cubes object to get access to positions of vertices
	cgv::media::mesh::streaming_mesh<float>* sm_ptr;
protected:
	/// interface for evaluation of the multivariate function
	float evaluate(const pnt_type& p) const;
	/// called when a new vertex is generated
	void new_vertex(unsigned int vertex_index);
	/// announces a new polygon defines by the vertex indices stored in the given vector
	void new_polygon(const std::vector<unsigned int>& vertex_indices);
	/// drop the currently first vertex that has the given global vertex index
	void before_drop_vertex(unsigned int vertex_index);

public:
	/// extract iso surface from volume with current iso_value and down_sampling_factor
	void extract_surface(int algorithm);
	/// remove all elements from surface
	void clear_surface();
	/// draw the surface with opengl
	void draw_surface(cgv::render::context& ctx);
	//@}
	/**@name slicing*/
	//@{
protected:
	/// slice normal in texture coordinates
	vec3 slice_normal_tex;
	/// distance to center of volume in texture coordinates
	float slice_center_distance_tex;
	//@}

	/**@name rendering*/
	//@{
protected:
	/// flag that tells whether to show the bounding box of the volume
	bool show_box;
	/// flag that tells whether to show the current voxel
	bool show_voxel;
	/// flag that tells whether to show the current block
	bool show_block;
	/// flag that tells whether to show the iso surface
	bool show_surface;
	/// 
	float texture_scale;
	/// 
	float texture_gamma;
	///
	surface_render_style iso_surface_style;
	/// material used to render the bounding box faces
	cgv::media::illum::phong_material box_mat;
private:
	/// volume slicing shader program not using blocks
	cgv::render::shader_program tex_slice_prog;
	/// used later to access block based volume data through indices
	cgv::render::shader_program indexed_tex_slice_prog;
	/// texture storing the volume data
	cgv::render::texture data_texture;
	/// whether to use reduced color depth to store texture on GPU
	bool use_reduced_color_depth;
	/// index texture will be used later to store for each block the index where it is stored in the data texture
	cgv::render::texture index_texture;
	///
	bool use_standard_prog;
	/// surface shading program 
	cgv::render::shader_program surface_prog;
	/// enables surface rendering with the given rendering style
	void enable_surface_rendering(cgv::render::context& ctx, const surface_render_style& style);
	/// disables surface rendering with the given rendering style
	void disable_surface_rendering(cgv::render::context& ctx, const surface_render_style& style);
public:
	/// called once per context for initialization of rendering stuff
	bool init(cgv::render::context& ctx);
	/// called once per frame before drawing
	void init_frame(cgv::render::context& ctx);
	/// draw a box specified in world coordinates in wire frame mode 
	void draw_box(cgv::render::context& ctx, const box3& B, const cgv::media::illum::phong_material& material);
	/// draw a box specified in world coordinates in wire frame mode 
	void draw_wire_box(cgv::render::context& ctx, const box3& B, const vec3& color);
	/// draw a voxel with a wire frame box
	void draw_voxel(cgv::render::context& ctx, const ivec3& voxel, const vec3& color);
	/// draw a block with a wire frame box
	void draw_block(cgv::render::context& ctx, const ivec3& block, const vec3& color, const vec3& overlap_color);
	/// draws all blocks that intersect the plane 
	void draw_blocks_in_plane(cgv::render::context& ctx, const vec3& color, const vec3& overlap_color);
	/// called to draw the frame
	void draw(cgv::render::context& ctx);
	/// draw textual information here
	void after_finish(cgv::render::context& ctx);
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
	/// voxel information to be drawn with font
	std::string voxel_value_info;
	/// last x and y position of mouse
	int last_x, last_y;
	/// store current voxel
	ivec3 current_voxel;
	/// store current block
	ivec3 current_block;
	/// return the point under the mouse pointer in world coordinates
	bool get_picked_point(int x, int y, vec3& p_pick_world);
	/// determine voxel location of mouse pointer
	void peek_voxel_values(int x, int y);
public:
	/// adjusts view to bounding box of all instances
	void auto_adjust_view();
	/// overload and implement this method to handle events
	bool handle(cgv::gui::event& e);
	/// overload to stream help information to the given output stream
	void stream_help(std::ostream& os);
	/// stream statistical information about volume
	void stream_stats(std::ostream& os);
	///
	void create_surface_gui(surface_render_style& style);
	/// called to create gui elements
	void create_gui();
	//@}

};

#include <cgv/config/lib_end.h>
