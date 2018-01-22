#include "volume_slicer.h"
#include <vol_data/volume_io.h>
#include <vol_data/sliced_volume_io.h>
#include <cgv/utils/scan.h>
#include <cgv/base/find_action.h>
#include <cgv/utils/file.h>
#include <cgv_gl/gl/gl.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/mouse_event.h>
#include <cgv/media/image/image_writer.h>


/// standard constructor
volume_slicer::volume_slicer() : cgv::base::node("volume_slicer")
{
	view_ptr = 0;
	rotate_sensitivity = 1.0;
	translation_sensitivity = 1.0;
	show_box = true;
	show_voxel = true;
	show_block = true;
	show_surface = true;
	show_surface_wireframe = false;
	slice_center_distance_tex = 0;
	slice_normal_tex = vec3(0, 0, 1);
	dimensions = ivec3(1, 1, 1);
	extent = vec3(1, 1, 1);
	interpolate = false;
	box_mat.set_ambient(cgv::media::illum::phong_material::color_type(0.8f, 0.8f, 0.8f, 1.0f));
	box_mat.set_diffuse(cgv::media::illum::phong_material::color_type(0.6f, 0.6f, 0.6f, 1.0f));
	box_mat.set_specular(cgv::media::illum::phong_material::color_type(0.0f, 0.0f, 0.0f, 1.0f));
	box_mat.set_shininess(20.0f);

	surface_mat.set_ambient(cgv::media::illum::phong_material::color_type(0.8f, 0.8f, 0.8f, 1.0f));
	surface_mat.set_diffuse(cgv::media::illum::phong_material::color_type(0.6f, 0.6f, 0.6f, 1.0f));
	surface_mat.set_specular(cgv::media::illum::phong_material::color_type(1.0f, 1.0f, 1.0f, 1.0f));
	surface_mat.set_shininess(60.0f);
#ifdef _DEBUG
	down_sampling_factor = 32;
#else
	down_sampling_factor = 4;
#endif
	negate_normals = true;
	iso_value = 0.5f;
	sm_ptr = 0;

	last_x = last_y = -1;
	current_voxel = current_block = ivec3(0, 0, 0);
	block_dimensions = ivec3(2, 2, 16);
	overlap = ivec3(1, 1, 1);

}

bool volume_slicer::ensure_view_ptr()
{
	if (view_ptr)
		return true;
	std::vector<cgv::render::gl::gl_view*> view_ptrs;
	cgv::base::find_interface(get_parent(), view_ptrs);
	if (view_ptrs.empty())
		return false;
	view_ptr = view_ptrs.front();
	return true;
}

/// adjusts view to bounding box of all instances
void volume_slicer::auto_adjust_view()
{
	if (ensure_view_ptr()) {
		view_ptr->set_scene_extent(box3(-extent, extent));
		view_ptr->set_default_view();
		view_ptr->set_z_near(extent.length() / dimensions.length());
		view_ptr->set_z_far(5.0f*extent.length());
		std::cout << "z_near = " << view_ptr->get_z_near() << ", z_far=" << view_ptr->get_z_far() << std::endl;
	}
}

/// self reflection is used to publish standard member variables to the set and get mechanism
bool volume_slicer::self_reflect(cgv::reflect::reflection_handler& rh)
{
	return
		rh.reflect_member("file_name", file_name) &&
		rh.reflect_member("show_box", show_box);
}

template <typename T>
cgv::math::fvec<T, 3> rotate(const cgv::math::fvec<T, 3>& v, const cgv::math::fvec<T, 3>& n, T a)
{
	cgv::math::fvec<T, 3> vn = dot(n, v)*n;
	return vn + cos(a)*(v - vn) + sin(a)*cross(n, v);
}

/// return the point under the mouse pointer in world coordinates
bool volume_slicer::get_picked_point(int x, int y, vec3& p_pick_world)
{
	// analyze the mouse location
	cgv::render::context& ctx = *get_context();
	cgv::render::context::mat_type* DPV_ptr, *DPV_other_ptr;
	int x_other, y_other, vp_col_idx, vp_row_idx, vp_width, vp_height;
	int eye_panel = view_ptr->get_DPVs(x, y, ctx.get_width(), ctx.get_height(), &DPV_ptr, &DPV_other_ptr, &x_other, &y_other, &vp_col_idx, &vp_row_idx, &vp_width, &vp_height);

	// get the possibly two (if stereo is enabled) different device z-values
	double z = ctx.get_z_D(x, y);
	double z_other = ctx.get_z_D(x_other, y_other);
	//  unproject to world coordinates with smaller (closer to eye) z-value one	
	if (z <= z_other) {
		if (DPV_ptr->ncols() != 4)
			return false;
		// use conversion to (double*) operator to map cgv::math::vec<double> to cgv::math::fvec<float,3>
		p_pick_world = (double*)ctx.get_point_W(x, y, z, *DPV_ptr);
	}
	else {
		if (DPV_other_ptr->ncols() != 4)
			return false;
		p_pick_world = (double*)ctx.get_point_W(x_other, y_other, z_other, *DPV_other_ptr);
	}
	return true;
}

///
void volume_slicer::peek_voxel_values(int x, int y)
{
	//  unproject pixel position of mouse
	volume::point_type p_pick_world;
	if (!get_picked_point(x, y, p_pick_world))
		return;

	// construct string stream
	std::stringstream ss;

	// convert point to texture coordinates of instance and add point
	box3 B(-0.5f*extent, 0.5f*extent);
	if (B.inside(p_pick_world)) {
		volume::point_type p_texture = texture_from_world_coordinates(p_pick_world);
		volume::point_type p_voxel = voxel_from_texture_coordinates(p_texture);
		current_voxel = p_voxel;
		current_block = block_from_voxel_coordinates(p_voxel);
		ss << p_texture << "\n" << p_voxel << "\n" << current_voxel;
	}
	voxel_value_info = ss.str();
	last_x = x;
	last_y = y;
	post_redraw();
}

bool volume_slicer::retrieve_block(const std::string& input_path, const std::string& output_path, ivec3& voxel_at)
{
	///REMOVE
	return true;

	ivec3 nr_blocks(
		unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), 
		unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), 
		unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	ivec3 block(
		unsigned(ceil(float(voxel_at(0)) / block_dimensions(0))),
		unsigned(ceil(float(voxel_at(1)) / block_dimensions(1))),
		unsigned(ceil(float(voxel_at(2)) / block_dimensions(2))));
	size_t block_size = (block_dimensions(1) + overlap(1))*((block_dimensions(0) + overlap(0))*cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB).get_entry_size()) * (block_dimensions(2) + overlap(2));

	std::cout << "Trying to get block at: " << block(0) <<  ", " <<  block(1) << ", " << block(2) << ", " << std::endl;

	std::stringstream ss;
	ss << input_path << "/level_00_blockslice_" << std::setw(3) << std::setfill('0') << block(2) << ".bvx";

	// compute block index and pointer to block data
	unsigned bi = block(1) * nr_blocks(0) + block(0);
	char* block_ptr = new char[block_size];

	// Not tested: throws exception at fread in some points 
	// Read block from dataset and write tiff from block slices
	// read only the block

	FILE* fp = fopen(ss.str().c_str(), "rb");
	fseek(fp, bi * block_size, SEEK_SET);
	

	if (fread(block_ptr, block_size, 1, fp) == 1)
		fclose(fp);
	else
		return false;

	return write_tiff_block(output_path + "/block_at_" + std::to_string(block(0)) + "_" + std::to_string(block(1)) + "_" + std::to_string(block(2)), block_ptr, "");
}

//! write a block located at \c data_ptr to a tiff file
bool volume_slicer::write_tiff_block(const std::string& file_name, const char* data_ptr, const std::string& options)
{
	// setup output format 
	cgv::data::data_format df;
	df.set_component_format(cgv::data::component_format(cgv::type::info::TI_UINT8, cgv::data::CF_RGB));
	df.set_width(block_dimensions(0) + 1);
	df.set_height(block_dimensions(1) + 1);
	cgv::media::image::image_writer iw(file_name + ".tif");
	iw.multi_set(options);
	size_t slice_size = df.get_nr_bytes();
	for (unsigned i = 0; i<(unsigned)(block_dimensions(2) + overlap(2)); ++i) {
		cgv::data::const_data_view dv(&df, data_ptr + slice_size*i);
		if (!iw.write_image(dv)) {
			std::cerr << "could not write slice " << i << " to tiff file " << file_name << std::endl;
			return false;
		}
	}
	return iw.close();
}

/// overload and implement this method to handle events
bool volume_slicer::handle(cgv::gui::event& e)
{
	// check for key events
	if (e.get_kind() == cgv::gui::EID_KEY) {
		cgv::gui::key_event& ke = static_cast<cgv::gui::key_event&>(e);
		switch (ke.get_key()) {
			// view all when space key pressed
		case cgv::gui::KEY_Space:
			auto_adjust_view();
			return true;
			// handle key press events for 'X'-key
		case 'X':
			// we only want to use event if user holds CTRL and SHIFT
			if (ke.get_modifiers() != cgv::gui::EventModifier(cgv::gui::EM_CTRL | cgv::gui::EM_SHIFT))
				// return false to allow other event_handlers to process event
				return false;
			// set normal direction to x-direction
			slice_normal_tex = vec3(1, 0, 0);
			// use on_set to notify of member value changes, such that gui is updated and 3D view redrawn
			// as there are controls for each components, on_set needs to be called for each vector component
			on_set(&slice_normal_tex[0]);
			on_set(&slice_normal_tex[1]);
			on_set(&slice_normal_tex[2]);
			return true;
		case 'Y':
			if (ke.get_modifiers() != cgv::gui::EventModifier(cgv::gui::EM_CTRL | cgv::gui::EM_SHIFT))
				return false;
			slice_normal_tex = vec3(0, 1, 0);
			on_set(&slice_normal_tex[0]);
			on_set(&slice_normal_tex[1]);
			on_set(&slice_normal_tex[2]);
			return true;
		case 'Z':
			if (ke.get_modifiers() != cgv::gui::EventModifier(cgv::gui::EM_CTRL | cgv::gui::EM_SHIFT))
				return false;
			slice_normal_tex = vec3(0, 0, 1);
			on_set(&slice_normal_tex[0]);
			on_set(&slice_normal_tex[1]);
			on_set(&slice_normal_tex[2]);
			return true;
		}
		return false;
	}
	// check for mouse event
	if (e.get_kind() == cgv::gui::EID_MOUSE) {
		cgv::gui::mouse_event& me = static_cast<cgv::gui::mouse_event&>(e);

		peek_voxel_values(me.get_x(), me.get_y());
		
		// we are interested only in mouse events with the Ctrl-modifier and no other modifier pressed
		if (me.get_modifiers() != cgv::gui::EM_CTRL)
			return false;
	
		switch (me.get_action()) {
			
		case cgv::gui::MA_RELEASE:
			// support drag and drop for left and right button
			if (me.get_button() == cgv::gui::MB_RIGHT_BUTTON) {
				
			}
			break;
		case cgv::gui::MA_PRESS:
			if (me.get_button() == cgv::gui::MB_LEFT_BUTTON || me.get_button() == cgv::gui::MB_RIGHT_BUTTON)
				return true;
			else
				return false;
		case cgv::gui::MA_DRAG:
			// check for non zero drag vector in pixel coordinates 
			// as well as availability of view_ptr and context_ptr 
			if ((me.get_dx() != 0 || me.get_dy() != 0) && ensure_view_ptr() && get_context() != 0) {
				// retrieve view coordinate system in world coordinates where z points towards eye and x/y to right/top of screen
				cgv::render::view::vec_type view_x_world, view_y_world, view_z_world;
				view_ptr->put_coordinate_system(view_x_world, view_y_world, view_z_world);
				// compute drag vector in view coordinates in range [-1,1]²
				// negate pixel y coordinates to make it point upward as the eye coordinate y-direction
				double drag_x = (double)me.get_dx() / get_context()->get_width();
				double drag_y = -(double)me.get_dy() / get_context()->get_height();
				// ctrl-left button rotates plane normal
				if (me.get_button_state() == cgv::gui::MB_LEFT_BUTTON) {
					// scale x/y-pixel motion to range [-1,1] and 
					// compute rotation axis orthogonal to drag vector in world coordinates
					vec3 axis_world = drag_x*view_y_world - drag_y*view_x_world;
					// normalize axis and use resulting length to scale rotation angle  
					double rotation_angle = 5.0 * axis_world.normalize() * rotate_sensitivity;
					// update slice_normal_tex by rotation around axis
					slice_normal_tex = rotate(slice_normal_tex, axis_world, (float)rotation_angle);
					// ensure that user interface has current values
					for (unsigned c = 0; c < 3; ++c)
						update_member(&slice_normal_tex[c]);
					// ensure that scene is redrawn with new slice_normal_tex
					post_redraw();
					return true;
				}
				// ctrl-right button moves plane
				if (me.get_button_state() == cgv::gui::MB_RIGHT_BUTTON) {
					// compute drag vector in world coordinates
					vec3 drag_world = drag_x*view_x_world + drag_y*view_y_world;
					// update distance by scaling world size (y_extent_at_focus) to 
					// texture coordinates (division by extent.length)
					slice_center_distance_tex += (float)(10.0f*dot(slice_normal_tex, drag_world)*
						view_ptr->get_y_extent_at_focus() / extent.length() / translation_sensitivity);
					// ensure that slice distance is in valid range
					if (slice_center_distance_tex > 0.5f)
						slice_center_distance_tex = 0.5f;
					else if (slice_center_distance_tex < -0.5f)
						slice_center_distance_tex = -0.5f;
					// ensure that user interface and rendering are up to date
					update_member(&slice_center_distance_tex);
					post_redraw();
					return true;
				}
				return false;

			}
		}
	}
	return false;
}

/// overload to stream help information to the given output stream
void volume_slicer::stream_help(std::ostream& os)
{
	os << "volume_slicer:\a\nSpace ... auto adjust view; Ctrl-Shift-X/Y/Z ... set slice normal\nCtrl-Left/Right Button ... rotate/move slice\b\n";
}
/// stream statistical information about volume
void volume_slicer::stream_stats(std::ostream& os)
{
	os << "volume_slicer: slice_normal=[" << slice_normal_tex << "], distance=" << slice_center_distance_tex << std::endl;
}
// read volume from sliced or single file representations
bool volume_slicer::open_block_volume(const std::string& directory_name)
{
	return false;
}

/// read regular volume file
bool volume_slicer::open_volume(const std::string& _file_name)
{
	// try to read volume from slice based or regular volume file 
	bool success = false;
	std::string ext = cgv::utils::to_upper(cgv::utils::file::get_extension(_file_name));
	if (ext == "SVX") {
		success = read_from_sliced_volume(file_name, V);
		if (success)
			extent = vec3(float(V.get_dimensions()(0)) / V.get_dimensions()(2), float(V.get_dimensions()(1)) / V.get_dimensions()(2), 1);
	}
	else {
		volume_info info;
		success = read_volume(file_name, V, &info);
		if (success)
			extent = info.extent;
	}

	// abort in case of failure 
	if (!success)
		return false;

	// copy dimensions
	dimensions = V.get_dimensions();

	// transfer volume to 3d texture 
	if (get_context() && get_context()->make_current()) {
		if (data_texture.is_created())
			data_texture.destruct(*get_context());
		if (!data_texture.create(*get_context(), V.get_data_view())) {
			cgv::data::data_format df(V.get_format());
			df.set_depth(df.get_depth() / 2);
			cgv::data::data_view dv(&df, V.get_data_ptr<char>());
			dimensions(2) = df.get_depth();
			if (!data_texture.create(*get_context(), dv)) {
				std::cerr << "error tex create" << std::endl;
			}
		}
	}

	// adjust view
	auto_adjust_view();
	return success;
}


/// convert world to texture coordinates
volume_slicer::vec3 volume_slicer::texture_from_world_coordinates(const vec3& p_world) const
{
	vec3 p_texture = (p_world + 0.5f*extent) / extent;
	return p_texture;
}
/// convert texture to voxel coordinates
volume_slicer::vec3 volume_slicer::voxel_from_texture_coordinates(const vec3& p_texture) const
{
	vec3 p_voxel = p_texture * dimensions;
	return p_voxel;
}

/// convert texture to world coordinates
volume_slicer::vec3 volume_slicer::world_from_texture_coordinates(const vec3& p_texture) const
{
	vec3 p_world = p_texture * extent - 0.5f*extent;
	return p_world;
}

/// convert texture to world coordinates
volume_slicer::vec3 volume_slicer::world_from_texture_normals(const vec3& n_texture) const
{
	vec3 n_world = n_texture / extent;
	n_world.normalize();
	return n_world;
}

/// convert voxel to texture coordinates
volume_slicer::vec3 volume_slicer::texture_from_voxel_coordinates(const vec3& p_voxel) const
{
	vec3 p_texture = p_voxel / dimensions;
	return p_texture;
}

/// convert voxel to block coordinates
volume_slicer::vec3 volume_slicer::block_from_voxel_coordinates(const vec3& p_voxel) const
{
	vec3 offset = overlap;
	offset *= 0.5f;
	vec3 p_block = (p_voxel - offset) / block_dimensions;
	return p_block;
}
/// convert voxel to block coordinates
volume_slicer::vec3 volume_slicer::voxel_from_block_coordinates(const vec3& p_block) const
{
	vec3 offset = overlap;
	offset *= 0.5f;
	vec3 p_voxel = p_block * block_dimensions + offset;
	return p_voxel;
}

///
std::string volume_slicer::get_type_name() const
{
	return "volume_slicer";
}
// extend method of volume_drawable to find and set view
bool volume_slicer::init(cgv::render::context& ctx)
{
	// generate a font for voxel labeling
	ctx.enable_font_face(cgv::media::font::find_font("Arial")->get_font_face(cgv::media::font::FFA_BOLD), 20);
	// ensure white background color
	ctx.set_bg_clr_idx(4);

	return true;
}

///
void volume_slicer::init_frame(cgv::render::context& ctx)
{
	if (!tex_slice_prog.is_created())
		tex_slice_prog.build_program(ctx, "tex_slicer.glpr", true);

	// upload more texture data for block based volumes

	/*
	if (data_texture.is_created()) {
	cgv::data::data_format df(170, 170, 170, cgv::type::info::TI_UINT8, cgv::data::CF_R);
	cgv::data::data_view dv(&df);
	cgv::type::uint8_type* ptr = dv.get_ptr<cgv::type::uint8_type>();
	std::fill(ptr, ptr + 170 * 170 * 170, 50);
	data_texture.replace(ctx, 100, 80, 60, dv, 0);
	}*/
}

void volume_slicer::extract_surface()
{
	clear_surface();
	box3 box(-0.5f*extent, 0.5f*extent);
	cgv::media::mesh::marching_cubes<float, float> mc(*this, this);
	sm_ptr = &mc;
	mc.extract(iso_value, box, dimensions(0) / down_sampling_factor, dimensions(1) / down_sampling_factor, dimensions(2) / down_sampling_factor,
		dimensions(0)*dimensions(1)*dimensions(2) > 10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor);
	sm_ptr = 0;
	std::cout << "iso_surface has " << positions.size() << " vertices and " << corner_indices.size() / 3 << " triangles" << std::endl;
	post_redraw();
}

/// remove all elements from surface
void volume_slicer::clear_surface()
{
	positions.clear();
	normals.clear();
	texture_coords.clear();
	corner_indices.clear();
	post_redraw();
}

/// return the value of a given voxel scaled to [0,1]
float volume_slicer::get_voxel_value(const ivec3& voxel) const
{
	switch (V.get_component_type()) {
	case cgv::type::info::TI_UINT8:
		return float(*V.get_voxel_ptr<cgv::type::uint8_type>(voxel(0), voxel(1), voxel(2))) / 255;
	case cgv::type::info::TI_UINT16:
		return float(*V.get_voxel_ptr<cgv::type::uint8_type>(voxel(0), voxel(1), voxel(2))) / 65535;
	}
	return 0;
}

/// interface for evaluation of the multivariate function
float volume_slicer::evaluate(const pnt_type& p) const
{
	// convert to voxel coordinates
	vec3 p_voxel = voxel_from_texture_coordinates(texture_from_world_coordinates(&p(0)));
	// quantize
	ivec3 voxel = p_voxel;
	// check if voxel is valid
	if (voxel(0) < 0 || voxel(1) < 0 || voxel(2) < 0)
		return 0.0f;
	if (voxel(0) >= dimensions(0) - 1 || voxel(1) >= dimensions(1) - 1 || voxel(2) >= dimensions(2) - 1)
		return 0.0f;
	// compute fractional part of voxel location needed for trilinear interpolation
	vec3 frac = voxel;
	frac = p_voxel - frac;

	// collect eight nearest voxels
	float v000 = get_voxel_value(voxel);
	float v100 = get_voxel_value(voxel + ivec3(1, 0, 0));
	float v010 = get_voxel_value(voxel + ivec3(0, 1, 0));
	float v110 = get_voxel_value(voxel + ivec3(1, 1, 0));
	float v001 = get_voxel_value(voxel + ivec3(0, 0, 1));
	float v101 = get_voxel_value(voxel + ivec3(1, 0, 1));
	float v011 = get_voxel_value(voxel + ivec3(0, 1, 1));
	float v111 = get_voxel_value(voxel + ivec3(1, 1, 1));

	// trilinear interpolation
	float v00 = (1 - frac(2))*v000 + frac(2)*v001;
	float v10 = (1 - frac(2))*v100 + frac(2)*v101;
	float v01 = (1 - frac(2))*v010 + frac(2)*v011;
	float v11 = (1 - frac(2))*v110 + frac(2)*v111;
	float v0 = (1 - frac(1))*v00 + frac(1)*v01;
	float v1 = (1 - frac(1))*v10 + frac(1)*v11;
	return (1 - frac(0))*v0 + frac(0)*v1;
}

/// called when a new vertex is generated
void volume_slicer::new_vertex(unsigned int vertex_index)
{
	positions.push_back(sm_ptr->vertex_location(vertex_index));
	normals.push_back(sm_ptr->vertex_normal(vertex_index));
	texture_coords.push_back(texture_from_world_coordinates(sm_ptr->vertex_location(vertex_index)));
}
/// announces a new polygon defines by the vertex indices stored in the given vector
void volume_slicer::new_polygon(const std::vector<unsigned int>& vertex_indices)
{
	// add corner indices and compute face normal
	vec3 face_nml(0, 0, 0);
	unsigned i;
	for (i = 0; i < vertex_indices.size(); ++i) {
		corner_indices.push_back(vertex_indices[i]);
		face_nml += cross(positions[vertex_indices[i]], positions[vertex_indices[(i + 1) % vertex_indices.size()]]);
	}
	face_nml.normalize();
	// distribute face normal to vertex normals
	for (i = 0; i < vertex_indices.size(); ++i)
		normals[vertex_indices[i]] += face_nml;
}

/// drop the currently first vertex that has the given global vertex index
void volume_slicer::before_drop_vertex(unsigned int vertex_index)
{
	normals[vertex_index].normalize();
	if (negate_normals)
		normals[vertex_index] = -normals[vertex_index];
}

///
void volume_slicer::draw_box(cgv::render::context& ctx, const box3& B, const cgv::media::illum::phong_material& material)
{
	// draw backfaces next
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	ctx.enable_material(material);
	ctx.tesselate_box(B, true);
	ctx.disable_material(material);
	glCullFace(GL_BACK);
	glDisable(GL_CULL_FACE);
}

///
void volume_slicer::draw_wire_box(cgv::render::context& ctx, const box3& B, const vec3& color)
{
	glColor3fv(color);
	glDisable(GL_CULL_FACE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	ctx.tesselate_box(B, false);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

/// draw a voxel with a wire frame box
void volume_slicer::draw_voxel(cgv::render::context& ctx, const ivec3& voxel, const vec3& color)
{
	box3 B(world_from_texture_coordinates(texture_from_voxel_coordinates(voxel)),
		world_from_texture_coordinates(texture_from_voxel_coordinates(voxel + ivec3(1, 1, 1))));
	draw_wire_box(ctx, B, color);
}

/// computes the distance taking into account the 
float volume_slicer::compute_distance_to_slice_tex(const vec3& p) const
{
	return cgv::math::dot(p - vec3(0.5f, 0.5f, 0.5f), slice_normal_tex) - slice_center_distance_tex;
}

bool volume_slicer::is_block_intersected(const box3& B_tex) 
{
	// intersection is calculated by analyzing the dot product of the corners and the normal of the slice
	// if all results share the same sign, it's not interesected. 

	float dist_corner;
	bool sign;
	bool first_time = true;

	// iterate over all the block's corners
	for (unsigned i = 0; i < 8; i++)
	{
		dist_corner = compute_distance_to_slice_tex(B_tex.get_corner(i));
		if (first_time) {
			sign = dist_corner < 0;
			first_time = false;
		} else {
			if ((dist_corner < 0) != sign) {
				return true;
			}
		}
	}
	return false;
}

float volume_slicer::block_distance(const box3& B_tex) 
{
	float dist_corner =0;
	// iterate over all the block's corners, adding the distance value
	for (unsigned i = 0; i < 8; i++){
		dist_corner += compute_distance_to_slice_tex(B_tex.get_corner(i));
	}
	return dist_corner;
}

/// draw a block with a wire frame box
void volume_slicer::draw_block(cgv::render::context& ctx, const ivec3& block, const vec3& color, const vec3& overlap_color, bool write_tiff)
{
	box3 B(world_from_texture_coordinates(texture_from_voxel_coordinates(voxel_from_block_coordinates(block))),
		world_from_texture_coordinates(texture_from_voxel_coordinates(voxel_from_block_coordinates(block + ivec3(1, 1, 1)))));
	draw_wire_box(ctx, B, color);

	vec3 offset = overlap;
	offset *= 0.5f;
	box3 B_overlap(
		world_from_texture_coordinates(texture_from_voxel_coordinates(voxel_from_block_coordinates(block) - offset)),
		world_from_texture_coordinates(texture_from_voxel_coordinates(voxel_from_block_coordinates(block + ivec3(1, 1, 1)) + offset)));
	draw_wire_box(ctx, B_overlap, overlap_color);

	//texture coordinates
	box3 B_tex(texture_from_voxel_coordinates(voxel_from_block_coordinates(block)),
		texture_from_voxel_coordinates(voxel_from_block_coordinates(block + ivec3(1, 1, 1))));

	if (is_block_intersected(B_tex)) {
		draw_wire_box(ctx, B, color);
		draw_wire_box(ctx, B_overlap, overlap_color);

		if (write_tiff) {
			std::string input_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices";
			std::string output_path = "D:/Users/JMendez/Documents/cgv-hdp-cr-local/data/visual_human/labeled/innerorgans/rgb_enlarged_slices/blocks";

			if (retrieve_block(input_path, output_path, current_voxel))
				std::cout << "successfully retrieved and wrote block in tiff format" << std::endl;
			else
				std::cout << "failed to retrieve block from block slices" << std::endl;
			
		}
	}
}

static bool abs_compare(float a, float b)
{
	return (std::abs(a) < std::abs(b));
}

void volume_slicer::update_intersected_blocks(cgv::render::context& ctx) {
	
	ivec3 nr_blocks(unsigned(ceil(float(dimensions(0)) / block_dimensions(0))), unsigned(ceil(float(dimensions(1)) / block_dimensions(1))), unsigned(ceil(float(dimensions(2)) / block_dimensions(2))));
	intersected_blocks.clear();

	int max_ind = std::distance(std::begin(slice_normal_tex), std::max_element(std::begin(slice_normal_tex), std::end(slice_normal_tex), abs_compare));
	
	std::vector<int> projected_dimensions_indices(2);
	std::size_t n(0);

	std::generate(projected_dimensions_indices.begin(), projected_dimensions_indices.end(),
	[&]{
		if (n == max_ind) { n++; }
		return n++; 
	});

	std::cout <<max_ind << std::endl;

	for (int projected_dim0 = 0; projected_dim0 < nr_blocks(projected_dimensions_indices.at(0)); projected_dim0++) {
		for (int projected_dim1 = 0; projected_dim1 < nr_blocks(projected_dimensions_indices.at(1)); projected_dim1++) {
			
			int l = 0, r = nr_blocks(max_ind);
			while (l <= r) {

				int m = l + (r - l) / 2;
			
				ivec3 block = get_box_indices_from_projection(max_ind, m, projected_dim0, projected_dim1);

				//texture coordinates
				box3 B_tex(texture_from_voxel_coordinates(voxel_from_block_coordinates(block)),
					texture_from_voxel_coordinates(voxel_from_block_coordinates(block + ivec3(1, 1, 1))));


				if (is_block_intersected(B_tex)) {
					intersected_blocks.push_back(block);
				}
				
				//if the deleted dimension is negative distance value should be inverted
				bool question = slice_normal_tex(max_ind) > 0 ? block_distance(B_tex) < 0 : block_distance(B_tex) > 0;

				if (question)
					l = m + 1;
				else
					r = m - 1;
					
			
			}
		}
	}
}

volume_slicer::ivec3 volume_slicer::get_box_indices_from_projection(int max_ind, int var_dim, int dim0, int dim1) {
	int x, y, z;
	switch (max_ind) {
		case 0: {
			x = var_dim;
			y = dim0;
			z = dim1;
			break;
		}
		case 1: {
			x = dim0;
			y = var_dim;
			z = dim1;
			break;
		}
		case 2: {
			x = dim0;
			y = dim1;
			z = var_dim;
			break;
		}
	}
	return ivec3(x,y,z);
}

void volume_slicer::draw_blocks_in_plane(cgv::render::context& ctx, const vec3& color, const vec3& overlap_color) {
	for (int i = 0; i < intersected_blocks.size(); i++) {
		draw_block(ctx, intersected_blocks[i], color, overlap_color, false);
	}
}

/// draw a block with a wire frame box
void volume_slicer::draw_surface(cgv::render::context& ctx)
{
	if (corner_indices.empty())
		return;
	if (show_surface_wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	glDisable(GL_COLOR_MATERIAL);
	ctx.enable_material(surface_mat);
	glVertexPointer(3, GL_FLOAT, 0, &positions[0]);
	glEnableClientState(GL_VERTEX_ARRAY);
	glNormalPointer(GL_FLOAT, 0, &normals[0]);
	glEnableClientState(GL_NORMAL_ARRAY);
	glDrawElements(GL_TRIANGLES, corner_indices.size(), GL_UNSIGNED_INT, &corner_indices[0]);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	ctx.disable_material(surface_mat);
	if (show_surface_wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

///
void volume_slicer::draw(cgv::render::context& ctx)
{
	// compute axis aligned 3d box of volume
	box3 B(-0.5f*extent, 0.5f*extent);
	if (show_box) {
		draw_wire_box(ctx, B, vec3(0, 1, 1));
		draw_box(ctx, B, box_mat);
	}
	if (show_voxel)
		draw_voxel(ctx, current_voxel, vec3(1, 0, 0));
	update_intersected_blocks(ctx);
	if (show_block) {
		draw_block(ctx, current_block, vec3(1, 0.5, 0), vec3(0, 1, 0), true);
		draw_blocks_in_plane(ctx, vec3(1, 0, 1), vec3(0, 0, 1));
	}
	if (show_surface)
		draw_surface(ctx);

	// render slice
	// compute two directions x/y_tex orthogonal to normal direction
	vec3 x_tex = cross(slice_normal_tex, vec3(1, 0, 0));
	if (x_tex.length() < 1e-6)
		x_tex = cross(slice_normal_tex, vec3(0, 1, 0));
	x_tex.normalize();
	vec3 y_tex = cross(slice_normal_tex, x_tex);
	y_tex.normalize();
	// scale direction vectors such to be longer than box extent
	x_tex *= 3.0f;
	y_tex *= 3.0f;
	// compute point on slice in texture coordinates 
	vec3 c_tex = slice_center_distance_tex*slice_normal_tex + vec3(0.5f, 0.5f, 0.5f);
	// define clipping planes in texture coordinates
	int i;
	for (i = 0; i < 6; ++i) {
		double cp[4] = { 0, 0, 0, 0 };
		if (i < 3) {
			cp[i] = 1;
			cp[3] = -B.get_min_pnt()(i);
		}
		else {
			cp[i - 3] = -1;
			cp[3] = B.get_max_pnt()(i - 3);
		}
		glClipPlane(GL_CLIP_PLANE0 + i, cp);
		glEnable(GL_CLIP_PLANE0 + i);
	}
	// disable culling
	glDisable(GL_CULL_FACE);
	// enable volume texture on texture unit 0
	if (data_texture.is_created()) {
		data_texture.set_mag_filter(interpolate ? cgv::render::TF_LINEAR : cgv::render::TF_NEAREST);
		data_texture.enable(ctx, 0);
	}
	// enable shader program and take volume texture from texture unit 0
	if (tex_slice_prog.is_linked()) {
		tex_slice_prog.enable(ctx);
		tex_slice_prog.set_uniform(ctx, "volume", 0);
	}
	// construct texture coordinates of quad that spans slice
	static std::vector<vec3> tex_coords(size_t(4));
	tex_coords[0] = c_tex + x_tex;
	tex_coords[1] = c_tex + y_tex;
	tex_coords[2] = c_tex - x_tex;
	tex_coords[3] = c_tex - y_tex;
	// map texture coordinates of quad corners to positions in world coordinates
	static std::vector<vec3> positions(size_t(4));
	for (unsigned i = 0; i<4; ++i)
		positions[i] = extent*(tex_coords[i] - vec3(0.5f, 0.5f, 0.5f));
	// draw quad using vertex array pointers which is a bit deprecated but still ok 
	glVertexPointer(3, GL_FLOAT, 0, &positions.front());
	glEnableClientState(GL_VERTEX_ARRAY);
	glTexCoordPointer(3, GL_FLOAT, 0, &tex_coords.front());
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glDrawArrays(GL_QUADS, 0, 4);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);

	// disable shader program
	if (tex_slice_prog.is_linked())
		tex_slice_prog.disable(ctx);
	// disable 3d texture
	if (data_texture.is_created())
		data_texture.disable(ctx);
	// disable clipping planes
	for (i = 0; i < 6; ++i)
		glDisable(GL_CLIP_PLANE0 + i);
}

/// draw textual information here
void volume_slicer::after_finish(cgv::render::context& ctx)
{
	if (!voxel_value_info.empty()) {
		glDisable(GL_LIGHTING);
		glDisable(GL_DEPTH_TEST);
		ctx.push_pixel_coords();
		float f = (float)last_x / ctx.get_width();
		glColor3f(1, 1, 0.3f);
		ctx.set_cursor(last_x - (int)(f*ctx.get_current_font_face()->measure_text_width(voxel_value_info, ctx.get_current_font_size())), last_y - 4);
		ctx.output_stream() << voxel_value_info.c_str();
		ctx.output_stream().flush();

		glColor3f(0.2f, 0.2f, 1);
		ctx.set_cursor(last_x + 2 - (int)(f*ctx.get_current_font_face()->measure_text_width(voxel_value_info, ctx.get_current_font_size())), last_y - 2);
		ctx.output_stream() << voxel_value_info.c_str();
		ctx.output_stream().flush();

		ctx.pop_pixel_coords();
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_LIGHTING);
	}
}

/// called to destruct the rendering objects
void volume_slicer::clear(cgv::render::context& ctx)
{
	tex_slice_prog.destruct(ctx);
	indexed_tex_slice_prog.destruct(ctx);
	data_texture.destruct(ctx);
	index_texture.destruct(ctx);
}

///
void volume_slicer::create_gui()
{
	// all add_member_control and add_gui functions will add a callback to the gui elements
	// that calls for every value change the on_set(member_ptr) method of the volume_slicer 
	// with member_ptr pointing to the member variable whose value changed


	// add a heading at the top of the volume slicer gui, increasing level decreases the font size, 
	// font size of level=2 is larger than other gui fonts
	add_decorator("volume slicer", "heading", "level=2");
	// begin a gui node that can be opened and closed with heading "Slicing" of next level (3)
	// use member slice_center_distance_tex to hash the node status and open node initially (true in 3rd argument)
	if (begin_tree_node("Slicing", slice_center_distance_tex, true, "level=3")) {
		// tabify gui
		align("\a");
		// add value_slider control for slice_center_distance_tex member specifying range [-1,1]
		add_member_control(this, "slice_center_distance_tex", slice_center_distance_tex, "value_slider", "min=-1;max=1;step=0.0001;ticks=true");
		// add a compount gui of type "direction" for the vector valued member slice_normal_tex
		// the direction gui ensures that the slice_normal_tex is of length one
		add_gui("slice_normal_tex", slice_normal_tex, "direction");
		// untabify gui
		align("\b");
		// end gui node definition specifying the same member used for hashing the node status 
		end_tree_node(slice_center_distance_tex);
	}
	// start tree node to interact with volume
	if (begin_tree_node("Volume", dimensions, false, "level=3")) {
		align("\a");
		// add gui of type "file_name" for member file_name
		// file_name gui adds an open button that opens a file dialog to query new file_name
		// the options 'title' and 'filter' configure the file dialog
		add_gui("file_name", file_name, "file_name", "title='open volume';filter='Volume Files(vox,qim,tif,avi) :*.vox;*.qim;*.tif;*.avi|All Files:*.*'");
		// add gui for the vector of pixel counts, where "dimensions" is used only in label of gui element for first vector component
		// using view as gui_type will only show the values but not allow modification
		// by align=' ' the component views are arranged with a small space in one row
		// by options='...' specify options for the component views
		// w=50 specifies the width of each component view in 50 pixels
		add_gui("dimensions", dimensions, "vector", "main_label='first';align=' ';gui_type='view';options='w=50'"); align("\n");
		// add gui for extent vector with component controls of type value_input
		add_gui("extent", extent, "vector", "main_label='first';gui_type='value_input';options='w=50;min=0;max=10';align=' '");
		// start new gui row with smaller Y-spacing by reducing height of gui rows by 6 pixels, going to next row and restoring row height
		align("%Y-=6\n%Y+=6");
		// add another gui for the extent vector with component guis of type wheel within one row
		add_gui("box_extent", extent, "vector", "gui_type='wheel';options='align=\"B\";w=50;h=10;step=0.01;min=0;max=10;log=true';align=' '");
		// start new row with regular spacing
		align("\n");
		add_member_control(this, "interpolate", interpolate, "toggle");
		align("\b");
		end_tree_node(dimensions);
	}
	// start tree node to interact with volume
	if (begin_tree_node("Blocks", block_dimensions, false, "level=3")) {
		align("\a");
		// w=50 specifies the width of each component view in 50 pixels
		add_gui("block_dimensions", block_dimensions, "vector", "main_label='first';align=' ';gui_type='value_input';options='w=50;min=1;max=1024;step=1'"); align("\n");
		// w=50 specifies the width of each component view in 50 pixels
		add_gui("overlap", overlap, "vector", "main_label='first';align=' ';gui_type='value_input';options='w=50;min=0;max=1;step=1'"); align("\n");
		align("\b");
		end_tree_node(block_dimensions);
	}
	// start tree node to interact with volume
	if (begin_tree_node("Iso Surface", down_sampling_factor, false, "level=3")) {
		align("\a");
		// w=50 specifies the width of each component view in 50 pixels
		connect_copy(add_button("extract")->click, cgv::signal::rebind(this, &volume_slicer::extract_surface));
		add_member_control(this, "iso_value", iso_value, "value_slider", "min=0;max=1;ticks=true");
		add_member_control(this, "negate_normals", negate_normals, "toggle");
		add_member_control(this, "down_sampling_factor", down_sampling_factor, "value_slider", "min=1;max=32;ticks=true;log=true");
		connect_copy(add_button("clear")->click, cgv::signal::rebind(this, &volume_slicer::clear_surface));
		align("\b");
		end_tree_node(down_sampling_factor);
	}
	if (begin_tree_node("Rendering", show_box, false, "level=3")) {
		align("\a");
		add_member_control(this, "show_box", show_box, "check");
		add_member_control(this, "show_voxel", show_voxel, "check");
		add_member_control(this, "show_block", show_block, "check");
		add_gui("box_mat", box_mat);
		add_member_control(this, "show_surface", show_surface, "check");
		add_member_control(this, "show_surface_wireframe", show_surface_wireframe, "check");
		add_gui("surface_mat", surface_mat);
		align("\b");
		end_tree_node(show_box);
	}
	if (begin_tree_node("Interaction", rotate_sensitivity, false, "level=3")) {
		align("\a");
		add_member_control(this, "rotate_sensitivity", rotate_sensitivity, "value_slider", "min=0.1;max=10;ticks=true;log=true");
		add_member_control(this, "translation_sensitivity", translation_sensitivity, "value_slider", "min=0.1;max=10;ticks=true;log=true");
		align("\b");
		end_tree_node(rotate_sensitivity);
	}
}


// 
void volume_slicer::on_set(void* member_ptr)
{
	// in case that file_name changed, read new volume
	if (member_ptr == &file_name) {
		open_volume(file_name);
		post_redraw();
		post_recreate_gui();
	}
	// update gui of member
	update_member(member_ptr);
	// for all value changes schedule a redraw of the 3D window
	post_redraw();
}

#include <cgv/base/register.h>

extern cgv::base::object_registration<volume_slicer> vol_slicer_reg("volume_slicer");
