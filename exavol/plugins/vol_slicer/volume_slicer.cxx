#include "volume_slicer.h"
#include <vol_data/volume_io.h>
#include <vol_data/sliced_volume_io.h>
#include <cgv/utils/advanced_scan.h>
#include <cgv/base/find_action.h>
#include <cgv/utils/scan.h>
#include <cgv/utils/file.h>
#include <cgv_gl/gl/gl.h>
#include <cgv_gl/gl/gl_tools.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/mouse_event.h>
#include <cgv/media/mesh/marching_cubes.h>
#include <cgv/media/mesh/dual_contouring.h>
#include <cgv/media/mesh/cuberille.h>
#include <cgv/media/image/image_writer.h>
#include <fstream>

/// set default values (no wireframe, no culling, no face normals, no negation, no two sided)
surface_render_style::surface_render_style()
{
	diffuse_mapping = DCM_MATERIAL;
	illumination = true;

	show_wireframe = false;
	wireframe_width = 1;
	culling_mode = cgv::render::CM_OFF;
	face_normals = false;
	negate_normals = false;
	two_sided = false;
	material.set_ambient(cgv::media::illum::phong_material::color_type(0.8f, 0.8f, 0.8f, 1.0f));
	material.set_diffuse(cgv::media::illum::phong_material::color_type(0.6f, 0.6f, 0.6f, 1.0f));
	material.set_specular(cgv::media::illum::phong_material::color_type(1.0f, 1.0f, 1.0f, 1.0f));
	material.set_shininess(60.0f);
}

/// standard constructor
volume_slicer::volume_slicer() : cgv::base::node("volume_slicer"), threaded_cache_manager(*this)
{
	extrusionLevel = 0.1f;
	draw_selected_cube_bool = false;

	selected_point_material.set_ambient(cgv::media::illum::phong_material::color_type(0.0f, 1.0f, 0.0f, 1.0f));
	selected_point_material.set_diffuse(cgv::media::illum::phong_material::color_type(0.0f, 1.0f, 0.0f, 1.0f));
	//selected_point_material.set_specular(cgv::media::illum::phong_material::color_type(0.0f, 0.0f, 0.0f, 1.0f));
	selected_point_material.set_shininess(20.0f);

	extrusion_material.set_ambient(cgv::media::illum::phong_material::color_type(1.0f, 1.0f, 0.0f, 1.0f));
	extrusion_material.set_diffuse(cgv::media::illum::phong_material::color_type(1.0f, 1.0f, 1.0f, 1.0f));
	//extrusion_material.set_specular(cgv::media::illum::phong_material::color_type(0.0f, 0.0f, 0.0f, 1.0f));
	extrusion_material.set_shininess(20.0f);

	texture_gamma = 1.0f;
	texture_scale = 1.0f;
	view_ptr = 0;
	value_name_map = 0;
	value_map = 0;
	use_reduced_color_depth = false;
	rotate_sensitivity = 1.0;
	translation_sensitivity = 1.0;
	show_box = true;
	show_voxel = true;
	show_block = true;
	show_surface = true;
	slice_center_distance_tex = 0;
	slice_normal_tex = vec3(0, 0, 1);
	dimensions = ivec3(1, 1, 1);
	extent = vec3(1, 1, 1);
	interpolate = false;
	box_mat.set_ambient(cgv::media::illum::phong_material::color_type(0.8f, 0.8f, 0.8f, 1.0f));
	box_mat.set_diffuse(cgv::media::illum::phong_material::color_type(0.6f, 0.6f, 0.6f, 1.0f));
	box_mat.set_specular(cgv::media::illum::phong_material::color_type(0.0f, 0.0f, 0.0f, 1.0f));
	box_mat.set_shininess(20.0f);

	selection_cube_style.diffuse_mapping = DCM_TEXTURE;
	selection_cube_style.illumination = false;

	use_standard_prog = false;
#ifdef _DEBUG
	down_sampling_factor = 32;
#else
	down_sampling_factor = 4;
#endif
	iso_value = 10;
	sm_ptr = 0;

	last_x = last_y = -1;
	current_voxel = current_block = ivec3(0, 0, 0);


	block_dimensions = ivec3(16,16,16);
	overlap = ivec3(1, 1, 1);
	threaded_cache_manager.init_listener();
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
		view_ptr->set_z_near(extent.length()/dimensions.length());
		view_ptr->set_z_far(5.0f*extent.length());
		std::cout << "z_near = " << view_ptr->get_z_near() << ", z_far=" << view_ptr->get_z_far() << std::endl;
	}
}

/// self reflection is used to publish standard member variables to the set and get mechanism
bool volume_slicer::self_reflect(cgv::reflect::reflection_handler& rh)
{
	return 
		rh.reflect_member("texture_scale", texture_scale) &&
		rh.reflect_member("texture_gamma", texture_gamma) &&
		rh.reflect_member("use_reduced_color_depth", use_reduced_color_depth) &&
		rh.reflect_member("value_name_map_file_name", value_name_map_file_name) &&
		rh.reflect_member("file_name", file_name) &&
		rh.reflect_member("iso_file_name", iso_file_name) &&
		rh.reflect_member("slices_path", block_config) &&
		rh.reflect_member("show_box", show_box);
}

bool volume_slicer::read_value_name_map(const std::string& file_name)
{
	std::string content;
	if (!cgv::utils::file::read(file_name, content, true)) {
		std::cerr << "ERROR: could not read value_name_map file <" << file_name << ">" << std::endl;
		return false;
	}
	std::vector<cgv::utils::line> lines;
	cgv::utils::split_to_lines(content, lines);
	value_map = new std::map<std::string, std::vector<int> >();
	value_name_map = new std::map<int, std::string>();
	for (auto l : lines) {
		if (l[0] != '"')
			continue;
		std::vector<cgv::utils::token> tokens;
		cgv::utils::split_to_tokens(l, tokens, "", false, "\"", "\"");
		std::string name = to_string(tokens[0]);
		for (int i = 1; i < (int)tokens.size(); ++i) {
			int idx;
			if (cgv::utils::is_integer(tokens[i].begin, tokens[i].end, idx)) {
				//std::cout << "map " << idx << " -> " << name << std::endl;
				(*value_name_map)[idx] = name;
				(*value_map)[name].push_back(idx);
			}
		}
	}
	return true;
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




double volume_slicer::get_points_distance(vec3 pointA, vec3 pointB) {
	return sqrt(pow((pointA[0] - pointB[0]), 2) + pow((pointA[1] - pointB[1]), 2) + pow((pointA[2] - pointB[2]), 2));
}

volume::point_type volume_slicer::get_selected_point_normal(volume::point_type point) {
	double currentDistance, smallestDistance;
	int currentIndex, smallestIndex;

	for (currentIndex = 0, smallestIndex = 0; currentIndex < positions.size(); currentIndex++) {
		currentDistance = get_points_distance(point, positions[currentIndex]);
		if (currentIndex == 0 || (currentIndex != 0 && currentDistance < smallestDistance)) {
			smallestIndex = currentIndex;
			smallestDistance = currentDistance;
		}
	}

	return normals[smallestIndex];
}

void volume_slicer::clear_selected_points() {
	selected_points.clear();
}

void volume_slicer::savePoint(int x, int y) {
	// unproject pixel position of mouse
	volume::point_type p_pick_world, normal;
	if (!get_picked_point(x, y, p_pick_world)) {
		return;
	}
	if (positions.size() > 0) {
		if (selected_points.size() < 4) {
			selected_points.push_back(p_pick_world);
			normal = get_selected_point_normal(p_pick_world);
			selected_normals.push_back(normal);
			std::cout << "Point saved at:\n" << "x: " << std::to_string(p_pick_world(0)) << " - y: " << std::to_string(p_pick_world(1)) << " - z: " << std::to_string(p_pick_world(2)) << std::endl;
			std::cout << "x: " << std::to_string(normal(0)) << " - y: " << std::to_string(normal(1)) << " - z: " << std::to_string(normal(2)) << std::endl << std::endl;
		}
		else
			std::cout << "Already selected 4 points" << std::endl;
	}
	else {
		std::cout << "NO SURFACE PRESNT!!!!!!!" << std::endl;
	}
}

void volume_slicer::draw_solid_block(cgv::render::context& ctx, const ivec3& voxel, const cgv::media::illum::phong_material& material) {
	box3 B(world_from_texture_coordinates(texture_from_voxel_coordinates(voxel)),
		world_from_texture_coordinates(texture_from_voxel_coordinates(voxel + ivec3(1, 1, 1))));
	draw_box(ctx, B, material);
}

void volume_slicer::bresenham3d(vec3 pointA, vec3 pointB, cgv::render::context& ctx) {
	volume::point_type p_texture, point1, point2, currentPoint;

	p_texture = texture_from_world_coordinates(pointA);
	point1 = voxel_from_texture_coordinates(p_texture);

	p_texture = texture_from_world_coordinates(pointB);
	point2 = voxel_from_texture_coordinates(p_texture);

	double pointAx = floor(point1[0]), pointAy = floor(point1[1]), pointAz = floor(point1[2]);
	double pointBx = floor(point2[0]), pointBy = floor(point2[1]), pointBz = floor(point2[2]);

	double sx = pointBx > pointAx ? 1 : pointBx < pointAx ? -1 : 0;
	double sy = pointBy > pointAy ? 1 : pointBy < pointAy ? -1 : 0;
	double sz = pointBz > pointAz ? 1 : pointBz < pointAz ? -1 : 0;

	double gx = pointAx;
	double gy = pointAy;
	double gz = pointAz;

	//Planes for each axis that will be crossed
	double gxp = pointAx + (pointBx > pointAx ? 1 : 0);
	double gyp = pointAy + (pointBy > pointAy ? 1 : 0);
	double gzp = pointAz + (pointBz > pointAz ? 1 : 0);

	//Used for multiplying up the error margins
	double vx = point2[0] == point1[0] ? 1 : point2[0] - point1[0];
	double vy = point2[1] == point1[1] ? 1 : point2[1] - point1[1];
	double vz = point2[2] == point1[2] ? 1 : point2[2] - point1[2];

	//Error is normalized to vx * vy * vz to only multiply up
	double vxvy = vx * vy;
	double vxvz = vx * vz;
	double vyvz = vy * vz;

	double errx = (gxp - point1[0]) * vyvz;
	double erry = (gyp - point1[1]) * vxvz;
	double errz = (gzp - point1[2]) * vxvy;

	double derrx = sx * vyvz;
	double derry = sy * vxvz;
	double derrz = sz * vxvy;

	do
	{
		box3 B(-extent, extent);
		vec3 point = vec3((float)gx, (float)gy, (float)gz);
		if (B.inside(point)) {
			draw_voxel(ctx, point, vec3(1.0, 1.0, 0.0));
		}

		if ((int)gx == (int)point2[0] && (int)gy == (int)point2[1] && (int)gz == (int)point2[2]) break;

		double xr = abs(errx);
		double yr = abs(erry);
		double zr = abs(errz);

		if (sx != 0 && (sy == 0 || xr < yr) && (sz == 0 || xr < zr)) {
			gx += sx;
			errx += derrx;
		}
		else if (sy != 0 && (sz == 0 || yr < zr)) {
			gy += sy;
			erry += derry;
		}
		else if (sz != 0) {
			gz += sz;
			errz += derrz;
		}
		else {
			std::cout << "the fiddle was a lie" << std::endl;
		}
	} while (true);
}

void volume_slicer::drawSelectedPoints(cgv::render::context& ctx) {
	std::vector<vec3> top_points, bottom_points;
	float scale = 0.5f*V.get_extent().length()*extrusionLevel;
	for (unsigned i = 0; i < selected_points.size(); ++i) {
		top_points.push_back(selected_points[i] + scale*selected_normal);
		bottom_points.push_back(selected_points[i] - scale*selected_normal);
	}

	volume::point_type p_voxel;// , p_voxel_top, p_voxel_bottom;
	volume::point_type p_texture;
	for (unsigned i = 0; i < selected_points.size(); ++i) {
		bresenham3d(bottom_points[i], top_points[i], ctx);
		p_texture = texture_from_world_coordinates(selected_points[i]);
		p_voxel = voxel_from_texture_coordinates(p_texture);
		draw_voxel(ctx, p_voxel, vec3(0, 1, 0));
	}
	for (auto &point : selected_points) {
		box3 B(-0.5f*extent, 0.5f*extent);
		if (B.inside(point)) {
			p_texture = texture_from_world_coordinates(point);
			p_voxel = voxel_from_texture_coordinates(p_texture);

			//p_voxel_top = p_voxel;
			//p_voxel_bottom = p_voxel;

			draw_voxel(ctx, p_voxel, vec3(0, 1, 0));
			//draw_solid_block(ctx, p_voxel, point_material);

			for (unsigned i = 0; i < selected_points.size(); ++i) {
				bresenham3d(bottom_points[i], top_points[i], ctx);
			}

			/*
			float counter = 0;
			while (counter < extrusionLevel) {
				p_voxel_top(1) += 1;
				p_voxel_bottom(1) -= 1;

				//draw_solid_block(ctx, p_voxel, extrusion_material);
				draw_voxel(ctx, p_voxel_top, vec3(1, 1, 0));
				draw_voxel(ctx, p_voxel_bottom, vec3(1, 1, 0));
				
				counter += 0.01;
			}
			*/
		}
	}
}

void volume_slicer::draw_selected_cube(cgv::render::context& ctx) {
	std::vector<vec3> top_points, bottom_points;
	float scale = 0.5f*V.get_extent().length()*extrusionLevel;
	for (unsigned i = 0; i < selected_points.size(); ++i) {
		top_points.push_back(selected_points[i] + scale*selected_normal);
		bottom_points.push_back(selected_points[i] - scale*selected_normal);
	}

	std::vector<vec3> positions, texcoords;
	for (unsigned i = 0; i < selected_points.size(); ++i) {
		positions.push_back(bottom_points[i]);
		positions.push_back(top_points[i]);
		texcoords.push_back(texture_from_world_coordinates(bottom_points[i]));
		texcoords.push_back(texture_from_world_coordinates(top_points[i]));
	}
	if (selected_points.size() == 4) {
		positions.push_back(bottom_points[0]);
		positions.push_back(top_points[0]);
		texcoords.push_back(texture_from_world_coordinates(bottom_points[0]));
		texcoords.push_back(texture_from_world_coordinates(top_points[0]));
	}

	if (selected_points.size() > 0) {
		enable_surface_rendering(ctx, selection_cube_style);
		glVertexPointer(3, GL_FLOAT, 0, &positions[0]);
		glEnableClientState(GL_VERTEX_ARRAY);
		glTexCoordPointer(3, GL_FLOAT, 0, &texcoords[0]);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, (GLsizei) positions.size());
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
		disable_surface_rendering(ctx, selection_cube_style);
	}
}

#include <cgv/gui/animate.h>
void volume_slicer::animate_view_to_selection(volume::point_type pointA, volume::point_type pointB, vec3 normal) {
	vec3 xDirection = pointB - pointA;
	vec3 upDirection = normal;
	if (iso_surface_style.negate_normals) {
		upDirection = -upDirection;
	}
	vec3 viewDirection = cross(xDirection, upDirection);
	viewDirection.normalize();
	double angle;
	cgv::render::view::pnt_type axis;
	/*view_ptr->compute_axis_and_angle(viewDirection, upDirection, axis, angle);
	cgv::gui::animate_with_axis_rotation(view_ptr->ref_view_dir(), axis, angle, 1);
	cgv::gui::animate_with_axis_rotation(view_ptr->ref_view_up_dir(), axis, angle, 1);
	cgv::gui::animate_with_geometric_blend(view_ptr->ref_y_extent_at_focus(), (double)1.4*V.get_extent().length()*extrusionLevel, 1);
	cgv::gui::animate_with_linear_blend(view_ptr->ref_focus(), cgv::render::view::pnt_type(0.5f*(pointA + pointB)), 1)->set_base_ptr(this);*/
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
		int value = (int)evaluate(p_pick_world.to_vec());
		ss << p_voxel << "\n" << value;
		if (value_name_map) {
			ss << "[" << (*value_name_map)[value] << "]";
			last_name = (*value_name_map)[value];
		}
		else
			last_name.clear();

	}
	voxel_value_info = ss.str();
	last_x = x;
	last_y = y;
	post_redraw();
}


/// overload and implement this method to handle events
bool volume_slicer::handle(cgv::gui::event& e)
{
	// check for key events
	if (e.get_kind() == cgv::gui::EID_KEY) {
		cgv::gui::key_event& ke = static_cast<cgv::gui::key_event&>(e);
		if (ke.get_action() == cgv::gui::KA_RELEASE)
			return false;
		switch (ke.get_key()) {
		// view all when space key pressed
		case cgv::gui::KEY_Space:
			auto_adjust_view();
			return true;
		case '1':
			if (selected_points.size() >= 2)
				animate_view_to_selection(selected_points[0], selected_points[1], selected_normal);
			return true;
		case '2':
			if (selected_points.size() >= 3)
				animate_view_to_selection(selected_points[1], selected_points[2], selected_normal);
			return true;
		case '3':
			if (selected_points.size() >= 4)
				animate_view_to_selection(selected_points[2], selected_points[3], selected_normal);
			return true;
		case '4':
			if (selected_points.size() >= 4)
				animate_view_to_selection(selected_points[3], selected_points[0], selected_normal);
			return true;
		case 'E' :
			if (!last_name.empty())
				extract_surface(4);
			return true;
		case 'M':
			extract_surface(0);
			return true;
		case 'S':
			show_surface = !show_surface;
			on_set(&show_surface);
			return true;
		case 'N':
			iso_surface_style.negate_normals = !iso_surface_style.negate_normals;
			on_set(&iso_surface_style.negate_normals);
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
		case cgv::gui::MA_PRESS:
			// support drag and drop for left and right button
			if (me.get_button() == cgv::gui::MB_LEFT_BUTTON)
				return true;
			else
				if (me.get_button() == cgv::gui::MB_RIGHT_BUTTON) {
					savePoint(me.get_x(), me.get_y());
					return true;
				}
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
				double drag_x =  (double)me.get_dx() / get_context()->get_width();
				double drag_y = -(double)me.get_dy() / get_context()->get_height();
				// ctrl-left button rotates plane normal
				if (me.get_button_state() == cgv::gui::MB_LEFT_BUTTON) {
					// scale x/y-pixel motion to range [-1,1] and 
					// compute rotation axis orthogonal to drag vector in world coordinates
					vec3 axis_world = drag_x*view_y_world - drag_y*view_x_world;
					// normalize axis and use resulting length to scale rotation angle  
					double rotation_angle = 5.0 * axis_world.normalize() * rotate_sensitivity;
					// update slice_normal_tex by rotation around axis
					//slice_normal_tex = rotate(slice_normal_tex, axis_world, (float)rotation_angle);
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
						view_ptr->get_y_extent_at_focus()/extent.length()/translation_sensitivity);
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
bool volume_slicer::open_volume(const std::string& _file_name, bool is_iso)
{
	volume V_local;
	volume& V_current = is_iso ? V : V_local;
	// try to read volume from slice based or regular volume file 
	bool success = false;
	std::string ext = cgv::utils::to_upper(cgv::utils::file::get_extension(_file_name));
	if (ext == "SVX") {
		success = read_from_sliced_volume(_file_name, V_current);
		if (success)
			extent = vec3(float(V_current.get_dimensions()(0))/ V_current.get_dimensions()(2), float(V_current.get_dimensions()(1)) / V_current.get_dimensions()(2), 1);
	}
	else {
		volume_info info;
		success = read_volume(_file_name, V_current, &info);
		if (success)
			extent = info.extent;
	}

	// abort in case of failure 
	if (!success)
		return false;

	// copy dimensions
	if (!is_iso) {
		dimensions = V_current.get_dimensions();

		// transfer volume to 3d texture 
		if (get_context() && get_context()->make_current()) {
			if (data_texture.is_created())
				data_texture.destruct(*get_context());

			static_cast<cgv::data::component_format&>(data_texture) = V_current.get_format();
			if (use_reduced_color_depth || !data_texture.create(*get_context(), V_current.get_data_view())) {
				data_texture.destruct(*get_context());
				data_texture.set_component_format("uint8[R:3,G:3,B:2]");
				if (!data_texture.create(*get_context(), V_current.get_data_view())) {
					std::cerr << "error tex create" << std::endl;
				}
			}
		}
	}
	if (is_iso && find_control(iso_value)) {
		if (V_current.get_component_type() == cgv::type::info::TI_UINT16)
			find_control(iso_value)->set("max", 65535);
		else
			find_control(iso_value)->set("max", 255);
	}
	// adjust view
	auto_adjust_view();
	return success;
}


/// Coordinates conversion

/// convert world to texture coordinates
volume_slicer::vec3 volume_slicer::texture_from_world_coordinates(const vec3& p_world) const
{
	vec3 p_texture = (p_world + 0.5f*extent)/extent;
	return p_texture;
}
/// convert texture to voxel coordinates
volume_slicer::vec3 volume_slicer::voxel_from_texture_coordinates(const vec3& p_texture) const
{
	vec3 p_voxel = p_texture * dimensions;
	return p_voxel;
}

volume_slicer::vec3 volume_slicer::voxel_from_texture_coordinates_blocks(const vec3& p_texture) const
{
	vec3 p_voxel = p_texture * slices_dimensions;
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

volume_slicer::vec3 volume_slicer::texture_from_voxel_coordinates_blocks(const vec3& p_voxel) const
{
	vec3 p_texture = p_voxel / slices_dimensions;
	return p_texture;
}

/// convert voxel to block coordinates
volume_slicer::vec3 volume_slicer::block_from_voxel_coordinates(const vec3& p_voxel) const
{
	vec3 offset = overlap;
	offset *= 0.5f;
	vec3 p_block = (p_voxel-offset) / block_dimensions;
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
	if (!surface_prog.is_created())
		surface_prog.build_program(ctx, "surface.glpr", true);

	/*
	if (data_texture.is_created()) {
		cgv::data::data_format df(170, 170, 170, cgv::type::info::TI_UINT8, cgv::data::CF_R);
		cgv::data::data_view dv(&df);
		cgv::type::uint8_type* ptr = dv.get_ptr<cgv::type::uint8_type>();
		std::fill(ptr, ptr + 170 * 170 * 170, 50);
		data_texture.replace(ctx, 100, 80, 60, dv, 0);
	}*/
}

template <typename I, typename T>
struct is_element
{
	std::vector<I> reference_values;
	is_element(const std::vector<I>& _values) : reference_values(_values) {}
	bool operator () (const T& _value) const { 
		for (auto v : reference_values)
			if (v == (I)_value)
				return true;
		return false;
	}
};

void volume_slicer::extract_surface(int algorithm)
{
	clear_surface();
	vec3 p = extent;
	p /= dimensions;
	box3 box(-0.5f*(extent - p), 0.5f*(extent-p));
	switch (algorithm) {
	case 0: {
			face_vertex_count = 3;
			interpolated_evaluation = true;
			cgv::media::mesh::marching_cubes<float, float> mc(*this, this);
			sm_ptr = &mc;
			mc.extract(
				iso_value, box,
				dimensions(0) / down_sampling_factor,
				dimensions(1) / down_sampling_factor,
				dimensions(2) / down_sampling_factor,
				dimensions(0)*dimensions(1)*dimensions(2) > int(10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor));
		}
		break;
	case 1: {
		face_vertex_count = 4;
		interpolated_evaluation = false;
		cgv::media::mesh::greater_equal<float> pred(iso_value);
		cgv::media::mesh::cuberille<float, float, cgv::media::mesh::greater_equal<float> > c(*this, this, pred);
		sm_ptr = &c;
		c.extract(box,
			dimensions(0) / down_sampling_factor,
			dimensions(1) / down_sampling_factor,
			dimensions(2) / down_sampling_factor,
			dimensions(0)*dimensions(1)*dimensions(2) > int(10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor));
	}
			break;
	case 2: {
		face_vertex_count = 4;
		interpolated_evaluation = false;
		cgv::media::mesh::equal<float> pred(iso_value);
		cgv::media::mesh::cuberille<float, float, cgv::media::mesh::equal<float> > c(*this, this, pred);
		sm_ptr = &c;
		c.extract(box,
			dimensions(0) / down_sampling_factor,
			dimensions(1) / down_sampling_factor,
			dimensions(2) / down_sampling_factor,
			dimensions(0)*dimensions(1)*dimensions(2) > int(10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor));
	}
			break;
	case 3: {
			face_vertex_count = 4;
			interpolated_evaluation = true;
			cgv::media::mesh::dual_contouring<float, float> dc(*this, this);
			sm_ptr = &dc;
			dc.extract(
				iso_value, box,
				dimensions(0) / down_sampling_factor,
				dimensions(1) / down_sampling_factor,
				dimensions(2) / down_sampling_factor,
				dimensions(0)*dimensions(1)*dimensions(2) > int(10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor));
		}
		break;
	case 4: {
		face_vertex_count = 4;
		interpolated_evaluation = false;
		std::cout << "extracting " << last_name << " [";
		for (auto v : (*value_map)[last_name])
			std::cout << " " << v;
		std::cout << " ]" << std::endl;
		is_element<int, float> pred((*value_map)[last_name]);
		cgv::media::mesh::cuberille<float, float, is_element<int, float> > c(*this, this, pred);
		sm_ptr = &c;
		c.extract(box,
			dimensions(0) / down_sampling_factor,
			dimensions(1) / down_sampling_factor,
			dimensions(2) / down_sampling_factor,
			dimensions(0)*dimensions(1)*dimensions(2) > int(10000 * down_sampling_factor*down_sampling_factor*down_sampling_factor));
	}
			break;
	}
	interpolated_evaluation = true;
	sm_ptr = 0;
	std::cout << "iso_surface has " << positions.size() << " vertices and " << corner_indices.size() / face_vertex_count << (face_vertex_count == 3 ? " triangles" : " quads") << std::endl;
	// triangulate quads
	if (face_vertex_count == 4) {
		unsigned n = (unsigned) (corner_indices.size() / 4);
		corner_indices.resize(6 * n);
		for (int i = n - 1; i >= 0; --i) {
			unsigned vi[4] = { corner_indices[4 * i], corner_indices[4 * i + 1], corner_indices[4 * i + 2], corner_indices[4 * i + 3] };
			corner_indices[6 * i]     = vi[0];
			corner_indices[6 * i + 1] = vi[1];
			corner_indices[6 * i + 2] = vi[2];
			corner_indices[6 * i + 3] = vi[0];
			corner_indices[6 * i + 4] = vi[2];
			corner_indices[6 * i + 5] = vi[3];
		}
		face_vertex_count = 3;
		std::cout << "iso_surface has " << positions.size() << " vertices and " << corner_indices.size() / face_vertex_count << (face_vertex_count == 3 ? " triangles" : " quads") << std::endl;
	}
	if (!show_surface) {
		show_surface = true;
		update_member(&show_surface);
	}
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
		return float(*V.get_voxel_ptr<cgv::type::uint8_type>(voxel(0), voxel(1), voxel(2)));
	case cgv::type::info::TI_UINT16:
		return float(*V.get_voxel_ptr<cgv::type::uint16_type>(voxel(0), voxel(1), voxel(2)));
	}
	return 0;
}

/// interface for evaluation of the multivariate function
float volume_slicer::evaluate(const pnt_type& p) const
{
	// convert to voxel coordinates
	vec3 p_voxel = voxel_from_texture_coordinates(texture_from_world_coordinates(&p(0))) - vec3(0.5f, 0.5f, 0.5f);
	// quantize
	ivec3 voxel = p_voxel;
	// check if voxel is valid
	if (voxel(0) < 0 || voxel(1) < 0 || voxel(2) < 0)
		return 0.0f;
	if (interpolated_evaluation) {
		if (voxel(0) >= dimensions(0)-1 || voxel(1) >= dimensions(1)-1 || voxel(2) >= dimensions(2)-1)
			return 0.0f;
	}
	else {
		if (voxel(0) >= dimensions(0) || voxel(1) >= dimensions(1) || voxel(2) >= dimensions(2))
			return 0.0f;
		return get_voxel_value(voxel);
	}
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
	float v00 = (1-frac(2))*v000 + frac(2)*v001;
	float v10 =	(1-frac(2))*v100 + frac(2)*v101;
	float v01 =	(1-frac(2))*v010 + frac(2)*v011;
	float v11 =	(1-frac(2))*v110 + frac(2)*v111;
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
	if (face_nml.length() > 1.0e-6f*extent.length()) {
		face_nml.normalize();
		// distribute face normal to vertex normals
		for (i = 0; i < vertex_indices.size(); ++i)
			normals[vertex_indices[i]] += face_nml;
	}
}

/// drop the currently first vertex that has the given global vertex index
void volume_slicer::before_drop_vertex(unsigned int vertex_index)
{

	if (normals[vertex_index].length() >= 1.0e-6f*extent.length()) {
		normals[vertex_index].normalize();
	}
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

static bool abs_compare(float a, float b){ return (std::abs(a) < std::abs(b)); }

void volume_slicer::update_intersected_blocks(cgv::render::context& ctx) {

	ivec3 nr_blocks(unsigned(ceil(float(slices_dimensions(0)) / block_dimensions(0))), 
		unsigned(ceil(float(slices_dimensions(1)) / block_dimensions(1))),
		unsigned(ceil(float(slices_dimensions(2)) / block_dimensions(2))));

	vec3 world_slice_normal = world_from_texture_normals(slice_normal_tex);
	int d_max = (int) std::distance(std::begin(world_slice_normal), std::max_element(std::begin(world_slice_normal), std::end(world_slice_normal), abs_compare));
	
	// set i0 to first other index than d_max and i1 to second other index
	int d0 = (d_max + 1) % 3;
	int d1 = (d_max + 2) % 3;

	previous_intersected_blocks.clear();
	previous_intersected_blocks.insert(intersected_blocks.begin(), intersected_blocks.end());
	intersected_blocks.clear();

	for (float i0 = 0; i0 < nr_blocks(d0); i0++) {
		for (float i1 = 0; i1 < nr_blocks(d1); i1++) {
			
			vec3 pi_tex[4]{
				texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(make_vec_d_max(d_max, 0, i0, i1))),
				texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(make_vec_d_max(d_max, 0, i0 + 1, i1))),
				texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(make_vec_d_max(d_max, 0, i0, i1 + 1))),
				texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(make_vec_d_max(d_max, 0, i0 + 1, i1 + 1)))
			};

			vec3 pi_block[4];
			
			float block_coord_min, block_coord_max;
			
			int k;

			for (k = 0; k < 4; k++) {

				pi_tex[k](d_max) = 
					1.0f / slice_normal_tex(d_max) * 
					(slice_center_distance_tex - 
					(pi_tex[k](d0) - 0.5f) * slice_normal_tex(d0) - 
					(pi_tex[k](d1) - 0.5f) * slice_normal_tex(d1)) + 0.5f;
				
				pi_block[k] = block_from_voxel_coordinates(voxel_from_texture_coordinates_blocks(pi_tex[k]));
				
				block_coord_max = (k == 0 ? pi_block[k][d_max] : std::max(block_coord_max, pi_block[k][d_max]));
				block_coord_min = (k == 0 ? pi_block[k][d_max] : std::min(block_coord_min, pi_block[k][d_max]));
			}

			int i_min = std::max(0,int(floor(block_coord_min)));
			int i_max = std::max(0,int(ceil(block_coord_max)));

			if (i_min > i_max) {
				std::cout << "this shouldn't happen"  << std::endl;
				int temp = i_min;
				i_min = i_max;
				i_max = i_min;
			}

			for (k = i_min; k < i_max; k++) {
				if (nr_blocks(d_max) > k && k >= 0) {
					ivec3 block = (ivec3) make_vec_d_max(d_max, (float) k, i0, i1);

					if (previous_intersected_blocks.find(block) == previous_intersected_blocks.end()) 
						new_batch = true;

					intersected_blocks.push_back(block);
				}
			}
		}
	}
}


volume_slicer::vec3 volume_slicer::make_vec_d_max(int d_max, float var_dim, float dim0, float dim1) {
	float x, y, z;
	switch (d_max) {
		case 0: {
			x = var_dim;
			y = dim0;
			z = dim1;
			break;
		}
		case 1: {
			x = dim1;
			y = var_dim;
			z = dim0;
			break;
		}
		case 2: {
			x = dim0;
			y = dim1;
			z = var_dim;
			break;
		}
	}
	return vec3(x,y,z);
}

/// draw a block with a wire frame box
void volume_slicer::draw_block(cgv::render::context& ctx, const ivec3& block, const vec3& color, const vec3& overlap_color)
{
	box3 B(world_from_texture_coordinates(texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(block))),
		world_from_texture_coordinates(texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(block + ivec3(1, 1, 1)))));
	draw_wire_box(ctx, B, color);
	vec3 offset = overlap;
	offset *= 0.5f;
	box3 B_overlap(
		world_from_texture_coordinates(texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(block)-offset)),
		world_from_texture_coordinates(texture_from_voxel_coordinates_blocks(voxel_from_block_coordinates(block + ivec3(1, 1, 1)) + offset)));
	draw_wire_box(ctx, B_overlap, overlap_color);
}

void volume_slicer::draw_blocks_in_plane(cgv::render::context& ctx, const vec3& color, const vec3& overlap_color) {
	
	for (auto it : intersected_blocks) {
		draw_block(ctx, it, color, overlap_color);
	}
}

void volume_slicer::enable_surface_rendering(cgv::render::context& ctx, const surface_render_style& style) 
{
	if (style.show_wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(style.wireframe_width);
	}
	switch (style.culling_mode) {
	case cgv::render::CM_OFF:
		glDisable(GL_CULL_FACE);
		break;
	case cgv::render::CM_FRONTFACE:
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		break;
	case cgv::render::CM_BACKFACE:
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		break;
	}
	if (use_standard_prog) {
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, style.two_sided ? 1 : 0);
		if (style.diffuse_mapping == DCM_COLOR)
			glEnable(GL_COLOR_MATERIAL);
		else
			glDisable(GL_COLOR_MATERIAL);
		ctx.enable_material(style.material);
	}
	else {
		ctx.enable_material(style.material);
		ctx.disable_material(style.material);
		data_texture.enable(ctx, 0);
		surface_prog.enable(ctx);
		cgv::render::gl::set_lighting_parameters(ctx, surface_prog);
		surface_prog.set_uniform(ctx, "texture_gamma", texture_gamma);
		surface_prog.set_uniform(ctx, "texture_scale", texture_scale);
		surface_prog.set_uniform(ctx, "diffuse_mapping", (int)style.diffuse_mapping);
		surface_prog.set_uniform(ctx, "two_sided", style.two_sided ? 1 : 0);
		surface_prog.set_uniform(ctx, "negate_normals", style.negate_normals);
		surface_prog.set_uniform(ctx, "illumination", style.illumination);
		surface_prog.set_uniform(ctx, "face_normals", style.face_normals);
		box3 B(-0.5f*extent, 0.5f*extent);
		surface_prog.set_uniform(ctx, "minp", B.get_min_pnt());
		surface_prog.set_uniform(ctx, "maxp", B.get_max_pnt());
		surface_prog.set_uniform(ctx, "diffuse_map", 0);
	}
}

/// disables surface rendering with the given rendering style
void volume_slicer::disable_surface_rendering(cgv::render::context& ctx, const surface_render_style& style) 
{
	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);
	if (use_standard_prog)
		ctx.disable_material(style.material);
	else {
		data_texture.disable(ctx);
		surface_prog.disable(ctx);
	}
	if (style.show_wireframe) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

/// draw a block with a wire frame box
void volume_slicer::draw_surface(cgv::render::context& ctx)
{
	if (corner_indices.empty())
		return;
	enable_surface_rendering(ctx, iso_surface_style);
		glVertexPointer(3, GL_FLOAT, 0, &positions[0]);
		glEnableClientState(GL_VERTEX_ARRAY);
		glNormalPointer(GL_FLOAT, 0, &normals[0]);
		glEnableClientState(GL_NORMAL_ARRAY);
			switch (face_vertex_count) {
			case 3: glDrawElements(GL_TRIANGLES, (GLsizei)corner_indices.size(), GL_UNSIGNED_INT, &corner_indices[0]); break;
			case 4: glDrawElements(GL_QUADS, (GLsizei) corner_indices.size(), GL_UNSIGNED_INT, &corner_indices[0]); break;
			default: std::cerr << "face vertex count " << face_vertex_count << "not supported!" << std::endl;
			}
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);
	disable_surface_rendering(ctx, iso_surface_style);
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
		draw_voxel(ctx, current_voxel, vec3(1,0,0));
	
	if (block_config != "") {
		update_intersected_blocks(ctx);
		// upload more texture data for block based volumes
		if (new_batch) {
			threaded_cache_manager.request_blocks(intersected_blocks);
			new_batch = false;
		}
	}
	
	if (show_block) {
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
	vec3 c_tex = slice_center_distance_tex*slice_normal_tex+vec3(0.5f, 0.5f, 0.5f);
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

	// Show selected voxels
	drawSelectedPoints(ctx);
	if (draw_selected_cube_bool)
		draw_selected_cube(ctx);

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
		tex_slice_prog.set_uniform(ctx, "texture_gamma", texture_gamma);
		tex_slice_prog.set_uniform(ctx, "texture_scale", texture_scale);
	}
	// construct texture coordinates of quad that spans slice
	static std::vector<vec3> tex_coords(size_t(4));
	tex_coords[0] = c_tex + x_tex;
	tex_coords[1] = c_tex + y_tex;
	tex_coords[2] = c_tex - x_tex;
	tex_coords[3] = c_tex - y_tex;
	// map texture coordinates of quad corners to positions in world coordinates
	static std::vector<vec3> positions(size_t(4));
	for (unsigned i=0; i<4; ++i)
		positions[i] = world_from_texture_coordinates(tex_coords[i]); // extent*(tex_coords[i] - vec3(0.5f, 0.5f, 0.5f));
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
	threaded_cache_manager.kill_listener();
}

///
void volume_slicer::create_surface_gui(surface_render_style& style)
{
	add_member_control(this, "diffuse_mapping", (cgv::type::DummyEnum&)style.diffuse_mapping, "dropdown", "enums='material,color,texture,normal,texcoord,position'");
	add_member_control(this, "show_wireframe", style.show_wireframe, "check");
	add_member_control(this, "wireframe_width", style.wireframe_width, "value_slider", "min=1;max=10;ticks=true");
	add_member_control(this, "culling_mode", (cgv::type::DummyEnum&)style.culling_mode, "dropdown", "enums='off,back,front'");
	add_member_control(this, "illumination", style.illumination, "check");
	add_member_control(this, "face_normals", style.face_normals, "check");
	add_member_control(this, "negate_normals", style.negate_normals, "check");
	add_member_control(this, "two_sided", style.two_sided, "check");
	add_gui("material", style.material);
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
		add_gui("iso_file_name", iso_file_name, "file_name", "title='open iso volume';filter='Volume Files(vox,qim,tif,avi) :*.vox;*.qim;*.tif;*.avi|All Files:*.*'");
		add_gui("slices_path", block_config, "file_name", "title='choose slices path';filter='Block Config Files:*.bsdc'");
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
	/*if (begin_tree_node("Blocks", block_dimensions, false, "level=3")) {
		align("\a");
		// w=50 specifies the width of each component view in 50 pixels
		add_gui("block_dimensions", block_dimensions, "vector", "main_label='first';align=' ';gui_type='value_input';options='w=50;min=1;max=1024;step=1'"); align("\n");
		// w=50 specifies the width of each component view in 50 pixels
		add_gui("overlap", overlap, "vector", "main_label='first';align=' ';gui_type='value_input';options='w=50;min=0;max=1;step=1'"); align("\n");
		align("\b");
		end_tree_node(block_dimensions);
	}*/
	// start tree node to interact with volume
	if (begin_tree_node("Iso Surface", down_sampling_factor, false, "level=3")) {
		align("\a");
		connect_copy(add_button("marching cubes")->click, cgv::signal::rebind(this, &volume_slicer::extract_surface, cgv::signal::_c<int>(0)));
		connect_copy(add_button("cuberille greaterequal")->click, cgv::signal::rebind(this, &volume_slicer::extract_surface, cgv::signal::_c<int>(1)));
		connect_copy(add_button("cuberille equal")->click, cgv::signal::rebind(this, &volume_slicer::extract_surface, cgv::signal::_c<int>(2)));
		connect_copy(add_button("dual contouring")->click, cgv::signal::rebind(this, &volume_slicer::extract_surface, cgv::signal::_c<int>(3)));
		add_member_control(this, "iso_value", iso_value, "value_slider", "min=0;max=255;step=1;ticks=true");
		if (!V.empty() && V.get_component_type() == cgv::type::info::TI_UINT16)
			find_control(iso_value)->set("max", 65535);

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
		add_member_control(this, "use_standard_prog", use_standard_prog, "check");
		add_member_control(this, "texture_gamma", texture_scale, "value_slider", "min=1;max=256;ticks=true");
		add_member_control(this, "texture_gamma", texture_gamma, "value_slider", "min=0.01;max=100;log=true;ticks=true");
		align("\a");
		add_decorator("iso surface style", "heading", "level=4");
		create_surface_gui(iso_surface_style);
		align("\b");
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
	if (begin_tree_node("Cutting", extrusionLevel, false, "level=3")) {
		align("\a");
		add_member_control(this, "extrusion", extrusionLevel, "value_slider", "min=0;max=1;ticks=true");
		add_member_control(this, "draw selected cube", draw_selected_cube_bool, "toggle");
		connect_copy(add_button("clear points")->click, cgv::signal::rebind(this, &volume_slicer::clear_selected_points));
		add_gui("selection normal", selected_normal, "direction");
		align("\b");
		end_tree_node(extrusionLevel);
	}
}

// 
void volume_slicer::on_set(void* member_ptr)
{
	if (member_ptr == &value_name_map_file_name) {
		read_value_name_map(value_name_map_file_name);
		return;
	}
	// in case that file_name changed, read new volume
	if (member_ptr == &file_name) {
		open_volume(file_name, false);
		block_config = "";
		threaded_cache_manager.set_block_folder(block_config);
		post_redraw();
		post_recreate_gui();
	}
	// in case that iso_file_name changed, read new iso volume
	if (member_ptr == &iso_file_name) {
		open_volume(iso_file_name, true);
		post_redraw();
		post_recreate_gui();
	}
	
	// in case that slices path changed, modify the folder in the cache
	if (member_ptr == &block_config) {
		set_block_structure(block_config);
		return;
	}

	// update gui of member
	update_member(member_ptr);
	// for all value changes schedule a redraw of the 3D window
	post_redraw();
}

void volume_slicer::set_block_structure(std::string) {
	threaded_cache_manager.set_block_folder("");
	std::cout << "\nsetting block dimensions..." << std::endl;
	std::string line_parse;
	std::ifstream infile;
	infile.open(block_config);

	//line 1: block structure dimensions
	//line 2: block dimensions
	//line 3: overlap
	//line 4: slices path

	bool line_checks[4] = { false, false, false, false };
	std::string slices_path;

	if (infile.is_open()) {
		int line_index = 0;
		while (std::getline(infile, line_parse)) {

			size_t pos = 0;
			std::string delimiter = "x";

			ivec3 variables[3];
			std::string variable_names[4] = { "block_structure_dimensions=", "block_dimensions=", "overlap=", "slices_folder=" };
			
			if (line_index < 3) {
				if ((pos = line_parse.find(variable_names[line_index])) != std::string::npos) {
					line_parse.erase(0, pos + variable_names[line_index].length());
					int i = 0;
					while ((pos = line_parse.find(delimiter)) != std::string::npos) {
						variables[line_index](i) = std::stoi(line_parse.substr(0, pos), NULL);
						line_parse.erase(0, pos + delimiter.length());
						i++;
					}
					variables[line_index](i) = std::stoi(line_parse.substr(0, pos), NULL);
					line_checks[line_index] = true;
				}

				switch (line_index) {
				case 0: { slices_dimensions = variables[line_index]; break; }
				case 1: { block_dimensions = variables[line_index]; break; }
				case 2: { overlap = variables[line_index]; break; }
				default: { break; }
				}
			}

			if (line_index == 3) {
				if ((pos = line_parse.find(variable_names[line_index])) != std::string::npos) {
					line_parse.erase(0, pos + variable_names[line_index].length());

					slices_path = line_parse.substr(0, std::string::npos);
					line_checks[line_index] = true;
				}
			}

			line_index += 1;
		}

	} else {
		block_config = "";
		std::cout << "couldn't read file \n" << std::endl;
	}

	infile.close();
	for (bool check : line_checks) {
		if (!check) {
			std::cout << "error reading file \n" << std::endl;
			return;
		}
	}

	threaded_cache_manager.set_block_folder(slices_path);
	std::cout << "success reading the block configuration file\n" << std::endl;
	threaded_cache_manager.request_blocks(intersected_blocks);

}

#include <cgv/base/register.h>

extern cgv::base::object_registration<volume_slicer> vol_slicer_reg("volume_slicer");
