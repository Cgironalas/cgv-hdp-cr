#include "volume_slicer.h"
#include <vol_data/volume_io.h>
#include <vol_data/sliced_volume_io.h>
#include <cgv/utils/scan.h>
#include <cgv/base/find_action.h>
#include <cgv/utils/file.h>
#include <cgv_gl/gl/gl.h>
#include <cgv/gui/key_event.h>
#include <cgv/gui/mouse_event.h>

/// standard constructor
volume_slicer::volume_slicer() : cgv::base::node("volume_slicer")
{
	view_ptr = 0;
	rotate_sensitivity = 1.0;
	translation_sensitivity = 1.0;
	show_box = true;
	slice_distance_tex = 0;
	slice_normal_tex = vec3(0, 0, 1);
	dimensions = ivec3(1, 1, 1);
	extent = vec3(1, 1, 1);
	interpolate = false;
	box_mat.set_ambient(cgv::media::illum::phong_material::color_type(0.8f, 0.8f, 0.8f, 1.0f));
	box_mat.set_diffuse(cgv::media::illum::phong_material::color_type(0.6f, 0.6f, 0.6f, 1.0f));
	box_mat.set_specular(cgv::media::illum::phong_material::color_type(0.0f, 0.0f, 0.0f, 1.0f));
	box_mat.set_shininess(20.0f);
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
		// we are interested only in mouse events with the Ctrl-modifier and no other modifier pressed
		if (me.get_modifiers() != cgv::gui::EM_CTRL)
			return false;

		switch (me.get_action()) {
		case cgv::gui::MA_PRESS:
			// support drag and drop for left and right button
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
					slice_distance_tex += (float)(10.0f*dot(slice_normal_tex, drag_world)*
						view_ptr->get_y_extent_at_focus() / extent.length() / translation_sensitivity);
					// ensure that slice distance is in valid range
					if (slice_distance_tex > 0.5f)
						slice_distance_tex = 0.5f;
					else if (slice_distance_tex < -0.5f)
						slice_distance_tex = -0.5f;
					// ensure that user interface and rendering are up to date
					update_member(&slice_distance_tex);
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
	os << "volume_slicer: slice_normal=[" << slice_normal_tex << "], distance=" << slice_distance_tex << std::endl;
}
// read volume from sliced or single file representations
bool volume_slicer::open_block_volume(const std::string& directory_name)
{
	return false;
}

/// read regular volume file
bool volume_slicer::open_volume(const std::string& _file_name)
{
	// volume data structure is only used to load data into 3d texture
	volume V;

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

///
std::string volume_slicer::get_type_name() const
{
	return "volume_slicer";
}
// extend method of volume_drawable to find and set view
bool volume_slicer::init(cgv::render::context& ctx)
{
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

///
void volume_slicer::draw(cgv::render::context& ctx)
{
	// compute axis aligned 3d box of volume
	box3 B(-0.5f*extent, 0.5f*extent);

	// render box
	if (show_box) {
		// first draw wireframe
		glColor4f(0, 1, 1, 1);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		ctx.tesselate_box(B, false);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// draw backfaces next
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
		ctx.enable_material(box_mat);
		ctx.tesselate_box(B, true);
		ctx.disable_material(box_mat);
		glCullFace(GL_BACK);
		glDisable(GL_CULL_FACE);
	}

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
	vec3 c_tex = slice_distance_tex*slice_normal_tex + vec3(0.5f, 0.5f, 0.5f);
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
	// use member slice_distance_tex to hash the node status and open node initially (true in 3rd argument)
	if (begin_tree_node("Slicing", slice_distance_tex, true, "level=3")) {
		// tabify gui
		align("\a");
		// add value_slider control for slice_distance_tex member specifying range [-1,1]
		add_member_control(this, "slice_distance_tex", slice_distance_tex, "value_slider", "min=-1;max=1;step=0.0001;ticks=true");
		// add a compount gui of type "direction" for the vector valued member slice_normal_tex
		// the direction gui ensures that the slice_normal_tex is of length one
		add_gui("slice_normal_tex", slice_normal_tex, "direction");
		// untabify gui
		align("\b");
		// end gui node definition specifying the same member used for hashing the node status 
		end_tree_node(slice_distance_tex);
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
	if (begin_tree_node("Rendering", show_box, false, "level=3")) {
		align("\a");
		add_member_control(this, "show_box", show_box, "check");
		add_gui("box_mat", box_mat);
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
