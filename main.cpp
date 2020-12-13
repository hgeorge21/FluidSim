#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>

#include <iostream>
#include <chrono>
#include <random>
#include <grid.h>
#include <particles.h>
#include <create_rectangle.h>

#include <advect_velocity.h>
#include <transfer_to_grid.h>
#include <add_gravity.h>

int main(int argc, char** argv) {
	igl::opengl::glfw::Viewer viewer;
	const int xid = viewer.selected_data_index;
	viewer.data().point_size = 2;
	viewer.append_mesh();


	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// igl::read_triangle_mesh((argc > 1 ? argv[1] : "../../../data/bunny.off"), V, F);
	
	double dt = 0.01;
	double height = 0.3;
	Eigen::Vector3d gravity = Eigen::Vector3d(0, -9.8, 0);
	Eigen::Vector3d h = Eigen::Vector3d(0.01, 0.01, 0.01);
	Eigen::Vector3d left_lower_corner = Eigen::Vector3d(-0.5, -0.5, -0.5);
	Eigen::Vector3d right_upper_corner = Eigen::Vector3d(0.5, 0.5, 0.5);

	// for visualizing boundary
	create_rectangle(left_lower_corner, right_upper_corner, V, F);

	Particle particles;
	Grid grid(left_lower_corner, right_upper_corner, h);

	// create grid	
	grid.setup(particles, height);

	

	// advection / move
	// advect_velocity(grid, particles, dt);
	// transfer to grid and save velocity
	transfer_to_grid(grid, particles);
	// add gravity
	add_gravity(grid, particles, gravity, dt);

	// do pressure






	std::cout << R"(
    )";
	std::cout << "\n";

	const auto add_points = [&]() {
		const Eigen::RowVector3d blue(0.1, 0.35, 0.75);
		viewer.data_list[xid].set_points(particles.q, (1. - (1. - blue.array()) * .8));
	};

	const auto update = [&]() {};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
			break;
		default:
			return false;
		}
		update();
		return true;
	};

	// start the viewer and set the mesh
	viewer.data().set_mesh(V, F);
	update();
	add_points();
	viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}