#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/slice_mask.h>

#include <iostream>
#include <chrono>
#include <random>

#include <timer.h>
#include <grid.h>
#include <particles.h>
#include <create_rectangle.h>

#include <advect_velocity.h>
#include <transfer_to_grid.h>
#include <add_gravity.h>
#include <grid_to_particle.h>

int main(int argc, char** argv) {
	igl::opengl::glfw::Viewer viewer;
	const int airID = viewer.selected_data_index;
	viewer.data().point_size = 2;
	viewer.append_mesh();
	const int fluidID = viewer.selected_data_index;
	viewer.data().point_size = 2;
	viewer.append_mesh();

	std::cout << R"(
    )";
	std::cout << "\n";


	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// igl::read_triangle_mesh((argc > 1 ? argv[1] : "../../../data/bunny.off"), V, F);
	
	double dt = 0.01;
	double height = 0.08;
	Eigen::Vector3d gravity = Eigen::Vector3d(0, -9.8, 0);
	Eigen::Vector3d h = Eigen::Vector3d(0.01, 0.01, 0.01);
	Eigen::Vector3d left_lower_corner = Eigen::Vector3d(-0.5, -0.5, -0.5);
	Eigen::Vector3d right_upper_corner = Eigen::Vector3d(0.5, 0.5, 0.5);

	// for visualizing boundary
	create_rectangle(left_lower_corner, right_upper_corner, V, F);

	Particle particles;
	Grid grid(left_lower_corner, right_upper_corner, h, FLIP);

	// create grid
	grid.init();
	grid.add_fluid(particles, height);

	// a complete time step
	const auto& step_update = [&]() {
		double t_before, t_after;

		// advection / move
		t_before = igl::get_seconds();
		advect_velocity(grid, particles, dt);
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";
		
		// transfer to grid and save velocity
		t_before = igl::get_seconds();
		transfer_to_grid(grid, particles);
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";
		
		// add gravity
		t_before = igl::get_seconds();
		add_gravity(grid, particles, gravity, dt);
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";

		// TODO: compute_distance_to_fluid
		// TODO: extend_velocity
		t_before = igl::get_seconds();
		grid.apply_boundary_condition();
		grid.pressure_projection();
		grid.get_divergence();
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";
		
		// TODO: make incompressible
		// TODO: extend_velocity
		t_before = igl::get_seconds();
		grid_to_particle_velocity_update(grid, particles);
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";
	};


	const auto show_particles = [&]() {
		Eigen::MatrixXd air_particles, fluid_particles;

		igl::slice_mask(particles.q, particles.type.array() == AIR_P, 1, air_particles);
		igl::slice_mask(particles.q, particles.type.array() == FLUID_P, 1, fluid_particles);

		const Eigen::RowVector3d blue(0.1, 0.35, 0.75);
		const Eigen::RowVector3d white(0.9, 0.9, 0.9);
		viewer.data_list[airID].set_points(air_particles, (1. - (1. - white.array()) * .8));
		viewer.data_list[fluidID].set_points(fluid_particles, (1. - (1. - blue.array()) * .8));
	};

	const auto update = [&]() {
		step_update();
		show_particles();
	};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
			for (int i = 0; i < 100; i++) {
				std::cerr << timer(update) << " seconds \n";
			}
			break;
		default:
			return false;
		}
		return true;
	};

	// start the viewer and set the mesh
	viewer.data().set_mesh(V, F);
	show_particles();
	viewer.data().show_lines = false;
	viewer.launch();
	return EXIT_SUCCESS;
}