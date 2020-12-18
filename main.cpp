#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/slice_mask.h>

#include <iostream>
#include <chrono>
#include <random>
#include <thread>
#include <omp.h>

#include <timer.h>
#include <grid.h>
#include <particles.h>
#include <create_rectangle.h>

#include <advect_velocity.h>
#include <transfer_to_grid.h>
#include <add_gravity.h>
#include <grid_to_particle.h>

int main(int argc, char** argv) {
	omp_set_num_threads(8);
	Eigen::setNbThreads(8);

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
	
	double dt = 1;
	double height = 0.3;
	Eigen::Vector3d h = Eigen::Vector3d(0.025, 0.08, 0.025);
	Eigen::Vector3d gravity = Eigen::Vector3d(0, -0.098, 0);
	Eigen::Vector3d left_lower_corner = Eigen::Vector3d(-0.2, -0.4, -0.2);
	Eigen::Vector3d right_upper_corner = Eigen::Vector3d(0.2, 0.4, 0.2);

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
		advect_velocity(grid, particles, dt);
		transfer_to_grid(grid, particles);
		add_gravity(grid, particles, gravity, dt);

		grid.apply_boundary_condition();

		t_before = igl::get_seconds();
		grid.pressure_projection();
		t_after = igl::get_seconds();
		std::cerr << t_after - t_before << "\n";
		
		grid_to_particle_velocity_update(grid, particles);
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

	const auto &update = [&]() {
		step_update();
		show_particles();
	};

	const auto& update_forever = [&]() {
		while (true) {
			std::cerr << "===========\n";
			std::cerr << timer(update) << " seconds\n";
		}
	};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
			for (int i = 0; i < 1; i++) {
				std::cerr << "===========\n";
				std::cerr << timer(update) << " seconds\n";
			}
			break;
		default:
			return false;
		}
		return true;
	};

	

	//std::thread simulation_thread(update_forever);
	//simulation_thread.detach();

	// start the viewer and set the mesh
	viewer.data().set_mesh(V, F);
	
	show_particles();
	viewer.data().show_lines = false;
	viewer.launch();
	return EXIT_SUCCESS;
}