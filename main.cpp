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
	Eigen::MatrixXd V, V1;
	Eigen::MatrixXi F, F1;
	
	double dt = 0.3;
	double height = 0.3;
	Eigen::Vector3d h = Eigen::Vector3d(0.1, 0.1, 0.1);
	Eigen::Vector3d gravity = Eigen::Vector3d(0, -0.098, 0);
	Eigen::Vector3d left_lower_corner = Eigen::Vector3d(-0.5, -0.5, -0.5);
	Eigen::Vector3d right_upper_corner = Eigen::Vector3d(0.5, 0.5, 0.5);

	// for visualizing boundary
	create_rectangle(left_lower_corner, right_upper_corner, V, F);
	create_rectangle(left_lower_corner + h, right_upper_corner - h, V1, F1);

	Particle particles;
	Grid grid(dt, left_lower_corner, right_upper_corner, h, PIC);

	// create grid
	grid.init();
	grid.add_fluid(particles, Eigen::Vector3d(-0.1, -0.1, -0.1), Eigen::Vector3d(0.2, 0.2, 0.2), Eigen::Vector3d(0.05, 0.05, 0.05), height);

	// a complete time step
	const auto& step_update = [&]() {
		double t_before, t_after;

		// advection / move

		std::cerr << "Now advect_velocity" << std::endl;
		advect_velocity(grid, particles, dt);
		std::cerr << "Now transfer_to_grid" << std::endl;
		transfer_to_grid(grid, particles);
		std::cerr << "Now add_gravity" << std::endl;
		add_gravity(grid, particles, gravity, dt);
		std::cerr << "Now apply_boundary_condition" << std::endl;
		grid.apply_boundary_condition();
		std::cerr << "Now pressure_projection" << std::endl;
		grid.pressure_projection();

		std::cerr << "Now grid_to_particle_velocity_update" << std::endl;
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

		std::cerr << "Now step update" << std::endl;
		step_update();

		std::cerr << "Now show particles" << std::endl;
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
	//viewer.data().set_mesh(V1, F1);
	
	show_particles();
	viewer.data().show_lines = false;
	viewer.launch();
	return EXIT_SUCCESS;
}
