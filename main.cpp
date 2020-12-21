#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>

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
#include <add_particles.h>
#include <grid_to_particle.h>
#include <pressure_projection.h>
#include <solve_pressure.h>
#include <solve_pressure_free_surface.h>

int main(int argc, char** argv) {
	igl::opengl::glfw::Viewer viewer;
	const int fluidID = viewer.selected_data_index;
	viewer.data().point_size = 5;
	viewer.append_mesh();
	igl::ColorMapType cmap_type = igl::COLOR_MAP_TYPE_MAGMA;

	int n_per_cell = 4;
	bool use_free_boundary = true;

	std::cout << R"(
  [space] Performs one time step
  +       Increase the number of particles per cell
  -       Decrease the number of particles per cell
  R,r     Resets the simulation
  M,m     Toggle between free surface or without
  F,f     Toggle between FLIP and PIC
    )";
	std::cout << "\n";

	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd CM;

	double dt = 0.1;
	double height = 0.3;
	Eigen::Vector3d h = Eigen::Vector3d(0.1, 0.1, 0.1);
	Eigen::Vector3d gravity = Eigen::Vector3d(0, -0.5, 0);
	Eigen::Vector3d left_lower_corner = Eigen::Vector3d(-1.0, -1.0, -1.0);
	Eigen::Vector3d right_upper_corner = Eigen::Vector3d(1.0, 1.0, 1.0);

	Eigen::Vector3d fluid_v1 = Eigen::Vector3d(-0.3, -0.9, -0.3);
	Eigen::Vector3d fluid_v2 = Eigen::Vector3d(0.3, 0.1, 0.3);
	Eigen::Vector3d fluid_hf = Eigen::Vector3d(0.04, 0.04, 0.04);

	// for visualizing boundary
	create_rectangle(left_lower_corner + h, right_upper_corner - h, V, F);
	Particle particles;
	Grid grid(dt, left_lower_corner, right_upper_corner, h, FLIP);


	const auto& setup = [&]() {
		// create grid and add particles
		grid.init();
		add_particles(grid, particles, fluid_v1, fluid_v2, fluid_hf, n_per_cell);
		double ymin = particles.q.col(1).minCoeff();
		double ymax = particles.q.col(1).maxCoeff();
		igl::colormap(cmap_type, particles.q.col(1), ymin, ymax, CM);
	};

	// a complete time step
	const auto& step_update = [&]() {
		advect_velocity(grid, particles, dt);
		transfer_to_grid(grid, particles);
		add_gravity(grid, particles, gravity, dt);
		grid.apply_boundary_condition();
		
		if (!use_free_boundary)
			solve_pressure(grid);
		else {
			grid.create_free_boundary();
			solve_pressure_free_surface(grid);
		}
		pressure_projection(grid);
		grid_to_particle_velocity_update(grid, particles);
	};


	const auto show_particles = [&]() {
		viewer.data_list[fluidID].set_points(particles.q, CM);
	};

	const auto& update = [&]() {
		step_update();
		show_particles();
	};

	const auto& update_forever = [&]() {
		int i = 0;
		while (true) {
			//std::cerr << "======================" << std::endl;
			//std::cout << "Step " << i << std::endl;
			//std::cerr << timer(update) << " seconds" << std::endl;
			update();
			i++;
		}
	};

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
			//std::cout << timer(update) << " seconds\n";
			update();
			break;
		case 'r':
		case 'R':
			setup();
			show_particles();
			break;
		case 'm':
		case 'M':
			use_free_boundary = !use_free_boundary;
			std::cout << "Changed to " << ((use_free_boundary) ? "free surface" : " without free surface") << std::endl;
			break;
		case '+':
			n_per_cell += 1;
			std::cout << "# Particles per cell: " << n_per_cell << std::endl;
			setup();
			show_particles();
			break;
		case '-':
			n_per_cell -= 1;
			std::cout << "# Particles per cell: " << n_per_cell << std::endl;
			setup();
			show_particles();
			break;
		case 'F':
		case 'f':
			grid.method = (grid.method == transfer_method::FLIP) ? transfer_method::PIC : transfer_method::FLIP;
			std::cout << "Change method to " << ((grid.method == transfer_method::FLIP) ? "FLIP" : "PIC") << std::endl;
			break;
		default:
			return false;
		}
		return true;
	};

	setup();
	show_particles();

	if (argc > 1) {
		std::thread simulation_thread(update_forever);
		simulation_thread.detach();
	}

	// start the viewer and set the mesh
	viewer.data().set_mesh(V, F);
	viewer.data().show_lines = true;
	viewer.data().show_faces = false;
	viewer.launch();
	return EXIT_SUCCESS;
}