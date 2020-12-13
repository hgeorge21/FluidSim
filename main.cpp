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

#include <transfer_to_grid.h>

int main(int argc, char** argv) {
	// load mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// igl::read_triangle_mesh((argc > 1 ? argv[1] : "../../../data/bunny.off"), V, F);
	create_rectangle(V, F);

	igl::opengl::glfw::Viewer viewer;
	const int xid = viewer.selected_data_index;
	viewer.data().point_size = 2;
	viewer.append_mesh();

	// create grid
	Grid grid(Eigen::Vector3d(-.5, -.3, -.5), Eigen::Vector3d(.5, -0.2, .5), 0.1);
	Particle particles;
	grid.setup(particles);


	transfer_to_grid(grid, particles);
	




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