#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <particles.h>

class Grid {

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

public:
	Eigen::Vector3d h;
	Eigen::Vector3d left_lower_corner;
	Eigen::Vector3d right_upper_corner;

	int nx, ny, nz;
	double density;
	Eigen::VectorXd Vx, Vy, Vz;
	Eigen::VectorXd pressure; // nx1 the pressure at each grid
	Eigen::VectorXi markers;  // nx1 marks the type of the cell
	Eigen::VectorXd divergence; // 3nx1 divergence of the each grid

	Eigen::SparseMatrix<double> Px, Py, Pz; // selection matrix
	Eigen::VectorXd Vx_, Vy_, Vz_;

	Grid(const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, double h_)
		:left_lower_corner(corner1), right_upper_corner(corner2), h(h_ * Eigen::Vector3d::Ones())
	{}

	Grid(const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, const Eigen::Vector3d &h_)
		:left_lower_corner(corner1), right_upper_corner(corner2), h(h_)
	{}
	
	int get_idx(const int& xi, const int& yi, const int& zi);

	void init();
	void add_fluid(Particle& particles, const double& height);
	void apply_boundary_condition();
	void pressure_projection();
	void get_divergence();
	void get_gradient();
	void save_grids();

	//// output:
	////grid - the grid index; ratio - dist from edge in each cell / length per cell
	//void ratio_from_edge(int& grid, int dim, float dist, float& ratio) {
	//	float h_per_cell = (dim == 0) ? h(0) : (dim == 1) ? h(1) : (dim == 2) ? h(2);
	//	float dist_from_side = (dim == 0) ? nx
	//		: (dim == 1) ? ny : (dim == 2) ? nz;
	//	dist_from_side /= h_per_cell;
	//	grid = floor(dist_from_side);
	//	ratio = dist_from_side - grid;
	//}

	//// output:
	////grid - the grid index; ratio - dist from center in each cell / length per cell
	//void ratio_from_center(int& grid, int dim, float dist, float& ratio) {
	//	int ncell = (dim == 0) ? nx
	//		: (dim == 1) ? ny
	//		: (dim == 2) ? nz;
	//	float h_per_cell = (dim == 0) ? h(0)
	//		: (dim == 1) ? h(1)
	//		: (dim == 2) ? h(2);
	//	float dist_from_side = ncell / h_per_cell - 0.5;
	//	grid = floor(dist_from_side);
	//	// if out of space - occurs when estimate
	//	if (grid < 0) {
	//		grid = 0;
	//		ratio = 0.0;
	//	}
	//	else if (grid > ncell - 2) { 
	//		grid = ncell - 2; 
	//		ratio = 1.0; }
	//	ratio = dist_from_side - grid;
	//}

	//void trilinear_interpolate(float x, float y, float z, float& out_x, float& out_y, float& out_z) {
	//	float rx, ry, rz;
	//	int grid_x, grid_y, grid_z;
	//	ratio_from_edge(grid_x, 0, x, rx);
	//	ratio_from_center(grid_y, 0, y, ry);
	//	pu = u.bilerp(i, j, fx, fy);
	//	dist_from_edge(px, i, fx);
	//	dist_from_center(py, j, fy);
	//	pv = v.bilerp(i, j, fx, fy);
	//}

private:
};



#endif // GRID_H
