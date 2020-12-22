#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <particles.h>

class Grid {

#define AIRCELL 1
#define FLUIDCELL -1
#define SOLIDCELL 0

public:
	double dt;
	double density = 1.0;
	int theta = 10; // for FLIP / PIC blend
	bool use_RK2_advection = true;

	Eigen::RowVector3d h;
	Eigen::RowVector3d g;
	Eigen::RowVector3d left_lower_corner; 
	Eigen::RowVector3d right_upper_corner;

	int num_fluid_cells;
	int nx, ny, nz, n_grids;
	Eigen::VectorXd Vx, Vy, Vz;
	Eigen::VectorXd pressure;   // nx1 the pressure at each cell
	Eigen::VectorXi markers;    // nx1 marks the type of the cell
	Eigen::SparseMatrix<double> A;

	Eigen::VectorXd phi; // distance function
	Eigen::VectorXi fluid_map; // mapping grid to fluid index
	Eigen::VectorXd Vx_, Vy_, Vz_; // saved grid velocity

	Grid(double dt_, const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2, const Eigen::Vector3d& gravity, double h_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), g(gravity), h(h_* Eigen::Vector3d::Ones())
	{}

	Grid(double dt_, const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2, const Eigen::Vector3d& gravity, const Eigen::Vector3d& h_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), g(gravity), h(h_)
	{}

	// initialize data structures
	void init();

	// gets the grid index give x, y, z positions
	int get_idx(const int& xi, const int& yi, const int& zi);

	// generates the free surface, phi, between liquid and air with kernel smoothing
	void create_free_boundary(Particle &particles);
	
	// apply boundary condition and filters out velocity on solid cells
	void apply_boundary_condition();

	// save the velocities for FLIP
	void save_grids();

private:
};

#endif // GRID_H
