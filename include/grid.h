#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <particles.h>

enum transfer_method {
	PIC,
	FLIP
};


class Grid {

#define AIRCELL 1
#define FLUIDCELL -1
#define SOLIDCELL 0

public:
	double dt;
	double density = 1.0;
	transfer_method method;

	Eigen::RowVector3d h;
	Eigen::RowVector3d left_lower_corner;
	Eigen::RowVector3d right_upper_corner;

	int num_fluid_cells;
	int nx, ny, nz, n_grids;
	Eigen::VectorXd Vx, Vy, Vz;
	Eigen::VectorXd pressure;   // nx1 the pressure at each grid
	Eigen::VectorXi markers;    // nx1 marks the type of the cell
	Eigen::SparseMatrix<double> A;

	Eigen::VectorXd phi; // distance function
	Eigen::VectorXi fluid_map; // mapping grid to fluid index
	Eigen::VectorXd Vx_, Vy_, Vz_; // saved grid velocity

	Grid(double dt_, const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2, double h_, transfer_method method_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), h(h_* Eigen::Vector3d::Ones()), method(method_)
	{}

	Grid(double dt_, const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2, const Eigen::Vector3d& h_, transfer_method method_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), h(h_), method(method_)
	{}

	int get_idx(const int& xi, const int& yi, const int& zi);

	void init();
	void create_free_boundary();
	void apply_boundary_condition();
	void save_grids();

private:
};

#endif // GRID_H
