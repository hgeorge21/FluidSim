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

#define AIRCELL 0
#define FLUIDCELL 1
#define SOLIDCELL 2

public:
	double dt;
	transfer_method method;

	Eigen::Vector3d h;
	Eigen::Vector3d left_lower_corner;
	Eigen::Vector3d right_upper_corner;

	int nx, ny, nz, n_grids;
	double density = 1.0;
	Eigen::VectorXd Vx, Vy, Vz;
	Eigen::VectorXd pressure;   // nx1 the pressure at each grid
	Eigen::VectorXi markers;    // nx1 marks the type of the cell
	Eigen::VectorXd divergence; // nx1 divergence of the each grid
	Eigen::VectorXd gradient;   // 7nx1 gradient;
	Eigen::SparseMatrix<double> A;

	Eigen::SparseMatrix<double> Px, Py, Pz; // selection matrix
	Eigen::SparseMatrix<double> Dx, Dy, Dz; // divergence operator
	Eigen::SparseMatrix<double> Gx, Gy, Gz; // gradient operator
	Eigen::VectorXd Vx_, Vy_, Vz_;

	Grid(double dt_, const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, double h_, transfer_method method_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), h(h_ * Eigen::Vector3d::Ones()), method(method_)
	{}

	Grid(double dt_, const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, const Eigen::Vector3d &h_, transfer_method method_)
		:dt(dt_), left_lower_corner(corner1), right_upper_corner(corner2), h(h_), method(method_)
	{}
	
	int get_idx(const int& xi, const int& yi, const int& zi);

	void init();
	void add_fluid(Particle& particles, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& hf, const double& height);
	void apply_boundary_condition();
	void pressure_projection();

	void get_divergence();
	void get_divergence2();

	void check_divergence();
	void check_laplacian();
	
	void solve_pressure();
	void update_velocity();
	void save_grids();

	void print_pressure();
	void check_pressure();
	void print_cell();

private:
	void get_divergence_operator();
	void get_laplacian_operator();
	void get_gradient_operator();
	void get_laplacian_operator2();
};



#endif // GRID_H
