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

	Eigen::RowVector3d h;
	Eigen::RowVector3d left_lower_corner;
	Eigen::RowVector3d right_upper_corner;

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
	//int get_idx2(const int &xi, const int &yi, const int &zi, const int &dim);

	void init();
	void add_fluid(Particle& particles, const Eigen::RowVector3d& v1, const Eigen::RowVector3d& v2, const Eigen::RowVector3d& hf, const double& height);
	void apply_boundary_condition();
	void pressure_projection();

	void get_divergence();
	
	void solve_pressure();
	void update_velocity();
	void save_grids();



	/* Free boundary related functions */
	Eigen::VectorXi cp; // closest points in each grid
	Eigen::VectorXd phi; // distance function
	
	void init_phi(Particle &particles);
	void signed_distance(const Eigen::RowVector3d& v1, const Eigen::RowVector3d& pt, double& d);
	void sweep_phi(Particle& particles);
	void distance_to_fluid(Particle& particles);

	void sweep_velocity();
	void sweep_dir(const int& dir);

private:
	void get_divergence_operator();
	void get_laplacian_operator();
	void get_gradient_operator();
};



#endif // GRID_H
