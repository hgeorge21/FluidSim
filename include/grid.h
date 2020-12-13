#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
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
	Eigen::VectorXd pressure; // the pressure at each grid
	Eigen::VectorXi markers;

	Eigen::VectorXd Vx_, Vy_, Vz_;

	Grid(const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, double h_)
		:left_lower_corner(corner1), right_upper_corner(corner2), h(h_ * Eigen::Vector3d::Ones())
	{}

	Grid(const Eigen::Vector3d &corner1, const Eigen::Vector3d &corner2, const Eigen::Vector3d &h_)
		:left_lower_corner(corner1), right_upper_corner(corner2), h(h_)
	{}
	

	void setup(Particle& particles);

	void save_grids();

private:
};



#endif // GRID_H
