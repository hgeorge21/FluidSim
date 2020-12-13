#ifndef GRID_H
#define GRID_H

#include <Eigen/Core>
#include <particles.h>

class Grid {
public:
	double h;
	Eigen::Vector3d left_lower_corner;
	Eigen::Vector3d right_upper_corner;

	int nx, ny, nz;
	double density;
	Eigen::VectorXd Vx, Vy, Vz;
	Eigen::VectorXd pressure; // the pressure at each grid

	Grid(Eigen::Vector3d corner1, Eigen::Vector3d corner2, double h_)
		:left_lower_corner(corner1), right_upper_corner(corner2), h(h_)
	{}

	void setup(Particle& particles);

private:
};



#endif // GRID_H
