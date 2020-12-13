#include <add_gravity.h>

void add_gravity(Grid& grid, Particle& particles, const Eigen::Vector3d& gravity, const double &dt) {
	Eigen::VectorXd g;
	g = gravity(0) * Eigen::VectorXd::Ones(grid.Vx.rows());
	grid.Vx = grid.Vx + dt * g;
	g = gravity(1) * Eigen::VectorXd::Ones(grid.Vy.rows());
	grid.Vy = grid.Vy + dt * g;
	g = gravity(2) * Eigen::VectorXd::Ones(grid.Vz.rows());
	grid.Vz = grid.Vz + dt * g;
}