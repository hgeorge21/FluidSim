#include <add_gravity.h>

void add_gravity(Grid& grid, Particle& particles, const Eigen::Vector3d& gravity, const double &dt) {
	grid.Vx = grid.Vx + dt * Eigen::VectorXd::Constant(grid.Vx.rows(), gravity(0));
	grid.Vy = grid.Vy + dt * Eigen::VectorXd::Constant(grid.Vy.rows(), gravity(1));
	grid.Vz = grid.Vz + dt * Eigen::VectorXd::Constant(grid.Vz.rows(), gravity(2));
}