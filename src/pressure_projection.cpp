#include <pressure_projection.h>

#include <finite_difference.h>

void pressure_projection(Grid &grid, Particle &particles, const double &dt) {
	/*Eigen::VectorXd dV = Eigen::VectorXd::Zero(grid.nx * grid.ny * grid.nz);
	Eigen::SparseMatrix<double> Dv;

	finite_difference(grid.nx, grid.ny, grid.nz, 0, grid.h(0), grid.Vx, Dv);
	dV = dV + Dv * grid.Vx;
	finite_difference(grid.nx, grid.ny, grid.nz, 1, grid.h(1), grid.Vy, Dv);
	dV = dV + Dv * grid.Vy;
	finite_difference(grid.nx, grid.ny, grid.nz, 2, grid.h(2), grid.Vz, Dv);
	dV = dV + Dv * grid.Vz;

	dV = grid.density / dt * dV;*/
}