#include <grid_to_particle.h>
#include <interpolate.h>

void grid_to_particle_velocity_update(Grid& grid, Particle& particles) {
	Eigen::SparseMatrix<double> Wx, Wy, Wz;
	Eigen::VectorXd Vx, Vy, Vz;
	if (grid.method == PIC) {
		Vx = grid.Vx;
		Vy = grid.Vy;
		Vz = grid.Vz;
	} else if (grid.method == FLIP) {
		Vx = grid.Vx - grid.Vx_;
		Vy = grid.Vy - grid.Vy_;
		Vz = grid.Vz - grid.Vz_;
	}

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	Eigen::RowVector3d h = grid.h.transpose();
	Eigen::RowVector3d corner = grid.left_lower_corner.transpose();

	Eigen::VectorXd sum;
	interpolate(grid.nx, grid.ny, grid.nz, 0, h, corner + h.cwiseProduct(Eigen::RowVector3d(0.5, 0, 0)), particles.q, sum, Wx);
	interpolate(grid.nx, grid.ny, grid.nz, 1, h, corner + h.cwiseProduct(Eigen::RowVector3d(0, 0.5, 0)), particles.q, sum, Wy);
	interpolate(grid.nx, grid.ny, grid.nz, 2, h, corner + h.cwiseProduct(Eigen::RowVector3d(0, 0, 0.5)), particles.q, sum, Wz);
	
	if (grid.method == PIC) {
		particles.v.col(0) = Wx * Vx;
		particles.v.col(1) = Wy * Vy;
		particles.v.col(2) = Wz * Vz;
	}
	else if (grid.method == FLIP) {
		particles.v.col(0) += Wx * Vx;
		particles.v.col(1) += Wy * Vy;
		particles.v.col(2) += Wz * Vz;
	}
}