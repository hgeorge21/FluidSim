#include <grid_to_particle.h>
#include <interpolate.h>

void grid_to_particle_velocity_update(Grid& grid, Particle& particles) {
	Eigen::SparseMatrix<double> Wx, Wy, Wz;
	Eigen::VectorXd Vx1, Vy1, Vz1;
	Eigen::VectorXd Vx2, Vy2, Vz2;
	
	if (grid.theta == 0) {
		Vx1 = grid.Vx;
		Vy1 = grid.Vy;
		Vz1 = grid.Vz;
	} 
	else if (grid.theta == 10) {
		Vx2 = grid.Vx - grid.Vx_;
		Vy2 = grid.Vy - grid.Vy_;
		Vz2 = grid.Vz - grid.Vz_;
	}
	else {
		Vx1 = grid.Vx;
		Vy1 = grid.Vy;
		Vz1 = grid.Vz;
		Vx2 = grid.Vx - grid.Vx_;
		Vy2 = grid.Vy - grid.Vy_;
		Vz2 = grid.Vz - grid.Vz_;
	}

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	Eigen::RowVector3d h = grid.h;
	Eigen::RowVector3d corner = grid.left_lower_corner;

	Eigen::VectorXd sum;
	interpolate(grid.nx, grid.ny, grid.nz, 0, h, corner, particles.q, sum, Wx);
	interpolate(grid.nx, grid.ny, grid.nz, 1, h, corner, particles.q, sum, Wy);
	interpolate(grid.nx, grid.ny, grid.nz, 2, h, corner, particles.q, sum, Wz);
	
	if (grid.theta == 0) {
		particles.v.col(0) = Wx * Vx1;
		particles.v.col(1) = Wy * Vy1;
		particles.v.col(2) = Wz * Vz1;
	}
	else if (grid.theta == 10) {
		particles.v.col(0) += Wx * Vx2;
		particles.v.col(1) += Wy * Vy2;
		particles.v.col(2) += Wz * Vz2;
	}
	else {
		Vx1 = Wx * Vx1;
		Vy1 = Wy * Vy1;
		Vz1 = Wz * Vz1;
		Vx2 = Wx * Vx2;
		Vy2 = Wy * Vy2;
		Vz2 = Wz * Vz2;
		Vx1 = Vx1 - particles.v.col(0);
		Vy1 = Vy1 - particles.v.col(1);
		Vz1 = Vz1 - particles.v.col(2);
		particles.v.col(0) += (1.0 - 0.1 * grid.theta) * Vx1 + 0.1 * grid.theta * Vx2;
		particles.v.col(1) += (1.0 - 0.1 * grid.theta) * Vy1 + 0.1 * grid.theta * Vy2;
		particles.v.col(2) += (1.0 - 0.1 * grid.theta) * Vz1 + 0.1 * grid.theta * Vz2;
	}
}