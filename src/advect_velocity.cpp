#include <advect_velocity.h>
#include <interpolate.h>


void advect_velocity(
	Grid& grid, 
	Particle& particles, 
	double dt) {

	// boundary	for clamping
	Eigen::Vector3d v1 = grid.left_lower_corner + 1.001 * grid.h;
	Eigen::Vector3d v2 = grid.right_upper_corner - 1.001 * grid.h;

	// Simple update
	if (!grid.use_RK2_advection) {
		particles.q = particles.q + dt * particles.v;
	}
	else {
		// Runge-kutta (RK2)
		int nx = grid.nx;
		int ny = grid.ny;
		int nz = grid.nz;
		Eigen::RowVector3d h = grid.h;
		Eigen::RowVector3d corner = grid.left_lower_corner;

		Eigen::VectorXd sum;
		Eigen::MatrixXd q = particles.q;
		Eigen::MatrixXd v = particles.v;
		Eigen::SparseMatrix<double> Wx, Wy, Wz;

		// first step
		interpolate(grid.nx, grid.ny, grid.nz, 0, h, corner, q, sum, Wx);
		interpolate(grid.nx, grid.ny, grid.nz, 1, h, corner, q, sum, Wy);
		interpolate(grid.nx, grid.ny, grid.nz, 2, h, corner, q, sum, Wz);
		if (grid.theta == 0) {
			q.col(0) = q.col(0) + 0.5 * dt * Wx * grid.Vx;
			q.col(1) = q.col(1) + 0.5 * dt * Wy * grid.Vy;
			q.col(2) = q.col(2) + 0.5 * dt * Wz * grid.Vz;
		}
		else {
			q.col(0) = q.col(0) + 0.5 * dt * (v.col(0) + Wx * (grid.Vx - grid.Vx_));
			q.col(1) = q.col(1) + 0.5 * dt * (v.col(1) + Wy * (grid.Vy - grid.Vy_));
			q.col(2) = q.col(2) + 0.5 * dt * (v.col(2) + Wz * (grid.Vz - grid.Vz_));
		}
		q.col(0) = q.col(0).unaryExpr([&](double x) { return std::max(v1(0), std::min(v2(0), x)); });
		q.col(1) = q.col(1).unaryExpr([&](double x) { return std::max(v1(1), std::min(v2(1), x)); });
		q.col(2) = q.col(2).unaryExpr([&](double x) { return std::max(v1(2), std::min(v2(2), x)); });

		// second step
		interpolate(grid.nx, grid.ny, grid.nz, 0, h, corner, q, sum, Wx);
		interpolate(grid.nx, grid.ny, grid.nz, 1, h, corner, q, sum, Wy);
		interpolate(grid.nx, grid.ny, grid.nz, 2, h, corner, q, sum, Wz);
		if (grid.theta == 0) {
			particles.q.col(0) = particles.q.col(0) + dt * Wx * grid.Vx;
			particles.q.col(1) = particles.q.col(1) + dt * Wy * grid.Vy;
			particles.q.col(2) = particles.q.col(2) + dt * Wz * grid.Vz;
		}
		else {
			particles.q.col(0) = particles.q.col(0) + dt * (v.col(0) + Wx * (grid.Vx - grid.Vx_));
			particles.q.col(1) = particles.q.col(1) + dt * (v.col(1) + Wy * (grid.Vy - grid.Vy_));
			particles.q.col(2) = particles.q.col(2) + dt * (v.col(2) + Wz * (grid.Vz - grid.Vz_));
		}
	}
	// Clamping the particles to be within boundary
	particles.q.col(0) = particles.q.col(0).unaryExpr([&](double x) { return std::max(v1(0), std::min(v2(0), x)); });
	particles.q.col(1) = particles.q.col(1).unaryExpr([&](double x) { return std::max(v1(1), std::min(v2(1), x)); });
	particles.q.col(2) = particles.q.col(2).unaryExpr([&](double x) { return std::max(v1(2), std::min(v2(2), x)); });
}