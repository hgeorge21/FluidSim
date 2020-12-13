#include <transfer_to_grid.h>
#include <interpolate.h>

void transfer_to_grid(Grid &grid, Particle &particles) {
	Eigen::SparseMatrix<double> Wx, Wy, Wz;
	
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	double h = grid.h;
	Eigen::RowVector3d corner = grid.left_lower_corner.transpose();

	interpolate(nx, ny, nz, 0, h, corner + h*Eigen::RowVector3d(0.5, 0, 0), particles.q, Wx);
	interpolate(nx, ny, nz, 1, h, corner + h*Eigen::RowVector3d(0, 0.5, 0), particles.q, Wy);
	interpolate(nx, ny, nz, 2,h, corner + h*Eigen::RowVector3d(0, 0, 0.5), particles.q, Wz);

	Eigen::VectorXd Vx, Vy, Vz;
	grid.Vx = Wx.transpose() * particles.v.col(0);
	grid.Vy = Wy.transpose() * particles.v.col(1);
	grid.Vz = Wz.transpose() * particles.v.col(2);
}