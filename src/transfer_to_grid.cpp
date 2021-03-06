#include <transfer_to_grid.h>
#include <interpolate.h>

void transfer_to_grid(Grid& grid, Particle& particles) {
	Eigen::SparseMatrix<double> Wx, Wy, Wz;

	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	Eigen::RowVector3d h = grid.h;
	Eigen::RowVector3d corner = grid.left_lower_corner;

	Eigen::VectorXd sum;
	interpolate(nx, ny, nz, 0, h, corner, particles.q, sum, Wx);
	grid.Vx = Wx.transpose() * particles.v.col(0);
	grid.Vx = (grid.Vx).cwiseQuotient(sum);

	interpolate(nx, ny, nz, 1, h, corner, particles.q, sum, Wy);
	grid.Vy = Wy.transpose() * particles.v.col(1);
	grid.Vy = (grid.Vy).cwiseQuotient(sum);

	interpolate(nx, ny, nz, 2, h, corner, particles.q, sum, Wz);
	grid.Vz = Wz.transpose() * particles.v.col(2);
	grid.Vz = (grid.Vz).cwiseQuotient(sum);

	// update the grid markers
	grid.markers.setConstant(AIRCELL);
	for (int n = 0; n < particles.q.rows(); n++) {
		Eigen::RowVector3d d = (particles.q.row(n) - corner).cwiseQuotient(h);
		int idx = int(d(0)) * (ny * nz) + int(d(1)) * nz + int(d(2));
		if (grid.markers(idx) != int(FLUIDCELL)) {
			grid.markers(idx) = int(FLUIDCELL);
		}
	}
	grid.save_grids();
}
