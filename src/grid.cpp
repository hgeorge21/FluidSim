#include <grid.h>
#include <chrono>
#include <random>


void Grid::setup(Particle& particles, const double& height) {
	Eigen::Vector3d ns = (right_upper_corner - left_lower_corner).cwiseQuotient(h);
	nx = ceil(ns(0));
	ny = ceil(ns(1));
	nz = ceil(ns(2));

	// initialization
	int n_grids = nx * ny * nz;
	pressure = Eigen::VectorXd::Zero(n_grids);
	markers = Eigen::VectorXi::Zero(n_grids);

	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	int ny_f = ceil(height / h(1));
	int n_fluid_grids = nx * ny_f * nz;
	particles.q = Eigen::MatrixXd(8 * n_fluid_grids, 3); // each cell has 8 particles
	particles.v = Eigen::MatrixXd::Zero(8 * n_fluid_grids, 3);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz; k++) {
				int idx = i * (ny_f * nz) + j * nz + k;
				
				// generate fluid particles
				Eigen::RowVector3d lower_corner = left_lower_corner + h.cwiseProduct(Eigen::Vector3d(i, j, k));
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(8, 3, [&]() { return dist(generator);  });
				q_.col(0) = h(0) * q_.col(0);
				q_.col(1) = h(1) * q_.col(1);
				q_.col(2) = h(2) * q_.col(2);

				q_ = q_.rowwise() + lower_corner;
				particles.q.block(8 * idx, 0, 8, 3) = q_;
			}
		}
	}
}


void Grid::save_grids() {
	Vx_ = Vx;
	Vy_ = Vy;
	Vz_ = Vz;
}