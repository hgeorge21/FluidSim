#include <grid.h>
#include <chrono>
#include <random>

#include <iostream>

void Grid::setup(Particle& particles) {
	Eigen::Vector3d ns = (right_upper_corner - left_lower_corner) / h;
	nx = ceil(ns(0));
	ny = ceil(ns(1));
	nz = ceil(ns(2));

	// initialization
	int n_grids = nx * ny * nz;
	particles.q = Eigen::MatrixXd(8 * n_grids, 3); // each cell has 8 particles
	particles.v = Eigen::MatrixXd::Zero(8 * n_grids, 3);
	Vx = Eigen::VectorXd::Zero(n_grids+ny*nz);
	Vy = Eigen::VectorXd::Zero(n_grids+nx*nz);
	Vz = Eigen::VectorXd::Zero(n_grids+nx*ny);
	pressure = Eigen::VectorXd::Zero(n_grids, 1);

	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				int idx = i * (ny * nz) + j * nz + k;
				
				// generate fluid particles
				Eigen::Vector3d lower_corner = left_lower_corner + h * Eigen::Vector3d(1. * i, 1. * j, 1. * k);
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(8, 3, [&]() { return dist(generator);  });
				particles.q.block(8 * idx, 0, 8, 3) = (lower_corner + h * q_).transpose();
			}
		}
	}
}