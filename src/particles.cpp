#include <particles.h>
#include <chrono>
#include <random>
#include <iostream>
Particle::Particle()
{
}

Particle::~Particle()
{
}

void Particle::generate_particles(Grid& grid) {
	double h = grid.h;
	Eigen::Vector3d ns = (grid.right_upper_corner - grid.left_lower_corner) / h;
	int nx = int(ns(0));
	int ny = int(ns(1));
	int nz = int(ns(2));

	int n_particles = 8 * nx*ny*nz;
	q = Eigen::MatrixXd(n_particles, 3); // each cell has 8 particles
	v = Eigen::MatrixXd::Zero(n_particles, 3);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				int idx = 8 * (i * (ny * nz) + j * nz + k);
				Eigen::Vector3d lower_corner = grid.left_lower_corner + h * Eigen::Vector3d(1.*i, 1.*j, 1.*k);
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(8, 3, [&]() { return dist(generator);  });
				q.block(idx, 0, 8, 3) = (lower_corner + h * q_).transpose();
			}
		}
	}
}
