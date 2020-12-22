#include <add_particles.h>
#include <sample_from_mesh.h>

#include <chrono>
#include <random>

void add_particles(
	Grid& grid, 
	Particle& particles, 
	const Eigen::RowVector3d& v1, 
	const Eigen::RowVector3d& v2, 
	const Eigen::RowVector3d& hf, 
	const int &n_per_cell) {

	// calculate for fluids
	Eigen::RowVector3d ns = (v2 - v1).cwiseQuotient(hf);
	int nx_ = ceil(ns(0));
	int ny_ = ceil(ns(1));
	int nz_ = ceil(ns(2));

	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	int ny_f = ny_;
	int nx_f = nx_;
	int nz_f = nz_;

	// loop over specified cells
	int n_fluid_grids = nx_f * ny_f * nz_f;
	particles.q = Eigen::MatrixXd(n_per_cell * n_fluid_grids, 3);
	particles.v = Eigen::MatrixXd::Zero(n_per_cell * n_fluid_grids, 3);
	particles.type = Eigen::VectorXi::Zero(n_per_cell * n_fluid_grids);
	for (int i = 0; i < nx_f; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz_f; k++) {
				int idx = i * (ny_f * nz_f) + j * nz_f + k;

				// generate fluid particles
				Eigen::RowVector3d lower_corner = v1 + hf.cwiseProduct(Eigen::RowVector3d(i, j, k));
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(n_per_cell, 3, [&]() { return dist(generator);  });
				q_.col(0) = hf(0) * q_.col(0);
				q_.col(1) = hf(1) * q_.col(1);
				q_.col(2) = hf(2) * q_.col(2);
				q_ = q_.rowwise() + lower_corner;
				particles.q.block(n_per_cell * idx, 0, n_per_cell, 3) = q_;
				particles.type.block(n_per_cell * idx, 0, n_per_cell, 1) = FLUID_P * Eigen::VectorXi::Ones(n_per_cell);
			}
		}
	}
}


void add_particles(
	Grid& grid, 
	Particle& particles, 
	const Eigen::RowVector3d& v1, 
	const Eigen::RowVector3d& v2, 
	const Eigen::RowVector3d& hf, 
	const int& n_per_cell,
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F) {

	Eigen::MatrixXd tmp;
	sample_from_mesh(V, F, 10000, tmp);

	Eigen::RowVector3d ns = (v2 - v1).cwiseQuotient(hf);
	int nx_ = ceil(ns(0));
	int ny_ = ceil(ns(1));
	int nz_ = ceil(ns(2));

	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	int ny_f = ny_;
	int nx_f = nx_;
	int nz_f = nz_;

	// loop over specified cells
	int n_fluid_grids = nx_f * ny_f * nz_f;
	particles.q = Eigen::MatrixXd(n_per_cell * n_fluid_grids + tmp.rows(), 3);
	particles.v = Eigen::MatrixXd::Zero(n_per_cell * n_fluid_grids + tmp.rows(), 3);
	particles.type = Eigen::VectorXi::Zero(n_per_cell * n_fluid_grids + tmp.rows());
	for (int i = 0; i < nx_f; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz_f; k++) {
				int idx = i * (ny_f * nz_f) + j * nz_f + k;

				// generate fluid particles
				Eigen::RowVector3d lower_corner = v1 + hf.cwiseProduct(Eigen::RowVector3d(i, j, k));
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(n_per_cell, 3, [&]() { return dist(generator);  });
				q_.col(0) = hf(0) * q_.col(0);
				q_.col(1) = hf(1) * q_.col(1);
				q_.col(2) = hf(2) * q_.col(2);
				q_ = q_.rowwise() + lower_corner;
				particles.q.block(n_per_cell * idx, 0, n_per_cell, 3) = q_;
				particles.type.block(n_per_cell * idx, 0, n_per_cell, 1) = FLUID_P * Eigen::VectorXi::Ones(n_per_cell);
			}
		}
	}
	
	// add the input mesh
	particles.q.block(n_per_cell * n_fluid_grids, 0, tmp.rows(), 3) = tmp;
	particles.type.block(n_per_cell * n_fluid_grids, 0, tmp.rows(), 1).setConstant(FLUID_P);
}