#include <grid.h>
#include <chrono>
#include <random>
#include <iostream>
#include <divergence.h>
#include <fstream>


void Grid::init() {
	// calculate number of cells in each dimension
	Eigen::Vector3d ns = (right_upper_corner - left_lower_corner).cwiseQuotient(h);
	nx = ceil(ns(0));
	ny = ceil(ns(1));
	nz = ceil(ns(2));

	// initialization
	n_grids = nx * ny * nz;
	pressure = Eigen::VectorXd::Zero(n_grids);
	markers = Eigen::VectorXi::Zero(n_grids);
	divergence = Eigen::VectorXd::Zero(n_grids);

	// sets up the selection matrix for the boundaries
	typedef Eigen::Triplet<double> T;
	std::vector<T> trip_x, trip_y, trip_z;

	Px.resize((nx + 1) * ny * nz, (nx + 1) * ny * nz);
	Py.resize(nx * (ny + 1) * nz, nx * (ny + 1) * nz);
	Pz.resize(nx * ny * (nz + 1), nx * ny * (nz + 1));
	Px.setZero();
	Py.setZero();
	Pz.setZero();

	// first & last two edges are 0
	const auto& set_sparse_vec = [&](const int& dim) {
		int ll0 = (dim == 0) ? 2 : 0;
		int ll1 = (dim == 1) ? 2 : 0;
		int ll2 = (dim == 2) ? 2 : 0;
		int ul0 = (dim == 0) ? nx-1 : nx;
		int ul1 = (dim == 1) ? ny-1 : ny;
		int ul2 = (dim == 2) ? nz-1 : nz;
		int dim0 = (dim == 0) ? nx+1 : nx;
		int dim1 = (dim == 1) ? ny+1 : ny;
		int dim2 = (dim == 2) ? nz+1 : nz;


		for (int i = ll0; i < ul0; i++)
			for (int j = ll1; j < ul1; j++)
				for (int k = ll2; k < ul2; k++) {
					int idx = i * dim1* dim2 + j * dim2 + k;

					if (dim == 0)
						trip_x.emplace_back(T(idx, idx, 1.0));
					if (dim == 1)
						trip_y.emplace_back(T(idx, idx, 1.0));
					if (dim == 2)
						trip_z.emplace_back(T(idx, idx, 1.0));
				}
	};

	set_sparse_vec(0);
	set_sparse_vec(1);
	set_sparse_vec(2);
	Px.setFromTriplets(trip_x.begin(), trip_x.end());
	Py.setFromTriplets(trip_y.begin(), trip_y.end());
	Pz.setFromTriplets(trip_z.begin(), trip_z.end());
}


void Grid::add_fluid(Particle& particles, const Eigen::RowVector3d &v1, const Eigen::RowVector3d &v2, const Eigen::RowVector3d &hf, const double& height) {
	// calculate for fluids
	Eigen::RowVector3d ns = (v2 - v1).cwiseQuotient(hf);
	int nx_ = ceil(ns(0));
	int ny_ = ceil(ns(1));
	int nz_ = ceil(ns(2));

	int n_per_cell = 8; // each cell has 8 particles
	
	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	int fluid_height = ceil(height / hf(1));
	int ny_f = ny_;
	int nx_f = nx_;
	int nz_f = nz_;

	int n_fluid_grids = nx_f * ny_f * nz_f;
	particles.q = Eigen::MatrixXd(n_per_cell * n_fluid_grids, 3);
	particles.v = Eigen::MatrixXd::Zero(n_per_cell * n_fluid_grids, 3);
	particles.type = Eigen::VectorXi::Zero(n_per_cell * n_fluid_grids);
	for (int i = 0; i < nx_f; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz_f; k++) {
				int idx = i * (ny_f * nz_f) + j * nz_f + k;

				// generate fluid particles
				Eigen::RowVector3d lower_corner = v1 + hf.cwiseProduct(Eigen::Vector3d(i, j, k));
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


void Grid::apply_boundary_condition() {
	int i, j, k;
	int id;

	// boundary cells
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			int k1 = 0; int k2 = nz-1;
			markers(get_idx(i, j, k1)) = SOLIDCELL;
			markers(get_idx(i, j, k2)) = SOLIDCELL;
		}
		for (k = 0; k < nz; k++) {
			int j1 = 0; int j2 = ny - 1;
			markers(get_idx(i, j1, k)) = SOLIDCELL;
			markers(get_idx(i, j2, k)) = SOLIDCELL;
		}
	}
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			int i1 = 0; int i2 = nx - 1;
			markers(get_idx(i1, j, k)) = SOLIDCELL;
			markers(get_idx(i2, j, k)) = SOLIDCELL;
		}
	}

	// filters out velocity at boundary grids
	Vx = Px * Vx;
	Vy = Py * Vy;
	Vz = Pz * Vz;
}


/*** 
// The operators, for loop would actually be faster...
***/
void Grid::get_divergence_operator() {
	divergence_op(nx, ny, nz, 0, h, markers, Dx);
	divergence_op(nx, ny, nz, 1, h, markers, Dy);
	divergence_op(nx, ny, nz, 2, h, markers, Dz);
}

void Grid::get_gradient_operator() {
	gradient_op(nx, ny, nz, 0, h, markers, Gx);
	gradient_op(nx, ny, nz, 1, h, markers, Gy);
	gradient_op(nx, ny, nz, 2, h, markers, Gz);
}

// Directly get laplacian operator - matrix A
void Grid::get_laplacian_operator() {
	Eigen::Vector3d inv_h = h.unaryExpr([](double x) { return 1. / pow(x, 2); });

	typedef Eigen::Triplet<double> T;
	std::vector<T> trip;

	// TODO: Apply Ghost Pressure

	A.resize(n_grids, n_grids);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				// omit air cell and solid cell
				int index = get_idx(i, j, k);
				int index2;
				if (markers[index] == FLUIDCELL) {
					trip.push_back(T(index, index, -2. * inv_h.sum()));

					index2 = get_idx(i - 1, j, k);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(0)));
					}
					else
						trip.push_back(T(index, index, inv_h(0))); // if solid, the entry for this in gradient matrix D would be 0

					index2 = get_idx(i + 1, j, k);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(0)));
					}
					else
						trip.push_back(T(index, index, inv_h(0)));

					index2 = get_idx(i, j - 1, k);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(1)));
					}
					else
						trip.push_back(T(index, index, inv_h(1)));

					index2 = get_idx(i, j + 1, k);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(1)));
					}
					else
						trip.push_back(T(index, index, inv_h(1)));

					index2 = get_idx(i, j, k - 1);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(2)));
					}
					else
						trip.push_back(T(index, index, inv_h(2)));

					index2 = get_idx(i, j, k + 1);
					if (markers[index2] != SOLIDCELL) {
						trip.push_back(T(index, index2, inv_h(2)));
					}
					else
						trip.push_back(T(index, index, inv_h(2)));
				}
			}
		}
	}
	A.setFromTriplets(trip.begin(), trip.end());
	// check self-adjoint
	//if (!A.transpose().conjugate().isApprox(A))
	//	std::cout << "Warning: Matrix A not self-adjoint" << std::endl;
	//else
	//	std::cout << "Matrix A IS self-adjoint, can switch ConjugateGradient to solve A instead" << std::endl;
}


int Grid::get_idx(const int& xi, const int& yi, const int& zi) {
	return xi * ny * nz + yi * nz + zi;
}


void Grid::pressure_projection() {
	get_gradient_operator();
	get_divergence();
	get_laplacian_operator();
	solve_pressure();
	update_velocity();
}


// Directly compute divergence
void Grid::get_divergence() {
	const auto& get_idx2 = [&](int xi, int yi, int zi, int dim) { 
		int dim0 = (dim == 0) ? nx + 1 : nx;
		int dim1 = (dim == 1) ? ny + 1 : ny;
		int dim2 = (dim == 2) ? nz + 1 : nz; 
		return xi * dim1 * dim2 + yi * dim2 + zi; 
	};
	divergence.resize(n_grids);
	divergence.setZero();
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nz; ++k) {
				int idx = get_idx(i, j, k);

				if (markers[idx] == FLUIDCELL) {
					divergence[idx] =
						(Vx[get_idx2(i + 1, j, k, 0)] - Vx[get_idx2(i, j, k, 0)]) / h(0) +
						(Vy[get_idx2(i, j + 1, k, 1)] - Vy[get_idx2(i, j, k, 1)]) / h(1) +
						(Vz[get_idx2(i, j, k + 1, 2)] - Vz[get_idx2(i, j, k, 2)]) / h(2);

				}
			}
		}
	}
	divergence = (density / dt) * divergence;
}


// Solve pressure by Conjugate Gradient Method
void Grid::solve_pressure() {
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-5);

	// not sure if A is self-adjoint - use AT*A instead
	cg.compute(A.transpose() * A);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver decomposition failed, given matrix is not self-adjoint" << std::endl;
	pressure = cg.solve(A.transpose() * divergence);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver solving failed. However decomposition seems work" << std::endl;
}


void Grid::update_velocity() {
	Vx -= dt / density * Gx * pressure;
	Vy -= dt / density * Gy * pressure;
	Vz -= dt / density * Gz * pressure;
}


void Grid::save_grids() {
	Vx_ = Vx;
	Vy_ = Vy;
	Vz_ = Vz;
}



void Grid::distance_to_fluid(Particle &particles) {
	init_phi(particles);
	sweep_phi(particles);
	sweep_phi(particles);
}


/* Signed Distance Stuff */

void Grid::init_phi(Particle &particles) {
	int n_grids = nx * ny * nz;
	phi = Eigen::VectorXd::Constant(n_grids, std::numeric_limits<double>::max());
	cp = Eigen::VectorXi::Constant(n_grids, -1);

	const auto& ks = [](double x) -> double { return std::max(0., pow(1 - x * x, 3)); };

	// initialize the arrays near the near geometry
	int i, j, k;
	Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n_grids, 3);
	Eigen::VectorXd w = Eigen::VectorXd::Zero(n_grids);
	Eigen::VectorXd r = Eigen::VectorXd::Zero(n_grids);

	Eigen::RowVector3d x;
	Eigen::RowVector3d ids;
	const auto& add_stuff = [&](const int& xi, const int& yi, const int& zi) {
		if (xi < 1 || xi > nx - 2)
			return;
		if (yi < 1 || yi > ny - 2)
			return;
		if (zi < 1 || zi > nz - 2)
			return;
		int grid_id = get_idx(xi, yi, zi);
		Eigen::RowVector3d p = ids.cwiseProduct(h) + 0.5 * h; // pressure point coordinates
		double wi = ks(0.5 * (x - p).norm() / h.maxCoeff());
		X.row(grid_id) += wi * x;
		r(grid_id) += wi * 0.5 * h.maxCoeff(); // should modify this
		w(grid_id) += wi;
	};

	for (int l = 0; l < particles.q.rows(); l++) {
		x = particles.q.row(l);
		ids = (x - left_lower_corner).cwiseQuotient(h);
		ids = ids.unaryExpr([](double x) { return floor(x); });

		for (int i0 = -1; i0 <= 1; i0++) {
			for (int j0 = -1; j0 <= 1; j0++) {
				for (int k0 = -1; k0 <= 1; k0++) {
					add_stuff(int(ids(0))+i0, int(ids(1))+j0, int(ids(2))+k0);
				}
			}
		}
	}

	for (int i = 1; i < nx-1; i++) {
		for (int j = 1; j < ny-1; j++) {
			for (int k = 1; k < nz-1; k++) {
				int grid_id = get_idx(i, j, k);
				if (w(grid_id) != 0) {
					Eigen::RowVector3d p = Eigen::RowVector3d(i+0.5, j+0.5, k+0.5).cwiseProduct(h);
					phi(grid_id) = (p - X.row(grid_id)).norm() - r(grid_id);

					if (phi(grid_id) < 0) std::cerr << "inside\n";
					if (phi(grid_id) > 0) std::cerr << "outside\n";
				}
			}
		}
	}
}



/*
void Grid::signed_distance(const Eigen::RowVector3d &v1, const Eigen::RowVector3d &pt, double& d) {
	Eigen::RowVector3d dist = (pt - v1).cwiseQuotient(h);
	if (dist.minCoeff() > 0 && dist.maxCoeff() < 1) {
		// inside the cell
		double d_1 = (pt - v1).maxCoeff();
		double d_2 = (v1 + h - pt).maxCoeff();
		d = std::max(d_1, d_2);
	}
	else {
		// outside the cell
		Eigen::RowVector3d x = pt;
		x = x.cwiseMax(v1);
		x = x.cwiseMin(v1 + h);
		d = (x - pt).norm();
	}
}

void Grid::sweep_phi(Particle& particles) {
	
	const auto& solve_distance = [&](const Eigen::RowVector3i &dir, const int& xi, const int& yi, const int& zi) {
		int idx = get_idx(xi, yi, zi);
		int idx0 = get_idx(xi + dir(0), yi, zi);
		int idx1 = get_idx(xi, yi + dir(1), zi);
		int idx2 = get_idx(xi, yi, zi + dir(2));
		
		if (markers(idx) != FLUIDCELL) {
			Eigen::RowVector3d grid_pt = Eigen::RowVector3d(xi, yi, zi).cwiseProduct(h);

			double d;
			if (cp(idx0) != -1) {
				signed_distance(grid_pt, particles.q.row(cp(idx0)), d);
				if (d < phi(idx)) {
					phi(idx) = d;
					cp(idx) = cp(idx0);
				}
			}
			if (cp(idx1) != -1) {
				signed_distance(grid_pt, particles.q.row(cp(idx1)), d);
				if (d < phi(idx)) {
					phi(idx) = d;
					cp(idx) = cp(idx1);
				}
			}
			if (cp(idx2) != -1) {
				signed_distance(grid_pt, particles.q.row(cp(idx2)), d);
				if (d < phi(idx)) {
					phi(idx) = d;
					cp(idx) = cp(idx2);
				}
			}
		}
	};

	// sweep in all eight directions
	for (int i = 1; i < nx; i++)
		for (int j = 1; j < ny; j++)
			for (int k = 1; k < nz; k++)
				solve_distance(Eigen::RowVector3i(-1, -1, -1), i, j, k);
	for (int i = 1; i < nx; i++)
		for (int j = 1; j < ny; j++)
			for (int k = nz-2; k >= 0; k--)
				solve_distance(Eigen::RowVector3i(-1, -1, 1), i, j, k);
	for (int i = 1; i < nx; i++)
		for (int j = ny-2; j >= 0; j--)
			for (int k = 1; k < nz; k++)
				solve_distance(Eigen::RowVector3i(-1, 1, -1), i, j, k);
	for (int i = 1; i < nx; i++)
		for (int j = ny - 2; j >= 0; j--)
			for (int k = nz - 2; k >= 0; k--)
				solve_distance(Eigen::RowVector3i(-1, 1, 1), i, j, k);
	for (int i = nx-2; i >= 0; i--)
		for (int j = 1; j < ny; j++)
			for (int k = 1; k < nz; k++)
				solve_distance(Eigen::RowVector3i(1, -1, -1), i, j, k);
	for (int i = nx - 2; i >= 0; i--)
		for (int j = 1; j < ny; j++)
			for (int k = nz - 2; k >= 0; k--)
				solve_distance(Eigen::RowVector3i(1, -1, 1), i, j, k);
	for (int i = nx - 2; i >= 0; i--)
		for (int j = ny - 2; j >= 0; j--)
			for (int k = 1; k < nz; k++)
				solve_distance(Eigen::RowVector3i(1, 1, -1), i, j, k);
	for (int i = nx - 2; i >= 0; i--)
		for (int j = ny - 2; j >= 0; j--)
			for (int k = nz - 2; k >= 0; k--)
				solve_distance(Eigen::RowVector3i(1, 1, 1), i, j, k);
}

void Grid::sweep_velocity() {
	int dim0, dim1, dim2;
	const auto& get_idx2 = [&](const int& xi, const int& yi, const int& zi) {
		return xi * dim1 * dim2 + yi * dim2 + zi;
	};
	
	// first for x
	sweep_dir(0);
	dim0 = nx + 1; dim1 = ny; dim2 = nz;
	for (int i = 0; i < dim0; i++) {
		for (int j = 0; j < dim1; j++) {
			Vx(get_idx2(i, j, 0)) = Vx(get_idx2(i, j, 1));
			Vx(get_idx2(i, j, dim2 - 1)) = Vx(get_idx2(i, j, dim2 - 2));
		}
		for (int k = 0; k < dim2; k++) {
			Vx(get_idx2(i, 0, k)) = Vx(get_idx2(i, 1, k));
			Vx(get_idx2(i, dim1 - 1, k)) = Vx(get_idx2(i, dim1 - 2, k));
		}
	}
	for (int j = 0; j < dim1; j++) {
		for (int k = 0; k < dim2; k++) {
			Vx(get_idx2(0, j, k)) = Vx(get_idx2(1, j, k));
			Vx(get_idx2(dim0, j, k)) = Vx(get_idx2(dim0 - 2, j, k));
		}
	}

	// now for y
	sweep_dir(1);
	dim0 = nx; dim1 = ny + 1; dim2 = nz;
	for (int i = 0; i < dim0; i++) {
		for (int j = 0; j < dim1; j++) {
			Vy(get_idx2(i, j, 0)) = Vy(get_idx2(i, j, 1));
			Vy(get_idx2(i, j, dim2 - 1)) = Vy(get_idx2(i, j, dim2 - 2));
		}
		for (int k = 0; k < dim2; k++) {
			Vy(get_idx2(i, 0, k)) = Vy(get_idx2(i, 1, k));
			Vy(get_idx2(i, dim1 - 1, k)) = Vy(get_idx2(i, dim1 - 2, k));
		}
	}
	for (int j = 0; j < dim1; j++) {
		for (int k = 0; k < dim2; k++) {
			Vy(get_idx2(0, j, k)) = Vy(get_idx2(1, j, k));
			Vy(get_idx2(dim0, j, k)) = Vy(get_idx2(dim0 - 2, j, k));
		}
	}

	// now for z
	sweep_dir(2);
	dim0 = nx; dim1 = ny; dim2 = nz + 1;
	for (int i = 0; i < dim0; i++) {
		for (int j = 0; j < dim1; j++) {
			Vz(get_idx2(i, j, 0)) = Vz(get_idx2(i, j, 1));
			Vz(get_idx2(i, j, dim2 - 1)) = Vz(get_idx2(i, j, dim2 - 2));
		}
		for (int k = 0; k < dim2; k++) {
			Vz(get_idx2(i, 0, k)) = Vz(get_idx2(i, 1, k));
			Vz(get_idx2(i, dim1 - 1, k)) = Vz(get_idx2(i, dim1 - 2, k));
		}
	}
	for (int j = 0; j < dim1; j++) {
		for (int k = 0; k < dim2; k++) {
			Vz(get_idx2(0, j, k)) = Vz(get_idx2(1, j, k));
			Vz(get_idx2(dim0, j, k)) = Vz(get_idx2(dim0 - 2, j, k));
		}
	}
}


void Grid::sweep_dir(const int& dir) {


	const auto& sweep_x = [&](int i0, int i1, int j0, int j1, int k0, int k1) {
		int di = (i0 < i1) ? 1 : -1;
		int dj = (j0 < j1) ? 1 : -1;
		int dk = (k0 < k1) ? 1 : -1;

		int grid_id;
		for (int i = i0; i != i1; i += di) {
			for (int j = j0; j != j1; j += dj) {
				for (int k = k0; k != k1; k += dk) {
					if (markers(get_idx(i - 1, j, k)) == AIRCELL && markers(get_idx(i, j, k) == AIRCELL)) {
						double dq = di * (phi(get_idx(i, j, k)) - phi(get_idx(i - 1, j, k)));
						if (dq < 0) continue;

					}

				}
			}
		}
	};

	int dim0 = (dir == 0) ? nx + 1 : nx;
	int dim1 = (dir == 1) ? ny + 1 : ny;
	int dim2 = (dir == 2) ? nz + 1 : nz;

	// sweep from all 8 directions
}
*/









