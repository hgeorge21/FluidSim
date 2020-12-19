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


void Grid::add_fluid(Particle& particles, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &hf, const double& height) {
	// calculate for fluids
	Eigen::Vector3d ns = (v2 - v1).cwiseQuotient(hf);
	int nx_ = ceil(ns(0));
	int ny_ = ceil(ns(1));
	int nz_ = ceil(ns(2));

	int n_per_cell = 8; // each cell has 8 particles
	
	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	// the first & last of each dimension is solid only
	int fluid_height = ceil(height / hf(1));
	int ny_f = ny_ - 2;
	int nx_f = nx_ - 2;
	int nz_f = nz_ - 2;

	int n_fluid_grids = nx_f * ny_f * nz_f;
	particles.q = Eigen::MatrixXd(n_per_cell * n_fluid_grids, 3);
	particles.v = Eigen::MatrixXd::Zero(n_per_cell * n_fluid_grids, 3);
	particles.type = Eigen::VectorXi::Zero(n_per_cell * n_fluid_grids);
	for (int i = 0; i < nx_f; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz_f; k++) {
				int idx = i * (ny_f * nz_f) + j * nz_f + k;

				// generate fluid particles
				Eigen::RowVector3d lower_corner = v1 + 2*h + hf.cwiseProduct(Eigen::Vector3d(i, j, k));
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(n_per_cell, 3, [&]() { return dist(generator);  });
				q_.col(0) = h(0) * q_.col(0);
				q_.col(1) = h(1) * q_.col(1);
				q_.col(2) = h(2) * q_.col(2);
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


void Grid::get_divergence_operator() {
	divergence_op(nx, ny, nz, 0, h, markers, Dx);
	divergence_op(nx, ny, nz, 1, h, markers, Dy);
	divergence_op(nx, ny, nz, 2, h, markers, Dz);
}


int Grid::get_idx(const int& xi, const int& yi, const int& zi) {
	return xi * ny * nz + yi * nz + zi;
}


void Grid::pressure_projection() {
	//print_cell();

	std::cout << "Now get gradient op" << std::endl;
	get_gradient_operator();

	std::cout << "Now get divergence op" << std::endl;
	get_divergence_operator();

	//get_divergence();
	std::cout << "Now get divergence2" << std::endl;
	get_divergence2();

	get_laplacian_operator();

	solve_pressure();
	update_velocity();
	//print_cell(); 
}


// Get divergence of v with operator
void Grid::get_divergence() {
	divergence = Dx * Vx + Dy * Vy + Dz * Vz;
	divergence = (density / dt) * divergence;
	std::cout << "Divergence norm: " << divergence.norm() << std::endl;
}


// Directly compute divergence
void Grid::get_divergence2() {
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

				//std::cout << "Now markers" << std::endl;

				if (markers[idx] == FLUIDCELL)

				//std::cout << "Now divergence and Vxyz" << std::endl;

					divergence[idx] =
					(Vx[get_idx2(i + 1, j, k, 0)] - Vx[get_idx2(i, j, k, 0)]) / h(0) +
					(Vy[get_idx2(i, j + 1, k, 1)] - Vy[get_idx2(i, j, k, 1)]) / h(1) +
					(Vz[get_idx2(i, j, k + 1, 2)] - Vz[get_idx2(i, j, k, 2)]) / h(2);
			}
		}
	}
	divergence = (density / dt) * divergence;
	std::cout << "Divergence norm: " << divergence.norm() << std::endl;

	check_divergence();
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
					if (markers[index2] != SOLIDCELL){
						trip.push_back(T(index, index2, inv_h(0)));
					}
					else
						trip.push_back(T(index, index, inv_h(0))); // if solid, the entry for this in gradient matrix D would be 0

					index2 = get_idx(i + 1, j, k);
					if (markers[index2] != SOLIDCELL){
						trip.push_back(T(index, index2, inv_h(0)));
					}
					else
						trip.push_back(T(index, index, inv_h(0))); 

					index2 = get_idx(i, j - 1, k);
					if (markers[index2] != SOLIDCELL){
						trip.push_back(T(index, index2, inv_h(1)));
					}
					else
						trip.push_back(T(index, index, inv_h(1)));

					index2 = get_idx(i, j + 1, k);
					if (markers[index2] != SOLIDCELL){
						trip.push_back(T(index, index2, inv_h(1)));
					}
					else
						trip.push_back(T(index, index, inv_h(1)));

					index2 = get_idx(i, j, k - 1);
					if (markers[index2] != SOLIDCELL){
						trip.push_back(T(index, index2, inv_h(2)));
					}
					else
						trip.push_back(T(index, index, inv_h(2)));

					index2 = get_idx(i, j, k + 1);
					if (markers[index2] != SOLIDCELL){
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
	if (!A.transpose().conjugate().isApprox(A))
		std::cout << "Warning: Matrix A not self-adjoint" << std::endl;
	else 
		std::cout << "Matrix A IS self-adjoint, can switch ConjugateGradient to solve A instead" << std::endl;

	check_laplacian();
}

// Get A = B * D
void Grid::get_laplacian_operator2() {
	get_gradient_operator();
	Eigen::MatrixXd B, Bx_, By_, Bz_, D, Dx_, Dy_, Dz_;

	Bx_ = Eigen::MatrixXd(Dx); // Note Dx is divergence instead of gradient matrix
	By_ = Eigen::MatrixXd(Dy);
	Bz_ = Eigen::MatrixXd(Dz);
	B.resize(n_grids, Dx.cols() + Dy.cols() + Dz.cols());
	B << Bx_, By_, Bz_;

	Dx_ = Eigen::MatrixXd(Gx);
	Dy_ = Eigen::MatrixXd(Gy);
	Dz_ = Eigen::MatrixXd(Gz);
	D.resize(Gx.rows() + Gy.rows() + Gz.rows(), n_grids);
	D << Dx_, Dy_, Dz_;

	A = (B * D).sparseView();
}

// Solve pressure by Conjugate Gradient Method
void Grid::solve_pressure() {
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-8);

	// not sure if A is self-adjoint - use AT*A instead
	cg.compute(A.transpose() * A);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver decomposition failed, given matrix is not self-adjoint" << std::endl;
	pressure = cg.solve(A.transpose() * divergence);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver solving failed. However decomposition seems work" << std::endl;
	print_pressure();
}


void Grid::update_velocity() {
	Vx -= dt / density * Gx * pressure;
	Vy -= dt / density * Gy * pressure;
	Vz -= dt / density * Gz * pressure;
}


// All the checking and printing
//
//

void Grid::print_pressure() {
	
	for (int level = ny - 1; level >= 0; level--) {
		double pressure_per_lv = 0.0;
		for (int x = 1; x < nx - 1; x++) {
			for (int z = 1; z < nz - 1; z++){
				int idx = get_idx(x, level, z);
				pressure_per_lv += pressure[idx];
			}
		}
		std::cout << "Avg pressure at level " << level << ": " << pressure_per_lv/(double)((nx-2)*(nz-2)) << std::endl;
	}
}

void Grid::print_cell() {
	for (int y = ny - 1; y >= 0 ; y--) {
		std::cout << "Level " << y << "-------------------------------\n";
		for (int x = 0; x < nx; x++) {
			for (int z = 0; z < nz; z++) {
				int idx = get_idx(x, y, z);
				char type = (markers[idx] == FLUIDCELL) ? 'F' : (markers[idx] == SOLIDCELL) ? 'S' : 'A';
				std::cout << type << ' ';
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}

// Assume to be called inside get_divergence2
void Grid::check_divergence() {
	Eigen::VectorXd divergence_ = Dx * Vx + Dy * Vy + Dz * Vz;
	divergence_ = (density / dt) * divergence_;
	if (!divergence_.isApprox(divergence)) // initially they are both 0, so not the same
		std::cerr << "Divergence is not the same" << std::endl;
	else
		std::cerr << "Divergence IS the same" << std::endl;
}

// check if gradient matrix A = B * D
void Grid::check_laplacian() {
	Eigen::MatrixXd B, Bx_, By_, Bz_, D, Dx_, Dy_, Dz_;

	Bx_ = Eigen::MatrixXd(Dx); // Note Dx is divergence instead of gradient matrix
	By_ = Eigen::MatrixXd(Dy);
	Bz_ = Eigen::MatrixXd(Dz);
	B.resize(n_grids, Dx.cols() + Dy.cols() + Dz.cols());
	B << Bx_ , By_ , Bz_;

	Dx_ = Eigen::MatrixXd(Gx);
	Dy_ = Eigen::MatrixXd(Gy);
	Dz_ = Eigen::MatrixXd(Gz);
	D.resize(Gx.rows() + Gy.rows() + Gz.rows(), n_grids);
	D << Dx_, Dy_, Dz_;

	if (!(B * D).sparseView().isApprox(A))
		std::cerr << "A != B * D" << std::endl;
	else
		std::cerr << "A = B * D" << std::endl;
}


void Grid::save_grids() {
	Vx_ = Vx;
	Vy_ = Vy;
	Vz_ = Vz;
}
