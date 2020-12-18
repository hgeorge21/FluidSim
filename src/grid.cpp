#include <grid.h>
#include <chrono>
#include <random>
#include <iostream>
#include <divergence.h>


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

	const auto& set_sparse_vec = [&](const int& dim) {
		int ll0 = (dim == 0) ? 1 : 0; // 1 or 2??
		int ll1 = (dim == 1) ? 1 : 0;
		int ll2 = (dim == 2) ? 1 : 0;
		int ul0 = (dim == 0) ? nx : nx; // nx-1 or nx ??
		int ul1 = (dim == 1) ? ny : ny;
		int ul2 = (dim == 2) ? nz : nz;
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


void Grid::add_fluid(Particle& particles, const double& height) {
	// set random seed for particle generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> dist(0., 1.);

	// the first & last of each dimension is solid only
	int fluid_height = ceil(height / h(1));
	int ny_f = ny - 2;
	int nx_f = nx - 2;
	int nz_f = nz - 2;

	int n_fluid_grids = nx_f * ny_f * nz_f;
	particles.q = Eigen::MatrixXd(8 * n_fluid_grids, 3); // each cell has 8 particles
	particles.v = Eigen::MatrixXd::Zero(8 * n_fluid_grids, 3);
	particles.type = Eigen::VectorXi::Zero(8 * n_fluid_grids);
	for (int i = 0; i < nx_f; i++) {
		for (int j = 0; j < ny_f; j++) {
			for (int k = 0; k < nz_f; k++) {
				int idx = i * (ny_f * nz_f) + j * nz_f + k;

				// claim the fluid cells
				int org_idx = (i + 1) * ny * nz + (j + 1) * nz + k + 1;
				markers(org_idx) = FLUIDCELL; // (j < fluid_height) ? FLUIDCELL : AIRCELL;

				// generate fluid particles
				Eigen::RowVector3d lower_corner = left_lower_corner + h + h.cwiseProduct(Eigen::Vector3d(i, j, k));
				Eigen::MatrixXd q_ = Eigen::MatrixXd::NullaryExpr(8, 3, [&]() { return dist(generator);  });
				q_.col(0) = h(0) * q_.col(0);
				q_.col(1) = h(1) * q_.col(1);
				q_.col(2) = h(2) * q_.col(2);

				q_ = q_.rowwise() + lower_corner;
				particles.q.block(8 * idx, 0, 8, 3) = q_;

				//if (j < fluid_height) {
					particles.type.segment<8>(8 * idx) = FLUID_P * Eigen::VectorXi::Ones(8);
				//}
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

			//id = i * ny * nz + j * nz + 0;
			//Vz(id) = 0;
			//id = i * ny * nz + j * nz + 1;
			//Vz(id) = 0;
			//id = i * ny * nz + j * nz + nz;
			//Vz(id) = 0;
			//id = i * ny * nz + j * nz + (nz-1);
			//Vz(id) = 0;
		}
		for (k = 0; k < nz; k++) {
			int j1 = 0; int j2 = ny - 1;
			markers(get_idx(i, j1, k)) = SOLIDCELL;
			markers(get_idx(i, j2, k)) = SOLIDCELL;

			//id = i * ny * nz + 0 * nz + k;
			//Vy(id) = 0;
			//id = i * ny * nz + 1 * nz + k;
			//Vy(id) = 0;
			//id = i * ny * nz + ny * nz + k;
			//Vy(id) = 0;
			//id = i * ny * nz + (ny-1) * nz + k;
			//Vy(id) = 0;
		}
	}
	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			int i1 = 0; int i2 = nx - 1;
			markers(get_idx(i1, j, k)) = SOLIDCELL;
			markers(get_idx(i2, j, k)) = SOLIDCELL;

			//id = 0 * ny * nz + j * nz + k;
			//Vx(id) = 0;
			//id = 1 * ny * nz + j * nz + k;
			//Vx(id) = 0;
			//id = nx * ny * nz + j * nz + k;
			//Vx(id) = 0;
			//id = (nx-1) * ny * nz + j * nz + k;
			//Vx(id) = 0;
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
	get_divergence_operator();
	get_laplacian_operator();

	get_divergence();
	solve_pressure();
	update_velocity();
}


// Get divergence of v
void Grid::get_divergence() {
	divergence = Dx * Vx + Dy * Vy + Dz * Vz;
	divergence = divergence;
}


// Get laplacian operator - matrix A
void Grid::get_laplacian_operator() {

	Eigen::Vector3d inv_h;
	inv_h << 1.0 / pow(h(0), 2), 1.0 / pow(h(1), 2), 1.0 / pow(h(2), 2);

	typedef Eigen::Triplet<double> T;
	std::vector<T> trip;

	// TODO: Apply Ghost Pressure

	A.resize(n_grids, n_grids);
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 1; k < nz - 1; k++) {
				// omit air cell and solid cell
				int index = get_idx(i, j, k);
				int index2;
				if (markers[index] == FLUIDCELL) {
					trip.push_back(T(index, index, -2. * inv_h.sum()));

					index2 = get_idx(i - 1, j, k);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(0)));

					index2 = get_idx(i + 1, j, k);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(0)));

					index2 = get_idx(i, j - 1, k);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(1)));

					index2 = get_idx(i, j + 1, k);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(1)));

					index2 = get_idx(i, j, k - 1);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(2)));

					index2 = get_idx(i, j, k + 1);
					if (markers[index2] != SOLIDCELL)
						trip.push_back(T(index, index2, inv_h(2)));
				}
			}
		}
	}
	A.setFromTriplets(trip.begin(), trip.end());
	// check self-adjoint
	// if (!A.transpose().conjugate().isApprox(A))
	//	std::cout << "Warning: Matrix A not self-adjoint" << std::endl;
	// else 
	//	std::cout << "Matrix A IS self-adjoint, can switch ConjugateGradient to solve A instead" << std::endl;
}


// Solve pressure by Conjugate Gradient Method
void Grid::solve_pressure() {
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-8);

	// not sure if A is self-adjoint - use AT*A instead
	cg.compute(A);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver decomposition failed, given matrix is not self-adjoint" << std::endl;
	pressure = cg.solve(divergence);
	if (cg.info() != Eigen::Success)
		std::cerr << "Warning: Conjugate Gradient Solver solving failed. However decomposition seems work" << std::endl;
}


void Grid::update_velocity() {
	get_gradient_operator();
	// TODO: is this plus or minus?
	Vx -= Gx * pressure;
	Vy -= Gy * pressure;
	Vz -= Gz * pressure;
}


void Grid::get_gradient_operator() {
	gradient_op(nx, ny, nz, 0, h, markers, Gx);
	gradient_op(nx, ny, nz, 1, h, markers, Gy);
	gradient_op(nx, ny, nz, 2, h, markers, Gz);
}


void Grid::save_grids() {
	Vx_ = Vx;
	Vy_ = Vy;
	Vz_ = Vz;
}
