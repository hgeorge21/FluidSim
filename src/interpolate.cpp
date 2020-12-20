#include <interpolate.h>
#include <omp.h>
#include <iostream>

#include <trilinear_interpolation.h>

void interpolate(
    const int &nx,
    const int &ny,
    const int &nz,
    const int &dim,
    const Eigen::RowVector3d& h,
    const Eigen::RowVector3d& corner,
    const Eigen::MatrixXd& q,
    Eigen::VectorXd& sum,
    Eigen::SparseMatrix<double>& W) {

    int tmp = (dim == 0) ? nx : (dim == 1) ? ny : nz;
    int dim0 = (dim == 0) ? nx + 1 : nx;
    int dim1 = (dim == 1) ? ny + 1 : ny;
    int dim2 = (dim == 2) ? nz + 1 : nz;

    int n_particles = q.rows();
    int i, j, k;
    double xd, yd, zd;
    Eigen::VectorXd w;

    W.resize(n_particles, dim0*dim1*dim2);
    W.reserve(8 * n_particles);
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> trip;


    const auto& bary = [](const double x, const double h, int& i, double& d) {
        double s = x / h;
        i = (int)s;
        d = s - floor(s);
    };

    const auto& bary_center = [&](const int dim, const double x, const double h, int& i, double& d) {
        int n = (dim == 0) ? nx : (dim == 1) ? ny : nz;
        
        double s = x / h - 0.5;
        i = (int)s;
        if (i < 0) {
            i = 0;
            d = 0.0;
        }
        else if (i > n - 2) {
            i = n - 2;
            d = 1.0;
        }
        else {
            d = s - floor(s);
        }
    };



    for (int l = 0; l < n_particles; l++) {
        Eigen::RowVector3d d;
        Eigen::RowVector3i ids;
        Eigen::RowVector3d x = q.row(l)-corner;
        
        int idx1 = (dim + 1) % 3;
        int idx2 = (dim + 2) % 3;
        bary(x(dim), h(dim), ids(dim), d(dim));
        bary_center(idx1, x(idx1), h(idx1), ids(idx1), d(idx1));
        bary_center(idx2, x(idx2), h(idx2), ids(idx2), d(idx2));

        i = int(ids(0)); j = int(ids(1)); k = int(ids(2));
        int ind = i * dim1 * dim2 + j * dim2 + k;

        trilinear_interpolation(d, w);

        trip.emplace_back(T(l, ind, w(0)));
        trip.emplace_back(T(l, ind + dim1 * dim2, w(1)));
        trip.emplace_back(T(l, ind + dim2, w(2)));
        trip.emplace_back(T(l, ind + dim1 * dim2 + dim2, w(3)));
        trip.emplace_back(T(l, ind + 1, w(4)));
        trip.emplace_back(T(l, ind + dim1 * dim2 + 1, w(5)));
        trip.emplace_back(T(l, ind + dim2 + 1, w(6)));
        trip.emplace_back(T(l, ind + dim1*dim2 + dim2 + 1, w(7)));
    
        if (ind + dim1 * dim2 + dim2 + 1 >= dim0 * dim1 * dim2) {
            std::cerr << "Dimension " << dim << "\n";
            std::cerr << x << "\n";
            std::cerr << d << "\n";
            std::cerr << ids << "\n";
            std::cerr << i << " " << j << " " << k << "\n";
            std::cerr << dim0 << " " << dim1 << " " << dim2 << "\n";
            std::cerr << "----------\n";
        }


    }

    W.setFromTriplets(trip.begin(), trip.end());
  
    sum = Eigen::VectorXd::Ones(n_particles);
    sum = W.transpose() * sum;
    sum = sum.unaryExpr([](double x) { return (x == 0.) ? 1.0 : x; });
}