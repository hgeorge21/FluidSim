#include <interpolate.h>

#include <iostream>

void interpolate(
    const int nx,
    const int ny,
    const int nz,
    const int dim,
    const double h,
    const Eigen::RowVector3d& corner,
    const Eigen::MatrixXd& q,
    Eigen::SparseMatrix<double>& W) {

    int n_particles = q.rows();
    int i, j, k;
    double xd, yd, zd;
    double c_000, c_001, c_010, c_011, c_100, c_101, c_110, c_111;

    typedef Eigen::Triplet<double> T;
    std::vector<T> trip;
    for (int l = 0; l < n_particles; l++) {
        Eigen::RowVector3d x = q.row(l) - corner;
        Eigen::RowVector3d d = x / h;
        Eigen::RowVector3d ids = d.unaryExpr([](double x){ return floor(std::max(x, 0.)); });
        d = d - ids;

        xd = d(0); yd = d(1); zd = d(2);
        i = int(ids(0)); j = int(ids(1)); k = int(ids(2));

        c_000 = (1 - xd) * (1 - yd) * (1 - zd);
        c_001 = xd * (1 - yd) * (1 - zd);
        c_010 = (1 - xd) * yd * (1 - zd);
        c_011 = xd * yd * (1 - zd);
        c_100 = (1 - xd) * (1 - yd) * zd;
        c_101 = xd * (1 - yd) * zd;
        c_110 = (1 - xd) * yd * zd;
        c_111 = xd * yd * zd;

        int ind = i * (ny+1) * (1+nz) + j * (1+nz) + k;

        trip.emplace_back(T(l, ind, c_000));
        trip.emplace_back(T(l, ind + (ny+1)*(nz+1), c_001));
        trip.emplace_back(T(l, ind + nz+1, c_010));
        trip.emplace_back(T(l, ind + (ny+1)*(nz+1) + (nz+1), c_011));
        trip.emplace_back(T(l, ind + 1, c_100));
        trip.emplace_back(T(l, ind + (nz+1) * (ny+1) + 1, c_101));
        trip.emplace_back(T(l, ind + (nz+1) + 1, c_110));
        trip.emplace_back(T(l, ind + (nz+1) + (nz+1) * (ny+1) + 1, c_111));
    }

    W.resize(n_particles, (nx+1)*(ny+1)*(nz+1));
    W.reserve(8 * n_particles);
    W.setFromTriplets(trip.begin(), trip.end());
}