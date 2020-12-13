#include <interpolate.h>

#include <iostream>

void interpolate(
    const int nx,
    const int ny,
    const int nz,
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
        Eigen::RowVector3d p = q.row(l) - corner;
        xd = p(0) / h;
        yd = p(1) / h;
        zd = p(2) / h;

        i = floor(xd);
        j = floor(yd);
        k = floor(zd);

        xd = xd - i;
        yd = yd - j;
        zd = zd - k;

        c_000 = (1 - xd) * (1 - yd) * (1 - zd);
        c_001 = xd * (1 - yd) * (1 - zd);
        c_010 = (1 - xd) * yd * (1 - zd);
        c_011 = xd * yd * (1 - zd);
        c_100 = (1 - xd) * (1 - yd) * zd;
        c_101 = xd * (1 - yd) * zd;
        c_110 = (1 - xd) * yd * zd;
        c_111 = xd * yd * zd;

        int ind = i + nx * (j + k * ny);
        trip.emplace_back(T(l, ind, c_000));
        trip.emplace_back(T(l, ind + 1, c_001));
        trip.emplace_back(T(l, ind + nx, c_010));
        trip.emplace_back(T(l, ind + nx + 1, c_011));
        trip.emplace_back(T(l, ind + nx * ny, c_100));
        trip.emplace_back(T(l, ind + nx * ny + 1, c_101));
        trip.emplace_back(T(l, ind + nx + nx * ny, c_110));
        trip.emplace_back(T(l, ind + nx + nx * ny + 1, c_111));
    }

    W.resize(n_particles, nx * ny * nz);
    std::cerr << W.cols() << "\n";
    std::cerr << trip.size();
    W.reserve(8 * n_particles);
    W.setFromTriplets(trip.begin(), trip.end());
}