#ifndef GRAVITY_H
#define GRAVITY_H

#include <Eigen/Core>

#include <grid.h>
#include <particles.h>

void add_gravity(Grid& grid, Particle& particles, const Eigen::Vector3d& gravity, const double& dt);

#endif //GRAVITY_H