#ifndef ADD_PARTICLES_H
#define ADD_PARTICLES_H

#include <Eigen/Core>

#include <grid.h>
#include <particles.h>

// Adds the particles in the specified domain
void add_particles(
	Grid& grid, 
	Particle& particles, 
	const Eigen::RowVector3d &v1, 
	const Eigen::RowVector3d &v2, 
	const Eigen::RowVector3d &hf, 
	const int &n_per_cell);

// Adds the particles in the specified domain along with an input mesh
void add_particles(
	Grid& grid, 
	Particle& particles, 
	const Eigen::RowVector3d& v1, 
	const Eigen::RowVector3d& v2, 
	const Eigen::RowVector3d& hf, 
	const int& n_per_cell,
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F);

#endif //ADD_PARTICLES_H