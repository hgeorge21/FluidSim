#ifndef GRID_TO_PARTICLE_H
#define GRID_TO_PARTICLE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <particles.h>
#include <grid.h>

// Transfers the velocity from grid to particles
// Apply either PIC, FLIP, or blend of PIC/FLIP
void grid_to_particle_velocity_update(
	Grid& grid, 
	Particle& particles);

#endif // GRID_TO_PARTICLE_H