#ifndef ADVECT_V_H
#define ADVECT_V_H

#include <Eigen/Core>

#include <grid.h>
#include <particles.h>

// Move the particles based on velocity
// Uses either one-step Euler or Runge-Kutta 2 (RK2)
void advect_velocity(
	Grid& grid, 
	Particle& particles, 
	double dt);

#endif //ADVECT_V_H