#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>


class Particle
{

#define AIR_P   0
#define FLUID_P 1

public:
	Eigen::MatrixXd q; // location
	Eigen::MatrixXd v; // velocity
	Eigen::VectorXi type; // type of particle (only fluids right now)

	Particle() {};
	~Particle() {};

private:

};

#endif // !PARTICLE_H
