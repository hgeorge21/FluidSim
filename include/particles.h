#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>

enum advection_method
{
	SIMPLE, RUNGE_KUTTA
};

class Particle
{

public:
	advection_method method;

	Eigen::MatrixXd q; // location
	Eigen::MatrixXd v; // velocity
	Eigen::VectorXi type;

	Particle();
	~Particle();

private:
	
};




#endif // !PARTICLE_H
