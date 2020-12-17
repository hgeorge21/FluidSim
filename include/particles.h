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

	Eigen::MatrixXd q;
	Eigen::MatrixXd v;
	Eigen::VectorXi type;

	Particle();
	~Particle();

private:
	
};




#endif // !PARTICLE_H
