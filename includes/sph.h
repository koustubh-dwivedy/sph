#include "sph.cpp"

class environment
{
private:
	//// Environment Parameters
	float gravity = 9.81;
	float time_step = 0.01;
	float temperature =  283.15;
	int atm_pressure = 101325;
	float density = 1000;
	float a  = 10;
	float b = 4;
	float d = 10;
	float c = 6;
	float w = 10;
	float particle_factor = 5;
	int number_of_particles = b*c*w*particle_factor*particle_factor*particle_factor; // 40*60*100 ~ b*c*w
	float mass = density*b*c*w/float(number_of_particles);
	int table_size = NextPrime(2*number_of_particles);

public:
	environment();
	~environment();

};