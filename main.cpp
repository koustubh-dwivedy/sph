//// 	WE ARE SIMULATING A DAM BREAK
	/*
	http://www.ijser.org/paper/Simulation-of-Dam-Break-Using-Modified-Incompressible-Smoothed/Image_034.jpg
	
		****w corresponds to depth****
	\						\
	\						\
	\ <-b->					\
	\\\\\\\\\				\
	\		\ 				\  |
	\		\ |				\  d
	\		\ c				\  |
	\		\ |				\
	\		\ 				\
	\		\				\
	\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
			<--a-->

			!!!!!This does not follow vector right hand rule!!!!!
			y     z
			|   /
			|  /
			| /
			--->x

	*/
////

#include "includes/sph.h"

int main(int argc, char **argv){
	//Create an array containing x, y, z coordinated of particles and pass it on to sph() function
	//CLASS -> environment
	//environment contains the following parameters
	
	// gravity = 9.81;
	// time_step = 0.01;
	// temperature =  283.15;
	// atm_pressure = 101325;
	// density = 1000;
	// a  = 10;
	// b = 4;
	// d = 10;
	// c = 6;
	// w = 10;
	// particle_factor = 5;
	// number_of_particles = b*c*w*particle_factor*particle_factor*particle_factor; // 40*60*100 ~ b*c*w
	// mass = density*b*c*w/float(number_of_particles);
	// table_size = NextPrime(2*number_of_particles);

	//environment contains the following methods

	//sph
	
}