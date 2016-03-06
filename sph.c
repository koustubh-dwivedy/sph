#include <stdio.h>
#include "CL/cl.h"

int main(int argv, char ** argc){
	/********
	Initializing SPH System
	********/
	//// Environment Parameters
	const float gravity = 9.81;
	const float time_step = 0.01;
	const float temperature =  283.15;
	const int atm_pressure = 101325;
	const int number_of_particles = 4000;
	const float a  = 10;
	const float b = 4;
	const float d = 10;
	const float c = 6;
	const float w = 10;
	////


	
	//// WE ARE SIMULATING A DAM BREAK
	/*
	http://www.ijser.org/paper/Simulation-of-Dam-Break-Using-Modified-Incompressible-Smoothed/Image_034.jpg


		****w corresponds to depth****
	\					\
	\					\
	\ <-b->					\
	\\\\\\\\\				\
	\	\ 				\  |
	\	\ |				\  d
	\	\ c				\  |
	\	\ |				\
	\	\ 				\
	\	\				\
	\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
			<--a-->
	*/
	////



	//// Creating n particles and setting their initial positions, velocities and mass
		
}
