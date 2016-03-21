//// 	WE ARE SIMULATING A DAM BREAK
	/*
	http://www.ijser.org/paper/Simulation-of-Dam-Break-Using-Modified-Incompressible-Smoothed/Image_034.jpg
	
		****w corresponds to depth****
	\										\
	\										\
	\ <-b->									\
	\\\\\\\\\								\
	\		\ 								\  |
	\		\ |								\  d
	\		\ c								\  |
	\		\ |								\
	\		\ 								\
	\		\								\
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
	int part[3] = {20, 30, 50};
	float**** particles;
	//Define your particles here

	environment environ1;
	environ1.gravity(9.81);
	environ1.time_step(0.01);
	environ1.temperature(283.15);
	environ1.atmosphericPressure(101325);
	environ1.density(1000);
	environ1.numParticles(part);
	environ1.fluidVolume(4*6*10);
	environ1.environmentLength(10);
	environ1.environmentHeight(10);
	environ1.environmentWidth(10);
	environ1.addParticles(particles);
	environ1.envInit();
	environ1.simulate();
	environ1.environmentFree();
	delete(part);//this is needed to be done manually	
}