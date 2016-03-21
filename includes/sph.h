#include <stdio.h>
#include <CL/cl.h>
#include <math.h>
#include <iostream>
#include "include.h"
#include "hash.h"

class environment
{
private:
	//// Environment Parameters
	float gravity;
	float time_step;
	float temperature;
	long int atm_pressure;
	float density;
	float a;
	float d;
	float w;
	int number_of_particles;
	int* number_of_particles_array;
	// this is an array of 3 int's which specify number of [articles in each dimension.
	// they may or may not be completely filled.
	// another tag of "isParticle" needs to be added in the ****particles to specify if that point is a particle or not.
	// enter number_of_particles_array as [x, y, z]
	/*
				!!!!!This does not follow vector right hand rule!!!!!
			y     z
			|   /
			|  /
			| /
			--->x
	*/
	float mass;
	int table_size;
	float**** particles;
	// construct particles as ->
	/*
				!!!!!This does not follow vector right hand rule!!!!!
			y     z
			|   /
			|  /
			| /
			--->x
	*/
	float volume_of_fluid;
	bool ready;

public:
	environment(){
		gravity = 9.81;
		time_step = 0.01;
		temperature = 283.15;
		atm_pressure = 101325;
		density = 0;
		a = 0;
		d = 0;
		w = 0;
		number_of_particles = 0;//b*c*w*particle_factor*particle_factor*particle_factor;
		mass = 0;
		table_size = 0;
	}
	void gravity(float gvt){
		if(gvt >= 0){
			gravity = gvt;
		}
		else{
			std::cout << "Please enter a valid gravity value. gravity has to be >=0.\n";
		}
	}
	void time_step(float stp){
		if(stp > 0){
			time_step = stp;
		}
		else{
			std::cout << "Please enter a valid time_step value. time_step has to be >0.\n";
		}
	}
	void temperature(float temp){
		if(temp >= 0){
			temperature = temp;
		}
		else{
			std::cout << "Please enter a valid temperature value. temperature has to be >=0.\n";
		}
	}
	void atmosphericPressure(long int atm_pr){
		if(atm_pr >= 0){
			atm_pressure = atm_pr;
		}
		else{
			std::cout << "Please enter a valid atmospheric_pressure value. atmospheric_pressure has to be >=0.\n";
		}
	}
	void density(float dens){
		if(dens > 0){
			density = dens;
		}
		else{
			std::cout << "Please enter a valid density value. density has to be >=0.\n";
		}
	}

		
	void numParticles(int* num){
		int temp;
		number_of_particles_array = num;
		temp = number_of_particles_array[0]*number_of_particles_array[1]*number_of_particles_array[2];

		if(temp > 0){
			number_of_particles = temp;
		}
		else{
			std::cout << "Please enter a valid numParticles value. numParticles has to be >0.\n";
		}
	}
	void fluidVolume(float vol){
		if(vol > 0){
			volume_of_fluid = vol;
		}
		else{
			std::cout << "Please enter a valid fluidVolume value. fluidVolume has to be >0.\n";
		}
	}
	void environmentLength(float length){
		if(length > 0){
			a = length;
		}
		else{
			std::cout << "Please enter a valid environmentLength value. environmentLength has to be >0.\n";
		}
	}
	void environmentHeight(float height){
		if(height > 0){
			d = height;
		}
		else{
			std::cout << "Please enter a valid environmentHeight value. environmentHeight has to be >0.\n";
		}
	}
	void environmentWidth(float width){
		if(width > 0){
			w = width;
		}
		else{
			std::cout << "Please enter a valid environmentWidth value. environmentWidth has to be >0.\n";
		}
	}
	void addParticles(float**** part){
		particles = part;
	}
	void envInit(void){
		if(density > 0){
			if(number_of_particles > 0){
				mass = volume_of_fluid*density/((float)number_of_particles);
				table_size = NextPrime(2*number_of_particles);
			}
			else if(number_of_particles == 0){
				std::cout << "PLEASE ENTER NUMBER OF PARTICLES. \n";
			}
		}

		std::cout << "INPUT PARAMS AS OF NOW: \n \n";
		std::cout << "density " << density << "\n";
		std::cout << "numParticles " << number_of_particles << "\n";
		std::cout << "fluidVolume " << volume_of_fluid << "\n";
		std::cout << "mass of 1 particle " << mass << "\n";
		std::cout << "gravity " << gravity << "\n";
		std::cout << "timeStep " << time_step << "\n";
		std::cout << "temperature " << temperature << "\n";
		std::cout << "pressure " << atm_pressure << "\n";
		std::cout << "hash table size " << table_size << "\n";
		std::cout << "environmentLength " << a << "\n";
		std::cout << "environmentHeight " << d << "\n";
		std::cout << "environmentWidth " << w << "\n \n \n";


		std::cout << "DO YOU WISH TO PROCEED WITH THE ABOVE VALUES? (1 or 0)" << "\n";
		bool check;
		std::cin >> check;
		if(check == 1){
			ready = 1;
		}
		else{
			ready = 0;
			std::cout << "THE PROGRAM WILL (SHOULD) NOW EXIT. PLEASE RECOMPILE WITH CORRECTED VALUES " << "\n";
		}
	}
	void simulate();
	void environmentFree(void){		
		/****freeing up memory*******/
		for(int i=0; i<number_of_particles_array[0]; i++){
			for(int j=0; j<number_of_particles_array[1]; j++){
				delete(particles[i][j]);
			}
			delete(particles[i]);
		}
		delete(particles);
		/****freeing up memory*******/
	}
};