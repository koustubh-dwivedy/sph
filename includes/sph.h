#include "hash.h"
#include "kernels.cpp"
#include "prime.cpp"
#include <iostream>
#include <CL/cl.h>

float default_kernel(float position[3], float particle[6], float radius);
void default_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float default_kernel_laplacian(float position[3], float particle[6], float radius);
float pressure_kernel(float position[3], float particle[6], float radius);
void pressure_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float pressure_kernel_laplacian(float position[3], float particle[6], float radius);
float viscosity_kernel(float position[3], float particle[6], float radius);
void viscosity_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float viscosity_kernel_laplacian(float position[3], float particle[6], float radius);

bool IsPrime(long int number);
long int NextPrime(long int a);








class environment
{
private:
	//// Environment Parameters
	float gravity;
	float time_step;
	float time_duration;
	float temperature;
	long int atm_pressure;
	float density;//rest density of fluid
	float a;
	float d;
	float w;
	float compact_support_radius;
	int estimatedNumNearestNeighbours;
	int number_of_particles;
	int* number_of_particles_array;
	// this is an array of 3 int's which specify number of [articles in each dimension.
	// they may or may not be completely filled.
	// another tag of "isParticle" needs to be added in the ****particles to specify if that point is a particle or not.
	// enter number_of_particles_array as [x, y, z]
	/*
		  z|   /y
		   |  /
		   | /
		   |/_____x
	*/
	float mass;
	long int table_size;
	float**** particles;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//particles_1 is the matrix which stores the content of (t-delta(t)) iteration
	//particles_2 is the matrix in which the compute takes place and later get flushed to particles_1
	float**** particles_1;//these are copies of particles with additional space for storing mass, density, force etc.
	//these will take part in the actual compute
	float**** particles_2;//these are copies of particles with additional space for storing mass, density, force etc.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// construct particles as ->
	/*
		  z|   /y
		   |  /
		   | /
		   |/_____x
	*/
	float volume_of_fluid;
	bool ready;




	float buoyancyDiffusionCoeff;
	float viscosityCoeff;
	float surfaceTensionCoeff;
	float thresholdCoeff;//What is threshold ???
	float gasStiffnessCoeff;
	float restitutionCoeff;


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
		number_of_particles = 0;
		mass = 0;
		table_size = 0;
		estimatedNumNearestNeighbours = 80;
		//THE ABOVE VALUE OF 2000 HAS BEEN GIVEN RANDOMLY.
	}
	void setGravity(float gvt){
		if(gvt >= 0){
			gravity = gvt;
		}
		else{
			std::cout << "Please enter a valid gravity value. gravity has to be >=0.\n";
		}
	}
	void setTimeStep(float stp){
		if(stp > 0){
			time_step = stp;
		}
		else{
			std::cout << "Please enter a valid time_step value. time_step has to be >0.\n";
		}
	}
	void setTimeDuration(float dur){
		if(dur > 0){
			time_duration = dur;
		}
		else{
			std::cout << "Please enter a valid time_duration value. time_duration has to be >0.\n";
		}
	}
	void setTemperature(float temp){
		if(temp >= 0){
			temperature = temp;
		}
		else{
			std::cout << "Please enter a valid temperature value. temperature has to be >=0.\n";
		}
	}
	void setAtmosphericPressure(long int atm_pr){
		if(atm_pr >= 0){
			atm_pressure = atm_pr;
		}
		else{
			std::cout << "Please enter a valid atmospheric_pressure value. atmospheric_pressure has to be >=0.\n";
		}
	}
	void setDensity(float dens){
		if(dens > 0){
			density = dens;
		}
		else{
			std::cout << "Please enter a valid density value. density has to be >=0.\n";
		}
	}		
	void numParticles(int num){
		if(num > 0){
			number_of_particles = num;
		}
		else{
			std::cout << "Please enter a valid numParticles value. numParticles has to be >0.\n";
		}
	}
	void particlesDimen(int* a){

		number_of_particles_array = a;
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
	void estimatedNumNearestNeighboursValue(int n){
		if(n > 0){
			estimatedNumNearestNeighbours = n;
		}
		else{
			std::cout << "Please enter a valid estimatedNumNearestNeighbours value. estimatedNumNearestNeighbours has to be >0.\n";
		}
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

		////copying contents of "particles" in "particles_1" and "particles_2"
		particles_1 = new float***[number_of_particles_array[0]];
		particles_2 = new float***[number_of_particles_array[0]];
		for(int i=0; i<number_of_particles_array[0]; i++){
			particles_1[i] = new float**[number_of_particles_array[1]];
			particles_2[i] = new float**[number_of_particles_array[1]];
			for(int j=0; j<number_of_particles_array[1]; j++){
				particles_1[i][j] = new float*[number_of_particles_array[2]];
				particles_2[i][j] = new float*[number_of_particles_array[2]];
				for(int k=0; k<number_of_particles_array[2]; k++){
					particles_1[i][j][k] = new float[10];
					particles_2[i][j][k] = new float[10];
					//0 -> is particle or not
					//1 -> x coord
					//2 -> y coord
					//3 -> z coord
					//4 -> Vx value
					//5 -> Vy value
					//6 -> Vz value
					//7 -> density
					//8 -> pressure
					//

					particles_1[i][j][k][0] = particles[i][j][k][0];
					particles_2[i][j][k][0] = particles[i][j][k][0];

					particles_1[i][j][k][1] = particles[i][j][k][1];
					particles_2[i][j][k][1] = particles[i][j][k][1];

					particles_1[i][j][k][2] = particles[i][j][k][2];
					particles_2[i][j][k][2] = particles[i][j][k][2];

					particles_1[i][j][k][3] = particles[i][j][k][3];
					particles_2[i][j][k][3] = particles[i][j][k][3];

					particles_1[i][j][k][4] = particles[i][j][k][4];
					particles_2[i][j][k][4] = particles[i][j][k][4];

					particles_1[i][j][k][5] = particles[i][j][k][5];
					particles_2[i][j][k][5] = particles[i][j][k][5];

					particles_1[i][j][k][6] = particles[i][j][k][6];
					particles_2[i][j][k][6] = particles[i][j][k][6];
				}
			}
		}

		// x has been hard-coded as 40 ~= 33 = (6)+(9*2+8)+(1) . IT CAN AND SHOULD BE TWEAKED
		//reference values for x can be found in the paper
		//////////////////////////////////////////////////////////////////////////////////////////
		compact_support_radius = pow(3*volume_of_fluid*40/(4*M_PI*number_of_particles),1/3);
		//////////////////////////////////////////////////////////////////////////////////////////

		std::cout << "INPUT PARAMS AS OF NOW: \n \n";
		std::cout << "density " << density << "\n";
		std::cout << "numParticles " << number_of_particles << "\n";
		std::cout << "fluidVolume " << volume_of_fluid << "\n";
		std::cout << "mass of 1 particle " << mass << "\n";
		std::cout << "gravity " << gravity << "\n";
		std::cout << "timeStep " << time_step << "\n";
		std::cout << "temperature " << temperature << "\n";
		std::cout << "pressure " << atm_pressure << "\n";
		std::cout << "compact_support_radius " << compact_support_radius << "\n";
		std::cout << "hash table size " << table_size << "\n";
		std::cout << "estimatedNumNearestNeighbours " << estimatedNumNearestNeighbours << "\n";
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
	void environmentFree(void){		
		/****freeing up memory*******/
		for(int i=0; i<number_of_particles_array[0]; i++){
			for(int j=0; j<number_of_particles_array[1]; j++){
				for(int k=0; k<number_of_particles_array[2]; k++){
					delete(particles[i][j][k]);
					delete(particles_1[i][j][k]);
					delete(particles_2[i][j][k]);
				}
				delete(particles[i][j]);
				delete(particles_1[i][j]);
				delete(particles_2[i][j]);
			}
			delete(particles[i]);
			delete(particles_1[i]);
			delete(particles_2[i]);
		}
		delete(particles);
		delete(particles_1);
		delete(particles_2);
		/****freeing up memory*******/
	}
	//this function defined for internal use in "simulate" only
	//it copies the contents of particles_2 to particles_1
	void copy(void){
		for(int i=0; i<number_of_particles_array[0]; i++){
			for(int j=0; j<number_of_particles_array[1]; j++){
				for(int k=0; k<number_of_particles_array[2]; k++){
					particles_1[i][j][k][0] = particles_2[i][j][k][0];
					particles_1[i][j][k][1] = particles_2[i][j][k][1];
					particles_1[i][j][k][2] = particles_2[i][j][k][2];
					particles_1[i][j][k][3] = particles_2[i][j][k][3];
					particles_1[i][j][k][4] = particles_2[i][j][k][4];
					particles_1[i][j][k][5] = particles_2[i][j][k][5];
					particles_1[i][j][k][6] = particles_2[i][j][k][6];
					particles_1[i][j][k][7] = particles_2[i][j][k][7];
					particles_1[i][j][k][8] = particles_2[i][j][k][8];
					particles_1[i][j][k][9] = particles_2[i][j][k][9];
				}
			}
		}
	}

	void simulate(void){
		if(ready == 1){

			hashTable table;
			table.insertParticles(particles_1, table_size, number_of_particles_array, compact_support_radius);
			int actualNearestNeighbours;//this varies from particle to particle
			float* nearestNeighbourAdresses[estimatedNumNearestNeighbours];
			//PASSING ON A CONSTANT SIZED ARRAY TO MAKE THINGS SIMPLE. (OTHERWISE DYNAMIC ALLOCATION WOULD NEED TO BE DONE)
			//THIS SIZE IS TO BE CHANGED DEPENDING ON PROBLEM TO PROBLEM. (or maybe not?)

			float temp_density;
			float temp_pressure;
			float temp_force_x, temp_force_y, temp_force_z;
			float receiver[3];

			
			//1. compute density and pressure
			//2. compute internal forces. internal <- pressure+viscosity
			//3. compute external forces
			//4. time integration and collision handling

			for(float t=0; t<time_duration; t = t + time_step){
				//1.
				for(int i=0; i<number_of_particles_array[0]; i++){
					for(int j=0; j<number_of_particles_array[1]; j++){
						for(int k=0; k<number_of_particles_array[2]; k++){
							table.particleQuery(particles_1[i][j][k][1], particles_1[i][j][k][2], particles_1[i][j][k][3], nearestNeighbourAdresses, estimatedNumNearestNeighbours, &actualNearestNeighbours, compact_support_radius, 1);
							temp_density = 0;
							temp_pressure = 0;
							for(int l=0; l<actualNearestNeighbours; l++){
								temp_density = temp_density + mass*default_kernel(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius);
							}
							temp_pressure = gasStiffnessCoeff*(temp_density - density);
							particles_1[i][j][k][7] = temp_density;
							particles_1[i][j][k][8] = temp_pressure;
						}
					}
				}
				//2.
				for(int i=0; i<number_of_particles_array[0]; i++){
					for(int j=0; j<number_of_particles_array[1]; j++){
						for(int k=0; k<number_of_particles_array[2]; k++){
							table.particleQuery(particles_1[i][j][k][1], particles_1[i][j][k][2], particles_1[i][j][k][3], nearestNeighbourAdresses, estimatedNumNearestNeighbours, &actualNearestNeighbours, compact_support_radius, 0);
							temp_force_x = 0;
							temp_force_y = 0;
							temp_force_z = 0;
							for(int l=0; l<actualNearestNeighbours; l++){
								//the following is pressure force
								pressure_kernel_gradient(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius, receiver);
								float temp_var = 0;
								temp_var = mass*((nearestNeighbourAdresses[l][8]/((float)((nearestNeighbourAdresses[l][7])*(nearestNeighbourAdresses[l][7])))) + (particles_1[i][j][k][8]/((float)(particles_1[i][j][k][7])*(particles_1[i][j][k][7]))));
								temp_force_x = (temp_force_x + receiver[0]*temp_var)*(-particles_1[i][j][k][7]);
								temp_force_y = (temp_force_y + receiver[1]*temp_var)*(-particles_1[i][j][k][7]);
								temp_force_z = (temp_force_z + receiver[2]*temp_var)*(-particles_1[i][j][k][7]);
								//
								//the following is viscous force
								temp_var = viscosityCoeff*(mass/((float)nearestNeighbourAdresses[l][7]))*viscosity_kernel_laplacian(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius);
								temp_force_x = temp_force_x + temp_var*(nearestNeighbourAdresses[l][4] - particles_1[i][j][k][4]);
								temp_force_y = temp_force_y + temp_var*(nearestNeighbourAdresses[l][5] - particles_1[i][j][k][5]);
								temp_force_z = temp_force_z + temp_var*(nearestNeighbourAdresses[l][6] - particles_1[i][j][k][6]);
								//
							}
						}
					}
				}











				//export data and hash table reinitialize
			}

			table.free();
		}
	}
};