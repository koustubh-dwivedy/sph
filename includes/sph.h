#include "hash.h"
#include "kernels.cpp"
#include "prime.cpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
//#include <CL/cl.h>

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
	//particles_1 is the matrix which stores the content of (time t) iteration. (it's a local copy)
	float**** particles_1;//this is the copy of particles with additional space for storing mass, density, force etc.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// construct particles as ->**this coordinate system is not valid for applying forces
	/*
		  z|   /y
		   |  /
		   | /
		   |/_____x
	*/
	float volume_of_fluid;
	bool is_liquid;//if isLiquid is 1 then buoyancy simulation is not performed
	//if isLiquid is 0 then it is assumed that fluid is gas
	bool ready;




	float buoyancyDiffusionCoeff;
	float viscosityCoeff;
	float surfaceTensionCoeff;
	float surfaceTensionThresholdCoeff;
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
		surfaceTensionThresholdCoeff = 7.065;//currently hardcoded it for water;
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
	void setBuoyancyDiffusionCoeff(float value){
		buoyancyDiffusionCoeff = value;
	}
	void setViscosityCoeff(float value){
		viscosityCoeff = value;
	}
	void setSurfaceTensionCoeff(float value){
		surfaceTensionCoeff = value;
	}
	void setSurfaceTensionThresholdCoeff(float value){
		surfaceTensionThresholdCoeff = value;
	}
	void setGasStiffnessCoeff(float value){
		gasStiffnessCoeff = value;
	}
	void setRestitutionCoeff(float value){
		restitutionCoeff = value;
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

		////copying contents of "particles" in "particles_1"
		particles_1 = new float***[number_of_particles_array[0]];
		for(int i=0; i<number_of_particles_array[0]; i++){
			particles_1[i] = new float**[number_of_particles_array[1]];
			for(int j=0; j<number_of_particles_array[1]; j++){
				particles_1[i][j] = new float*[number_of_particles_array[2]];
				for(int k=0; k<number_of_particles_array[2]; k++){
					particles_1[i][j][k] = new float[12];
					//0 -> is particle or not
					//1 -> x coord
					//2 -> y coord
					//3 -> z coord
					//4 -> Vx value THIS IS NOT STORING VELOCITY VALUE AT TIME t BUT AT t (+-) 0.5*delta(t)
					//5 -> Vy value THIS IS NOT STORING VELOCITY VALUE AT TIME t BUT AT t (+-) 0.5*delta(t)
					//6 -> Vz value THIS IS NOT STORING VELOCITY VALUE AT TIME t BUT AT t (+-) 0.5*delta(t)
					//7 -> density
					//8 -> pressure
					//9 -> force_x
					//10-> force_y
					//11-> force_z

					particles_1[i][j][k][0] = particles[i][j][k][0];

					particles_1[i][j][k][1] = particles[i][j][k][1];

					particles_1[i][j][k][2] = particles[i][j][k][2];

					particles_1[i][j][k][3] = particles[i][j][k][3];

					particles_1[i][j][k][4] = particles[i][j][k][4];

					particles_1[i][j][k][5] = particles[i][j][k][5];

					particles_1[i][j][k][6] = particles[i][j][k][6];
				}
			}
		}

		// x has been hard-coded as 40 ~= 33 = (6)+(9*2+8)+(1) . IT CAN AND SHOULD BE TWEAKED
		//reference values for x can be found in the paper
		//////////////////////////////////////////////////////////////////////////////////////////
		compact_support_radius = pow(3*volume_of_fluid*40/(4*M_PI*number_of_particles),1/3.0);
		//////////////////////////////////////////////////////////////////////////////////////////

		std::cout << "INPUT PARAMS AS OF NOW: \n \n";
		std::cout << "density " << density << "\n";
		std::cout << "isLiquid " << is_liquid << "\n";
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
		std::cout << "BuoyancyDiffusionCoeff " << buoyancyDiffusionCoeff << "\n";
		std::cout << "ViscosityCoeff " << 	viscosityCoeff << "\n";
		std::cout << "SurfaceTensionCoeff " << 	surfaceTensionCoeff << "\n";
		std::cout << "SurfaceTensionThresholdCoeff " << 	surfaceTensionThresholdCoeff << "\n";
		std::cout << "GasStiffnessCoeff " << 	gasStiffnessCoeff << "\n";
		std::cout << "RestitutionCoeff " << 	restitutionCoeff << "\n";
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
	void isLiquid(bool a){

		is_liquid = a;
	}
	void environmentFree(void){
		/****freeing up memory*******/
		for(int i=0; i<number_of_particles_array[0]; i++){
			for(int j=0; j<number_of_particles_array[1]; j++){
				for(int k=0; k<number_of_particles_array[2]; k++){
					delete(particles[i][j][k]);
					delete(particles_1[i][j][k]);
				}
				delete(particles[i][j]);
				delete(particles_1[i][j]);
			}
			delete(particles[i]);
			delete(particles_1[i]);
		}
		delete(particles);
		delete(particles_1);
		/****freeing up memory*******/
	}
	void simulate(void){
		if(ready == 1){

			//1. compute density and pressure
			//2. compute internal forces. internal <- pressure+viscosity
			//3. compute external forces
			//4. time integration and collision handling
			//5. export data and hash table reinitialize


			//a thing to note-> the if condition which checks if the material is liquid or gas can be taken outside the ~30000*100 loops
			//this can result in significant speedup
			for(float t=0; t<time_duration; t = t + time_step){

				hashTable table;
				table.insertParticles(particles_1, table_size, number_of_particles_array, compact_support_radius);


				//1.
				#pragma omp parallel for
				for(int i=0; i<number_of_particles_array[0]; i++){
					for(int j=0; j<number_of_particles_array[1]; j++){
						for(int k=0; k<number_of_particles_array[2]; k++){
							if(particles_1[i][j][k][0] == 1){
								int actualNearestNeighbours;//this varies from particle to particle
								float* nearestNeighbourAdresses[estimatedNumNearestNeighbours];
								//PASSING ON A CONSTANT SIZED ARRAY TO MAKE THINGS SIMPLE. (OTHERWISE DYNAMIC ALLOCATION WOULD NEED TO BE DONE)
								//THIS SIZE IS TO BE CHANGED DEPENDING ON PROBLEM TO PROBLEM. (or maybe not?)

								float temp_density;
								float temp_pressure;

								table.particleQuery(particles_1[i][j][k][1], particles_1[i][j][k][2], particles_1[i][j][k][3], nearestNeighbourAdresses, estimatedNumNearestNeighbours, &actualNearestNeighbours, compact_support_radius, 0);//changed ij from 1 to 0 while testing
								temp_density = density;
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
				}
				//2. and 3.
				#pragma omp parallel for
				for(int i=0; i<number_of_particles_array[0]; i++){
					for(int j=0; j<number_of_particles_array[1]; j++){
						for(int k=0; k<number_of_particles_array[2]; k++){
							if(particles_1[i][j][k][0] == 1){
								int actualNearestNeighbours;//this varies from particle to particle
								float* nearestNeighbourAdresses[estimatedNumNearestNeighbours];
								//PASSING ON A CONSTANT SIZED ARRAY TO MAKE THINGS SIMPLE. (OTHERWISE DYNAMIC ALLOCATION WOULD NEED TO BE DONE)
								//THIS SIZE IS TO BE CHANGED DEPENDING ON PROBLEM TO PROBLEM. (or maybe not?)

								table.particleQuery(particles_1[i][j][k][1], particles_1[i][j][k][2], particles_1[i][j][k][3], nearestNeighbourAdresses, estimatedNumNearestNeighbours, &actualNearestNeighbours, compact_support_radius, 0);

								float temp_force_x=0, temp_force_y=0, temp_force_z=0;
								float surf_norm_x=0, surf_norm_y=0, surf_norm_z=0;

								for(int l=0; l<actualNearestNeighbours; l++){
									//the following is pressure force
									float receiver[3];
									pressure_kernel_gradient(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius, receiver);
									float temp_var = 0;
									temp_var = mass*((nearestNeighbourAdresses[l][8]/((float)((nearestNeighbourAdresses[l][7])*(nearestNeighbourAdresses[l][7])))) + (particles_1[i][j][k][8]/((float)(particles_1[i][j][k][7])*(particles_1[i][j][k][7]))));
									temp_force_x = temp_force_x + (receiver[0]*temp_var)*(-particles_1[i][j][k][7]);
									temp_force_y = temp_force_y + (receiver[1]*temp_var)*(-particles_1[i][j][k][7]);
									temp_force_z = temp_force_z + (receiver[2]*temp_var)*(-particles_1[i][j][k][7]);
									//
									//the following is viscous force
									temp_var = viscosityCoeff*(mass/((float)nearestNeighbourAdresses[l][7]))*viscosity_kernel_laplacian(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius);
									temp_force_x = temp_force_x + temp_var*(nearestNeighbourAdresses[l][4] - particles_1[i][j][k][4]);
									temp_force_y = temp_force_y + temp_var*(nearestNeighbourAdresses[l][5] - particles_1[i][j][k][5]);
									temp_force_z = temp_force_z + temp_var*(nearestNeighbourAdresses[l][6] - particles_1[i][j][k][6]);
									//
									if(is_liquid){//finding normal to fluid for finding force due to surface tension
										//the following is surface normal
										default_kernel_gradient(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius, receiver);
										surf_norm_x = surf_norm_x + receiver[0]*mass/((float)nearestNeighbourAdresses[l][8]);
										surf_norm_y = surf_norm_y + receiver[1]*mass/((float)nearestNeighbourAdresses[l][8]);
										surf_norm_z = surf_norm_z + receiver[2]*mass/((float)nearestNeighbourAdresses[l][8]);
										//
									}
								}
								//finding surface tension force from surface normal if the material is liquid
								if(is_liquid){
									//the following is gravity
									temp_force_z = temp_force_z + particles_1[i][j][k][7]*gravity;
									//
									float surf_norm_magn = sqrt(surf_norm_x*surf_norm_x + surf_norm_y*surf_norm_y + surf_norm_z*surf_norm_z);//this was another hell of an error "who writes a + in place of a *????????"
									for(int l=0; l<actualNearestNeighbours; l++){
										float temp_var = (mass/(nearestNeighbourAdresses[l][7]))*default_kernel_laplacian(nearestNeighbourAdresses[l], particles_1[i][j][k], compact_support_radius);
										temp_force_x = temp_force_x + (-surfaceTensionCoeff)*(surf_norm_x/((float)surf_norm_magn))*temp_var;
										temp_force_y = temp_force_y + (-surfaceTensionCoeff)*(surf_norm_y/((float)surf_norm_magn))*temp_var;
										temp_force_z = temp_force_z + (-surfaceTensionCoeff)*(surf_norm_z/((float)surf_norm_magn))*temp_var;
									}
								}
								else if(!is_liquid){//forces specific to gases
									//for buoyancy
									temp_force_z = temp_force_z + buoyancyDiffusionCoeff*gravity*(particles_1[i][j][k][7] - density);
									//
								}

								particles_1[i][j][k][9] = temp_force_x;
								particles_1[i][j][k][10] = temp_force_y;
								particles_1[i][j][k][11] = temp_force_z;
							}
						}
					}
				}

				//4.
				//the following if condition is to take care of initial velocity offset
				if(t == 0){
					for(int i=0; i<number_of_particles_array[0]; i++){
						for(int j=0; j<number_of_particles_array[1]; j++){
							for(int k=0; k<number_of_particles_array[2]; k++){
								if(particles_1[i][j][k][0] == 1){

									float temp_accn_x, temp_accn_y, temp_accn_z;

									temp_accn_x = particles_1[i][j][k][9]/((float)particles_1[i][j][k][7]);
									temp_accn_y = particles_1[i][j][k][10]/((float)particles_1[i][j][k][7]);
									temp_accn_z = particles_1[i][j][k][11]/((float)particles_1[i][j][k][7]);

									particles_1[i][j][k][4] = particles_1[i][j][k][4] - 0.5*time_step*temp_accn_x;
									particles_1[i][j][k][5] = particles_1[i][j][k][5] - 0.5*time_step*temp_accn_y;
									particles_1[i][j][k][6] = particles_1[i][j][k][6] - 0.5*time_step*temp_accn_z;
								}
							}
						}
					}
				}
				else{
					for(int i=0; i<number_of_particles_array[0]; i++){
						for(int j=0; j<number_of_particles_array[1]; j++){
							for(int k=0; k<number_of_particles_array[2]; k++){
								if(particles_1[i][j][k][0] == 1){

									float temp_accn_x, temp_accn_y, temp_accn_z;

									temp_accn_x = particles_1[i][j][k][9]/((float)particles_1[i][j][k][7]);
									temp_accn_y = particles_1[i][j][k][10]/((float)particles_1[i][j][k][7]);
									temp_accn_z = particles_1[i][j][k][11]/((float)particles_1[i][j][k][7]);

									particles_1[i][j][k][4] = particles_1[i][j][k][4] + time_step*temp_accn_x;
									particles_1[i][j][k][5] = particles_1[i][j][k][5] + time_step*temp_accn_y;
									particles_1[i][j][k][6] = particles_1[i][j][k][6] + time_step*temp_accn_z;

									particles_1[i][j][k][1] = particles_1[i][j][k][1] + time_step*particles_1[i][j][k][4];
									particles_1[i][j][k][2] = particles_1[i][j][k][2] + time_step*particles_1[i][j][k][5];
									particles_1[i][j][k][3] = particles_1[i][j][k][3] + time_step*particles_1[i][j][k][6];


									//detecting collision
									if(particles_1[i][j][k][1] < 0){
										float depth = -particles_1[i][j][k][1];
										particles_1[i][j][k][1] = 0;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][4] = -(restitutionCoeff*depth*particles_1[i][j][k][4])/((float)time_step*vel_magn);
									}
									else if(particles_1[i][j][k][1] > a){
										float depth = particles_1[i][j][k][1] - a;
										particles_1[i][j][k][1] = a;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][4] = -(restitutionCoeff*depth*particles_1[i][j][k][4])/((float)time_step*vel_magn);
									}
									if(particles_1[i][j][k][2] < 0){
										float depth = -particles_1[i][j][k][2];
										particles_1[i][j][k][2] = 0;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][5] = -(restitutionCoeff*depth*particles_1[i][j][k][5])/((float)time_step*vel_magn);
									}
									else if(particles_1[i][j][k][2] > w){
										float depth = particles_1[i][j][k][2] - w;
										particles_1[i][j][k][2] = w;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][5] = -(restitutionCoeff*depth*particles_1[i][j][k][5])/((float)time_step*vel_magn);
									}
									if(particles_1[i][j][k][3] < 0){
										float depth = -particles_1[i][j][k][3];
										particles_1[i][j][k][3] = 0;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][6] = -(restitutionCoeff*depth*particles_1[i][j][k][6])/((float)time_step*vel_magn);
									}
									else if(particles_1[i][j][k][3] > d){
										float depth = particles_1[i][j][k][3] - d;
										particles_1[i][j][k][3] = d;
										float vel_magn = sqrt(particles_1[i][j][k][4]*particles_1[i][j][k][4] + particles_1[i][j][k][5]*particles_1[i][j][k][5] + particles_1[i][j][k][6]*particles_1[i][j][k][6]);
										particles_1[i][j][k][6] = -(restitutionCoeff*depth*particles_1[i][j][k][6])/((float)time_step*vel_magn);
									}
								}
							}
						}
					}
				}
				//5.
				table.free();
				//method for exporting data
				
				if(t == 0){
					std::cout << "Writing Data\n";
				}
				std::ofstream outfile;
				std::stringstream sstm;
				sstm << "output/sphResult_step" << (int)(t/time_step) << ".vtk";
				std::string name;
				name = sstm.str();
				std::cout << name << std::endl;
				char *cstr = new char[name.length() + 1];
				std::strcpy(cstr, name.c_str());
				outfile.open(cstr);
				outfile << "# vtk DataFile Version 2.0\n";
				outfile << "SPH Data\n";
				outfile << "ASCII\n";
				outfile << "DATASET POLYDATA\n";
				outfile << "POINTS " << number_of_particles << " float\n";
				for(int i=0; i<number_of_particles_array[0]; i++){
					for (int j=0; j<number_of_particles_array[1]; j++){
						for(int k=0; k<number_of_particles_array[2]; k++){
							outfile << (float)1.0*particles_1[i][j][k][1] << " " << (float)1.0*particles_1[i][j][k][2] << " " << (float)1.0*particles_1[i][j][k][3] << "\n";
						}
					}
				}

			  	outfile.close();
				

				/*
				//the following is only for testing
				for(int i=0; i<number_of_particles_array[0]; i++){
					for(int j=0; j<number_of_particles_array[1]; j++){
						std::cout << particles_1[i][j][3][1] << " " << particles_1[i][j][3][2] << " " << particles_1[i][j][3][3] << std::endl;
					}
					std::cout << std::endl;
					std::cout << std::endl;
					std::cout << std::endl;
				}
				*/
			}
		}
	}
};
