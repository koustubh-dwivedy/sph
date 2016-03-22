#include "sph.h"

void environment::simulate(){
	if(ready == 1){

		//// Initialize the smoothing kernels using (5.14) to compute the compact support radius
		// x has been hard-coded as 13
		float compact_support_radius = pow(3*b*c*w*13/(4*M_PI*number_of_particles),1/3.);



















	}
}