#include "kernels.cpp"

float default_kernel(float position[3], float particle[6], float radius);
void default_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float default_kernel_laplacian(float position[3], float particle[6], float radius);
float pressure_kernel(float position[3], float particle[6], float radius);
void pressure_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float pressure_kernel_laplacian(float position[3], float particle[6], float radius);
float viscosity_kernel(float position[3], float particle[6], float radius);
void viscosity_kernel_gradient(float position[3], float particle[6], float radius, float* receiver);
float viscosity_kernel_laplacian(float position[3], float particle[6], float radius);
