# SPH
Smoothed Particle Hydrodynamics simulation using OpenCL <br />
Implementation of following paper -> http://www.glowinggoo.com/sph/bin/kelager.06.pdf <br />
Part of ME 766 High Performance Scientific Computing course project <br />
Code built and tested on Ubuntu 14.04 <br />
As of now, code cannot handle number of particles greater than ~30,000. To handle them, variables will need to be converted from int to long int. <br /><br />
*NOTE that in the VTK file generated, the velocity values present are not at that moment. For getting the correct velocity values, you need to take the average of that velocity and the next time step's velocity.<br />