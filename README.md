CERN
====

2D advection-diffusion Peaceman-Rachford algorithm implementation for charge carrier dynamics in Medipix-3 silicon pixel detectors (no aluminum layer on bottom surface).


============================= DESCRIPTION ==================================

The source code contained in this folder numerically propagates a linear
advective-diffusive system of charge carriers in pure silicon with an 
inhomogeneous applied electric field in either one or two dimensions, as 
specified by the user. In one dimension, the Crank-Nicolson implicit scheme
is used. More specifically, the Laplacian is treated implicitly in time
and the advection term is treated explicitly. In two dimensions, a 
Alternating-Direction Implicit scheme is used; in particular, the Peaceman-Rachford algorithm. Again, the Laplacian is
treated impilcitly and the advection terms are treated expilcitly with
central difference stencils. The simulation is capable of propagating either
electrons or holes through the medium.


============================ DEPENDENCIES ===================================

The vector and matrix objects used in the one and two dimensional algorithms
are taken from the Class Library for High Energy Physics, which can be
obtained at:

http://proj-clhep.web.cern.ch/proj-clhep/

Some of the plotting functions defined in the header file incorporate objects
from the ROOT library. Information concerning the installation of the ROOT
library from source can be found at:
 
http://root.cern.ch/drupal/content/installing-root-source


=========================== IMPLEMENTATION ==================================

1. Create an executable from the source code and the makefile.

2. Navigate to the directory containing the executable.

3. Run " ./(name_of_executable) (dimensions) (species) (run_type) (bias_type) (voltage) (temperature) (time_steps) (time_step_size) "

where: (dimensions) = integer (either 1 or 2) specifying the dimensionality of the computational domain.
       (species) = string (either 'electron' or 'hole') specifying the type of charge carrier to propagate.
       (run_type) = char (either 'a','n', or 'd') specifying either an analytic, numeric, or difference run, respectively.
       (bias_type) = char (either 'f' or 'b') specifying the applied voltage to be either forward or reverse biased, respectively.  
       (voltage) = double specifying the magnitude of the applied voltage in volts.
       (temperature) = integer specifying the ambient temperature in kelvin.
       (time_steps) = integer specifying the number of time steps.
       (time_step_size) = double specifying the temporal step size.

4. EXAMPLE: Let the name of an example executable be "Diffusion". Then, if we want to run the 2D numerical algorithm with electrons using
a forward bias applied voltage of 87.9 V at 300 K with 100 times steps of step size 1E-11, we would run:

      ./Diffusion 2 electron n f 87.9 300 100 1E-11


Running the executable produces two output files. The 1D algorithm produces files "1DDiffusionINITIAL.txt" and "1DDiffusionFINAL.txt",
which provide the charge density distribution at the initial and final time steps, respectively. Likewise, the 2D algorithm produces 
files "2DDiffusionINITIAL.txt" and "2DDiffusionFINAL.txt", which provide the coordinates of the charge density distribution at the 
initial and final time steps, respectively. 

=========================== STABILITY ISSUES/PRECAUTIONS ======================

Although the purely diffusive 1D Crank-Nicolson and 2D alternating-direction 
implicit schemes are unconditionally stable, the incorporation of the advective
term (i.e. the electric field) ruins this property. Thus, both algorithms 
break down if the time step used in the propagation is too large. For the 1D
algorithm, the largest time step that can be used before break down is dt = 1E-6.
For the 2D algorithm, the largest time step that can be used before break down
is dt = 1E-11. 




/////////////////////////////////////////////////////////////////////////////


CONTACT INFO:

Mason Biamonte
masonbiamonte@gmail.com

John Idarraga
CERN
idarraga@cern.ch
