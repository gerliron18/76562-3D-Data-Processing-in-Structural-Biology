# Immune Synapse Dynamics Simulation  

This is a basic simulation of the dynamics in an immune synapse, 
specifically simulating different particles in a T-Cell bound to 
an APC (antigen-presenting cell).
For the full GUI experience run gui\GUI.py. Otherwise, run the code using the API in 
run_simulation.py under the main folder.
The program allows the user to either run a MCMC simulation with Brownian
Dynamics, based on input images\frames, or just run the 
Brownian Dynamics simulation.  

##Start screen

![Alt text](https://github.cs.huji.ac.il/drordod/StructuralBiologyHackathon/blob/master/gui/readme_images/start_screen.JPG?raw=true)

##Brownian Dynamics Screen

![Alt text](https://github.cs.huji.ac.il/drordod/StructuralBiologyHackathon/blob/master/gui/readme_images/bd_screen.JPG?raw=true)



##MCMC

![Alt text](https://github.cs.huji.ac.il/drordod/StructuralBiologyHackathon/blob/master/gui/readme_images/mcmc_screen.JPG?raw=true)


Biological background:
- 
Our goal in this project was to develop and apply a simulation of the
 dynamics in an immunological synapse. As a start, we focused on two types of
 molecules – T Cell Receptor (TCR) and CD45. 
The T cell receptor is a protein complex found on the membrane surface of T-Cells, which plays a major role in detecting antigen peptides. Antigen
 peptides are represented by a Major Histocompatibility Complex (MHC) and
  the T cell receptor recognizes those peptides according to its specificity
  . These situation is referred to as the immunological synapse. After the
   MHC-TCR complexs bind together, other proteins are recruited and a series of
    events occurrs to activate the complex. An important protein in this
     process is the CD45, a tyrosine phosphatase expressed on all nucleated
      hematopoietic cells. The CD45 can phosphorylate or dephosphorylate the
       TCRs, thus creating an interaction between those proteins.

In the project we focused on the dynamics of these two molecules and their behavior as a function of time. 
There are many challenges and questions we needed to consider, such as the
 interactions and forces between the proteins, how these change as a function
  of time, how each protein interacts with itself, the density of the
   proteins in the cell and a lot more. The research focused on the questions mentioned above. 
- Our major assumptions are as following: 

	- The number TCR proteins in an area of squared micron is 1000 and the number of CD45 is 300. 
	- For the Brownian dynamics simulation we used diffusion coefficients of 0.01 (μm^2)/s for TCR proteins and 0.037 (μm^2)/s for CD45. 
	- There are two types of forces we took into consideration: A force due
	 to the position of a particle in the grid (that represents the cell
	 ), and a force due to the location of the particle and its interactions
	  with, or affects due, to other proteins.
	  
	- There is almost no force between a TCR and a CD45, as the CD45 is in charge of phosphorylation and then it’s more stable at another place. 
	- Two proteins which have the same type pull each other, thus applying a
	 positive force. The force between two TCR proteins is stronger than
	  between two CD45 proteins, as it is more likely that TCRs will cluster
	   together. 

There are a lot of things to take into consideration, and understanding this
 deeply can affect our model. However, we believe this simulation is a good
  start and models the central biological process that occurs in an immune
   synapse. 



Internals:
- 
- MCMC Flow: 
    - A user should provide two frames, one for the start scenario and one for the end scenario. We use
     simple image processing to detect the different particles in the cell (under 
     image_processing\extract_parameters_from_frame)
    - Besides the frames, the simulation gets a **force_function** from the user (defined in 
    Forces\choose_force_function) and other relevant parameters for the run. 
    - Run a MCMC (Markov Chain Monte Carlo) simulation on the configuration space relevant to the 
    requested force function. The user can define the number of random start places for the MCMC.
    - Every MCMC iteration, the MCMC runs a Brownian Dynamics simulation from starting at the first frames time and 
    ending in the second frames time. 
    - The MCMC will calculate the loss function between the given
     information (second frame) and the 
    simulation output and will keep the information on the configurations that provided a good loss score.
    - Finally the MCMC will return the configurations that had the best value and a full mp4 video 
    of the simulation. 

- Brownian Dynamics Flow:
    - This part of the program can either receive an initial frame and
     extract the information from the image, or given a frame size the
      program will randomly generate a simulation, where the number and
       types of particles are proportional to what's accepted in the
        literature. 
    - The user inputs the number of iterations the simulation should run, and 
    the program saves an mp4 file with the particles movements. 
    
    Although we aspired to a writing generic program, if a user wants to
     add more particles or change the force 
    function or MCMC parameters, this will required working with the
     actual code and making the following changes:
     
- Using a Different Force Function:
    - A user can input his own force function to use in the Brownian
     Dynamics simulation. To do so he must add the name of the function to the
      ALL_FORCE_FUNCTIONS list in the Forces/force_calculation.py file.
    - Add a condition to the get_force_function for this function, should
     return the function and a list of the coefficients needed to run the function.

- Adding More Particles:
    - Particles
        - To add additional types of particles to the simulation, add the
         new class (inherited from Particle) in the Particle.py file, and
          implement the get_particle_type function. 
    - Initial configuration definition
        - From picture - need to add support in the get_particles_from_frame function in extract_parameters_from_frame.py
        - From frame size - to run a _random_ Brownian Dynamics simulation with the new particle type, you need to add
            the particle type to the random_BD() function (more detailed
             information and an example implementation are in the function)
    - GUI 
        - Add condition to the get_color function in the gen_png_mp4.py file
    - Forces: 
        - Add conditions to the _get_interaction_mat for every new particle pair 
        - Add a condition to _location_force with the location dependant force function for this particle
        - If the conditions aren't added the forces for each of the new particles will be zero

- Notes:
    - The image processing assume that each pixel represent one particle
     with no radius.
    