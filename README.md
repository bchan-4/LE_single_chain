# Active loop extrusion simulation on a single polymer chain
## Brian Chan, Duke University, 2023

This document will instruct you how to compile and run a hybrid Molecular Dynamics - Monte Carlo (MD-MC) simulation code for steady state active loop extrusion on a single polymer chain.

--------
Prerequisites
----
The main piece of code is ```ExtrSim.c```, which uses LAMMPS (23Jun2022 version) as a C library and also uses MPI for parallelization. To compile this code, you must first build LAMMPS as a shared library (see LAMMPS documentation https://docs.lammps.org/Build_link.html). You must also ensure that the caller code can find the library, which may require editing your environemntal variables. Note that this code was tested with OpenMPI 4.1.1.

------
Compilation
----
Once LAMMPS is built as a shared library, make sure that the following files are in your working directory:
- ```Makefile``` will compile the code, after appropriate modifications.
- ```ExtrSim.c``` is the main code with active loop extrusion. 

Before compiling the code, load your MPI implementation, for example with ``` module load OpenMPI/4.1.1 ``` or something similar. **Make sure to you use the same MPI implementation used for building LAMMPS.** 

Then, compile with ```make ExtrSim```. Note that the Makefile will need to be modified to fit your specific machine and directory structure. Specifically, the following snippets could be changed:
- ```CC=mpicc``` to fit your C compiler
- ```-I/path/to/LAMMPS/src/``` to specify the path of the LAMMPS src directory
- ```-llammps_mpi``` to specify the name and/or path of the LAMMPS shared library. In this case, the LAMMPS library was called ```liblammps_mpi.so``` and its path was in the ```LD_LIBRARY_PATH``` environment variable.

There should now be an executable called ```ExtrSim``` in your working directory.

----------------
Running a simulation
---
Ensure that your working directory has the following files:
- ```in.FirstSteps``` is a LAMMPS input file used for initial MD setup.
- ```data.InitConformation``` is the LAMMPS data file with the initial chain conformation
- ```ExampleSubmission.joblst``` is a barebones example submission script. This will depend on the machine on which you run the simulation.
  
You will need to load the MPI implementation. Then, run a simulation with the following command:

```mpirun -n <NTASKS> ./ExtrSim <MaxCoh> <p_jump> <NumBeads> <NumMCSteps> <outputfileSuffix> <K_coh> <R0_coh> <MDtauperMC> <MDtauperdump> <WorkPath> <InitialData> <BoxSize> <InitialMDTau> <t_int> <bindsite> <r_c> <calcR2contfrequency> <printR2contfrequency> <ReeRgfrequency> <ring> <RandSeed> <p_on> <p_off> <stopfilename> <NumberOfStops>```

All but the last two arguments are required.

- ```<NTASKS>``` is the number of MPI tasks that you are using.
- ```<MaxCoh>``` is the maximum number of extruders allowed to be simultaneously bound to the polymer chain. Usually set this to a number 4-5 times the expected average number of bound cohesins.
- ```<p_jump>``` this is the probability of cohesin domains passing/traversing each other. 1 means unimpeded. 0 means fully impeded.
- ```<NumBeads>``` is the number of beads in the simulation
- ```<NumMCSteps>``` is the number of MC steps you want to run for
- ```<outputfileSuffix>``` is the suffix you want to be appended to all of your output files
- ```<K_coh>``` is the spring constant of the cohesin bond. Set it to 1.
- ```<R0_coh>``` is the max length of the cohesin bond. Set it to 4.
- ```<MDtauperMC>``` is the ratio of MD tau (Lennard Jones time) to MC steps. We usually use 5. This sets the average extrusion velocity.
- ```<MDtauperdump>``` is how frequently you want to dump out bead coordinates in units of MD Tau. Usually we use 1, 5, or 10.
- ```<WorkPath>``` is the path to the directory you want all of the output files to go
- ```<InitialData>``` is the initial data file. Use data.InitConformation
- ```<InitialMDTau>``` is the number of MD tau you want to run before starting extrusion.
- ```<t_int>``` is the integration time step for MD. Typically use 0.01.
- ```<bindsite>``` is the bead index of cohesin binding, if needed. Set it to 0 for nonspecific binding.
- ```<r_c>``` is the minimum distance between beads considered for contact, in units of MD simulation distance units (Lennard Jones sigma). Typically use something close to 1.5.
- ```<printR2contfrequency>``` is how frequently you want to print the contact matrix and internal mean square distances for the chain (units of MC steps).
- ```<ring>``` should be 0 if you are simulating a linear chain, 1 if you are simulationg a ring. Use 0.
- ```<RandSeed>``` is the seed used for the random number generator.
- ```<p_on>``` is the probability that a cohesin binds a bead during one MC timestep. This is the probability **per bead**
- ```<p_off>``` is the probability that a cohesin unbinds during one MC timestep. This is the probability **per bound cohesin**
- ```<stopfilename>``` is the name of the file that has the TAD anchors/stop signs. Not necessary if there are no anchors. The file ```stopfile.txt``` has an example with TAD anchor pairs at beads (1,200) and (500,700).
- ```<NumberOfStops>``` is the number of TAD anchor/stop sign pairs. Not necessary if there are no anchors.
--------
Outputs
----
This code will have several output files that will be in ```<WorkPath>```. Each will end in the ```<outputfileSuffix>``` of your choosing.
- ```Status....txt``` has some input information and provides status updates (when the simulation is 25%, 50%, 75% done).
- ```NumCoh....txt``` prints the number of cohesins bound to the chain every step
- ```Acc....txt``` has four columns. The first column is the number of attempted moves made by a left cohesin domain. The second column is the number of accepted moves. The 3rd and 4th are the same, but for a right cohesin domain. Each column is a **running sum**.
- ```Bind....txt``` has the MC step at which a cohesin bound to the chain
- ```FinalLife....txt``` has two columns. The first column is the residence time of a cohesin when it unbinds (in MC steps). The second column is the loop length when it unbinds
- ```coordsLEUnwrap....txt``` has the coordinates of the chain formatted as a LAMMPS dump file with unwrapped coordinates.
- ```Contacts....txt``` has the number of conformations in which two beads were within the cutoff distance r_c. First column is bead i, second column bead j, third column is number of contacts.
- ```R2....txt``` has the mean squared distance between beads i and j, organized like ```Contacts....txt```.
- ```CohLocs...txt``` has 2*T rows, where T is the number of MC steps. Each column shows the bead indices to which an extruder is bound. Ex., [1 10 50; 100 32 52] means three extruders were bound, one at beads 1-100, one at 10-32, and one at 50-52.
