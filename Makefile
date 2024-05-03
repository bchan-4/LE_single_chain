CC=mpicc

%.o: %.c
  $(CC) -O3 -I/path/to/LAMMPS/src/ -c $< -lm

ExtrSim: ExtrSim.o
  $(CC) $< -o $@ -lm -L/path/to/LAMMPS/src/ -llammps_mpi

ExtrSim_MoreOutputs: ExtrSim_MoreOutputs.o
  $(CC) $< -o $@ -lm -L/path/to/LAMMPS/src/ -llammps_mpi
