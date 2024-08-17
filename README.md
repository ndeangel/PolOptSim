# POLAR-2 Polarimeter Optical Simulation

If you find this code/work, please cite [De Angelis N. et al. Optimizing the light output of a plastic scintillator and SiPM based detector through optical characterization and simulation: A case study for POLAR-2. JINST 2024](https://arxiv.org/abs/2407.10741)

## How to compile and run the G4 simulation

* cd path/to/PolOptSim/build
  - Go to the 'build' directory of the 'PolOptSim' simulations

* ../g4cmake.sh ..
  - Inside the 'build' directory, you should run 'cmake ../'. This will automatically setup everything for you regarding Geant4, ROOT, scripts, etc. It will produce a Makefile.
    It will also create a 'bin' (for all binaries) and a 'data/G4out' (for simulation files) directory.

* make -jN
  - Compile by replacing N by the number of cores you want to use.

* ./OptSim gps.mac run.mac -c config.mac -o output.root
  - Run the simulations. To run with visualization, use the '-g vis.mac' option.


## Version note

All the work described in [De Angelis Nicolas. Development of the Next Generation Space-based Compton Polarimeter and Energy Resolved Polarization Analysis of Gamma-Ray Bursts Prompt Emission. PhD Thesis 2023.](https://doi.org/10.13097/archive-ouverte/unige:173869) and [De Angelis N. et al. Optimizing the light output of a plastic scintillator and SiPM based detector through optical characterization and simulation: A case study for POLAR-2. JINST 2024](https://arxiv.org/abs/2407.10741) has been performed with Geant4 version 10.6.1 and ROOT version 6.18/05.


