Compton Scattering and Pi0 Photoproduction Event Generator
Designed for use with the MAMI A2 Geant4 Simulation

Author - P. Martel
Version - 21 June 2011

Relevant files:
   BaseGen.h
   BasePart.h
   ComptonGen.h
   EventGen.cxx
   Makefile
   out/
   par/Compton*
   par/Pi0_maid*
   par/Pi0_said*
   physics.h
   Pi0PhotGen.h
   readme.txt

1) To compile, type 'make' from within this main directory (EventGen).

2) To run, type './EventGen' from within this main directory, and follow the prompts to select the desired conditions.

3) Possible selections are:
   a) Process - Compton, Pi0

   b) Type - Normal, Incoherent, Coherent

   c) Database:
      i) For Compton - Dispersion (Pasquini)
      ii) For Pi0 - MAID, SAID

   d) Target:
      i) For Normal - proton
      ii) For Incoherent - 12C (16O in progress)
      iii) For Coherent - 3He, 4He, 12C, 16O

   e) Polarization Settings
      i) Beam Pol -1.0 to 1.0 (0.7 default)
      ii) Target Pol -1.0 to 1.0 (0.8 default for normal, 0 for incoherent/coherent)

   f) Beam energy range
      i) For Compton - 100-450 MeV, in 5 MeV steps
      ii) For Pi0 - 145-450 MeV, in 5 MeV steps

   g) Number of events to run

4) For Normal and Incoherent reactions, parameter files are read in from the 'par' sub-directory, which is what causes the restriction on the energy step.
   (NOTE: If a finer step size is required, please contact me and I can provide the appropriate files)

   a) Compton files layout:
      i) Angle (lab)
      ii) Unpolarized CS (nb)
      iii) Positive CS (in direction of target pol, right hel beam)
                       (opposite dir of target pol, left hel beam)
      iv) Negative CS (in direction of target pol, left hel beam)
                      (opposite dir of target pol, right hel beam)

   b) Pi0 files layout (Observables nomenclature):
      i) Angle (cm)
      ii) DSG (ub)
      iii) T
      iv) E
      v) F

5) For Coherent reactions, cross sections produced using a method provided by R. Miskimen.

6) Two root files are produced as output in the 'out' sub-directory. 'hist.root' contains histograms for testing, energy and angular distributions, etc. 'ntpl.root' is then the ntuple for use in the A2 Geant4 simulation.


Please feel free to let me know of any errors, concerns, suggestions, or requests. Thanks!

-Phil
