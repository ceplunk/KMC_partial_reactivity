# KMC_partial_reactivity

All files for the kinetic Monte Carlo simulations for the work on diffusion above a plane with partially reactive patches. The three primary simulation routines are:
1. KMC_robin_shapes.m

   Calculate the splitting probability for being absorbed or not by a single patch
   
2. KMC_robin_full.m

   Calculate the time of absorption by a plane with an infinite number of partially-reactive circular patches on a square unit grid
   
3. KMC_local_patch.m

   Calculate the local time accumulated at a circular patch
