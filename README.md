InjectionModelWithGelation contains MATLAB code for a research-driven needle injection model containing a drug molecule of interest and a phase-separating agent that gels in the simulation domain as the simulation progresses. The model consists of the main function and three underlying functions:

main.m is the main file that handles parameter initialization, execution of all physics, and plotting of results.

continuumMechanics.m uses the theory provided by Netti et al., 2003, to calculate key continuum mechanics (CM) outputs, including the dilatation, displacement, fluid velocity and fluid velocity gradient. It does so by inputting the analytical results of CM theory in Laplace space and using an inverse Laplace function to compute the results in time. The analytical results in Laplace space are given by the following functions: getDilatationLaplace.m, getDisplacementLaplace.m, getVelocityLaplace.m and getVelocityGradLaplace.m.

cahnHilliard.m computes the order, or the composition, of the phase-separating agent in space and time using Cahn Hilliard (CH) theory and finite differences in spherical coordinates. It is coupled to the CM section via the radius of the cavity that forms during the injection, which is calculated as the first spatial point in the solid displacement field plus the initial cavity radius. This cavity forms a boundary condition that is key to the CH simulations. 

massTransport.m calculates the concentration of the drug using the convection-diffusion equation and mass transport (MT) theory and finite differences in spherical coordinates. It takes as input the outputs of the CM and CH simulations. The radius of the cavity (from CM) sets the moving boundary for MT, and CH affects the diffusion rate of drug in space. We assume that the diffusion coefficient of the drug changes linearly with respect to the order that CH outputs, in other words, given a maximum and minimum diffusion coefficient corresponding to order = 0 and order = 1, the local diffusion coefficient is given by a linear relationship with the local composition of the phase-separating or gelating agent. 

In solving and using these equations, we often pass the outputs of each as splines. Due to the adaptive temporal sampling for performance reasons, MT, CH and CM may have time vectors of different sizes, so we pass through spline instead that we can use to interpolate the value of any metric of interest (displacement, velocity, order, concentration, etc) at any relevant point in space and time.

For the spatial vector, we use r_mass for CH and MT since they require simulation of what happens inside the needle, but use the subset r for MT, since this theory does not care about what happens inside the needle. We also use relevant_r_bound to increase resolution within a small area around the needle tip. We do not increase it throughout the domain to save on computational power.

Daniel Adrianzen, 10/07/2024