# MicroWorld
In this repository we brought together the most popular singularity solutions of Stokes’
equation and created a program that can plot those. 
We presented the fundamental solution for free space and incompressible fluid, called Stokeslet, that is Green’s Function. 
Later on we took a closer look at its derivatives, that also satisfies our differential equation. 
Those derivatives ale called singularity solutions and in this work we focused at dipole, Rotlet,
Stresslet, source, sink and source doublet (source dipole). We also put in biological context to describe how bifalgellated
algae and a spermatozoon swims.
We also described flow for a point force next to interface of two liquids. As we relay on the
ratio of two fluids we can describe limiting cases as one of the fluids is acting as a free surface
or a solid wall. To solve this kind of problem with boundary conditions we used method of
images.
Software MicroWorld was written in Python and its most useful libraries are; NumPy, SciPy
and Matplotlib. In the program are two possible ways to plot vector field, the first one plot
is fully based on Matplotlib and function streamplot that creates streamlines automatically.
The second one is using integration with SciPy to locally calculate streamlines positions. This
data can be exported and type of its extension whas choosen with thought of future development.
We started creating visualisation in 3D for Stokeslet with library k3d and using Jupyter Notebook.

We highly encourage to write with any questions about the current version or future development to the following e-mail address g.niechwiado@student.uw.edu.pl .

As examples are shown singularities but also film of swimming algae.

Work was created at University of Warsaw.
