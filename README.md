# Fiber Bragg Grating Monte Carlo

A Python implementation of the Monte Carlo simulation for fiber Bragg gratings from the paper "Modeling Light Propagation Through Bragg Reflection Gratings With an Adapted Monte Carlo Method" by Ben Pelleg and Adam Fontecchio.

The simulation is only 1D and currently uses a simplified Monte Carlo method where photons "jump" between interfaces rather than simulating over small step size. It was written this way since it does not include absorption or scattering, but the code could be modified fairly easily to include these.
