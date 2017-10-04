# Fiber Bragg Grating Monte Carlo

A Python implementation of the Monte Carlo simulation for fiber Bragg gratings from the paper "Modeling Light Propagation Through Bragg Reflection Gratings With an Adapted Monte Carlo Method" by Ben Pelleg and Adam Fontecchio.

The simulation is only 1D and currently uses a simplified Monte Carlo method where photons "jump" between interfaces rather than simulating over small step size. It was written this way since it does not include absorption or scattering, but the code could be modified fairly easily to include these.

**Note**: This code was written in a very short amount of time for a course project, so it's only a simple working example and isn't modular. But since the paper didn't release any code afaik, hopefully this can be used as a starting point. Any pull requests to clean up the code are welcome :)
