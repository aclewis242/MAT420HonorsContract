### *Written report included*

This project simulates the trajectory of a beam of light inside a rotating cylinder, where the index of refraction varies according to the behavior of the compressible fluid inside it.

Each of the three simulations -- rotating cylinder, stationary cylinder, and simple medium-to-medium refraction -- has a designated section in `main.py`. Computational parameters (i.e. maximum time, step size) are defined at the beginning of the first section (for the rotating cylinder). Physical parameters (defining the cylinder's properties) can be provided as arguments in the cylinder's initialisation, but it defaults to a set of base values defined in `cyl.py` itself.

Required packages:
- `numpy`
- `matplotlib`
- `math`
- `PIL`