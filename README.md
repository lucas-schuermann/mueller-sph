`mueller-sph` is a concise 2D implementation of Müller's "Particle-Based Fluid Simulation for Interactive Applications" (SPH) [paper](https://matthias-research.github.io/pages/publications/sca03.pdf) in C++. See [here](https://github.com/cerrno/mueller-sph-rs) for a more recent implementation in Rust including a parallel solver.

Please see the accompanying [tutorial](https://lucasschuermann.com/writing/implementing-sph-in-2d) for more information.

## Running
```bash
# install dependencies (debian/ubuntu)
apt install libopengl-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev libeigen3-dev

# uncomment header in `Makefile` depending on platform
make

# launch demo
./sph
```

## License
This project is distributed under the [MIT license](LICENSE.md).