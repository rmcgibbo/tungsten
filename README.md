Tungsten
========

Installation
------------
Tungsten is a single binary, built with cmake build system. If you don't already
have it, [get cmake](http://www.cmake.org/cmake/resources/software.html). Then,
configure the build using `ccmake .`, generate the makefile, and run `make`. You
should now have the `tungsten` binary, which you can put whereever you like.

Dependencies
------------
Tungsten has three dependencies.
- [OpenMM](https://simtk.org/home/openmm) for running simulations.
- [netCDF](http://www.unidata.ucar.edu/software/netcdf/docs/index.html) for saving and
loading trajectories to disk (tungsten writes AMBER NetCDF-compatible trajectories).
- [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface) for message passing and
intra-processor coordination.

Tungsten's build system should be pretty good about finding OpenMM and netcdf on your system,
if they're installed. For MPI, it uses the `mpicc` and `mpicxx` compiler wrappers, which should
be installed by any competent MPI distribution.

If you don't already have them, you can install netCDF and OpenMPI on an ubuntu box with 

```
sudo apt-get install ibopenmpi-dev openmpi-bin libnetcdf-dev
```
