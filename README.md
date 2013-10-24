Tungsten
========

Overview
--------

Usage
-----

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

Note: If you're building netCDF from source for this application, note that tungsten does not
require netCDF4/HDF5, so you can pass `--disable-netcdf-4` to netCDF's `configure` script
if you like. This can make the build a little easier.

License
-------
Tungsten is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser
General Public License as published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

Tungsten is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
Public License for more details. You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software Foundation, Inc., 51 Franklin St,
Fifth Floor, Boston, MA 02110-1301 USA.
