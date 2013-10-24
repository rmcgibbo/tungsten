Tungsten
========

Overview
--------
Tungsten is an application, based on the OpenMM molecular dynamics library,
for Markov state model accelerated molecular simulation. Tungsten runs a set
of independent MD simulations in parallel over MPI. With a certain frequency,
the simulations pause and the nodes collaborativly build a Markov state model
(see, for example Bowman et. al, *J. Chem. Theory Comput.*, **2010**, 6 (3),
787). The conformational states that have fewest incomming transitions are
then scattered back to the nodes for further simulation, ala. Weber and
Pande's "min-counts sampling" (*J. Chem. Theory Comput.*, **2011**, 7 (10),
3405.).

Usage
-----
`tungsten` requires 4 inputs, specified on the command line. The first thee
are OpenMM xml files, `system.xml`, `integrator.xml`, and `state.xml`, which
specify the complete Hamiltonian, integrator, and initial conditions for
the simulation. These files can be constructed from an RCSB PDB, amber
incprd/prmtop, gromacs gro/top, or desmond dms file using the OpenMM python
app.

`tungsten` also requires a short configuration file, in a simple key=value
format to specify, among other things, the frequency of clustering, the
clustering RMSD cutoff, and the OpenMM platform to use (e.g. CUDA, OpenCL,
etc). With these inputs, you can run the suite of simulations over `N`
parallel processors using

```mpirun -np <N> tungsten system.xml integrator.xml state.xml config.ini```

The configuration options for the ini file (last command line option)
are available:

```
numRounds:int                   number of adaptive sampling rounds
numStepsPerRound:int            number of MD steps in each round, between clusterings
numStepsPerWrite:int            trajectory save frequency, in steps
outputRootPath:path             directory where output trajectories are saved
kcentersRmsdIndicesFile:path    path to a space or newline separated file listing the
                                zero-based indices of the atoms to use in the rmsd
                                computation, for clustering.
kcentersRmsdCutoff:float        RMSD distance cutoff for clustering, in units of
                                nanometers. A lower cutoff will result in more
                                numerous, smaller conformational states.
openmmPlatform:string           OpenMM platform to run the simulations. May be one of
                                CUDA, OpenCL, Reference, etc.
```

For example, the following is a valid tungsten config

```
; tungsten config file.
# a line is considered a comment if it starts with a
# semicolon or pound sign / hash.
numRounds = 2
numStepsPerRound = 500
numStepsPerWrite = 50
outputRootPath = data/
kcentersRmsdIndicesFile = data/AtomIndices.indx
kcentersRmsdCutoff = 0.005
openmmPlatform = OpenCL
```

Installation
------------
Tungsten is a single binary, built with cmake build system. If you don't
already have it, [get cmake](http://www.cmake.org/cmake/resources/software.html).
Then, configure the build using `ccmake .`, generate the makefile, and run
`make`. You should now have the `tungsten` binary, which you can put where
ever you like.

Dependencies
------------
Tungsten has three dependencies.
- [OpenMM](https://simtk.org/home/openmm) for running simulations.
- [netCDF](http://www.unidata.ucar.edu/software/netcdf/docs/index.html) for
saving and loading trajectories to disk (tungsten writes AMBER NetCDF-compatible
trajectories).
- [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface) for message
passing and intra-processor coordination.

Tungsten's build system should be pretty good about finding OpenMM and netcdf
on your system, if they're installed. For MPI, it uses the `mpicc` and `mpicxx`
compiler wrappers, which should be installed by any competent MPI distribution.

If you don't already have them, you can install netCDF and OpenMPI on an ubuntu
box with 

```
sudo apt-get install ibopenmpi-dev openmpi-bin libnetcdf-dev
```

Note: If you're building netCDF from source for this application, note that
tungsten does not require netCDF4 / HDF5, so you can pass `--disable-netcdf-4`
to netCDF's `configure` script if you like. This can make the netCDF build a
little easier.

License
-------
```
Tungsten is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

Tungsten is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details. You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
```
