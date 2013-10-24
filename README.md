Tungsten
========

Overview
--------
Tungsten is an application, based on the OpenMM molecular dynamics library,
for high performance, GPU-based, Markov state model accelerated molecular
simulation. Tungsten runs a set of independent MD simulations in parallel
over MPI. With a certain frequency, the simulations pause and the nodes
collaboratively build a Markov state model (see, for example Bowman et. al,
*J. Chem. Theory Comput.*, **2010**, 6 (3), 787). The conformational states
that have fewest incomming transitions are then scattered back to the nodes
for further simulation, ala. Weber and Pande's "min-counts sampling" (*J.
Chem. Theory Comput.*, **2011**, 7 (10), 3405.).

Usage
-----
`tungsten` requires 4 inputs files, specified on the command line. The first
thee are OpenMM xml files, `system.xml`, `integrator.xml`, and `state.xml`,
which specify the complete Hamiltonian, integrator, and initial conditions
for the simulation. These files can be constructed from an RCSB PDB, amber
.incprd/.prmtop, gromacs .gro/.top, or desmond .dms file using the OpenMM
python app. You can modify and use the script `examples/buildTungstenXML.py`,
included with the tungsten source.

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
kcentersRmsdIndicesFile = data/AtomIndices.dat
kcentersRmsdCutoff = 0.005
openmmPlatform = OpenCL
```

Execution
---------
When you run a command under MPI using `mpirun` (or `aprun`, `mpiexec`, or
another vendor specific MPI process manager), multiple copies of the
application are instantiated which are capable of communicating with each
other via the exchange of messages. The processes may be on the same physical
machine, or more likely, distributed across a cluster of nodes. Each MPI
process is referred to by its "rank", an index starting from zero. In
tungsten, each rank caries it own independent OpenMM context -- the force
computation and integration is not shared over MPI, and each rank propagates
an independent copy of the simulation (the ranks do engage in collective
communication during the clustering and Markov state model parameterization
stage, but this occupies a small fraction of the overall number of clock
cycles).

For optimal performance with tungsten, each rank should to correspond to a
single GPU. In general, you don't get a net speedup from putting two
simulations on the same GPU, since the GPU must timeshare between the two
contexts. There are generally a variety of flags you can set with `mpirun` (or
your alternative MPI process manager) to control how the ranks are bound to
specific physical nodes. These are importance, because inadvertently
scheduling all of your tungsten ranks on a single node/gpu and leaving your
other nodes unoccupied will (obviously) lead to a dramatic slowdown. When
tungsten starts up, each rank will print information about the node and device
its Context is bound to, which can help confirm that your run is set up
correctly.

In tungsten, each rank associates itself with a single output file that it
reads and writes to exclusively. These output files are AMBER compatible
NetCDF trajectory files (you can load them with VMD, etc), and numbered by the
rank in a `trj-%05d.nc` pattern (i.e. `trj-00001.nc`,` trj-00002.nc`, etc).
The trajectory file stores the positions, periodic box dimensions, and
simulation time. Velocities and forces are not stored. When a round of
adaptive sampling is completed and the min-counts clusters are identified,
each rank obtains a new conformation (generally from another rank), the
simulation time is reset to zero, new velocities are chosen from the Boltzmann
distribution, and simulation continues. The new coordinates are written to the
same trajectory file. When loading trajectory files, you can determine where
the "breaks" are by looking at the time dimension in the trajectory file, and
marking where the simulation clock is reset to zero. If your analysis code
doesn't show the time variable, you can always use the `ncdump` command line
application which is installed with the NetCDF C library, and can dump the
trajectory content to stdout (e.g. `ncdump -v time trj-00000.nc`).

Note that tungsten is not suitable for distributed or heterogeneous hardware;
each of the ranks will effectively slow down and wait for the slowest member
of the team, because each rank must finish the requested amount of MD
simulation per round. If one of the ranks is much slower than the rest, the
clustering sub-step will not really start until that rank has finished its
simulation. (it will start, the process will hang waiting for messages
from the straggler).

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
saving and loading trajectories to disk (tungsten writes AMBER
NetCDF-compatible trajectories).
- [MPI](http://en.wikipedia.org/wiki/Message_Passing_Interface) for message
passing and intra-processor coordination.

Tungsten's build system should be pretty good about finding OpenMM and netcdf
on your system, if they're installed. For MPI, it uses the `mpicc` and
`mpicxx` compiler wrappers, which should be installed by any competent MPI
distribution.

If you don't already have them, you can install netCDF and OpenMPI on an
Ubuntu box with `sudo apt-get install ibopenmpi-dev openmpi-bin
libnetcdf-dev`. It should be basically the same on any other platform with a
package manager.

Note: If you're building netCDF from source for this application, note that
tungsten does not require netCDF4 / HDF5, so you can pass `--disable-netcdf-4`
to netCDF's `configure` script if you like. This can make the netCDF build a
little easier.

License
-------
```
Tungsten is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

Tungsten is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details. You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software Foundation,
Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
```
