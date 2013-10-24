from __future__ import print_function
from os.path import exists, join, dirname
try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
except ImportError as err:
    print("Failed to import OpenMM packages:", err.message)
    print("Make sure OpenMM is installed and the library path is set correctly.")
    exit()


#*****************************************************************************#
#  Customize these lines to change the OpenMM simulation setup                #
#  See the documentation at https://simtk.org/api_docs/openmm/api5_2/python/  #
#  for details on the available options                                       #
#*****************************************************************************#

input_pdb = join(dirname(__file__), 'input.pdb')
pdb = PDBFile(input_pdb)
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 2.0*femtosecond)


# Now that the system is setup, write out all of the files to disk for tungsten

with open('system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))
    print('saved system.xml')
with open('integrator.xml', 'w') as f:
    f.write(XmlSerializer.serialize(integrator))
    print('saved integrator.xml')
context = Context(system, integrator)
context.setPositions(pdb.positions)
context.setVelocitiesToTemperature(300*kelvin)
state = context.getState(getPositions=True, getVelocities=True)
with open('state.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
    print('saved state.xml')
with open('AtomIndices.dat', 'w') as f:
    for atom in pdb.topology.atoms():
        if atom.name == 'CA':
            f.write('%d\n' % atom.index)
    print('saved AtomIndices.dat')

print("All done.")
