from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm

pdb = app.PDBFile('input.pdb')
ff = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, constraints=app.HBonds)
integrator = mm.LangevinIntegrator(300, 1, 0.002)

with open('system.xml', 'w') as f:
    print('Saving system.xml')
    f.write(mm.XmlSerializer.serialize(system))

with open('integrator.xml', 'w') as f:
    print('Saving integrator.xml')
    f.write(mm.XmlSerializer.serialize(integrator))

context = mm.Context(system, integrator, mm.Platform.getPlatformByName('Reference'))
context.setPositions(pdb.positions)
state = context.getState(getPositions=True, getVelocities=True)

with open('state.xml', 'w') as f:
    print('Saving state.xml')
    f.write(mm.XmlSerializer.serialize(state))
