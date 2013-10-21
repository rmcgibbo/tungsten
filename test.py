import mdtraj as md
import numpy as np

atom_indices = np.loadtxt('data/AtomIndices.dat', int)
t1 = md.load('data/trj-00000.nc', top='data/input.pdb', atom_indices=atom_indices)
t2 = md.load('data/trj-00001.nc', top='data/input.pdb', atom_indices=atom_indices)
t = t1 + t2

r = md.rmsd_cache(t, 'axis')

print "Coords"
print r.cords[0]
print 'RMSDs\n', r.rmsds_to(r, 0)
print 'G\n', r._g[:len(t1)]
print r._g[len(t1):]

