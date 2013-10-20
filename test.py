import mdtraj as md
t = md.load('trj-00000.nc', top='input.pdb')
t2 = md.load('trj-00001.nc', top='input.pdb')
t = t + t2

r = md.rmsd_cache(t, 'atom')

print 'RMSDs\n', r.rmsds_to(r, 0)
print 'G\n', r._g
print 'Coords0\n', r.cords[0, 0]
print 'Coords1\n', r.cords[1, 0]
