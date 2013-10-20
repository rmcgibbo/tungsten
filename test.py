import mdtraj as md
t1 = md.load('data/trj-00000.nc', top='data/input.pdb')
t2 = md.load('data/trj-00001.nc', top='data/input.pdb')
t = t1 + t2

r = md.rmsd_cache(t, 'atom')

print 'RMSDs\n', r.rmsds_to(r, 0)
print 'G\n', r._g[:len(t1)]
print r._g[len(t1):]

