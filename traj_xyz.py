from ase.io import read
from ase.io import write
import ase
# wrap and save traj in .xyz --- the .traj is a non human readable database file, xyz is much better
out_traj = ase.io.read('waterdyn.traj', ':')
for at in out_traj:
    at.wrap()
    if 'momenta' in at.arrays: del at.arrays['momenta']
ase.io.write('waterdyn_fromtraj.xyz', out_traj, 'xyz')
