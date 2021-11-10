"""Demonstrates molecular dynamics with constant energy."""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from deepmd.calculator import DP
from ase import Atoms
#from ase.calculators.psi4 import Psi4
import numpy as np
from ase.build import molecule
from xtb.ase.calculator import XTB

from ase.io import Trajectory
from ase.io import write
from ase.io.xyz import write_xyz
from ase.io.extxyz import write_extxyz

#traj = Trajectory('waterdyn.traj')

# Use Asap for a huge performance increase if it is installed
use_asap = False

#loading data
"""traj_load = np.load('atom_traj.npy')
vel_load = np.load('atom_vel.npy')"""

def printenergy(a,label=""):
    """Function to print the potential, kinetic and total energy"""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print(label,'Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

if use_asap:
    from asap3 import EMT
    size = 10
else:
    from ase.calculators.emt import EMT
    size = 3

# Set up a crystal
atoms = Atoms('H2O',positions=[(0.7601, 1.9270, 1),(1.9575, 1, 1),(1., 1., 1.)],cell=[100,100,100])

# Describe the interatomic interactions with the Effective Medium Theory
#atoms.calc = EMT()
#calc_slow = Psi4(method='HF', memory='500MB',basis='sto-3g')
#calc_fast = Psi4(method='b3lyp', memory='500MB',basis='6-311g_d_p_')
calc_fast = DP(model="/global/homes/y/yzamora/deepmd-kit/examples/water/se_e2_a/graph.pb")
calc_slow = XTB(method="GFN2-xTB")

atoms.calc = calc_slow

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 1 * units.fs, trajectory='waterdyn.traj')  # 5 fs time step.
#dyn = VelocityVerlet(atoms, 1 * units.fs)  # 5 fs time step.
"""for a in traj:
    printenergy(a,"Reference")
    save_calc = a.calc
    a.calc = calc_fast
    printenergy(a, "ML")
"""

"""nsteps = 100
x0 = atoms.get_positions()
traj_x = np.empty((nsteps, *x0.shape))
traj_x[0] = x0

v0 = atoms.get_velocities()
vel_x = np.empty((nsteps,*x0.shape))
vel_x[0] = v0"""


nsteps = 2
# Now run the dynamics
printenergy(atoms)
"""for i in range(nsteps-1):
    dyn.run(100)
    #write_extxyz('dyn.extxyz',atoms)
    #traj_x[i+1] = atoms.get_positions()
    #vel_x[i+1] = atoms.get_velocities()
    #import pdb;pdb.set_trace()
    if i%2 ==0:
        atoms.calc = calc_fast
    else:
        atoms.calc = calc_slow
    if i%10:
        printenergy(atoms)"""
#print(atoms.get_forces())
#import pdb; pdb.set_trace()
dyn.run(100)
#write('dyn_write.xyz',atoms)
#write_xyz('dyn.xyz',atoms)
"""np.save('atom_traj.npy',traj_x)
np.save('atom_vel.npy',vel_x)"""
