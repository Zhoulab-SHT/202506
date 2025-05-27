# load packages
import numpy as np
import pandas as pd
import sys
import os
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit

import mdtraj
try:
    import nglview
except ImportError:
    print('Please install nglview to visualize molecules in the jupyter notebooks.')

sys.path.append('../../')
from openabc.forcefields.parsers import HPSParser
from openabc.forcefields import HPSModel
from openabc.utils.helper_functions import build_straight_CA_chain, write_pdb
from openabc.utils.insert import insert_molecules

# set simulation platform
platform_name = 'CUDA'

sequence = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS'
ca_lc_pdb = 'init_LC_CA.pdb'
ca_lc_atoms = build_straight_CA_chain(sequence, r0=0.38)
write_pdb(ca_lc_atoms, ca_lc_pdb)
sequence2 = 'peptide sequence'
ca_ara_pdb = 'init_ARA_CA.pdb'
ca_ara_atoms = build_straight_CA_chain(sequence2, r0=0.38)
write_pdb(ca_ara_atoms, ca_ara_pdb)
lc_parser = HPSParser(ca_lc_pdb)
ara_parser = HPSParser(ca_ara_pdb)

# insert molecules into the simulation box randomly
n1 = 30
n2 = 230
if not os.path.exists('start.pdb'):
    insert_molecules('init_LC_CA.pdb', 'tmp.pdb', n_mol=n1, box=[15, 15, 280])
    insert_molecules('init_ARA_CA.pdb', 'start.pdb',n_mol=n2, existing_pdb='tmp.pdb', box=[15,15,280])
protein = HPSModel()
for i in range(n1):
    protein.append_mol(lc_parser)
for i in range(n2):
    protein.append_mol(ara_parser)
top = app.PDBFile('start.pdb').getTopology()
init_coord = app.PDBFile('start.pdb').getPositions()
protein.create_system(top, box_a=15, box_b=15, box_c=280)

temperature = 150
end_temperature = 310
salt_conc = 10*unit.millimolar
protein.add_protein_bonds(force_group=1)
protein.add_contacts('Urry', mu=1, delta=0.08, force_group=2)
protein.add_dh_elec(force_group=3)
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
protein.set_simulation(integrator, platform_name, init_coord=init_coord)
simulation = protein.simulation
protein.simulation.minimizeEnergy()
output_interval = 200000
output_dcd = 'LC_slab.dcd'
protein.add_reporters(output_interval, output_dcd)
protein.simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
n_iterations = 100
for temperature_i in np.linspace(temperature, end_temperature, n_iterations):
    integrator.setTemperature(temperature_i*unit.kelvin)
    simulation.step(1000000)
protein.simulation.step(500000000) # 5us total time
integrator.setTemperature(end_temperature*unit.kelvin)
protein.simulation.context.setVelocitiesToTemperature(end_temperature*unit.kelvin)
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True,
                                    getParameters=True, enforcePeriodicBox=True)

with open('noNPT_simulation__state.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(state))

