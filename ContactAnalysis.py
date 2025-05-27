import numpy as np
from MDAnalysis.analysis.distances import distance_array
import mdtraj
import MDAnalysis as mdana
from sklearn.metrics.pairwise import euclidean_distances
from tqdm import tqdm
import matplotlib.pyplot  as plt
import seaborn as sns
traj = mdtraj.load_dcd('LC_slab.dcd',top='start.pdb')
n_frames = traj.n_frames
n_atoms = traj.n_atoms
n_monomers_Mol1 =
n_monomers_Mol2 = 
n_atoms_per_monomer_Mol1 = 163
n_atoms_per_monomer_Mol2 = 40

contactmap_1 = np.zeros(163*163).reshape(163,163)
contactmap_2 = np.zeros(163*40).reshape(163,40)
interact_num_11,interact_num_12 = 0,0
intact_coord = np.zeros((n_frames, n_atoms, 3))
for i in tqdm(n_frames/2,n_frames)):
    for j in range(n_monomers_Mol1):
        for n in range(j+1,n_monomers_Mol1):
            coord_LC1 = traj.xyz[i,j*n_atoms_per_monomer_Mol1:(j+1)*n_atoms_per_monomer_Mol1].copy()
            coord_LC2 = traj.xyz[i,n*n_atoms_per_monomer_Mol1:(n+1)*n_atoms_per_monomer_Mol1].copy()
            pairwise_distance = euclidean_distances(coord_LC1, coord_LC2)
            normalized_distance = np.where(pairwise_distance < 2.0, 1, 0)
            interact_num_11 += np.count_nonzero(normalized_distance)
            contactmap_1 = np.add(contactmap_1,normalized_distance)

for i in tqdm(n_frames/2,n_frames)):
    for j in range(n_monomers_Mol1):
        for n in range(n_monomers_Mol2):
            coord_Mol1 = traj.xyz[i,j*n_atoms_per_monomer_Mol1:(j+1)*n_atoms_per_monomer_Mol1].copy()
            coord_Mol2 = traj.xyz[i,(n_monomers_Mol1*n_atoms_per_monomer_Mol1+n*n_atoms_per_monomer_Mol2):(n_monomers_Mol1*n_atoms_per_monomer_Mol1+(n+1)*n_atoms_per_monomer_Mol2)].copy()
            pairwise_distance = euclidean_distances(coord_Mol1, coord_Mol2)
            normalized_distance2 = np.where(pairwise_distance < 2.0, 1, 0)
            interact_num_12 += np.count_nonzero(normalized_distance2)
            contactmap_2 = np.add(contactmap_2,normalized_distance2)
Norm_ContactMap_11 = np.divide(contactmap_1,float(n_frames))
Norm_ContactMap_12 = np.divide(contactmap_2,float(n_frames))
Lij = interact_num_12*(n_monomers_Mol1+n_monomers_Mol2)/(interact_num_12+interact_num_11)/n_monomers_Mol2
print('interaction number of i-i is ',interact_num_11)
print('interaction number of i-j is ',interact_num_12)
print('interaction parameter is',Lij)

plt.figure(figsize=(20, 20))
sns.heatmap(Norm_ContactMap_11,cmap="YlGnBu")
np.save('contactmap11_2.0',Norm_ContactMap_11)
plt.savefig( "ContactMap_11-2.0.png")
sns.heatmap(Norm_ContactMap_12,cmap="YlGnBu")
plt.savefig("ContactMap_12-2.0.png")
np.save('contactmap12_2.0',Norm_ContactMap_12)

