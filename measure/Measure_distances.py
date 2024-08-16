#Original name Measurements.py
import glob
import os
import argparse
import mdtraj as md
import numpy as np
from tqdm import tqdm
import pickle

class DistanceCalculator:
    def __init__(self, input_file):
        self.file = input_file
        self.features = pickle.load(open('/home/pdb3/SMO/APO/Analysis/pkl/features.pkl', 'rb'))
        self.dist_1, self.dist_2 = np.array_split(self.features['Distances'], 2)
        self.topology_file = self._select_topology_file()
        self.traj = md.load(self.file, top=self.topology_file)
        self.n_frames = self.traj.n_frames

    def _select_topology_file(self):
        return '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7' if '6XBL' not in self.file else '/home/pdb3/SMO/APO/Analysis/6XBL_Apo_HMass_stripped.parm7'

    def calculate_distances(self):
        distI = np.empty([self.n_frames, len(self.dist_1)])
        distInorm = np.empty([self.n_frames, len(self.dist_1)])
        
        for k, l, p in zip(self.dist_1, self.dist_2, range(len(self.dist_1))):
            atom_name_1, atom_no_1 = k.split('_')
            atom_name_2, atom_no_2 = l.split('_')
            
            atom_1 = self.traj.topology.select(f'name {atom_name_1} and resid {atom_no_1}')
            atom_2 = self.traj.topology.select(f'name {atom_name_2} and resid {atom_no_2}')
            
            dist_atoms = md.compute_distances(self.traj, [[atom_1[0], atom_2[0]]], periodic=False)
            print(dist_atoms, dist_atoms.shape if p < 1 else '')
            
            distI[:, [p]] = dist_atoms
            distInorm[:, [p]] = (dist_atoms - np.mean(np.concatenate(dist_atoms))) / (np.std(np.concatenate(dist_atoms)))
        
        return distI, distInorm

    def save_results(self, distI, distInorm):
        filename = os.path.basename(self.file).split(".dcd")[0]
        np.save(f"/home/pdb3/SMO/APO/Analysis/npy/{filename}_extra.npy", distI)
        np.save(f"/home/pdb3/SMO/APO/Analysis/npy_normalized/{filename}_norm_extra.npy", distInorm)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', type=str, required=True)
    args = parser.parse_args()

    distance_calculator = DistanceCalculator(args.input)
    distI, distInorm = distance_calculator.calculate_distances()
    distance_calculator.save_results(distI, distInorm)

if __name__ == '__main__':
    main()

