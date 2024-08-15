import argparse
import mdtraj as md
import numpy as np
import pickle
import os

class DistanceCalculator:
    def __init__(self, input_file, features_path, output_dir):
        self.file = input_file
        self.features = pickle.load(open(features_path, 'rb'))
        self.dists = self.features['Distancestunnel']
        self.dist_1, self.dist_2 = np.array_split(self.dists, 2)
        self.topology_file = self._select_topology()
        self.traj = md.load(self.file, top=self.topology_file)
        self.nFrames = self.traj.n_frames
        self.output_dir = output_dir

    def _select_topology(self):
        if '6XBL' not in self.file:
            return '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7'
        else:
            return '/home/pdb3/SMO/APO/Analysis/6XBL_Apo_HMass_stripped.parm7'

    def calculate_distances(self):
        distI = np.empty([self.nFrames, len(self.dist_1)])
        distInorm = np.empty([self.nFrames, len(self.dist_1)])

        for k, l, p in zip(self.dist_1, self.dist_2, range(len(self.dist_1))):
            atom_name_1, atom_no_1 = k.split('_')
            atom_name_2, atom_no_2 = l.split('_')
            atom_1 = self.traj.topology.select(f'name {atom_name_1} and resid {atom_no_1}')
            atom_2 = self.traj.topology.select(f'name {atom_name_2} and resid {atom_no_2}')
            dist_atoms = md.compute_distances(self.traj, [[atom_1[0], atom_2[0]]], periodic=False)
            distI[:, p] = dist_atoms[:, 0]
            distInorm[:, p] = (dist_atoms[:, 0] - np.mean(dist_atoms[:, 0])) / np.std(dist_atoms[:, 0])

        self._save_results(distI, distInorm)

    def _save_results(self, distI, distInorm):
        filename = os.path.basename(self.file).replace('.dcd', '')
        np.save(os.path.join(self.output_dir, f'{filename}_extra_tunnel.npy'), distI)
        np.save(os.path.join(self.output_dir, f'{filename}_norm_extra_tunnel.npy'), distInorm)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', type=str, required=True)
    args = parser.parse_args()

    calculator = DistanceCalculator(
        input_file=args.input,
        features_path='/home/pdb3/SMO/APO/Analysis/pkl/features.pkl',
        output_dir='/home/pdb3/SMO/APO/Analysis/npy/'
    )
    calculator.calculate_distances()

if __name__ == '__main__':
    main()

