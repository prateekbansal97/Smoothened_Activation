import argparse
import mdtraj as md
import numpy as np
import itertools
import os

NAME = 'IonicLock'

class TrajectoryAnalysis:
    def __init__(self, input_file):
        self.file = input_file
        self.topI = self._select_topology()
        self.t = md.load(self.file, top=self.topI)
        self.c = self.t.xyz
        self.nFrames = self.t.n_frames
        self.distI = np.empty([self.nFrames, 7])

    def _select_topology(self):
        if '6XBL' not in self.file:
            return '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7'
        else:
            return '/home/pdb3/SMO/APO/Analysis/6XBL_Apo_HMass_stripped.parm7'

    def compute_classA_network(self):
        contact_resis = [281, 364, 391, 395, 477]
        contact_combs = list(itertools.combinations(contact_resis, 2))
        g = md.compute_contacts(self.t, contacts=contact_combs, scheme='closest-heavy', periodic=False)
        return g[0]

    def compute_ionic_lock(self):
        Wa = ['CD1', 'CZ2', 'CE3']
        Wn = 477
        W_coords = [np.concatenate(self.c[:, self.t.topology.select(f'name {atom} and resid {Wn}'), :]) for atom in Wa]

        Ra = ['NH1', 'NH2', 'NE', 'CZ']
        Rn = 393
        R_coords = [np.concatenate(self.c[:, self.t.topology.select(f'name {atom} and resid {Rn}'), :]) for atom in Ra]

        dist, angle = self._calculate_distances_and_angles(W_coords, R_coords)

        self._update_distI(dist, angle)
    
    def _calculate_distances_and_angles(self, W_coords, R_coords):
        dist, angle = [], []
        for i in range(R_coords[0].shape[0]):
            Rp = plane(R_coords[0][i], R_coords[1][i], R_coords[2][i])
            Wp = plane(W_coords[0][i], W_coords[1][i], W_coords[2][i])
            cent = np.mean(W_coords, axis=0)[i]
            dist.append(la.norm(cent - R_coords[3][i]) * 10)
            angle.append(acute_plane_angle(Rp, Wp) * 180 / np.pi)
        return np.array(dist), np.array(angle)

    def _update_distI(self, dist, angle):
        self.distI[:, 0] = dist
        self.distI[:, 1] = angle
        p = self.t.topology.select('backbone and resid 267')
        di = md.compute_dihedrals(self.t, [p])
        self.distI[:, 2] = np.concatenate(di)
        dihed = md.compute_phi(self.t)
        dihed2 = md.compute_psi(self.t)
        self.distI[:, 3] = dihed[1][:, 266]
        self.distI[:, 4] = dihed2[1][:, 266]
        self.distI[:, 5] = dihed2[1][:, 267]
        self.distI[:, 6] = dihed[1][:, 267]

    def save_metrics(self):
        filename = os.path.basename(self.file).replace('.dcd', '')
        np.save(f"/home/pdb3/SMO/APO/Analysis/npy/{filename}_{NAME}_metrics.npy", self.distI)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', type=str, required=True)
    args = parser.parse_args()

    analysis = TrajectoryAnalysis(args.input)
    analysis.compute_ionic_lock()
    analysis.save_metrics()

if __name__ == '__main__':
    main()

