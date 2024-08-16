import numpy as np
import pyemma
import pickle
import random
from tqdm import tqdm

class MSMAnalysis:
    def __init__(self, msm_path, traj_metrics_path, traj_path, totdist_path):
        self.msm = pickle.load(open(msm_path, 'rb'))
        self.weights = self.msm.trajectory_weights()
        self.traj = np.load(traj_metrics_path)
        self.traj2 = np.load(traj_path)
        self.totdist = np.load(totdist_path, allow_pickle=True)
        self.subslen = int(0.8 * len(self.traj))

    def compute_correlation_coefficients(self, iterations=200):
        co = []
        for _ in tqdm(range(iterations)):
            subs = random.sample(self.traj.tolist(), self.subslen)
            w, xs, ys = [], [], []
            for t in subs:
                loc = np.where(self.traj2 == t)[0][0]
                loc2 = np.where(self.traj == t)[0][0]
                if len(self.weights[loc]) == len(self.totdist[loc2][:, 0]):
                    w.append(self.weights[loc])
                    xs.append(self.totdist[loc2][:, 0])
                    ys.append(self.totdist[loc2][:, 1])
            w = np.concatenate(w)
            xs = np.concatenate(xs)
            ys = np.concatenate(ys)
            mxw = np.sum(np.multiply(w, xs)) / np.sum(w)
            myw = np.sum(np.multiply(w, ys)) / np.sum(w)
            xy = np.sum(w * (xs - mxw) * (ys - myw))
            xx = np.sum(w * (xs - mxw) ** 2)
            yy = np.sum(w * (ys - myw) ** 2)
            co.append(xy / np.sqrt(xx * yy))
        np.save('Corrcoeff_pi_cation_lock.npy', co)
        return co

    def compute_rho(self):
        txx = np.concatenate(self.totdist)
        x = txx[:, 0]
        y = txx[:, 1]
        return np.corrcoef(x, y)

def main():
    analysis = MSMAnalysis(
        msm_path='/home/pdb3/SMO/APO/Analysis/MSM/MSMobjs/MSMobj_9D_1200.py',
        traj_metrics_path='./dist_I_metrics.npy',
        traj_path='./dist_I.npy',
        totdist_path='/home/pdb3/SMO/APO/Analysis/totdistcombined_metrics.npy'
    )
    
    co = analysis.compute_correlation_coefficients()
    print(co)
    
    rho = analysis.compute_rho()
    print(rho)

if __name__ == '__main__':
    main()

