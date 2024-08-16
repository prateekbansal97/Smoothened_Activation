import numpy as np
import matplotlib.pyplot as plt
import pyemma
import matplotlib as mpl
import pandas as pd
from tqdm import tqdm
import pickle

class FreeEnergyPlotter:
    def __init__(self, msm_path, df_path, data_path, features1, features2):
        self.msm = pickle.load(open(msm_path, 'rb'))
        self.weights = self.msm.trajectory_weights()
        self.df = pd.read_csv(df_path)
        self.features1 = features1
        self.features2 = features2
        self.yc = np.concatenate(np.load(data_path, allow_pickle=True))

    def BWfrommdt(self, resi):
        BWe = []
        for res in resi:
            r, rn = res.split('_')
            loc = np.where(self.df.mdtraj_num_active == int(rn))[0][0]
            BWe.append(f'{self.df.iloc[loc].Residue}{self.df.iloc[loc].Resno}_{self.df.iloc[loc].Location}.{int(self.df.iloc[loc].BW)}')
        return BWe

    def free_energy_plot(self, x_data, y_data, rs, bws):
        r1, r2, r3, r4 = rs
        bw1, bw2, bw3, bw4 = bws
        R, T = 0.001987, 310
        fig_wid, fig_hig = 9, 6
        cmap, Max_energy = mpl.cm.jet, 3.5
        x_bins, y_bins = 300, 300
        x_data_min, y_data_min = np.min(x_data), np.min(y_data)
        x_data_max, y_data_max = np.max(x_data), np.max(y_data)
        x_hist_lim_low, y_hist_lim_low = x_data_min - 0.5, y_data_min - 0.5
        x_hist_lim_high, y_hist_lim_high = x_data_max + 0.5, y_data_max + 0.5
        xspace, yspace = abs(x_data_min - x_data_max) / 10, abs(y_data_min - y_data_max) / 10
        x_lim_low, x_lim_high = x_data_min - xspace, x_data_max + xspace
        y_lim_low, y_lim_high = y_data_min - yspace, y_data_max + yspace

        hist = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins],
                              range=[[x_hist_lim_low, x_hist_lim_high], [y_hist_lim_low, y_hist_lim_high]],
                              density=True, weights=None)
        prob_density = hist[0]
        xedge, yedge = hist[1], hist[2]
        x_bin_size, y_bin_size = xedge[1] - xedge[0], yedge[1] - yedge[0]
        free_energy = -R * T * np.log(prob_density * x_bin_size * y_bin_size)
        delta_free_energy = free_energy - np.min(free_energy)
        xx = [(xedge[i] + xedge[i + 1]) / 2 for i in range(len(xedge) - 1)]
        yy = [(yedge[i] + yedge[i + 1]) / 2 for i in range(len(yedge) - 1)]

        fig, axs = plt.subplots(1, 1, figsize=(fig_wid, fig_hig))
        axs.contour(xx, yy, delta_free_energy.T, levels=6, colors='k', linewidths=0.2)
        cd = axs.contourf(xx, yy, delta_free_energy.T, np.linspace(0, Max_energy, 30), vmin=0.0, vmax=Max_energy, cmap=cmap)
        fig.colorbar(cd, ticks=range(int(Max_energy) + 1))
        axs.set_xlim(x_lim_low, x_lim_high)
        axs.set_ylim(y_lim_low, y_lim_high)
        axs.set_xticks(range(int(x_lim_low), int(x_lim_high), 2))
        axs.set_xticklabels(range(int(x_lim_low), int(x_lim_high), 2))
        axs.set_yticks(np.around(np.arange(y_lim_low, y_lim_high, 2), 2))
        axs.set_yticklabels(np.around(np.arange(y_lim_low, y_lim_high, 2), 2))
        plt.ylabel(rf'Distance {r1}$^{{{bw1}}}$ - - {r2}$^{{{bw2}}}$ (Å)', fontsize=18)
        plt.xlabel(rf'Distance {r3}$^{{{bw3}}}$ - - {r4}$^{{{bw4}}}$ (Å)', fontsize=18)
        plt.grid(ls='-', lw=0.2)
        plt.savefig(f'Distance_{r3}--{r4}_vs_{r1}--{r2}.png', dpi=500)
        plt.close('all')

    def plot_all(self):
        t1 = self.BWfrommdt(self.features1)
        t2 = self.BWfrommdt(self.features2)
        for k in tqdm(range(0, len(t1), 2)):
            r1, bw1 = t1[k].split('_')
            r2, bw2 = t2[k].split('_')
            r3, bw3 = t1[k + 1].split('_')
            r4, bw4 = t2[k + 1].split('_')
            x_data = self.yc[:, k] * 10
            y_data = self.yc[:, k + 1] * 10
            self.free_energy_plot(x_data, y_data, [r1, r2, r3, r4], [bw1, bw2, bw3, bw4])
            j = pyemma.plots.plot_free_energy(x_data, y_data, cmap='jet', nbins=200, ncontours=300, kT=2.479 / 4.184,
                                              cbar_label='kcal/mol')
            j[1].set_ylabel(f'Distance {r1} - - {r2} (Å)', fontsize=18)
            j[1].set_xlabel(f'Distance {r3} - - {r4} (Å)', fontsize=18)
            j[0].savefig(f'Distance_{r3}--{r4}_vs_{r1}--{r2}_pyemma.png')
            plt.close('all')

def main():
    plotter = FreeEnergyPlotter(
        msm_path='/home/pdb3/SMO/APO/Analysis/MSM/MSMobjs/MSMobj_9D_1200.py',
        df_path='/home/pdb3/SMO/APO/Analysis/SMO_BW_Numbering.csv',
        data_path='totdistcombined_extra_tunnel.npy',
        features1=['NH2_342', 'NH2_342', 'OD1_415', ' CG_205', 'CD2_394', 'NH1_393', 'CD1_209', 'CD1_277', 'CD1_277', 'CD1_397',
                   'CE1_397', 'CG1_263', 'CD2_267', ' CZ_274', 'CG_281', 'CG2_208', ' N_398', 'NH2_393', 'CE2_404', 'ND1_412', 'CD2_464'],
        features2=['OD2_415', 'OE1_460', 'OE1_460', 'CD1_394', 'CZ3_477', 'CZ3_477', 'CE3_477', 'CD1_397', 'CD2_209', 'CD2_209',
                   'CE3_477', ' OH_264', ' CE_268', 'CZ2_273', 'CA_280', 'CD2_209', 'CZ_397', 'CD2_394', 'CA_405', 'CD1_413', 'ND2_463']
    )
    plotter.plot_all()

if __name__ == "__main__":
    main()

