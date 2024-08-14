import numpy as np
import matplotlib.pyplot as plt
import pyemma
import matplotlib as mpl
from matplotlib import font_manager
import pandas as pd

class BWNumbering:
    def __init__(self, csv_path):
        self.df = pd.read_csv(csv_path)

    def get_BW_from_mdt(self, resi):
        BWe = []
        for res in resi:
            r, rn = res.split('_')
            loc = np.where(self.df.mdtraj_num_active == int(rn))[0][0]
            BWe.append(f'{self.df.iloc[loc].Residue}{self.df.iloc[loc].Resno}_{self.df.iloc[loc].Location}.{int(self.df.iloc[loc].BW)}')
        return BWe

class FreeEnergyPlotter:
    def __init__(self, y_data_file, temperature=310):
        self.R = 0.001987
        self.T = temperature
        self.yc = np.concatenate(np.load(y_data_file, allow_pickle=True))

    def calculate_histogram(self, x_data, y_data, x_bins=300, y_bins=300):
        hist = np.histogram2d(
            x_data, y_data, 
            bins=[x_bins, y_bins],
            range=[[np.min(x_data) - 0.5, np.max(x_data) + 0.5], [np.min(y_data) - 0.1, np.max(y_data) + 0.1]],
            density=True
        )
        return hist

    def calculate_free_energy(self, prob_density, x_bin_size, y_bin_size):
        free_energy = -self.R * self.T * np.log(prob_density * x_bin_size * y_bin_size)
        return free_energy - np.min(free_energy)

    def plot_free_energy(self, rs, bws, max_energy=3.5, x_bins=300, y_bins=300):
        r1, r2 = rs
        bw1, bw2 = bws
        
        x_data = self.yc[:, 0]
        y_data = self.yc[:, 1] * np.pi / 180
        
        hist = self.calculate_histogram(x_data, y_data, x_bins, y_bins)
        prob_density = hist[0]
        xedge = hist[1]
        yedge = hist[2]
        
        x_bin_size = xedge[1] - xedge[0]
        y_bin_size = yedge[1] - yedge[0]
        delta_free_energy = self.calculate_free_energy(prob_density, x_bin_size, y_bin_size)
        
        xx = [(xedge[i] + xedge[i + 1]) / 2 for i in range(len(xedge) - 1)]
        yy = [(yedge[i] + yedge[i + 1]) / 2 for i in range(len(yedge) - 1)]
        
        fig, axs = plt.subplots(1, 1, figsize=(10, 7))
        contours = np.linspace(0, max_energy, 5)
        axs.contour(xx, yy, delta_free_energy.T, levels=6, colors='k', linewidths=0.2)
        cd = axs.contourf(xx, yy, delta_free_energy.T, np.linspace(0, max_energy, 30), vmin=0.0, vmax=max_energy, cmap=mpl.cm.jet)
        cbar = fig.colorbar(cd, ticks=range(int(max_energy) + 1))
        cbar.ax.set_yticklabels(range(int(max_energy) + 1))
        
        axs.set_xlim(3, 15)
        axs.set_ylim(0, np.pi / 2)
        axs.set_xticks(range(int(np.min(x_data)) - 1, int(np.max(x_data)) + 2, 2))
        axs.set_xticklabels(range(int(np.min(x_data)) - 1, int(np.max(x_data)) + 2, 2))
        axs.set_yticks(np.around(np.arange(0, np.pi / 2 + 0.1, 0.2), 2))
        axs.set_yticklabels(np.around(np.arange(0, np.pi / 2 + 0.1, 0.2), 2))
        
        plt.ylabel(f'{r1}-{r2} ({bw1}-{bw2})', fontsize=18)
        plt.xlabel('Distance (Å)', fontsize=18)    
        plt.rc('xtick', labelsize=10)
        plt.rc('ytick', labelsize=10)
        plt.grid(ls='-', lw=0.2)
        plt.savefig(f'{r1}_{r2}_Free_Energy.png', dpi=500)

def main():
    # BW Numbering setup
    bw_numbering = BWNumbering('/home/pdb3/SMO/APO/Analysis/SMO_BW_Numbering.csv')
    
    # Free Energy Plotting setup
    plotter = FreeEnergyPlotter('totdistcombined_metrics.npy')
    resids = ['Res1', 'Res2']  # Replace with actual residues
    bws = bw_numbering.get_BW_from_mdt(resids)
    
    # Plot the free energy
    plotter.plot_free_energy(rs=resids, bws=bws)

if __name__ == "__main__":
    main()

