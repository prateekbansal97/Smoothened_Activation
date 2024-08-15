import pyemma
import matplotlib.pyplot as plt
import numpy as np

class TICAAnalysis:
    def __init__(self, data_path, lag=300, var_cutoff=0.95):
        self.data = np.load(data_path, allow_pickle=True).tolist()
        self.lag = lag
        self.var_cutoff = var_cutoff
        self.tica = None
        self.tica_concatenated = None

    def perform_tica(self):
        self.tica = pyemma.coordinates.tica(self.data, lag=self.lag, var_cutoff=self.var_cutoff)
        tica_output = self.tica.get_output()
        self.tica_concatenated = np.concatenate(tica_output)

    def extract_feature(self, index=34):
        feat = [traj[:, index] for traj in self.data]
        return np.array(feat)

    def plot_free_energy(self):
        if self.tica_concatenated is not None:
            plt.figure()
            pyemma.plots.plot_free_energy(
                *self.tica_concatenated[:, :2].T,
                vmax=5,
                kT=2.479 / 4.184,
                cbar_label='Free Energy (kcal/mol)',
                cmap='jet',
                nbins=150
            )
            plt.show()

def main():
    analysis = TICAAnalysis(data_path='../RRCS_MSM/totdist_RRCSg40.npy')
    analysis.perform_tica()
    feature = analysis.extract_feature(index=34)
    analysis.plot_free_energy()
    print(feature.shape)

if __name__ == "__main__":
    main()

