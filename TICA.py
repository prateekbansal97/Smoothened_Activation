import numpy as np
import pyemma
import matplotlib.pyplot as plt
import pickle

class TICAAnalysis:
    def __init__(self, data_path, trajnames_path, dims=20, cluster_centers=3000):
        self.data = np.load(data_path, allow_pickle=True).tolist()
        self.dims = dims
        self.cluster_centers = cluster_centers
        self.trajnames = np.load(trajnames_path, allow_pickle=True)
        self.feature_labels = ['IC{}'.format(i + 1) for i in range(dims)]

    def perform_tica(self):
        tica = pyemma.coordinates.tica(self.data, dim=self.dims)
        tica_output = tica.get_output()
        pickle.dump(tica, open(f'./pkl/tica_{self.dims}D_{self.cluster_centers}_ticaobj.pkl', 'wb'))
        return tica, tica_output

    def cluster_tica_output(self, tica_output):
        cluster_tica = pyemma.coordinates.cluster_kmeans(tica_output, k=self.cluster_centers, max_iter=300, stride=1)
        pickle.dump(cluster_tica, open(f'./pkl/cluster_{self.cluster_centers}_tica_{self.dims}D_clusobj.pkl', 'wb'))
        return cluster_tica

    def save_tica_correlation(self, tica):
        lp = tica.feature_TIC_correlation
        pickle.dump(lp, open(f'./pkl/TICA_coordinates_breakdown_{self.dims}D_{self.cluster_centers}.npy', 'wb'))

    def save_dtrajs(self, cluster_tica):
        cluster_tica.save_dtrajs(trajfiles=self.trajnames, prefix='tica_dtraj', output_dir='./tica_dtrajs/')

    def plot_its(self, cluster_tica):
        its_tica = pyemma.msm.its(dtrajs=cluster_tica.dtrajs, lags=50, nits=10)
        pyemma.plots.plot_implied_timescales(its_tica, units='ns', dt=0.1)
        plt.title('ITS-TICA')
        plt.savefig('its-tica.png', dpi=300)
        pyemma.plots.plot_implied_timescales(its_tica, units='ns')

def main():
    analysis = TICAAnalysis(data_path='SMO_APO_totdist.npy', trajnames_path='./File_I.npy')
    
    tica, tica_output = analysis.perform_tica()
    cluster_tica = analysis.cluster_tica_output(tica_output)
    
    analysis.save_tica_correlation(tica)
    analysis.save_dtrajs(cluster_tica)
    analysis.plot_its(cluster_tica)

if __name__ == "__main__":
    main()

