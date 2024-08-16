import numpy as np
import pyemma
import pickle
import matplotlib.pyplot as plt

class TICAAnalysis:
    def __init__(self, file_path='./File_I.npy', data_path='totdistcombined_extra.npy'):
        self.y = np.load(file_path)
        self.z = [i for i, j in enumerate(self.y) if '6XBL' in j]
        self.tip = self.z[-1] + 1
        self.data = np.load(data_path, allow_pickle=True)

    def split_data(self):
        inac = self.data[self.tip:]
        ac = self.data[:self.tip]
        np.save('./totdist5L7D_extra.npy', inac)
        np.save('./totdist6XBL_extra.npy', ac)
        return inac, ac

    def perform_tica(self):
        data_list = self.data.tolist()
        tica = pyemma.coordinates.tica(data_list, dim=20)
        tica_output = tica.get_output()
        tica_concatenated = np.concatenate(tica_output)
        return tica, tica_concatenated

    def transform_data(self, tica, inac, ac):
        inac_trans = tica.transform(np.concatenate(inac))
        ac_trans = tica.transform(np.concatenate(ac))
        pickle.dump(inac_trans, open("./inac_transformed.pkl", "wb"))
        pickle.dump(ac_trans, open("./ac_transformed.pkl", "wb"))
        return inac_trans, ac_trans

    def plot_free_energy(self, tica_concatenated, inac_trans, ac_trans):
        fig, ax = plt.subplots(1, 3, figsize=(15, 3.75))
        pyemma.plots.plot_free_energy(*tica_concatenated[:, :2].T, cmap='jet', nbins=250, ax=ax[2], vmin=0, vmax=6)
        ax[2].set_xlabel('IC 1')
        ax[2].set_ylabel('IC 2')
        xmin, xmax, ymin, ymax = ax[2].axis()

        pyemma.plots.plot_free_energy(*inac_trans[:, :2].T, cmap='Blues', nbins=250, ax=ax[1], cbar=False, vmin=0, vmax=6)
        pyemma.plots.plot_free_energy(*ac_trans[:, :2].T, cmap='Reds', nbins=250, ax=ax[1], cbar=False, vmin=0, vmax=6)
        ax[1].set_xlabel('IC 1')
        ax[1].set_ylabel('IC 2')
        ax[1].set_xlim(xmin, xmax)
        ax[1].set_ylim(ymin, ymax)

        pyemma.plots.plot_free_energy(*ac_trans[:, :2].T, cmap='Reds', nbins=250, ax=ax[0], cbar=False, vmin=0, vmax=6)
        pyemma.plots.plot_free_energy(*inac_trans[:, :2].T, cmap='Blues', nbins=250, ax=ax[0], cbar=False, vmin=0, vmax=6)
        ax[0].set_xlabel('IC 1')
        ax[0].set_ylabel('IC 2')
        ax[0].set_xlim(xmin, xmax)
        ax[0].set_ylim(ymin, ymax)

        fig.savefig('./TICA.png', dpi=500)

def main():
    analysis = TICAAnalysis()
    
    inac, ac = analysis.split_data()
    tica, tica_concatenated = analysis.perform_tica()
    
    inac_trans, ac_trans = analysis.transform_data(tica, inac, ac)
    analysis.plot_free_energy(tica_concatenated, inac_trans, ac_trans)

if __name__ == "__main__":
    main()

