import pandas as pd
import numpy as np

class BWNumbering:
    def __init__(self, csv_path):
        self.df = pd.read_csv(csv_path)

    def get_mdtraj_num(self, resids, BW=False):
        results = []
        for res in resids:
            r, rn = res[0], int(res[1:])
            loc = np.where(self.df.Resno == rn)
            mdt = int(self.df.iloc[loc].mdtraj_num_active)
            result = {'mdtraj_num': mdt}
            if BW:
                result['BW'] = f'{r}{rn}_{self.df.iloc[rn-1].Location}.{int(self.df.iloc[rn-1].BW)}'
            results.append(result)
        return results

    def get_BW_from_mdt(self, resi):
        BWe = []
        for res in resi:
            r, rn = res.split('_')
            loc = np.where(self.df.mdtraj_num_active == int(rn))[0][0]
            BWe.append(f'{self.df.iloc[loc].Residue}{self.df.iloc[loc].Resno}_{self.df.iloc[loc].Location}.{int(self.df.iloc[loc].BW)}')
        return BWe

class TunnelFeatureComparison:
    def __init__(self, bw_numbering):
        self.bw_numbering = bw_numbering

    def compare_features(self, features1, features2):
        t1 = self.bw_numbering.get_BW_from_mdt(features1)
        t2 = self.bw_numbering.get_BW_from_mdt(features2)
        return list(zip(t1, t2))

bw_numbering = BWNumbering('/home/pdb3/SMO/APO/Analysis/SMO_BW_Numbering.csv')

resids = ['W535', 'R451', 'L325']
mdtraj_results = bw_numbering.get_mdtraj_num(resids, BW=True)
for result in mdtraj_results:
    print(result)

features_tunnel_extra1 = ['NH2_342','NH2_342','OD1_415','CG_205','CD2_394','NH1_393','CD1_209','CD1_277','CD1_277','CD1_397','CE1_397','CG1_263','CD2_267','CZ_274','CG_281','CG2_208','N_398','NH2_393','CE2_404','ND1_412','CD2_464']
features_tunnel_extra2 = ['OD2_415','OE1_460','OE1_460','CD1_394','CZ3_477','CZ3_477','CE3_477','CD1_397','CD2_209','CD2_209','CE3_477','OH_264','CE_268','CZ2_273','CA_280','CD2_209','CZ_397','CD2_394','CA_405','CD1_413','ND2_463']

tunnel_comparison = TunnelFeatureComparison(bw_numbering)
comparison_results = tunnel_comparison.compare_features(features_tunnel_extra1, features_tunnel_extra2)
for i, j in comparison_results:
    print(f'{i} -- {j}')
