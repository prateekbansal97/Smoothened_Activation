import os
import pickle
import random
import numpy as np
from tqdm import tqdm
from pdb3 import lsext

class ClusterGenerator:
    def __init__(self, traj_dir='./tica_dtrajs/', ext='tica_dtraj'):
        self.dtrajs, self.ndtrajs = lsext(traj_dir, ext, sort=True, preapp=True)
        self.clusters = {j: [] for j in range(3000)}

    def generate_clusters(self):
        for i in tqdm(self.dtrajs, total=len(self.dtrajs)):
            with open(i) as f:
                lines = f.readlines()
            for j, k in enumerate(lines):
                key = int(k.rstrip())
                traj_file = '_'.join(i.split('/')[-1].split('_')[2:-1]).strip('-strip') + '.nc'
                self.clusters[key].append(f'{traj_file} {j + 1}')
        self.save_clusters()

    def save_clusters(self, filepath='./pkl/traj_to_clusters.pkl'):
        with open(filepath, 'wb') as f:
            pickle.dump(self.clusters, f)

class CpptrajGenerator:
    def __init__(self, clusters, choice, output_dir='/home/pdb3/ds02/SMO/DeltaCRD/System_Building/startrstgen/'):
        self.clusters = clusters
        self.choice = choice
        self.output_dir = output_dir
        self.startpath = '/home/pdb3/SMO/DeltaCRD/System_Building/startrst/'
        self.parm_6XBL_APO = "path_to_6XBL_APO_parm"  # Replace with actual path
        self.parm_5L7D_APO = "path_to_5L7D_APO_parm"  # Replace with actual path

    def generate_cpptraj_files(self):
        c, d = 1, 1
        for j, i in tqdm(enumerate(self.choice), total=len(self.choice)):
            file = i.split()[0]
            if '6XBL' in file:
                parm = self.parm_6XBL_APO
                pdb = self.startpath + f'6XBL_Delta_CRD_{c}'
                c += 1
            else:
                parm = self.parm_5L7D_APO
                pdb = self.startpath + f'5L7D_Delta_CRD_{d}'
                d += 1
            self._write_cpptraj_file(j, parm, file, pdb, i.split()[1])

    def _write_cpptraj_file(self, index, parm, trajin_file, pdb_output, frame_num):
        filename = os.path.join(self.output_dir, f'SMO_DeltaCRDAPO_R_01_pdbgen_{index + 1}')
        with open(filename, 'w+') as f:
            f.write(f'''parm {parm}
trajin /home/pdb3/SMO/APO/Analysis/nc/{trajin_file}
center origin
autoimage origin
strip :WAT
trajout {pdb_output}.pdb onlyframes {frame_num}
run
quit''')

def main():
    cluster_gen = ClusterGenerator()
    cluster_gen.generate_clusters()

    clusters = pickle.load(open('./pkl/traj_to_clusters.pkl', 'rb'))
    choice = [random.choice(clusters[i]) for i in np.random.randint(0, 3000, 1000)]
    
    with open('./pkl/frames_for_SMO_DeltaCRD.pkl', 'wb') as f:
        pickle.dump(choice, f)
    
    cpptraj_gen = CpptrajGenerator(clusters, choice)
    cpptraj_gen.generate_cpptraj_files()

if __name__ == "__main__":
    main()

