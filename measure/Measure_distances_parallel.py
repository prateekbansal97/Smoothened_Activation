import subprocess
import os
import numpy as np

class NodeManager:
    def __init__(self):
        self.nodes = []
        self.act_nodes = []
        self.act_nodes_free = []
        self.percent_free_nodes = []
        self.n_nodes = 0
        self.n_traj = 0
        self.ratio_dist = []
        self.chunks = []

    def find_active_nodes(self):
        subprocess.call('qstat -f > act_nodes', shell=True)
        with open('act_nodes') as f:
            lines = f.readlines()
        self.nodes = [j.split()[0] for j in lines if 'local' in j]
        self.act_nodes = [j.split()[0] for j in lines if 'local' in j and j.split()[-1] == 'linux-x64']
        self.percent_free_nodes = [
            (int(j.split()[2].split('/')[2]) - int(j.split()[2].split('/')[1]))
            for j in lines if 'local' in j and j.split()[-1] == 'linux-x64'
        ]
        self.act_nodes_free = [
            j for i, j in enumerate(self.act_nodes) if self.percent_free_nodes[i] != 0
        ]
        self.n_nodes = len(self.act_nodes_free)

    def generate_traj_list(self):
        subprocess.call('source /home/pdb3/SMO/APO/Analysis/trajlistgen.sh', shell=True)
        self.n_traj = range(len(open('/home/pdb3/SMO/APO/Analysis/trajlist').readlines()))

    def proper_round(self, num, dec=0):
        num_str = str(num)[:str(num).index('.') + dec + 2]
        if num_str[-1] >= '5':
            a = num_str[:-2 - (not dec)]
            b = int(num_str[-2 - (not dec)]) + 1
            return float(a) + b**(-dec + 1) if a and b == 10 else float(a + str(b))
        return float(num_str[:-1])

    def calculate_ratio_dist(self):
        self.ratio_dist = [
            int(self.proper_round(self.n_traj[-1] * j / sum(self.percent_free_nodes)))
            for j in self.percent_free_nodes if j != 0
        ]
        if sum(self.ratio_dist) != int(self.n_traj[-1] + 1):
            diff = int(self.n_traj[-1] + 1) - sum(self.ratio_dist)
            position = np.random.randint(len(self.ratio_dist))
            self.ratio_dist[position] += diff

    def split_chunks(self):
        chunk = np.cumsum(self.ratio_dist)
        chunk = np.array([0] + list(chunk))
        self.chunks = [[chunk[i] + 1, chunk[i + 1]] for i in range(len(chunk) - 1)]

    def create_job_scripts(self):
        for n, chunk in enumerate(self.chunks):
            with open(f'./Dist_node_{self.act_nodes_free[n].split("@")[1]}', 'w+') as f:
                f.write(f'''### Use this script for CPU jobs
### Lines starting with "#$" are options for the SGE scheduling system
#$ -S /bin/bash    # Set shell to run job
#$ -q {self.act_nodes_free[n]}        # Choose queue to run job in
#$ -pe onenode 1      # Request processors, other options include onenode, distribute, orte
#$ -cwd            # Run job from my current working directory
#$ -t {chunk[0]}-{chunk[-1]}
TRAJ=`awk "NR==$SGE_TASK_ID" /home/pdb3/SMO/APO/Analysis/trajlist`
python /home/pdb3/SMO/APO/Analysis/Metrics.py -i "${{TRAJ}}"''')

    def create_submission_script(self):
        with open('./Submission_Script.sh', 'w+') as f:
            f.write(f'qsub -N job0 Dist_node_{self.act_nodes_free[0].split("@")[1]}\n')
            for i in range(1, len(self.chunks)):
                if i % 2 != 0:
                    f.write(f'qsub -N job{i} Dist_node_{self.act_nodes_free[i].split("@")[1]}\n')
                else:
                    f.write(f'qsub -N job{i} -hold_jid job{i-1} Dist_node_{self.act_nodes_free[i].split("@")[1]}\n')

def main():
    manager = NodeManager()
    manager.find_active_nodes()
    manager.generate_traj_list()
    manager.calculate_ratio_dist()
    manager.split_chunks()
    manager.create_job_scripts()
    manager.create_submission_script()

if __name__ == '__main__':
    main()

