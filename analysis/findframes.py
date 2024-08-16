'''
Written by Prateek Bansal
University of Illinois
November 23 2021
pdb3@illinois.edu
'''
'''
Example usage:

CRD_TMD_Angle.csv:
140.04192872117397, 24.305882352941175
144.73794549266245, 23.670588235294115
148.84696016771488, 24.270588235294117
144.44444444444446, 24.97647058823529


Command:
python findframes.py -i ./totdistcombined_BP_6XBL_APO.npy -x 1 -y 0 -xl 145 150 -yl 22.5 25 -td ../dcd/ -tn ./dist_I_BP_6XBL_APO.npy --ymul -1 --framename CRD_Angle_APO -parp strip -s 6XBL 5L7D -pa CRD_TMD_Angle.csv

'''

import numpy as np
import argparse
from pdb3 import lsext
import subprocess
import os
from tqdm import tqdm
import matplotlib.path as mpltPath
from numpy import genfromtxt

class find_frames():

    def __init__(self, inp, x, y, xl, yl, t, tn, name, parm=None,path=None,sys=None,ext='rst7', xm=1, ym=1, trajtype='.dcd', parmprompt=None):
        self.input = np.load(inp,allow_pickle=True)        
        self.x = x
        self.y = y
        self.xl = xl
        self.yl = yl
        self.t = t
        self.tn = np.load(tn,allow_pickle=True)
        self.parm = parm
        self.path = path
        self.sys = sys
        self.name = name
        self.ext = ext
        self.xm = xm
        self.ym = ym
        self.trajtype = trajtype
        self.parmprompt = parmprompt
    
    def find_parm(self): 
        ''' 
            Automatically searches the parent and superparent directory for parm files
        '''
        try:
            if not self.parm:
                self.parm = lsext('../../','.parm7',preapp=True,abs=True,extra=self.parmprompt)[0]
                if not self.parm:
                    self.parm = lsext('../','.parm7',preapp=True,abs=True,extra=self.parmprompt)[0]
            self.parm = os.path.abspath(self.parm) if isinstance(self.parm,str) else [os.path.abspath(j) for j in self.parm]
        except: 
            print('No parm file supplied or found')
        return self.parm


    
    def traj_get(self):
        '''
            Gets Trajectory files from the specified traj folder
        '''
        trajs = lsext(self.t,self.trajtype,nat=True,abs=True,preapp=True)[0]
        return trajs

    
    def dirsgen(self):
        '''
       Makes directories
        '''
        dirs = [name for name in os.listdir(".") if os.path.isdir(name)]
        if 'frames' not in dirs:
            subprocess.call('mkdir frames',shell=True)
        elif 'framesgen' not in dirs:
            subprocess.call('mkdir framesgen',shell=True)


    def get_frames(self):
        '''
       Generates frames given path in csv file or xlims and ylims.
        '''
        self.parm = self.find_parm()
        self.dirsgen()
        if not self.path:
            c=1
            totframes=0
            trajs_frames_selected={}
            for i in tqdm(range(len(self.input))):
                
                if self.sys:
                    sys = [f'{j}' for j in self.sys if f'{j}' in self.tn[i]][0] 
                    try:
                        parmu = [f'{j}' for j in self.parm if sys in j][0]
                    except IndexError:
                        parmu = '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7' if 'strip' in self.parmprompt else '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass.parm7' 
                else:
                    parmu = self.parm
                
                full_trajs = self.traj_get()
                htraj = self.tn[i]
                traj = [j for j in full_trajs if htraj in j][0]
                points = self.input[i]
                xf = points[:,self.x]*self.xm
                yf = points[:,self.y]*self.ym
                pos = np.where((xf<self.xl[1]) & (xf>self.xl[0]) & (yf>self.yl[0]) & (yf<self.yl[1]))[0]
                totframes+=len(pos)  
                
                if pos.size!=0:
                    trajs_frames_selected[traj]=pos
                    for j in range(len(pos)):
                        rst = f'../frames/{self.name}_{c}.{self.ext}'
                        f=open(f'./framesgen/{self.name}_rstgen_{c}_cpp','w+')
                        f.write(f'''parm {parmu}
trajin {traj}
trajout {rst} onlyframes {pos[j]+1}
run
quit''')
                        f.close()
                        c+=1
            #print(trajs_frames_selected)
            print('Total no. of frames in the specified region: {}'.format(totframes))
        
        else:
            c=1
            totframes=0
            polygon = genfromtxt(self.path, delimiter=',')
            print(polygon)
            path = mpltPath.Path(polygon)
            for i in tqdm(range(len(self.input))):
                xpoints = self.input[i][:,[self.x]]*self.xm
                ypoints = self.input[i][:,[self.y]]*self.ym
                points=np.hstack((xpoints,ypoints))
                inside2 = path.contains_points(points)
                frames=inside2.nonzero()[0]
                totframes+=len(frames) 
                #print(totframes)
                if self.sys:
                    try:
                        sys = [f'{j}' for j in self.sys if f'{j}' in self.tn[i]][0]
                        parmu = [f'{j}' for j in self.parm if sys in j][0]
                    except IndexError:
                        parmu = '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass_stripped.parm7' if 'strip' in self.parmprompt else '/home/pdb3/SMO/APO/Analysis/5L7D_Apo_HMass.parm7'
                else:
                    parmu = self.parm

                full_trajs = self.traj_get()
                htraj = self.tn[i]
                traj = [j for j in full_trajs if htraj in j][0]

                if len(inside2.nonzero()[0])!=0:
                     #print(frames)
                     for j in range(len(frames)):
                         rst = f'../frames/{self.name}_{c}.{self.ext}'
                         f=open(f'./framesgen/{self.name}_rstgen_{c}_cpp','w+')
                         f.write(f'''parm {parmu}
trajin {traj}
trajout {rst} onlyframes {frames[j]+1}
run
quit''')
                         f.close()
                         c+=1
            print('Total no. of frames in the specified region: {}'.format(totframes))

def get_args():
    my_parser = argparse.ArgumentParser(description='Extract Frames')
    my_parser.add_argument('-i','--input',action='store',type=str,required=True)
    my_parser.add_argument('-x','--xind',action='store',type=int,required=True)
    my_parser.add_argument('-y','--yind',action='store',type=int,required=True)
    my_parser.add_argument('-xl','--xlims',action='store',type=float,required=True,nargs=2)
    my_parser.add_argument('-yl','--ylims',action='store',type=float,required=True,nargs=2)
    my_parser.add_argument('-td','--trajdir',action='store',type=str,required=True)
    my_parser.add_argument('-tn','--trajnames',action='store',type=str,required=True)
    my_parser.add_argument('-n','--framename',action='store',type=str,required=True)
    my_parser.add_argument('-p','--parm',action='store',type=str,required=False)
    my_parser.add_argument('-pa','--path',action='store',type=str,required=False)
    my_parser.add_argument('-e','--ext',action='store',type=str,required=False,default='rst7')
    my_parser.add_argument('-s','--sys',action='store',type=str,required=False,nargs='*')
    my_parser.add_argument('--xmul',action='store',type=float,required=False,default=1)
    my_parser.add_argument('--ymul',action='store',type=float,required=False,default=1)
    my_parser.add_argument('-parp','--parmprompt',action='store',type=str,required=False)
    args = my_parser.parse_args()
    return args

def print_args(args):
    for j in args.args():
        print(j)
if __name__=='__main__':
    args = get_args()
    #print_args(args)
    frames = find_frames(args.input, args.xind, args.yind, args.xlims, args.ylims, args.trajdir, args.trajnames, args.framename, parm=args.parm, path=args.path, sys=args.sys, ext=args.ext, xm=args.xmul, ym=args.ymul, parmprompt = args.parmprompt)
    frames.get_frames()

