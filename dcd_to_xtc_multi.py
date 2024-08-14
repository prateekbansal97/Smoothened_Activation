from multiprocessing import Pool
import subprocess
from pdb3 import lsext

class DCDConverter:
    def __init__(self, directory, extension):
        self.xtc_list = self._get_xtc_list(directory, extension)

    def _get_xtc_list(self, directory, extension):
        return lsext(directory, extension)[0]

    def _convert_dcd_to_xtc(self, cpp):
        subprocess.call(f'cpptraj -i {cpp}', shell=True)

    def convert_all(self, pool_size=64):
        with Pool(pool_size) as pool:
            pool.map(self._convert_dcd_to_xtc, self.xtc_list)

def main():
    converter = DCDConverter('./dcd_to_xtc/', 'cpp')
    converter.convert_all()

if __name__ == '__main__':
    main()

