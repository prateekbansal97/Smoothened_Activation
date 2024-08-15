import pandas as pd
import argparse

class FrameExtractor:
    def __init__(self, input_pdb, csv_path='../SMO_BW_Numbering.csv'):
        self.input_pdb = input_pdb
        self.csv_path = csv_path
        self.lines = self._read_pdb()
        self.smobw = pd.read_csv(self.csv_path, delimiter=',')

    def _read_pdb(self):
        with open(self.input_pdb, 'r') as f:
            return f.readlines()

    def update_residue_numbers(self):
        for j, line in enumerate(self.lines):
            if line.startswith('ATOM'):
                pymol_num = int(line.split()[4])
                actual_resno = int(self.smobw.loc[self.smobw['pymol_num_active'] == pymol_num].Resno)
                actual_resnostr = f"{actual_resno:>3}"
                self.lines[j] = f'{line[:23]}{actual_resnostr}{line[26:]}'

    def save_pdb(self):
        output_pdb = f"{self.input_pdb[:-4]}_conv.pdb"
        with open(output_pdb, 'w') as g:
            g.write(''.join(self.lines))

def main():
    parser = argparse.ArgumentParser(description='Extract Frames')
    parser.add_argument('-i', '--input', action='store', type=str, required=True)
    args = parser.parse_args()

    extractor = FrameExtractor(args.input)
    extractor.update_residue_numbers()
    extractor.save_pdb()

if __name__ == '__main__':
    main()

