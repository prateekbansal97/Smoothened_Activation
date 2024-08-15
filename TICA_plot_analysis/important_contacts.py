import numpy as np
import pandas as pd
import os
from pdb3 import lsext

class ContactAnalyzer:
    def __init__(self, table_dir, bw_file):
        self.tables = lsext(table_dir, 'table', preapp=True, abs=True, nat=True)[0]
        self.bw = pd.read_csv(bw_file, delimiter=',')

    def analyze_contacts(self):
        for table in self.tables:
            df = self._load_table(table)
            df = self._process_table(df)
            self._save_results(df, table)

    def _load_table(self, table):
        df = pd.read_csv(table, delimiter='\t')
        df.columns = ['Residue1', 'Residue2', 'high_SAG', 'low_SANT1']
        df['Ago_Ant_diff'] = np.abs(df['high_SAG'] - df['low_SANT1'])
        return df

    def _process_table(self, df):
        intrahelical_or_not = []
        residue1actual, residue2actual, location1, location2 = [], [], [], []
        residue1_vmd, residue2_vmd = [], []

        for j in range(len(df)):
            resi1, resi2 = self._get_residue_indices(df, j)
            residue1_vmd.append(resi1)
            residue2_vmd.append(resi2)
            resiactual1, resiactual2 = self._get_actual_residues(resi1, resi2)
            residue1actual.append(resiactual1)
            residue2actual.append(resiactual2)
            loc1, loc2 = self._get_locations(resi1, resi2)
            location1.append(loc1)
            location2.append(loc2)
            intrahelical_or_not.append(self._determine_type(resi1, resi2))

        df['Location1'] = location1
        df['Location2'] = location2
        df['Type'] = intrahelical_or_not
        df['Residue1_vmd'] = residue1_vmd
        df['Residue2_vmd'] = residue2_vmd
        df['Residue1'] = residue1actual
        df['Residue2'] = residue2actual
        return df

    def _get_residue_indices(self, df, index):
        resi1 = int(df.iloc[index]['Residue1'].split(':')[-1])
        resi2 = int(df.iloc[index]['Residue2'].split(':')[-1])
        return resi1, resi2

    def _get_actual_residues(self, resi1, resi2):
        resiactual1 = int(self.bw.loc[self.bw['pymol_num_active'] == resi1].Resno)
        resiactual2 = int(self.bw.loc[self.bw['pymol_num_active'] == resi2].Resno)
        resiactual1_str = f"{resiactual1:>3}"
        resiactual2_str = f"{resiactual2:>3}"
        return resiactual1_str, resiactual2_str

    def _get_locations(self, resi1, resi2):
        bwnum1 = self._format_bw_number(resi1)
        bwnum2 = self._format_bw_number(resi2)
        loc1 = str(self.bw.loc[self.bw['pymol_num_active'] == resi1].Location).split('\n')[0].split()[1] + bwnum1
        loc2 = str(self.bw.loc[self.bw['pymol_num_active'] == resi2].Location).split('\n')[0].split()[1] + bwnum2
        return loc1, loc2

    def _format_bw_number(self, resi):
        bw_value = self.bw.loc[self.bw['pymol_num_active'] == resi].BW
        return f".{str(float(bw_value)).split('.')[0]}" if not pd.isna(bw_value) else ''

    def _determine_type(self, resi1, resi2):
        return 'Intra-Helical' if abs(resi1 - resi2) <= 5 else 'Extra-Helical'

    def _save_results(self, df, table):
        table_name = table.split('/')[-1].split('_')[-1].upper()
        important_contacts = df[df['Ago_Ant_diff'] > 0.5]
        important_contacts.to_excel(f"./contacts_{table_name}_SAG_SANT1.xlsx")
        important_contacts.to_csv(f"./contacts_{table_name}_SAG_SANT1.csv", sep=',')
        self._print_vmd_commands(important_contacts)

    def _print_vmd_commands(self, df):
        r1 = df['Residue1'].to_list()
        r2 = df['Residue2'].to_list()
        print(' or resid '.join([str(r.split(':')[-1]) for r in r1]))
        print('\n\n\n')
        print(' or resid '.join([str(r.split(':')[-1]) for r in r2]))
        print('\n\n\n\n\n\n\n\n')

def main():
    analyzer = ContactAnalyzer(table_dir='./', bw_file='../SMO_BW_Numbering.csv')
    analyzer.analyze_contacts()

if __name__ == "__main__":
    main()

