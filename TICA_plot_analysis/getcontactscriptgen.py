from pdb3 import lsext

class ContactScriptGenerator:
    def __init__(self, output_file, parm_files, contact_types, cutoffs, names):
        self.output_file = output_file
        self.parm_files = parm_files
        self.contact_types = contact_types
        self.cutoffs = cutoffs
        self.names = names

    def generate_script(self):
        with open(self.output_file, 'w+') as f:
            self._write_dynamic_contacts(f)
            self._write_contact_frequencies(f)
            self._write_contact_fingerprints(f)

    def _write_dynamic_contacts(self, f):
        for j in self.contact_types:
            for i, p in zip(self.names, self.parm_files):
                f.write(f'get_dynamic_contacts.py --topology {p} --trajectory ./W339_3.50_{i}_combined.dcd '
                        f'--itypes {j} --output tsv/contacts_{i}_{j}.tsv\n')

    def _write_contact_frequencies(self, f):
        for j in self.contact_types:
            for i in self.names:
                f.write(f'get_contact_frequencies.py --input_files tsv/contacts_{i}_{j}.tsv '
                        f'--output_file tsv/resfrequencies_{i}_{j}.tsv\n')

    def _write_contact_fingerprints(self, f):
        for j, cutoff in zip(self.contact_types, self.cutoffs):
            inputs = ' '.join([f'tsv/resfrequencies_{i}_{j}.tsv' for i in self.names])
            columns = ' '.join(self.names)
            f.write(f'get_contact_fingerprints.py --input_frequencies {inputs} '
                    f'--column_headers {columns} --frequency_cutoff {cutoff} '
                    f'--plot_output fingerprint_{j}.png --pymol_output pymol_{j}.pml '
                    f'--table_output table_op_{j}\n')

def main():
    parm_files = [
        '../../../6XBL_SAG/Analysis/6XBL_SAG_HMass_stripped_vmd2.psf',
        '../../../5L7D_SANT1/Analysis/5L7D_SANT1_HMass_stripped.parm7'
    ]
    contact_types = ["hp", "sb", "pc", "ps", "ts", "vdw", "hb"]
    cutoffs = [0.1, 0.6, 0.5, 0.1, 0.2, 0.5, 0.6]
    names = ['high_SAG', 'low_SANT1']
    
    generator = ContactScriptGenerator(
        output_file='./source_this_script_for_getcontacts',
        parm_files=parm_files,
        contact_types=contact_types,
        cutoffs=cutoffs,
        names=names
    )
    generator.generate_script()

if __name__ == "__main__":
    main()

