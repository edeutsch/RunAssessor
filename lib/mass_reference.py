#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import copy
import re
import json


####################################################################################################
#### MassReference class
class MassReference:

    ####################################################################################################
    #### Constructor
    def __init__(self, verbose=0):

        # Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        self.atomic_masses = None
        self.ion_series_attributes = None

        self.aa_formulas = None
        self.aa_masses = None
        self.double_aa_masses = None
        self.nterm_modifications = None
        self.nterm_aa_masses = None
        self.nterm_double_aa_masses = None
        self.reporter_ions = None
        self.terminal_modifications = None

        self.neutral_losses = {}
        self.neutral_losses_by_residue = {}
        self.neutral_losses_by_formula = {}

        self.prepare_mass_tables()


    ####################################################################################################
    def prepare_mass_tables(self):

        # Define a subset of useful atomic masses and the proton
        self.atomic_masses = {
            'proton': 1.00727646688,
            'hydrogen': 1.007825035,
            'carbon': 12.0000000,
            'nitrogen': 14.0030740,
            'oxygen': 15.99491463,
            'phosphorus': 30.973762,
            'sulfur': 31.9720707,
        }
        self.atomic_masses_by_letter = {
            '+': 1.00727646688,
            'H': 1.007825035,
            'C': 12.0000000,
            'N': 14.0030740,
            'O': 15.99491463,
            'P': 30.973762,
            'S': 31.9720707,
        }

        # Define the set of peptide ions series along with some attributes of each
        self.ion_series_attributes = {
            'a': { 'terminus_type': 'nterm', 'mass': -1 * self.atomic_masses['carbon'] - self.atomic_masses['oxygen'], 'complement': 'a' },
            'b': { 'terminus_type': 'nterm', 'mass': 0.0, 'complement': 'y' },
            'c': { 'terminus_type': 'nterm', 'mass': 3 *  self.atomic_masses['hydrogen'] + self.atomic_masses['nitrogen'], 'complement': 'z' },
            'x': { 'terminus_type': 'cterm', 'mass': self.atomic_masses['carbon'] + 2 * self.atomic_masses['oxygen'], 'complement': 'a' },
            'y': { 'terminus_type': 'cterm', 'mass': 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'], 'complement': 'b'},
            'z': { 'terminus_type': 'cterm', 'mass': self.atomic_masses['oxygen'] - self.atomic_masses['nitrogen'], 'complement': 'c' },
            'm': { 'terminus_type': 'nterm', 'mass': 0.0, 'complement': 'none' }
        }

        # Define the amino acid formulas and transform
        self.aa_formulas = {
            'G': 'C2H3ON',
            'A': 'C3H5ON',
            'S': 'C3H5O2N',
            'P': 'C5H7ON',
            'V': 'C5H9ON',
            'T': 'C4H7O2N',
            'C': 'C3H5ONS',
            'L': 'C6H11ON',
            'I': 'C6H11ON',
            'N': 'C4H6O2N2',
            'D': 'C4H5O3N',
            'Q': 'C5H8O2N2',
            'K': 'C6H12ON2',
            'E': 'C5H7O3N',
            'M': 'C5H9ONS',
            'H': 'C6H7ON3',
            'F': 'C9H9ON',
            'R': 'C6H12ON4',
            'Y': 'C9H9O2N',
            'W': 'C11H10ON2',
        }
        for aa,formula in self.aa_formulas.items():
            atoms = self.parse_atomic_formula(formula)
            self.aa_formulas[aa] = { 'formula_string': formula, 'atoms': atoms }

        # Define the basic set of amino acids
        # From https://proteomicsresource.washington.edu/protocols06/masses.php
        self.aa_masses = {
            'G': 57.021463735,
            'A': 71.037113805,
            'S': 87.032028435,
            'P': 97.052763875,
            'V': 99.068413945,
            'T': 101.047678505,
            'C': 103.009184505,
            'L': 113.084064015,
            'I': 113.084064015,
            'N': 114.042927470,
            'D': 115.026943065,
            'Q': 128.058577540,
            'K': 128.094963050,
            'E': 129.042593135,
            'M': 131.040484645,
            #'O': 132.089877680,  # probably not useful to have this here. But consider UNIMOD:372 Arg->Orn
            'H': 137.058911875,
            'F': 147.068413945,
            'U': 150.953633405,  # selenocysteine
            'R': 156.101111050,
            'Y': 163.063328575,
            'W': 186.079312980,
            #'X': 113.0840636,
            #'B': 114.534935,
            #'Z': 128.550585,
        }

        # Create a set of non-redundant amino-acid masses without the I
        self.nr_aa_masses = copy.copy(self.aa_masses)
        del(self.nr_aa_masses['I'])

        # Create a set of mass_equivalencies
        self.mass_equivalencies = { 'by_mass': {}, 'by_residue': {} }
        for residue,mass in self.nr_aa_masses.items():
            mass_4digit = '{:.4f}'.format(mass)
            self.mass_equivalencies['by_mass'][mass_4digit] = [ residue ]
            self.mass_equivalencies['by_residue'][residue] = [ residue ]


        # Should make a catalog of equivalencies and not search them at the beginning, but add them in later during final evaluation
        # Would be more efficient to compute these rather than curate by hand
        # Q[Deamidated] = E
        # E[Glu->pyro-Glu] = Q[Gln->pyro-Glu]

        # Add some modifications
        modifications = {
            'Carbamidomethyl': { 'delta_mass': 57.021464, 'formula': '+H3C2NO', 'residues': [ 'C' ] },
            'Oxidation': { 'delta_mass': 15.994915, 'formula': '+O', 'residues': [ 'M', 'W', 'P', 'H' ] },
            'Phospho': { 'delta_mass': 79.966331, 'formula': '+HPO3', 'residues': [ 'S', 'T', 'Y' ], 'neutral_losses': [ 97.976896 ] },
            'Deamidated': { 'delta_mass': 0.984016, 'formula': '-HN+O', 'residues': [ 'N' ] }, # Q[Deamidated] is the same as E, so don't bother here
            'Acetyl': { 'delta_mass': 42.010565, 'formula': '+C2H2O', 'residues': [ 'K' ] },
            #'Pyrophospho': { 'delta_mass': 159.932662, 'residues': [ 'S', 'T', 'Y' ], 'neutral_losses': [ 177.943227, 79.966331 ] },
            #'Beta-methythiolation': { 'delta_mass': 45.987721, 'residues': [ 'C' ] },
            #'Acetylation': { 'delta_mass': 42.010565, 'residues': [ 'K' ] },
            #'Methylthio': { 'delta_mass': 45.987721, 'residues': [ 'K', 'D', 'N', 'C' ] },
            #'Methylation': { 'delta_mass': 14.015650, 'residues': [ 'K', 'R' ] },
            #'Ubiquitination': { 'delta_mass': 114.042927, 'residues': [ 'K' ] },
            #'Carbamylation': { 'delta_mass': 43.005814, 'residues': [ 'K' ] },

            #'TMT': { 'delta_mass': 224.152478, 'residues': [ 'K' ] },
            #'TMT6': { 'delta_mass': 229.162932, 'residues': [ 'K' ] },
            #'TMT6plex': { 'delta_mass': 229.162932, 'residues': [ 'K' ] },
            #'iTRAQ4': { 'delta_mass': 144.102063, 'residues': [ 'K' ] },
            'iTRAQ4plex': { 'delta_mass': 144.102063, 'residues': [ 'K' ] },
        }
        for modification in modifications:
            for residue in modifications[modification]['residues']:
                mod_residue = f"{residue}[{modification}]"
                mass = self.aa_masses[residue] + modifications[modification]['delta_mass']
                self.aa_masses[mod_residue] = mass
                mass_4digit = '{:.4f}'.format(mass)

                if mass_4digit in self.mass_equivalencies['by_mass']:
                    self.mass_equivalencies['by_mass'][mass_4digit].append(mod_residue)
                    reference_residue = self.mass_equivalencies['by_mass'][mass_4digit][0]
                    self.mass_equivalencies['by_residue'][reference_residue].append(mod_residue)
                else:
                    self.mass_equivalencies['by_mass'][mass_4digit] = [mod_residue]
                    self.mass_equivalencies['by_residue'][mod_residue] = [mod_residue]
                    self.nr_aa_masses[mod_residue] = mass



        # Add n-terminal modifications
        self.nterm_aa_modifications = {
            'Glu->pyro-Glu': { 'delta_mass': -18.010565, 'residues': [ 'E' ] },
            #'Gln->pyro-Glu': { 'delta_mass': -17.026549, 'residues': [ 'Q' ] } # only need to do one, because they end up the same mass
        }

        self.nterm_aa_masses = copy.deepcopy(self.aa_masses)
        nterm_residues = {}
        for modification in self.nterm_aa_modifications:
            for residue in self.nterm_aa_modifications[modification]['residues']:
                self.nterm_aa_masses[f"{residue}[{modification}]"] = self.aa_masses[residue] + self.nterm_aa_modifications[modification]['delta_mass']
                nterm_residues[f"{residue}[{modification}]"] = 1




        # Create all the amino acid pairs for de novo
        self.double_aa_masses = {}
        self.nr_double_aa_masses = {}
        nomod_double_aa_masses = {}
        for aa1 in sorted(self.aa_masses.keys()):
            for aa2 in self.aa_masses:

                #### Manually skip I here. We can't use nr_masses above because it might skip modified residues
                if aa1 == 'I' or aa2 == 'I':
                    continue

                residue_pair = aa1 + aa2
                mass = self.aa_masses[aa1] + self.aa_masses[aa2]
                mass_4digit = '{:.4f}'.format(mass)

                if mass_4digit in self.mass_equivalencies['by_mass']:
                    self.mass_equivalencies['by_mass'][mass_4digit].append(residue_pair)
                    reference_residue = self.mass_equivalencies['by_mass'][mass_4digit][0]
                    self.mass_equivalencies['by_residue'][reference_residue].append(residue_pair)
                else:
                    self.mass_equivalencies['by_mass'][mass_4digit] = [residue_pair]
                    self.mass_equivalencies['by_residue'][residue_pair] = [residue_pair]
                    self.nr_double_aa_masses[residue_pair] = { 'mass': mass, 'residues': [ aa1, aa2 ] }

                self.double_aa_masses[residue_pair] = mass
                if '[' not in residue_pair:
                    nomod_double_aa_masses[residue_pair] = mass
        #print(json.dumps(self.mass_equivalencies, indent=2, sort_keys=True))
        #exit()

        # Create all the n-terminal amino acid pairs
        self.nterm_double_aa_masses = {}
        for aa1 in self.nterm_aa_masses:
            for aa2 in self.nterm_aa_masses:
                if aa2 in nterm_residues:
                    continue
                if aa2 + aa1 in self.nterm_double_aa_masses:
                    continue
                self.nterm_double_aa_masses[aa1 + aa2] = self.nterm_aa_masses[aa1] + self.nterm_aa_masses[aa2]



        # Define some potential neutral losses
        self.neutral_losses = {
            #'water': { 'formula': 'H2O', 'residues': [ 'S', 'T', 'E', 'D' ],              # canonical
            'water': { 'formula': 'H2O', 'residues': [ 'S', 'T', 'E', 'D', 'K', 'A', 'Y', 'C[Carbamidomethyl]' ],          # observed
                'delta_mass': self.atomic_masses['hydrogen'] * 2 + self.atomic_masses['oxygen'] },
            #'ammonia': { 'formula': 'NH3', 'residues': [ 'R', 'K', 'N', 'Q' ],              # canonical
            'ammonia': { 'formula': 'NH3', 'residues': [ 'R', 'K', 'N', 'Q', 'G' ],          # observed
                'delta_mass': self.atomic_masses['nitrogen'] + self.atomic_masses['hydrogen'] * 3 },
            'carbon monoxide': { 'formula': 'CO', 'residues': [ 'b-ion' ],
                'delta_mass': self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
            #'phosphoric acid': { 'formula': 'H3PO4', 'residues': [ 'S[Phospho]', 'T[Phospho]', 'Y[Phospho]' ], # Removed this in favor of letting HPO3 and H2O work together. Otherwise, you can get both phosphoric acid and metaphosphoric acid on the same residue. You'd somehow need to encode which losses are exclusive and which are combinable
            #    'delta_mass': self.atomic_masses['hydrogen'] * 3 + self.atomic_masses['phosphorus'] + self.atomic_masses['oxygen'] * 4 },
            'metaphosphoric acid': { 'formula': 'HPO3', 'residues': [ 'S[Phospho]', 'T[Phospho]', 'Y[Phospho]' ],
                'delta_mass': self.atomic_masses['hydrogen'] * 1 + self.atomic_masses['phosphorus'] + self.atomic_masses['oxygen'] * 3 },
            'methanesulfenic acid': { 'formula': 'CH4OS', 'residues': [ 'M[Oxidation]' ],
                'delta_mass': self.atomic_masses['carbon'] * 1 + self.atomic_masses['hydrogen'] * 4 + self.atomic_masses['oxygen'] * 1 + self.atomic_masses['sulfur'] * 1 },
        }

        # Also key neutral losses by residue and by formula
        for neutral_loss_name, neutral_loss in self.neutral_losses.items():
            for residue in neutral_loss['residues']:
                if residue not in self.neutral_losses_by_residue:
                    self.neutral_losses_by_residue[residue] = []
                self.neutral_losses_by_residue[residue].append(neutral_loss)
            self.neutral_losses_by_formula[neutral_loss['formula']] = neutral_loss

        self.immonium_ions = {}
        for residue,residue_mass in self.aa_masses.items():
            immonium_name = f"I{residue}"
            immonium_mass = residue_mass - ( self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] ) + self.atomic_masses['proton']
            self.immonium_ions[immonium_name] = immonium_mass
            # FIXME: What should the rest of the immonium losses be? See Immonium email from Tytus
            #immonium_name = f"I{residue}-C3H6"
            #immonium_mass = immonium_mass - ( self.atomic_masses['carbon'] * 3 + self.atomic_masses['hydrogen'] * 6 )
            #self.immonium_ions[immonium_name] = immonium_mass

        self.low_mass_ions = copy.deepcopy(self.immonium_ions)
        for residue in [ 'R', 'K' ]:
            ion_name = f"y1({residue})"
            ion_mass = self.aa_masses[residue] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = f"y1({residue})-H20"
            ion_mass = self.aa_masses[residue] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

        # Add to the low mass ions, all the b2 and a2 ions, which are very common and make for a nice calibrator
        # But exclude a few that seem to generate some side noise
        #double_aas_to_skip = { 'QH', 'YY', 'LW' }
        for residue_pair,mass in nomod_double_aa_masses.items():
            #if residue_pair in double_aas_to_skip or 'M' in residue_pair:
            #    continue
            ion_name = f"b2({residue_pair})"
            ion_mass = mass + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = f"a2({residue_pair})"
            ion_mass = mass - self.atomic_masses['carbon'] - self.atomic_masses['oxygen'] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass


        # Define reporter ions to look for
        self.reporter_ions = {
            'TMT6/10_126': { 'type': 'TMT', 'mz': 126.127726 },     # clean in TMT6
            'TMT6/10_127N': { 'type': 'TMT', 'mz': 127.124761 },    # clean in TMT6
            'TMT10_127C': { 'type': 'TMT', 'mz': 127.131081 },      # Need to first find and calibrate on TMT6 and only then do a narrow look here
            'TMT10_128N': { 'type': 'TMT', 'mz': 128.128116 },
            'TMT6/10_128C': { 'type': 'TMT', 'mz': 128.134436 },   # clean in TMT6
            'TMT6/10_129N': { 'type': 'TMT', 'mz': 129.131471 },   # clean in TMT6
            'TMT10_129C': { 'type': 'TMT', 'mz': 129.137790 },
            'TMT10_130N': { 'type': 'TMT', 'mz': 130.134825 },
            'TMT6/10_130C': { 'type': 'TMT', 'mz': 130.141145 },
            'TMT6/10_131': { 'type': 'TMT', 'mz': 131.138180 },
            'TMT11_131C': { 'type': 'TMT', 'mz': 131.138180 },
            'TMT_nterm': { 'type': 'TMT', 'mz': 229.162932 + self.atomic_masses['proton'] },
            'TMTzero_nterm': { 'type': 'TMT', 'mz': 224.152478 + self.atomic_masses['proton'] },

            # My old numbers from somewhere
            'iTRAQ4_114': { 'type': 'iTRAQ', 'mz': 114.11068 },
            'iTRAQ4_115': { 'type': 'iTRAQ', 'mz': 115.107715 },
            'iTRAQ4_116': { 'type': 'iTRAQ', 'mz': 116.111069 },
            'iTRAQ4_117': { 'type': 'iTRAQ', 'mz': 117.114424 },
            'iTRAQ4_nterm_114': { 'type': 'iTRAQ4', 'mz': 144.105918 + self.atomic_masses['proton'] },
            'iTRAQ4_nterm_115': { 'type': 'iTRAQ4', 'mz': 144.099599 + self.atomic_masses['proton'] },
            'iTRAQ4_nterm_1167': { 'type': 'iTRAQ4', 'mz': 144.102063 + self.atomic_masses['proton'] },
            'iTRAQ4_nterm_H2O': { 'type': 'iTRAQ4', 'mz': 144.102063 + self.atomic_masses['proton'] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] },

            # Jimmy's numbers from https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php
            #'iTRAQ4_114': { 'type': 'iTRAQ4', 'mz': 114.1112 },
            #'iTRAQ4_115': { 'type': 'iTRAQ4', 'mz': 115.1083 },
            #'iTRAQ4_116': { 'type': 'iTRAQ4', 'mz': 116.1116 },
            #'iTRAQ4_117': { 'type': 'iTRAQ4', 'mz': 117.1150 },
            #'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': (144.105918 + 144.099599 + 144.102063 + 144.102063)/4 + self.atomic_masses['proton'] },

            # My old numbers from somewhere
            'iTRAQ8_118': { 'type': 'iTRAQ8', 'mz': 118.111459 },   #confounder?
            'iTRAQ8_119': { 'type': 'iTRAQ8', 'mz': 119.114814 },
            'iTRAQ8_121': { 'type': 'iTRAQ8', 'mz': 121.121524 },
            'iTRAQ8_113': { 'type': 'iTRAQ8', 'mz': 113.107325 },   #confounder??

            # Numbers from Jimmy: 113.1078	114.1112	115.1082	116.1116	117.1149	118.1120	119.1153	121.1220
        }

        self.terminal_modifications = {
            'Acetyl': { 'terminus': 'nterm', 'name:': 'Acetyl', 'mass': 42.010565, 'formula': '+C2H2O', 'frequency': 'common', 'type': 'natural', 'is_labile': True },
            'Carbamyl': { 'terminus': 'nterm', 'name:': 'Carbamyl', 'mass': 43.0058, 'formula': '+HCNO', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
            'Methyl': { 'terminus': 'cterm', 'name:': 'Methyl', 'mass': 14.0157, 'formula': '+H2C', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
            #'TMT': { 'terminus': 'nterm', 'name:': 'TMT6plex', 'mass': 224.152478, 'frequency': 'high', 'type': 'label', 'is_labile': True },
            #'TMT6plex': { 'terminus': 'nterm', 'name:': 'TMT6plex', 'mass': 229.162932, 'frequency': 'high', 'type': 'label', 'is_labile': True },
            'iTRAQ4plex': { 'terminus': 'nterm', 'name:': 'iTRAQ4plex', 'mass': 144.102063, 'frequency': 'high', 'type': 'label', 'is_labile': True },
            #'iTRAQ8plex': { 'terminus': 'nterm', 'name:': 'iTRAQ8plex', 'mass': 304.205360, 'frequency': 'high', 'type': 'label', 'is_labile': True },
        }




        # Add in other misc low mass known ions: Use this: https://www.chemcalc.org/mf-finder
        other_ions = {
            '_Cytosine': 112.050538,
            '_Uracil': 113.034554,
            '_Thymine': 127.050204,
            '_Adenine': 136.061772,
            '_Guanine': 152.056686,
            '_UrocanicAcid': 139.050204,
            '_Xanthine': 153.040702,
            '_Dimethylglycine': self.calc_mass(self.parse_atomic_formula('C4H9NO2')) + self.atomic_masses['proton'],
            '_Methyladenine': self.calc_mass(self.parse_atomic_formula('C6H7N5')) + self.atomic_masses['proton'],
            '_C5H9NO': 100.07567,
            '_C8H12NO4': 186.07608,
            #'_C8H14NO5': 204.08665, HexNAc
            '_C7H8NO2': 138.05495,
            '_Hexose': 163.0601,
            '_HexNAc': 204.0869,
            '_Deoxyhexose': 147.0652,
            '_Pentose': 133.0496,
            '_NeuNAc': 292.1027,
            '_NeuGc': 308.0976,
            '_Oxonium-N': 204.0872,
            '_Oxonium-NH': 366.1397,
            '_Oxonium-NHH': 528.1935,
            '_Oxonium-NHS': 657.2357,
            '_C6H6O2': self.calc_mass(self.parse_atomic_formula('C6H6O2')) + self.atomic_masses['proton']
        }
        for other_ion_name,other_ion_mz in other_ions.items():
            self.low_mass_ions[other_ion_name] = other_ion_mz

        # Define a set of additional mass modifications of immonium ions from NIST
        self.aa_immonium_losses = {
            'G': [],
            'A': [],
            'S': [],
            'P': [],
            'V': [ '-CH2-NH3', '-NH3', '+CO-NH3-CH2' ],
            'T': [ '+CO-NH3'],
            'C': [],
            'L': [ '-C3H6', '-CH2' ],
            'I': [ '-C3H6', '-CH2' ],
            'N': [ '-NH3' ],
            'D': [ '-H2O' ],
            'Q': [ '-CO-NH3', '-NH3', '+CO'],
            'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '+CO+H2ON2', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO+H2O'],
            'E': [],
            'M': [ '-C2H2-NH3'],
            'H': [ '-CH2N', '+CO-NH2', '+CO-NH3', '+CO-NH', '+CO+H2O' ],
            'F': [ '-CH3N'],
            'R': [ '-C3H6N2', '-CH5N3', '-CH6N2', '-C2H4N2', '-CH2N2', '-CH3N', '-NH3', '-C4H7N', '+H2O+H2O-N3H7', '+CO+H2O' ],
            'Y': [ '-CO-NH3', '-CH3N' ],
            'W': [ '+CO', '-C4H6N2', '-C2H4N', '-CH3N', '-CHN', '+CO-NH3', '-NH3'],
        }

        #### Add them to the immonium ion masses
        for aa,immonium_loss_list in self.aa_immonium_losses.items():
            for loss_formula in immonium_loss_list:
                diff_atoms = self.parse_atomic_formula(loss_formula)
                diff_mass = self.calc_mass(diff_atoms)
                mz = self.immonium_ions[f"I{aa}"] + diff_mass
                ion_name = f"I{aa}{loss_formula}"
                self.low_mass_ions[ion_name] = mz
                #print(f"{ion_name}={mz}")



    def calc_mass(self,molecule):
        mass = 0.0
        for atom,copies in molecule.items():
            mass += self.atomic_masses_by_letter[atom] * copies
        return mass


    def parse_atomic_formula(self,input_string):
        atoms = {}
        atom = ''
        number_buffer = ''
        direction = 1
        for character in input_string:
            #print(f"**{character}")
            match = re.match(r'[A-Z]',character)
            if match:
                if atom == '':
                    atoms[character] = direction
                else:
                    if number_buffer > '':
                        atoms[atom] += direction * (int(number_buffer) - 1)
                    if character not in atoms:
                        atoms[character] = 0
                    atoms[character] += direction
                atom = character
                number_buffer = ''
            elif character == '+':
                direction = 1
            elif character == '-':
                direction = -1
            else:
                match = re.match(r'\d',character)
                if match:
                    number_buffer += character
                else:
                    print(f"ERROR: Unable to parse character {character}")

        if number_buffer > '':
            atoms[atom] += direction * (int(number_buffer) - 1)

        return atoms


    def subtract_atomic_formula(self,molecule, loss):
        result = copy.copy(molecule)
        for atom,copies in loss.items():
            if atom not in result:
                result[atom] = 0
            result[atom] -= copies
        return result



####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Creates a set of data structures of reference masses')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # Define all the reference mass information
    mass_reference = MassReference()

    # Print out a few items
    print(f"Mass of a proton = {mass_reference.atomic_masses['proton']}")


#### For command line usage
if __name__ == "__main__": main()
