#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

import argparse
import copy
import re
import json

# Common background ions:
# https://www.waters.com/webassets/cms/support/docs/bkgrnd_ion_mstr_list.pdf
# 445 is a common one that Rob and Mike mention

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
        self.special_label_losses = None
        self.terminal_modifications = None

        self.neutral_losses = {}
        self.neutral_losses_by_residue = {}
        self.neutral_losses_by_formula = {}

        self.prepare_mass_tables()


    ####################################################################################################
    def prepare_mass_tables(self, labels=None):

        # Define a subset of useful atomic masses and the proton
        self.atomic_masses = {
            'proton': 1.00727646688,
            'electron': 0.000548,
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
        for aa, formula in self.aa_formulas.items():
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

        # Define a basic set of mass modifications that we would encounter
        modifications = {
            'Carbamidomethyl': { 'delta_mass': 57.021464, 'formula': '+H3C2NO', 'residues': [ 'C', 'K' ] },
            'Oxidation': { 'delta_mass': 15.994915, 'formula': '+O', 'residues': [ 'M', 'W', 'P', 'H', 'Y' ] },
            'Dioxidation': { 'delta_mass': 2 * 15.994915, 'formula': '+O2', 'residues': [ 'W' ] },
            'Trp->Kynurenin': { 'delta_mass': 3.994915, 'formula': '+O-C', 'residues': [ 'W' ] },
            'Trioxidation': { 'delta_mass': 3 * 15.994915, 'formula': '+O3', 'residues': [ 'C' ] },
            'Phospho': { 'delta_mass': 79.966331, 'formula': '+HPO3', 'residues': [ 'S', 'T', 'Y', 'H' ], 'neutral_losses': [ 97.976896 ] },
            'Deamidated': { 'delta_mass': 0.984016, 'formula': '-HN+O', 'residues': [ 'N' ] }, # Q[Deamidated] is the same as E, so don't bother here
            'Acetyl': { 'delta_mass': 42.010565, 'formula': '+C2H2O', 'residues': [ 'K' ] },
            'Dimethyl': { 'delta_mass': 28.031300, 'formula': '+C2H4', 'residues': [ 'K' ] },
            'Formyl': { 'delta_mass': 27.9949, 'formula': '+CO', 'residues': [ 'S', 'T', 'K' ] },   #is K[Formyl] the same as R???
            'Carbonyl': { 'delta_mass': 13.979265, 'formula': '+O-H2', 'residues': [ 'H' ] },
            'Methylation': { 'delta_mass': 14.015650, 'residues': [ 'K', 'R', 'H' ] },
            'Carbamyl': { 'delta_mass': 43.005814, 'residues': [ 'K' ] },
            'Sulfide': { 'delta_mass': 31.972071, 'residues': [ 'S' ] },
            #'Glycerinyl': { 'delta_mass': 88.016, 'formula': '+C', 'residues': [ 'K' ] },
            #'Pyrophospho': { 'delta_mass': 159.932662, 'residues': [ 'S', 'T', 'Y' ], 'neutral_losses': [ 177.943227, 79.966331 ] },
            #'Beta-methythiolation': { 'delta_mass': 45.987721, 'residues': [ 'C' ] },
            #'Acetylation': { 'delta_mass': 42.010565, 'residues': [ 'K' ] },
            #'Methylthio': { 'delta_mass': 45.987721, 'residues': [ 'K', 'D', 'N', 'C' ] },
            #'Ubiquitination': { 'delta_mass': 114.042927, 'residues': [ 'K' ] },
        }

        #### If a labeling type was specified, then add modifications for that label
        if labels is not None:
            known_labeling_modifications = {
                #'TMT': { 'delta_mass': 224.152478, 'residues': [ 'K' ] },
                'TMT6': { 'delta_mass': 229.162932, 'residues': [ 'K' ] },
                'TMT6plex': { 'delta_mass': 229.162932, 'residues': [ 'K' ] },
                'TMTpro': { 'delta_mass': 304.207146, 'residues': [ 'K', 'S', 'T', 'R' ] },
                'iTRAQ4': { 'delta_mass': 144.102063, 'residues': [ 'K' ] },
                'iTRAQ4plex': { 'delta_mass': 144.102063, 'residues': [ 'K' ] },
                #'mtraq': { 'delta_mass': 140.094963, 'residues': [ 'K' ] },
            }
            for label in labels:
                if label in known_labeling_modifications:
                    modifications[label] = known_labeling_modifications[label]

        #### Then for all of the modifications, expand to a full list of mass deltas
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

                self.double_aa_masses[residue_pair] = { 'mass': mass, 'residues': [ aa1, aa2 ] }
                if '[' not in residue_pair:
                    nomod_double_aa_masses[residue_pair] = { 'mass': mass, 'residues': [ aa1, aa2 ] }
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
            #'lysine +CO-H2O': { 'formula': 'H2O+CO', 'residues': [ 'K' ],
            #    'delta_mass': self.atomic_masses['hydrogen'] * 2 - self.atomic_masses['carbon'] * 1 },
            #'water': { 'formula': 'H2O', 'residues': [ 'S', 'T', 'E', 'D' ],              # canonical
            'water': { 'formula': 'H2O', 'residues': [ 'S', 'T', 'E', 'D', 'K', 'A', 'Y', 'C[Carbamidomethyl]', 'S[Phospho]', 'T[Phospho]', 'H' ],          # observed
                'delta_mass': self.atomic_masses['hydrogen'] * 2 + self.atomic_masses['oxygen'] },
            #'weird': { 'formula': 'H2O+CO', 'residues': [ 'S', 'T', 'E', 'D', 'K', 'A', 'Y', 'C[Carbamidomethyl]' ],          # observed
            #    'delta_mass': self.atomic_masses['hydrogen'] * 2 - self.atomic_masses['carbon'] },
            #'ammonia': { 'formula': 'NH3', 'residues': [ 'R', 'K', 'N', 'Q' ],              # canonical
            'ammonia': { 'formula': 'NH3', 'residues': [ 'R', 'K', 'N', 'Q', 'G' ],          # observed
                'delta_mass': self.atomic_masses['nitrogen'] + self.atomic_masses['hydrogen'] * 3 },
            'carbon monoxide': { 'formula': 'CO', 'residues': [ 'b-ion' ],
                'delta_mass': self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
            #'phosphoric acid': { 'formula': 'H3PO4', 'residues': [ 'S[Phospho]', 'T[Phospho]', 'Y[Phospho]' ], # Removed this in favor of letting HPO3 and H2O work together. Otherwise, you can get both phosphoric acid and metaphosphoric acid on the same residue. You'd somehow need to encode which losses are exclusive and which are combinable
            #    'delta_mass': self.atomic_masses['hydrogen'] * 3 + self.atomic_masses['phosphorus'] + self.atomic_masses['oxygen'] * 4 },
            'metaphosphoric acid': { 'formula': 'HPO3', 'residues': [ 'S[Phospho]', 'T[Phospho]', 'Y[Phospho]', 'H[Phospho]' ],
                'delta_mass': self.atomic_masses['hydrogen'] * 1 + self.atomic_masses['phosphorus'] + self.atomic_masses['oxygen'] * 3 },
            'methanesulfenic acid': { 'formula': 'CH4OS', 'residues': [ 'M[Oxidation]' ],
                'delta_mass': self.atomic_masses['carbon'] * 1 + self.atomic_masses['hydrogen'] * 4 + self.atomic_masses['oxygen'] * 1 + self.atomic_masses['sulfur'] * 1 },
            'Tryp Mannose for Kristian': { 'formula': 'C4H8O4', 'residues': [ 'W[Hex]' ],
                'delta_mass': self.atomic_masses['carbon'] * 4 + self.atomic_masses['hydrogen'] * 8 + self.atomic_masses['oxygen'] * 4 },
            'Observed C[Trioxidation] loss': { 'formula': 'H2O3S', 'residues': [ 'C[Trioxidation]' ],
                'delta_mass': self.atomic_masses['hydrogen'] * 2 + self.atomic_masses['oxygen'] * 3 + self.atomic_masses['sulfur'] * 1},
            'Potential C[Trioxidation] loss': { 'formula': 'CO4S', 'residues': [ 'C[Trioxidation]' ],
                'delta_mass': self.atomic_masses['carbon'] * 1 + self.atomic_masses['oxygen'] * 4 + self.atomic_masses['sulfur'] * 1},
            'Crazy potential C[Trioxidation] loss': { 'formula': 'C5H5NO2', 'residues': [ 'C[Trioxidation]' ],
                'delta_mass': self.atomic_masses['carbon'] * 5 + self.atomic_masses['hydrogen'] * 5 + self.atomic_masses['nitrogen'] * 1 + self.atomic_masses['oxygen'] * 2 },
            'Observed C[Dioxidation] loss': { 'formula': 'H2O2S', 'residues': [ 'C[Dioxidation]' ],
                'delta_mass': self.atomic_masses['hydrogen'] * 2 + self.atomic_masses['oxygen'] * 2 + self.atomic_masses['sulfur'] * 1},
            'Observed K[Carbamyl] loss': { 'formula': 'HCNO', 'residues': [ 'K[Carbamyl]' ],
                'delta_mass': self.atomic_masses['hydrogen'] * 1 + self.atomic_masses['carbon'] * 1 + self.atomic_masses['nitrogen'] * 1 + self.atomic_masses['oxygen'] * 1},
            'Crazy K[Carbamyl] loss. All but NH3': { 'formula': 'C7H10N2O2', 'residues': [ 'K[Carbamyl]' ],
                'delta_mass': self.atomic_masses['carbon'] * 7 + self.atomic_masses['hydrogen'] * 10 + self.atomic_masses['nitrogen'] * 2 + self.atomic_masses['oxygen'] * 2 },

            # ADP-Ribosyl
            'ADP-Ribosyl-H2O loss': { 'formula': 'ADP-Ribose-H2O', 'residues': [ 'S[ADP-Ribosyl]', 'T[ADP-Ribosyl]', 'R[ADP-Ribosyl]' ],
                'delta_mass': 542.067 - self.atomic_masses['proton'] },
            'ADP-Ribosyl-H2O-H2O loss': { 'formula': 'ADP-Ribose-H2O-H2O', 'residues': [ 'S[ADP-Ribosyl]', 'T[ADP-Ribosyl]', 'R[ADP-Ribosyl]' ],
                'delta_mass': 524.057 - self.atomic_masses['proton'] },
            'AMP loss': { 'formula': 'AMP', 'residues': [ 'S[ADP-Ribosyl]', 'T[ADP-Ribosyl]', 'R[ADP-Ribosyl]' ],
                'delta_mass': 348.069 - self.atomic_masses['proton'] },

            # ETD?
            #'oxygen': { 'formula': 'O', 'residues': [ 'M[Oxidation]' ],
            #    'delta_mass': self.atomic_masses['oxygen'] },

        }

        # Also key neutral losses by residue and by formula
        for neutral_loss_name, neutral_loss in self.neutral_losses.items():
            for residue in neutral_loss['residues']:
                if residue not in self.neutral_losses_by_residue:
                    self.neutral_losses_by_residue[residue] = []
                self.neutral_losses_by_residue[residue].append(neutral_loss)
            self.neutral_losses_by_formula[neutral_loss['formula']] = neutral_loss


        #### A few immonium ions are isobaric with another, so have a list of ones that we should just skip
        immonium_ions_to_exclude = {
            'IK[Formyl]': True,     # This interferes with IK+CO
            'IK[Dimethyl]': True,   # This probably isn't a thing, and could conflict with r[TMT129C]
            'IK[Carbamidomethyl]': True,   # This probably isn't a thing, and is isobaric with a2{GK} which is probably more likely
            'IT[Formyl]': True,             # This probably isn't a thing, and is isobaric with IE which is more likely
            }

        self.immonium_ions = {}
        for residue, residue_mass in self.aa_masses.items():
            immonium_name = f"I{residue}"
            if immonium_name in immonium_ions_to_exclude:
                continue
            immonium_mass = residue_mass - ( self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] ) + self.atomic_masses['proton']
            self.immonium_ions[immonium_name] = immonium_mass


        self.low_mass_ions = copy.deepcopy(self.immonium_ions)
        for residue in [ 'R', 'K' ]:
            ion_name = "y1" + "{" + residue +"}"
            ion_mass = self.aa_masses[residue] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = "y1" + "{" + residue +"}-H2O"
            ion_mass = self.aa_masses[residue] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

        # Add to the low mass ions all the b2 and a2 ions, which are very common
        for residue_pair_str,residue_pair in nomod_double_aa_masses.items():
            mass = residue_pair['mass']
            #if residue_pair in double_aas_to_skip or 'M' in residue_pair:
            #    continue
            ion_name = "b2" + "{" + residue_pair_str +"}"
            ion_mass = mass + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass

            ion_name = "a2" + "{" + residue_pair_str +"}"
            ion_mass = mass - self.atomic_masses['carbon'] - self.atomic_masses['oxygen'] + self.atomic_masses['proton']
            self.low_mass_ions[ion_name] = ion_mass


        # Define reporter ions to look for
        self.reporter_ions = {
            # Jimmy's numbers from https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php
            'TMT126': { 'type': 'TMT', 'mz': 126.127726 },
            'TMT127N': { 'type': 'TMT', 'mz': 127.124761 },
            'TMT127C': { 'type': 'TMT', 'mz': 127.131081 },
            'TMT128N': { 'type': 'TMT', 'mz': 128.128116 },
            'TMT128C': { 'type': 'TMT', 'mz': 128.134436 },
            'TMT129N': { 'type': 'TMT', 'mz': 129.131471 },
            'TMT129C': { 'type': 'TMT', 'mz': 129.137790 },
            'TMT130N': { 'type': 'TMT', 'mz': 130.134825 },
            'TMT130C': { 'type': 'TMT', 'mz': 130.141145 },
            'TMT131N': { 'type': 'TMT', 'mz': 131.138180 },
            'TMT131C': { 'type': 'TMT', 'mz': 131.1445 },
            'TMT132N': { 'type': 'TMT', 'mz': 132.141535 },
            'TMT132C': { 'type': 'TMT', 'mz': 132.147855 },
            'TMT133N': { 'type': 'TMT', 'mz': 133.14489 },
            'TMT133C': { 'type': 'TMT', 'mz': 133.15121 },
            'TMT134N': { 'type': 'TMT', 'mz': 134.148245 },
            'TMT134C': { 'type': 'TMT', 'mz': 134.154565 },
            'TMT135N': { 'type': 'TMT', 'mz': 135.151600 },

            #### TMT-35 peaks
            'TMT127D': { 'type': 'TMT', 'mz': 127.134003 },
            'TMT128ND': { 'type': 'TMT', 'mz': 128.131038 },
            'TMT128CD': { 'type': 'TMT', 'mz': 128.137358 },
            'TMT129ND': { 'type': 'TMT', 'mz': 129.134393 },
            'TMT129CD': { 'type': 'TMT', 'mz': 129.140713 },
            'TMT130ND': { 'type': 'TMT', 'mz': 130.137748 },
            'TMT130CD': { 'type': 'TMT', 'mz': 130.144068 },
            'TMT131ND': { 'type': 'TMT', 'mz': 131.141103 },
            'TMT131CD': { 'type': 'TMT', 'mz': 131.147423 },
            'TMT132ND': { 'type': 'TMT', 'mz': 132.144458 },
            'TMT132CD': { 'type': 'TMT', 'mz': 132.150778 },
            'TMT133ND': { 'type': 'TMT', 'mz': 133.147813 },
            'TMT133CD': { 'type': 'TMT', 'mz': 133.154133 },
            'TMT134ND': { 'type': 'TMT', 'mz': 134.151171 },
            'TMT134CD': { 'type': 'TMT', 'mz': 134.157491 },
            'TMT135ND': { 'type': 'TMT', 'mz': 135.154526 },
            'TMT135CD': { 'type': 'TMT', 'mz': 135.160846 },

            #### When dimethyl hangs on the TMTpro tag!
            #'TMT126+C2H4': { 'type': 'TMT', 'mz': 126.127726 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT127N+C2H4': { 'type': 'TMT', 'mz': 127.124761 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT127C+C2H4': { 'type': 'TMT', 'mz': 127.131081 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT128N+C2H4': { 'type': 'TMT', 'mz': 128.128116 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT128C+C2H4': { 'type': 'TMT', 'mz': 128.134436 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT129N+C2H4': { 'type': 'TMT', 'mz': 129.131471 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT129C+C2H4': { 'type': 'TMT', 'mz': 129.137790 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT130N+C2H4': { 'type': 'TMT', 'mz': 130.134825 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT130C+C2H4': { 'type': 'TMT', 'mz': 130.141145 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT131N+C2H4': { 'type': 'TMT', 'mz': 131.138180 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT131C+C2H4': { 'type': 'TMT', 'mz': 131.1445 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT132N+C2H4': { 'type': 'TMT', 'mz': 132.141535 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT132C+C2H4': { 'type': 'TMT', 'mz': 132.147855 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT133N+C2H4': { 'type': 'TMT', 'mz': 133.14489 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT133C+C2H4': { 'type': 'TMT', 'mz': 133.15121 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT134N+C2H4': { 'type': 'TMT', 'mz': 134.148245 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT134C+C2H4': { 'type': 'TMT', 'mz': 134.154565 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMT135N+C2H4': { 'type': 'TMT', 'mz': 135.151600 + 2 * self.atomic_masses['carbon'] + 4 * self.atomic_masses['hydrogen'] },
            #'TMTpro+C2H4': { 'type': 'TMT', 'mz': 304.207146 + self.atomic_masses['proton'] + 4 * self.atomic_masses['hydrogen'] + 2 * self.atomic_masses['carbon'] },

            'TMTzero': { 'type': 'TMT', 'mz': 224.152478 + self.atomic_masses['proton'] },
            'TMT2': { 'type': 'TMT', 'mz': 225.155833 + self.atomic_masses['proton'] },
            'TMT6plex': { 'type': 'TMT', 'mz': 229.162932 + self.atomic_masses['proton'] },
            'TMT6plex+H2O': { 'type': 'TMT', 'mz': 229.162932 + self.atomic_masses['proton'] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] },
            'TMTprozero': { 'type': 'TMT', 'mz': 295.189592 + self.atomic_masses['proton'] },
            'TMTpro': { 'type': 'TMT', 'mz': 304.207146 + self.atomic_masses['proton'] },
            'TMTpro+H20': { 'type': 'TMT', 'mz': 304.207146 + self.atomic_masses['proton'] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] },
            'TMTpro+NH3': { 'type': 'TMT', 'mz': 304.207146 + self.atomic_masses['proton'] + 3 * self.atomic_masses['hydrogen'] + self.atomic_masses['nitrogen'] },
            'TMTpro+C2NH5': { 'type': 'TMT', 'mz': 304.207146 + self.atomic_masses['proton'] + 5 * self.atomic_masses['hydrogen'] + self.atomic_masses['nitrogen'] + 2 * self.atomic_masses['carbon'] },

            # My old numbers from somewhere
            'iTRAQ114': { 'type': 'iTRAQ', 'mz': 114.11068 },
            'iTRAQ115': { 'type': 'iTRAQ', 'mz': 115.107715 },
            'iTRAQ116': { 'type': 'iTRAQ', 'mz': 116.111069 },
            'iTRAQ117': { 'type': 'iTRAQ', 'mz': 117.114424 },
            'iTRAQ4plex_114': { 'type': 'iTRAQ4', 'mz': 144.105918 + self.atomic_masses['proton'] },
            'iTRAQ4plex_115': { 'type': 'iTRAQ4', 'mz': 144.099599 + self.atomic_masses['proton'] },
            'iTRAQ4plex_1167': { 'type': 'iTRAQ4', 'mz': 144.102063 + self.atomic_masses['proton'] },
            'iTRAQ4plex+H2O': { 'type': 'iTRAQ4', 'mz': 144.102063 + self.atomic_masses['proton'] + 2 * self.atomic_masses['hydrogen'] + self.atomic_masses['oxygen'] },


            # Jimmy's numbers from https://proteomicsresource.washington.edu/protocols03/isotopic_labeling.php
            #'iTRAQ4_114': { 'type': 'iTRAQ4', 'mz': 114.1112 },
            #'iTRAQ4_115': { 'type': 'iTRAQ4', 'mz': 115.1083 },
            #'iTRAQ4_116': { 'type': 'iTRAQ4', 'mz': 116.1116 },
            #'iTRAQ4_117': { 'type': 'iTRAQ4', 'mz': 117.1150 },
            #'iTRAQ4_nterm': { 'type': 'iTRAQ4', 'mz': (144.105918 + 144.099599 + 144.102063 + 144.102063)/4 + self.atomic_masses['proton'] },

            # My old numbers from somewhere
            'iTRAQ113': { 'type': 'iTRAQ8', 'mz': 113.107325 },
            'iTRAQ118': { 'type': 'iTRAQ8', 'mz': 118.111459 },   #confounder?
            'iTRAQ119': { 'type': 'iTRAQ8', 'mz': 119.114814 },
            'iTRAQ121': { 'type': 'iTRAQ8', 'mz': 121.121524 },
            # Numbers from Jimmy: 113.1078	114.1112	115.1082	116.1116	117.1149	118.1120	119.1153	121.1220

            #'mTRAQNterm': { 'type': 'mTRAQ', 'mz': 140.094963 + self.atomic_masses['proton'] }, # this is the same as 0@_"a2(AP)" and it will rarely be mTRAQ. FIXME
        }


        self.special_label_losses = {
            'TMT6plex': {
                '[TMT126]-CO': { 'type': 'TMT', 'mz': 126.127726 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT127N]-CO': { 'type': 'TMT', 'mz': 127.124761 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT127C]-CO': { 'type': 'TMT', 'mz': 127.131081 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT128N]-CO': { 'type': 'TMT', 'mz': 128.128116 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT128C]-CO': { 'type': 'TMT', 'mz': 128.134436 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT129N]-CO': { 'type': 'TMT', 'mz': 129.131471 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT129C]-CO': { 'type': 'TMT', 'mz': 129.137790 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT130N]-CO': { 'type': 'TMT', 'mz': 130.134825 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT130C]-CO': { 'type': 'TMT', 'mz': 130.141145 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT131N]-CO': { 'type': 'TMT', 'mz': 131.138180 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT131C]-CO': { 'type': 'TMT', 'mz': 131.1445 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT132N]-CO': { 'type': 'TMT', 'mz': 132.141535 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT132C]-CO': { 'type': 'TMT', 'mz': 132.147855 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT133N]-CO': { 'type': 'TMT', 'mz': 133.14489 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT133C]-CO': { 'type': 'TMT', 'mz': 133.15121 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT134N]-CO': { 'type': 'TMT', 'mz': 134.148245 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT134C]-CO': { 'type': 'TMT', 'mz': 134.154565 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT135N]-CO': { 'type': 'TMT', 'mz': 135.151600 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
           },
            'TMTpro': {
                '[TMT126]-CO': { 'type': 'TMT', 'mz': 126.127726 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT127N]-CO': { 'type': 'TMT', 'mz': 127.124761 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT127C]-CO': { 'type': 'TMT', 'mz': 127.131081 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT128N]-CO': { 'type': 'TMT', 'mz': 128.128116 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT128C]-CO': { 'type': 'TMT', 'mz': 128.134436 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT129N]-CO': { 'type': 'TMT', 'mz': 129.131471 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT129C]-CO': { 'type': 'TMT', 'mz': 129.137790 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT130N]-CO': { 'type': 'TMT', 'mz': 130.134825 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT130C]-CO': { 'type': 'TMT', 'mz': 130.141145 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT131N]-CO': { 'type': 'TMT', 'mz': 131.138180 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT131C]-CO': { 'type': 'TMT', 'mz': 131.1445 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT132N]-CO': { 'type': 'TMT', 'mz': 132.141535 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT132C]-CO': { 'type': 'TMT', 'mz': 132.147855 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT133N]-CO': { 'type': 'TMT', 'mz': 133.14489 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT133C]-CO': { 'type': 'TMT', 'mz': 133.15121 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT134N]-CO': { 'type': 'TMT', 'mz': 134.148245 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT134C]-CO': { 'type': 'TMT', 'mz': 134.154565 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                '[TMT135N]-CO': { 'type': 'TMT', 'mz': 135.151600 - self.atomic_masses['proton'] + self.atomic_masses['carbon'] + self.atomic_masses['oxygen'] },
                #### Seems like most precursor-reporter losses have -CO, but not always, see mzspec:PXD999001:20231003-TMTpool167-Fr05-01:scan:480:[TMTpro]-SGQSEDR/2
                '[TMT126]': { 'type': 'TMT', 'mz': 126.127726 - self.atomic_masses['proton'] },
                '[TMT127N]': { 'type': 'TMT', 'mz': 127.124761 - self.atomic_masses['proton'] },
                '[TMT127C]': { 'type': 'TMT', 'mz': 127.131081 - self.atomic_masses['proton'] },
                '[TMT128N]': { 'type': 'TMT', 'mz': 128.128116 - self.atomic_masses['proton'] },
                '[TMT128C]': { 'type': 'TMT', 'mz': 128.134436 - self.atomic_masses['proton'] },
                '[TMT129N]': { 'type': 'TMT', 'mz': 129.131471 - self.atomic_masses['proton'] },
                '[TMT129C]': { 'type': 'TMT', 'mz': 129.137790 - self.atomic_masses['proton'] },
                '[TMT130N]': { 'type': 'TMT', 'mz': 130.134825 - self.atomic_masses['proton'] },
                '[TMT130C]': { 'type': 'TMT', 'mz': 130.141145 - self.atomic_masses['proton'] },
                '[TMT131N]': { 'type': 'TMT', 'mz': 131.138180 - self.atomic_masses['proton'] },
                '[TMT131C]': { 'type': 'TMT', 'mz': 131.1445 - self.atomic_masses['proton'] },
                '[TMT132N]': { 'type': 'TMT', 'mz': 132.141535 - self.atomic_masses['proton'] },
                '[TMT132C]': { 'type': 'TMT', 'mz': 132.147855 - self.atomic_masses['proton'] },
                '[TMT133N]': { 'type': 'TMT', 'mz': 133.14489 - self.atomic_masses['proton'] },
                '[TMT133C]': { 'type': 'TMT', 'mz': 133.15121 - self.atomic_masses['proton'] },
                '[TMT134N]': { 'type': 'TMT', 'mz': 134.148245 - self.atomic_masses['proton'] },
                '[TMT134C]': { 'type': 'TMT', 'mz': 134.154565 - self.atomic_masses['proton'] },
                '[TMT135N]': { 'type': 'TMT', 'mz': 135.151600 - self.atomic_masses['proton'] },

           }
        }


        #### Set up some basic terminal modifications
        self.terminal_modifications = {
            'Acetyl': { 'terminus': 'nterm', 'name:': 'Acetyl', 'mass': 42.010565, 'formula': '+C2H2O', 'frequency': 'common', 'type': 'natural', 'is_labile': True },
            'Carbamyl': { 'terminus': 'nterm', 'name:': 'Carbamyl', 'mass': 43.0058, 'formula': '+HCNO', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
            'Methyl': { 'terminus': 'cterm', 'name:': 'Methyl', 'mass': 14.0157, 'formula': '+H2C', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
            'Carbamidomethyl': { 'terminus': 'nterm', 'name:': 'Carbamidomethyl', 'mass': 57.021464, 'formula': '+H3C2NO', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
            'Formyl': { 'terminus': 'nterm', 'name:': 'Formyl', 'mass': 27.9949, 'formula': '+CO', 'frequency': 'common', 'type': 'natural', 'is_labile': False },
        }

       #### If a labeling type was specified, then add modifications for that label
        if labels is not None:
            known_labeling_modifications = {
                #'TMT': { 'terminus': 'nterm', 'name:': 'TMT6plex', 'mass': 224.152478, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                'TMT6plex': { 'terminus': 'nterm', 'name:': 'TMT6plex', 'mass': 229.162932, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                'TMTpro': { 'terminus': 'nterm', 'name:': 'TMTpro', 'mass': 304.207146, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                'iTRAQ4plex': { 'terminus': 'nterm', 'name:': 'iTRAQ4plex', 'mass': 144.102063, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                'iTRAQ8plex': { 'terminus': 'nterm', 'name:': 'iTRAQ8plex', 'mass': 304.205360, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                #'mtraq': { 'terminus': 'nterm', 'name:': 'mtraq', 'mass': 140.094963, 'frequency': 'high', 'type': 'label', 'is_labile': True },
                }
            for label in labels:
                if label in known_labeling_modifications:
                    self.terminal_modifications[label] = known_labeling_modifications[label]


        # Add in other misc low mass known ions: Use this: https://www.chemcalc.org/mf-finder
        other_ions = {
            'Cytosine': 112.050538,
            'Uracil': 113.034554,
            'Thymine': 127.050204,
            'Adenine': 136.061772,
            'Guanine': 152.056686,
            'UrocanicAcid': 139.050204,
            'Xanthine': 153.040702,

            'Cytosine+C5H7O6P': 112.050538 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')),
            'Cytosine+C5H7O6P-H2O': 112.050538 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')) - self.calc_mass(self.parse_atomic_formula('H2O')),
            'Adenine+C5H7O6P': 136.061772 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')),
            'Adenine+C5H7O6P-H2O': 136.061772 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')) - self.calc_mass(self.parse_atomic_formula('H2O')),
            'Guanine+C5H7O6P': 152.056686 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')),
            'Guanine+C5H7O6P-H2O': 152.056686 + self.calc_mass(self.parse_atomic_formula('C5H7O6P')) - self.calc_mass(self.parse_atomic_formula('H2O')),

            # From ADP-Ribosyl fragmentation as per https://pubs.acs.org/doi/suppl/10.1021/jasms.0c00040/suppl_file/js0c00040_si_008.pdf
            'Adenosine-H2O': 250.093,
            'AMP': 348.069,
            'ADP': 428.035,
            'ADP-Ribose-H2O-H2O': 524.057,
            'ADP-Ribose-H2O': 542.067,
            'ADPr-Carbodiimide': 584.089,

            'Dimethylglycine': self.calc_mass(self.parse_atomic_formula('C4H9NO2')) + self.atomic_masses['proton'],
            'Methyladenine': self.calc_mass(self.parse_atomic_formula('C6H7N5')) + self.atomic_masses['proton'],

            'f{C5H9NO}': self.calc_mass(self.parse_atomic_formula('C5H9NO')) + self.atomic_masses['proton'],
            'f{C8H12NO4}': self.calc_mass(self.parse_atomic_formula('C8H12NO4')) + self.atomic_masses['proton'],
            'f{C7H8NO2}': self.calc_mass(self.parse_atomic_formula('C7H8NO2')) + self.atomic_masses['proton'],
            'f{C6H6O2}': self.calc_mass(self.parse_atomic_formula('C6H6O2')) + self.atomic_masses['proton'],

            #Interesting monosaccharide web pages:
            #https://www.ionsource.com/Card/carbo/sugar.htm
            #https://www.shimadzu.co.jp/aboutus/ms_r/archive/files/OlgScchrdTable_11.pdf
            'HexNAc-frag1': 126.05495,
            'Pentose': 133.0496,
            'HexNAc-frag2': 138.05495,
            'HexNAc-frag3': 144.0661,
            'Hexose-H2O': 145.04953,
            'Deoxyhexose': 147.0652,
            'Hexose': 163.0601,
            'HexNAc-2H2O': 168.0655,
            'HexNAc-H2O': 186.0761,
            'HexNAc': 204.08666,
            #'Oxonium-N': 204.08666, # Same as HexNAc
            'Neu5Ac-2H2O': 256.08156,
            'Neu5Ac-H2O': 274.0921,
            'Neu5Ac': 292.1026915,   # Also NeuNAc
            'NeuGc': 308.0976,
            'HexHex': 325.1129,
            'Oxonium-NH': 366.1397, # HexHexNAc
            'Oxonium-NHH': 528.1935,
            'Oxonium-NHS': 657.2357,
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
            #'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '+CO+H2ON2', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO+H2O', '+CO-H2O'],
            'K': [ '+CO-NH3', '-NH3', '+CO', '-C2H4-NH3', '-NH3', '-C4H7N', '+CO+CO-C2H3N3', '+CO-H2O'],
            'E': [],
            'M': [ '-C2H2-NH3'],
            'H': [ '-CH2N', '+CO-NH2', '+CO-NH3', '+CO-NH', '+CO+H2O' ],
            'F': [ '-CH3N'],
            'R': [ '-C3H6N2', '-CH5N3', '-CH6N2', '-C2H4N2', '-CH2N2', '-CH3N', '-NH3', '-C4H7N', '+H2O+H2O-N3H7' ],
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

    argparser = argparse.ArgumentParser(description='Class that creates and holds data structures of reference masses')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('--show_all_independent_ions', action='count', help='If set, print a table of all independent ions' )
    argparser.add_argument('--show_reporter_ions', action='count', help='If set, print the reporter ion information' )
    argparser.add_argument('--show_low_mass_ions', action='count', help='If set, print the other low mass ion information' )
    argparser.add_argument('--show_neutral_losses', action='count', help='If set, print the neutral losses information' )
    argparser.add_argument('--show_atomic_masses', action='count', help='If set, print the atomics masses' )
    params = argparser.parse_args()

    #### Set verbose
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # Define all the reference mass information
    mass_reference = MassReference()


    # Show reporter ion information
    if params.show_all_independent_ions:
        all_ions = {}
        for key in mass_reference.reporter_ions.keys():
            all_ions[key] = mass_reference.reporter_ions[key]['mz']
        for key in mass_reference.low_mass_ions.keys():
            all_ions[key] = mass_reference.low_mass_ions[key]
        for key in mass_reference.immonium_ions.keys():
            all_ions[key] = mass_reference.immonium_ions[key]
        all_ions_list = []
        for key, value in all_ions.items():
            all_ions_list.append([key, value])
        all_ions_list.sort(key=lambda x: x[1])
        print(f"ion_name\tmz")
        for item in all_ions_list:
            print(f"{item[0]}\t{item[1]}")
        return

    # Show reporter ion information
    if params.show_reporter_ions:
        for key in sorted(mass_reference.reporter_ions.keys()):
            print(f"{key}\t{mass_reference.reporter_ions[key]['mz']}")
        return

    # Show other low mass ion information
    if params.show_low_mass_ions:
        for key in sorted(mass_reference.low_mass_ions.keys()):
            print(f"{key}\t{mass_reference.low_mass_ions[key]}")
        return

    # Show neutral losses information
    if params.show_neutral_losses:
        print(json.dumps(mass_reference.neutral_losses, indent=2, sort_keys=True))
        return

    # Show atomic mass information
    if params.show_atomic_masses:
        print("Atom masses:")
        for item,mass in mass_reference.atomic_masses.items():
            print(f"{item:15s}\t{mass}")
        return

    argparser.print_help()

#### For command line usage
if __name__ == "__main__": main()
