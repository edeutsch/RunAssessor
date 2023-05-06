#!/usr/bin/env python3
import sys
def eprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)

#### Import some standard modules
import os
import argparse
import os.path
import re


####################################################################################################
#### Peptidoform class
class Peptidoform:

    ####################################################################################################
    #### Constructor
    def __init__(self, peptidoform_string=None, verbose=0):

        # Set verbosity
        if verbose is None:
            verbose = 0
        self.verbose = verbose

        self.nterm = { 'string': '' }
        self.cterm = { 'string': '' }
        self.unlocalized_mods = []
        self.residues = []
        self.stripped_sequence = ''

        self.peptidoform_string = peptidoform_string

        if peptidoform_string is not None:
            self.parse_peptidoform_string()


    ####################################################################################################
    #### Decompose a supplied peptidoform string into individual elements
    def parse_peptidoform_string(self, peptidoform_string=None):

        self.nterm = { 'string': '', 'name': '' }
        self.cterm = { 'string': '', 'name': '' }
        self.unlocalized_mods = []
        self.residues = []

        # Use the passed peptidoform_string or else the previously supplied one or return None
        if peptidoform_string is None:
            if self.peptidoform_string is None:
                return
            else:
                simplified_string = self.peptidoform_string
        else:
            simplified_string = peptidoform_string
            self.peptidoform_string = peptidoform_string

        # Identify and strip labile or unlocalized components
        done = False
        while not done:
            match = re.match(r'({.+?})',simplified_string)
            if match:
                unlocalized_mod = match.group(1)
                self.unlocalized_mods.append( { 'string': unlocalized_mod, 'name': unlocalized_mod[1:-1] } )
                simplified_string = simplified_string[len(unlocalized_mod):]
            else:
                done = True


        # Identify and strip N terminal components
        match = re.match(r'n-',simplified_string)
        if match:
            simplified_string = simplified_string[2:]
        else:
            match = re.match(r'(\[.+?\])-(.+)',simplified_string)
            if match:
                self.nterm['string'] = match.group(1)
                self.nterm['name'] = self.nterm['string'][1:-1]
                simplified_string = match.group(2)
                #print(f"{self.nterm['string']}, {self.nterm['name']}, {simplified_string}")

        # Identify and strip C terminal components
        match = re.search(r'-c$',simplified_string)
        if match:
            simplified_string = simplified_string[0:len(simplified_string)-2]
        else:
            match = re.match(r'(.+)-(\[.+?\])$',simplified_string)
            if match:
                self.cterm[ 'string' ] = match.group(2)
                simplified_string = match.group(1)

        # Parse the internal part of the peptide
        characters = list(simplified_string)
        index = 1
        brace_level = 0
        for character in characters:
            if character == '[':
                self.residues[-1]['string'] += character
                brace_level += 1
            elif character == ']':
                self.residues[-1]['string'] += character
                brace_level -= 1
                if brace_level < 0:
                    self.is_valid = False
                    return
            else:
                if brace_level > 0:
                    self.residues[-1]['string'] += character
                else:
                    self.residues.append( { 'index': index, 'string': character } )
                    self.stripped_sequence += character
                    index += 1



    ####################################################################################################
    #### Return a printable buffer string of the details of the peptidoform
    def show(self):

        buf = ''
        buf += f"Peptidoform_string={self.peptidoform_string}\n"
        for mod in self.unlocalized_mods:
            buf += f"unlocalized: {mod['string']}\n"
        buf += f"   n: {self.nterm['string']}\n"
        for residue in self.residues:
            buf += '  ' + '{:2d}'.format(residue['index']) + f": {residue['string']}\n"
        buf += f"   c: {self.cterm['string']}"
        return buf


####################################################################################################
#### For command-line usage
def main():

    argparser = argparse.ArgumentParser(description='Class representing a peptidoform')
    argparser.add_argument('--verbose', action='count', help='If set, print more information about ongoing processing' )
    argparser.add_argument('--version', action='version', version='%(prog)s 0.5')
    argparser.add_argument('peptidoform_string', type=str, nargs='*', help='Optional peptidoform strings to parse')
    params = argparser.parse_args()

    # Set verbose mode
    verbose = params.verbose
    if verbose is None:
        verbose = 1

    # If there are supplied peptidoform strings, use those, else set a default example
    peptidoform_strings_list = params.peptidoform_string
    if len(peptidoform_strings_list) == 0:
        peptidoform_strings_list = [ 'SNACK', 'n-GQEY[phospho]LLLEK-c', '[Acetyl]-GQEY[phospho]LLSWLEK-c', '{Phospho}HANGRYISHM[Oxidation]ACK',
                                     'DINNE[+27.9949]R' ]

    #### Loop over all the peptidoform strings and decompose them
    for peptidoform_string in peptidoform_strings_list:
        print('**************************************************************************')
        peptidoform = Peptidoform(peptidoform_string, verbose=verbose)
        #peptidoform.parse_peptidoform_string(peptidoform_string)
        print(peptidoform.show())

#### For command line usage
if __name__ == "__main__": main()
