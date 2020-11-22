#!/usr/bin/env python3



    ####################################################################################################
    # This is a simplistic deprecated method to read a spectrum denovo from the termini. It sometimes
    # works, but not well enough.
    def read_de_novo_from_termini(self, spectrum, sequencing_parameters):

        # Get the number of peaks and ensure that there is an index
        n_peaks = spectrum.attributes['number of peaks']
        if len(spectrum.peak_index) == 0:
            self.index_peaks(spectrum)

        #### Extract some sequence parameters
        precursor_mz = sequencing_parameters['precursor_mz']
        precursor_charge = sequencing_parameters['precursor_charge']
        tolerance = sequencing_parameters['tolerance']
        fragmentation_type = sequencing_parameters['fragmentation_type']

        #XXXXXXXXXX Crufty stuff to remove
        #spectrum['proc_interpretations'] = {}
        verbose = 1

        #### Set up a structure to hold intermediate sequencing results for later stitching together
        partial_sequencing_index = {}

        # Create a list for all possible sequences to emerge from here
        sequences_list = []
        finished_sequences_list = []

        # Get handles for some needed reference masses
        masses = self.mass_reference.atomic_masses
        residue_masses = self.mass_reference.aa_masses
        ion_series_attr = self.mass_reference.ion_series_attributes
        neutral_losses = self.mass_reference.neutral_losses
        neutral_losses_by_residue = self.mass_reference.neutral_losses_by_residue
        neutral_losses_by_formula = self.mass_reference.neutral_losses_by_formula
        #terminal_modifications = self.mass_reference.terminal_modifications

        # Compute the full precursor molecular weight
        precursor_mw = ( precursor_mz - masses['proton'] ) * precursor_charge

        ion_series_list = [ 'b', 'y']
        terminal_modifications = {
            'nterm': { 'n-': 0.0, '[Acetyl]-': 42.010565 },
            #'nterm': { 'n-': 0.0, '[Acetyl]-': 42.010565, '[Carbamidomethyl]-': 57.021464 },   # [Carbamidomethyl]- is the same as having a G!
            #'nterm': { 'n-': 0.0, '[TMT6]-': 229.162932, '[iTRAQ4]-': 144.102063 },
            'cterm': { '-c': 0.0 },
            #'cterm': { '-c': 0.0, '[AlphaAmidation]':-58.005479, '[Methylthio]-': 45.987721 },
        }

        # Set up a trivial terminus reverser
        terminal_reversal = { 'nterm': 'cterm', 'cterm': 'nterm' }

        # Set up the list of sequences to start with
        for ion_series in ion_series_list:
            terminus_type = ion_series_attr[ion_series]['terminus_type']
            for terminal_modification in terminal_modifications[terminus_type]:
                sequence = {
                    'terminus_type': terminus_type,
                    'working_position': 0,
                    'working_mass': ion_series_attr[ion_series]['mass'] + terminal_modifications[terminus_type][terminal_modification] + masses['proton'],
                    'complement_mass': precursor_mw - ion_series_attr[ion_series]['mass'] - ion_series_attr[ion_series_attr[ion_series]['complement']]['mass'],
                    'residues': [ terminal_modification ],
                    'preliminary_score': 0.0,
                    'status': 'GO',
                }
                sequences_list.append(sequence)

        # Main loop to try to find the next possible residues
        done = False
        while not done:
            new_sequences_list = []
            for sequence in sequences_list:

                found_match = False
                root = ''.join(sequence['residues'])
                if verbose: print(f"Looking to extend {root} with single residues")

                # If we're starting from the cterm and we might be nearing the nterm, look for completion with nterm residues
                aa_masses = residue_masses
                if sequence['terminus_type'] == 'cterm' and sequence['working_position'] != 0:
                     aa_masses = self.mass_reference.nterm_aa_masses

                #### For each of the potential terminal modifications on the opposite end that we're starting from
                for terminal_modification_label,terminal_modification_mass in terminal_modifications[terminal_reversal[sequence['terminus_type']]].items():

                    #### And for each potential residue that might be the last one, see if we can complete the sequence
                    for residue,residue_mass in aa_masses.items():
                        # If this residue is the missing piece, declare victory
                        unexplained_mass_ppm = (sequence['complement_mass'] - residue_mass - terminal_modification_mass) / precursor_charge / precursor_mz * 1e6
                        #if verbose: print(f"For residue {residue}, we still have {unexplained_mass_ppm} ppm to go")
                        if abs(unexplained_mass_ppm) < tolerance:
                            if verbose: print("Complete!")
                            new_sequence = copy.deepcopy(sequence)
                            new_sequence['working_position'] += 1
                            new_sequence['working_mass'] += residue_mass + terminal_modification_mass
                            new_sequence['complement_mass'] -= residue_mass + terminal_modification_mass
                            new_sequence['preliminary_score'] += 1.0
                            new_sequence['residues'].append(terminal_modification_label + residue)
                            new_sequence['status'] = 'COMPLETE'
                            finished_sequences_list.append(new_sequence)
                            found_match = True
                            continue

                #### If we didn't manange to complete all the way to the other terminus, see which residues might take us another step
                if not found_match:

                    #### Set the usual residues, except if we're just starting out from the nterm, then use possible nterm residues
                    aa_masses = residue_masses
                    if sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0:
                        aa_masses = self.mass_reference.nterm_aa_masses

                    #### For each potential residue that might be the next one, see if we can take one more step
                    for residue,residue_mass in aa_masses.items():
                        lookahead_mz = sequence['working_mass'] + residue_mass
                        #if verbose: print(f"  Look for {residue} at {lookahead_mz}...")
                        matches = self.find_close_ions(spectrum,lookahead_mz,tolerance)
                        for match in matches:
                            i_match_peak = match[INT_REFERENCE_PEAK]
                            if verbose: print(f"  Found match for {root}{residue} at {match[0]} with delta_ppm {match[3]} ppm and score {match[4]}")
                            new_sequence = copy.deepcopy(sequence)
                            new_sequence['working_position'] += 1
                            new_sequence['working_mass'] = match[0]
                            new_sequence['complement_mass'] -= residue_mass
                            new_sequence['preliminary_score'] += match[4]
                            new_sequence['residues'].append(residue)
                            new_sequences_list.append(new_sequence)
                            found_match = True

                            #### Store this partial sequencing result in the partial_sequencing_index for possible later stitching
                            integer_index_mass = int(new_sequence['working_mass'])
                            if integer_index_mass not in partial_sequencing_index:
                                partial_sequencing_index[integer_index_mass] = []
                            partial_sequencing_index[integer_index_mass].append(new_sequence)

                            #### FIXME For now, we're just taking the first match. not great
                            break

                #### If we didn't find a next step to take, or always at the nterm, try looking ahead two steps
                #### There's a danger here that if a single step leads us astray, then we don't have the opportunity to recover the correct chain with a double hope if we never try it
                if not found_match or ( sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0 ):

                    if verbose: print(f"Looking to extend {root} with double residues")
                    found_match = False

                    # If we're starting from the cterm and we might be nearing the nterm, look for completion with nterm residues
                    double_aa_masses = self.mass_reference.double_aa_masses
                    if sequence['terminus_type'] == 'cterm' and sequence['working_position'] != 0:
                        double_aa_masses = self.mass_reference.nterm_double_aa_masses

                    #### For each of the potential terminal modifications on the opposite end that we're starting from
                    for terminal_modification_label,terminal_modification_mass in terminal_modifications[terminal_reversal[sequence['terminus_type']]].items():

                        #### And for each potential residue pair that might be the last one, see if we can complete the sequence
                        #### FIXME. Maybe this should be done for a triplet eventually as well
                        for residue,residue_mass in double_aa_masses.items():
                            # If this residue is the missing piece, declare victory
                            unexplained_mass_ppm = (sequence['complement_mass'] - residue_mass - terminal_modification_mass) / precursor_charge / precursor_mz * 1e6
                            #if verbose: print(f"For residue {terminal_modification_label}({residue}), we still have {unexplained_mass_ppm} ppm to go")
                            if abs(unexplained_mass_ppm) < tolerance:
                                if verbose: print("Complete!")
                                new_sequence = copy.deepcopy(sequence)
                                new_sequence['working_position'] += 2
                                new_sequence['working_mass'] += residue_mass + terminal_modification_mass
                                new_sequence['complement_mass'] -= residue_mass + terminal_modification_mass
                                new_sequence['preliminary_score'] += 1.0
                                new_sequence['residues'].append(f"{terminal_modification_label}({residue})")
                                new_sequence['status'] = 'COMPLETE'
                                finished_sequences_list.append(new_sequence)
                                found_match = True

                    #### If we didn't manange to complete all the way to the other terminus, see which residue pair might take us another step
                    if not found_match:

                        #### Set the usual residue pairs, except if we're just starting out from the nterm, then use possible nterm residues
                        double_aa_masses = self.mass_reference.double_aa_masses
                        if sequence['terminus_type'] == 'nterm' and sequence['working_position'] == 0:
                            double_aa_masses = self.mass_reference.nterm_double_aa_masses

                        #### For each potential residue that might be the next one, see if we can take one more step
                        for residue,residue_mass in double_aa_masses.items():
                            lookahead_mz = sequence['working_mass'] + residue_mass
                            #if verbose: print(f"  Look for {root}{residue} at {lookahead_mz}...")
                            matches = self.find_close_ions(spectrum,lookahead_mz,tolerance)
                            for match in matches:
                                i_match_peak = match[INT_REFERENCE_PEAK]
                                if verbose: print(f"  Found match for {root}{residue} at {match[0]} with delta_ppm {match[3]} ppm and score {match[4]}")
                                new_sequence = copy.deepcopy(sequence)
                                new_sequence['working_position'] += 2
                                new_sequence['working_mass'] = match[0]
                                new_sequence['complement_mass'] -= residue_mass
                                new_sequence['preliminary_score'] += match[4]
                                new_sequence['residues'].append(f"({residue})")
                                new_sequences_list.append(new_sequence)
                                found_match = True

                #### If we haven't found anywhere to continue, then we're just as a dead end here with nowhere farther to go
                if not found_match:
                    sequence['status'] = 'STOPPED'
                    finished_sequences_list.append(sequence)

            #### If there are sequences that need more work, copy them back to the primary sequences list
            if len(new_sequences_list) > 0:
                if verbose: print("Copying new_sequences_list to sequences_list")
                sequences_list = new_sequences_list
            #### Or else if there are no new sequences then we're done
            else:
                done = True

            # Relief valve for now
            if len(sequences_list) > 1000:
                print("WARNING: Exceeding 1000 test sequences. Giving up for now..")
                return


        #return

        # Print out the results
        for terminus_type in [ 'cterm', 'nterm' ]:
            print(f"{terminus_type} sequencing results:")
            sequences_list = []
            preliminary_scores_list = []
            for sequence in finished_sequences_list:
                if sequence['terminus_type'] == terminus_type:
                    sequence['preliminary_score'] /= ( precursor_mw / 1000 )
                    if terminus_type == 'cterm':
                        residues = sequence['residues']
                        residues.reverse()
                        complement_mass = '{:.4f}'.format(sequence['complement_mass'])
                        if sequence['status'] == 'COMPLETE':
                            sequences_list.append(''.join(residues))
                            preliminary_scores_list.append(sequence['preliminary_score'])
                        else:
                            sequences_list.append(f"n-(+{complement_mass})" + ''.join(residues))
                            preliminary_scores_list.append(sequence['preliminary_score'])
                            #### See if there's something with the complement mass
                            integer_index_mass = int(sequence['complement_mass'])
                            if integer_index_mass in partial_sequencing_index:
                                for potential_sequence_match in partial_sequencing_index[integer_index_mass]:
                                    complement_mass = '{:.4f}'.format(potential_sequence_match['complement_mass'])
                                    potential_sequence_residues = sequence['residues']
                                    potential_sequence_residues.reverse()
                                    print("    " + f"n-(+{complement_mass})" + ''.join(potential_sequence_residues))
                    else:
                        residues = sequence['residues']
                        complement_mass = '{:.4f}'.format(sequence['complement_mass'])
                        if sequence['status'] == 'COMPLETE':
                            sequences_list.append(''.join(residues) + "-c")
                            preliminary_scores_list.append(sequence['preliminary_score'])
                        else:
                            sequences_list.append(''.join(residues) + f"(+{complement_mass})-c")
                            preliminary_scores_list.append(sequence['preliminary_score'])

            #### Sort the final list by reverse score
            column_name_list = [ 'score', 'sequence' ]
            zipped_list =  list(zip(preliminary_scores_list, sequences_list))
            table = pd.DataFrame(zipped_list, columns = column_name_list)
            table.sort_values(by=['score'], ascending=False, inplace=True, ignore_index=True )

            scores_list = table['score']
            sequences_list = table['sequence']
            i_sequence = 0
            for score in scores_list:
                print("\t".join( [ '{:.3f}'.format(score),sequences_list[i_sequence] ] ))
                i_sequence += 1
                if i_sequence > 100:
                    break








