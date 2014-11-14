
# coding: utf-8

## Define functions

# In[5]:

import sys
import pickle
import copy             # !! Not used
import json

## !! DEFINE MODULE LEVEL CONSTANTS AT THE TOP.

def load_json(file_name, period=32):
    """
    Loads a .json file in the correct way, so it can still be read by cadnano2.

    Parameters
    ----------
    file_name : str
        The name and location of the *.json file
    period : int
        The period of repeating segments of the staples, in the square lattice this is by default 32.

    Returns
    -------
    data : list
        all data that is loaded from the *.json file
    vstrands : list
        vstrands data
    num_helices : int
        number of helices in the design
    num_bases : int
        number of bases in the design
    idx : dict
    polarity
    period

    Example
    -------
    data, vstrands, num_helices, num_bases, idx, polarity, per = load_json('ruler_design_12nov_1353.json')
    """
    with open(file_name) as f:
        data = json.load(f)

    vstrands = data['vstrands']     # !! vstrands also in outer script.
    num_helices = len(vstrands)
    num_bases = len(vstrands[0]['scaf'])
    idx = {} #Generate dictionary for translating helix_num to vstrand_num
    polarity = {}
    for helix_num in range(num_helices):
        idx[vstrands[helix_num]['num']] = helix_num
        polarity[helix_num] = vstrands[idx[helix_num]]['col'] + vstrands[idx[helix_num]]['row']
    return data, vstrands, num_helices, num_bases, idx, polarity, period, file_name

def save_json(file_name):
    """
    Saves the .json file in the correct way, so it can still be read by cadnano2.

    Parameters
    ----------
    file_name : str
        The name and location of the *.json file
    """
    data['vstrands'] = vstrands     # !! Accessing global variables
    with open(file_name, 'wb') as outfile:
        json.dump(data, outfile)


def comp_seq_FN(raw_sequence):
    """
    Returns the complementary sequence and makes all characters in uppercase.
    """
    uppercase = {'a':'A', 'A':'A', 'c':'C', 'C':'C', 'g':'G', 'G':'G', 't':'T', 'T':'T'}
    complement = {'a':'T', 'A':'T', 'c':'G', 'C':'G', 'g':'C', 'G':'C', 't':'A', 'T':'A'}
    antisense_seq = ''
    for letter in raw_sequence:
        if letter in uppercase:
            antisense_seq = complement[letter] + antisense_seq
    return antisense_seq

def stap_color_string_FN(stap_color_int):
    return color_dc[stap_color_int]          # !! Accessing global variables

def initVars():
    """
    Initialization of variables, gives the colors a name and initiates the null_bp

    Example
    -------
    stap_color_dc, null_bp = initVars()
    """
    null_bp = [-1, -1]       # !! Redefining global variables
    # Do this in a single statement:
    stap_color_dc  = {}
    stap_color_dc[13369344] = 'red'
    stap_color_dc[16204552] = 'red orange'
    stap_color_dc[16225054] = 'light orange'
    stap_color_dc[11184640] = 'olive'
    stap_color_dc[5749504]  = 'light green'
    stap_color_dc[29184]    = 'dark green'
    stap_color_dc[243362]   = 'cyan'
    stap_color_dc[1507550]  = 'blue'
    stap_color_dc[7536862]  = 'purple'
    stap_color_dc[12060012] = 'magenta'
    stap_color_dc[3355443]  = 'dark gray'
    stap_color_dc[8947848]  = 'light gray'
    return stap_color_dc, null_bp

def openPickledFile(f):
    """ Loads a pickled file """
    input_file = file(f, 'r')       # DONT USE file() !!
    loaded_txt = pickle.load(input_file)
    input_file.close()
    return loaded_txt

def openFile(f):
    """" Opens a txt file in standard format. """
    input_file = file(f, 'r')       # DONT USE file() !!
    loaded_txt = input_file.read()
    input_file.close()
    return loaded_txt

def give_sequences(cadnano_file):
    data, vstrands, num_vstrands, num_helix_bases, idx, polarity, per, json_f = load_json(cadnano_file)
    seq_counter_dc = {}
    for seq_length in seq_dc.keys():
            seq_counter_dc[seq_length] = 0

    #Generate arrays for stap paths and scaf paths
    scaf_path_ra = []
    stap_path_ra = []
    for vstrand_num in range(num_vstrands):
        for helix_base_num in range(num_helix_bases):
            for [parity, path_ra] in [['scaf', scaf_path_ra], ['stap', stap_path_ra]]:
                [prev_helix_num, prev_helix_base_num, next_helix_num, next_helix_base_num] = vstrands[vstrand_num][parity][helix_base_num]
                if ([prev_helix_num, prev_helix_base_num] == null_bp) and ([next_helix_num, next_helix_base_num] != null_bp):
                    sub_ra = []
                    curr_helix_num, curr_helix_base_num = vstrands[vstrand_num]['num'], helix_base_num
                    end_of_strand = False
                    while not end_of_strand:
                        sub_ra.append([curr_helix_num, curr_helix_base_num])
                        [curr_helix_num, curr_helix_base_num] = vstrands[idx[curr_helix_num]][parity][curr_helix_base_num][2:]
                        end_of_strand = [curr_helix_num, curr_helix_base_num] == null_bp
                    path_ra.append(sub_ra)

    #Determine length of maxiscaf
    maxiscaf_path_ra = []
    for sub_ra in scaf_path_ra:
        if len(sub_ra) > len(maxiscaf_path_ra):
            maxiscaf_path_ra = sub_ra
    maxiscaf_length = 0
    for [helix_num, helix_base_num] in maxiscaf_path_ra:
        maxiscaf_length += 1 + vstrands[idx[helix_num]]['skip'][helix_base_num] + vstrands[idx[helix_num]]['loop'][helix_base_num]
        if vstrands[idx[helix_num]]['scaf'][helix_base_num][2] not in [helix_num, -1]:
            if square_lattice and [helix_num, helix_base_num] not in no_loop_exception_ra:
                if (helix_base_num - helix_num%2)%8 == 7:
                    maxiscaf_length += on_lattice_xover_scaf_loop_length
                else:
                    maxiscaf_length += off_lattice_xover_scaf_loop_length

    if len(maxiscaf_seq) < maxiscaf_length:
        print "Error: ", maxiscaf_seq_filename, "only has ", len(maxiscaf_seq), "bases, whereas", maxiscaf_length, "bases are required for", json_filename, "."
        sys.exit()

    print "According to", json_f, ", the maxiscaf length should be", maxiscaf_length
    print "The maxiscaf sequence length in ", maxiscaf_seq_filename, "is", len(maxiscaf_seq)
    seq_dc[maxiscaf_length] = [maxiscaf_seq]
    seq_counter_dc[maxiscaf_length] = 0
    print

    #Assign scaf base sequences
    sub_ra = ['.' for i in range(num_helix_bases)]
    scaf_base_seq_dc = {}
    for vstrand_num in range(num_vstrands):
        scaf_base_seq_dc[vstrands[vstrand_num]['num']] = sub_ra[:]
    for sub_ra in scaf_path_ra:
        seq_length = 0
        for [helix_num, helix_base_num] in sub_ra:
            seq_length += 1 + vstrands[idx[helix_num]]['skip'][helix_base_num] + vstrands[idx[helix_num]]['loop'][helix_base_num]
            if vstrands[idx[helix_num]]['scaf'][helix_base_num][2] not in [helix_num, -1]:
                if square_lattice and [helix_num, helix_base_num] not in no_loop_exception_ra:
                    if (helix_base_num - helix_num%2)%8 == 7:
                        seq_length += on_lattice_xover_scaf_loop_length
                    else:
                        seq_length += off_lattice_xover_scaf_loop_length
        seq = seq_dc[seq_length][seq_counter_dc[seq_length]]
        seq_counter_dc[seq_length] += 1
        seq_pointer = 0
        for [helix_num, helix_base_num] in sub_ra:
            base_length = 1 + vstrands[idx[helix_num]]['skip'][helix_base_num] + vstrands[idx[helix_num]]['loop'][helix_base_num]
            scaf_base_seq_dc[helix_num][helix_base_num] = seq[seq_pointer:seq_pointer + base_length]
            seq_pointer += base_length
            if vstrands[idx[helix_num]]['scaf'][helix_base_num][2] not in [helix_num, -1]:
                if square_lattice and [helix_num, helix_base_num] not in no_loop_exception_ra:
                    if (helix_base_num - helix_num%2)%8 == 7:
                        seq_pointer += on_lattice_xover_scaf_loop_length
                    else:
                        seq_pointer += off_lattice_xover_scaf_loop_length

    #Print vstrand sequences
    for vstrand_num in range(num_vstrands):
        helix_num = vstrands[vstrand_num]['num']
        if helix_num < 10:
            helix_num_string = '0' + str(helix_num)
        else:
            helix_num_string = str(helix_num)
        print helix_num_string, ''.join(scaf_base_seq_dc[helix_num])

    #Set up stap_color_ra to help with matching caDNAno colors to staple strands
    sub_ra = [-1 for i in range(num_helix_bases)]
    stap_color_ra = [sub_ra[:] for vstrand_num in range(num_vstrands)]
    for vstrand_num in range(num_vstrands):
        for [helix_base_num, stap_color_int] in vstrands[vstrand_num]['stap_colors']:
            stap_color_ra[vstrand_num][helix_base_num] = stap_color_int

    #Generate staple strand output
    stap_output_ra = []
    for sub_ra in stap_path_ra:
        seq = ''
        for [helix_num, helix_base_num] in sub_ra:
            seq += comp_seq_FN(scaf_base_seq_dc[helix_num][helix_base_num])
        start_pointer = sub_ra[0]
        end_pointer = sub_ra[-1]
        [first_helix_num, first_helix_base_num] = start_pointer
        stap_color_int = stap_color_ra[idx[first_helix_num]][first_helix_base_num]
        stap_output_ra.append([stap_color_dc[stap_color_int], len(seq), start_pointer, end_pointer, seq])

    #Generate scaffold strand output
    scaf_output_ra = []
    for sub_ra in scaf_path_ra:
        seq = ''
        for [helix_num, helix_base_num] in sub_ra:
            seq += scaf_base_seq_dc[helix_num][helix_base_num]
        start_pointer = sub_ra[0]
        end_pointer = sub_ra[-1]
        scaf_output_ra.append([len(seq), start_pointer, end_pointer, seq])


    #Sort staple strands according to caDNAno color
    sorted_stap_output_ra = sorted(stap_output_ra, key = lambda stap_output:stap_output[0])
    return scaf_output_ra, sorted_stap_output_ra, vstrands

def print_sequences():
    print
    print
    #Print sorted staple strand sequences with annotations
    for sub_ra in sorted_stap_output_ra:
        seq = sub_ra[-1]
        note = 'stap strand\t' + str(sub_ra[0]) + '\t' + str(sub_ra[1]) + 'mer\t' + str(sub_ra[2]) + '\t' + 'start\t' + str(sub_ra[3]) + '\t' + 'end'
        seq_note = seq + '\t' + note
        print seq_note

    #Print short scaffold-parity strand sequences with annotations
    for sub_ra in scaf_output_ra:
        seq = sub_ra[-1]
        if len(seq) < 100:
            note = ' short scaf strand\t' + str(sub_ra[0]) + 'mer\t' + str(sub_ra[1]) + '\t' + 'start\t' + str(sub_ra[2]) + '\t' + 'end'
            seq_note = seq + '\t' + note
            print seq_note

    #Print additional miscellaneous information
    num_stap_bases = 0
    for sub_ra in sorted_stap_output_ra:
        seq = sub_ra[-1]
        num_stap_bases += len(seq)
    print "number of stap bases is", num_stap_bases


    num_short_scaf_bases = 0
    for sub_ra in scaf_output_ra[1:]:
        seq = sub_ra[-1]
        num_short_scaf_bases += len(seq)
    print "number of short scaf bases is", num_short_scaf_bases


## Settings and files

# In[6]:

#################################################################
#adjustable parameters
square_lattice = True #Set to False if honeycomb lattice is used
on_lattice_xover_scaf_loop_length = 2 #ssDNA loops at scaffold crossovers placed on the square lattice points?
off_lattice_xover_scaf_loop_length = 0 #ssDNA loops at scaffold crossovers placed off the square lattice points?
no_loop_exception_ra = [[17, 24], [16, 247]] #no scaffold loops allowed after the indicated base pointer positions
prefix = 'outputfile\t'    #Annotation prefix

#Read files
##Need to read in designed sequences for every length of scaffold strand in the caDNAno json file,
##otherwise an error will occur
seq_dc = openPickledFile('130126_1127_designed_seqs_05-69.txt')
stap_color_dc, null_bp = initVars() #initialize variables
maxiscaf_seq_filename = 'p7308_cadnanoversion.txt'
maxiscaf_seq = openFile(maxiscaf_seq_filename) #Import maxiscaf seq
maxiscaf_seq = comp_seq_FN(comp_seq_FN(maxiscaf_seq))
maxiscaf_seq = maxiscaf_seq[30:] + maxiscaf_seq[:30]


## Ruler + handles

# In[8]:

scaf_output_ra, sorted_stap_output_ra, vstrands = give_sequences('ruler_output.json')
eight_mers = openPickledFile('141110_1628_ortho_8mers_2303.txt')
handle_color = {'cyan' : [0, 5], 'blue' : [1, 6], 'red orange' : [2, 7], 'light gray' : [3, 8], 'magenta' : [4, 9]}

for color, seqidx in handle_color.items():
    for row in sorted_stap_output_ra:
        if row[0] == color and row[3][0] in (0, 1):
            # Add 'TT' to staple seq:
            row[-1] += 'TT' + eight_mers[seqidx[row[3][0]]]

print_sequences()


## Plate

# In[10]:

scaf_output_ra, sorted_stap_output_ra, vstrands = give_sequences('7x3_output.json')
eight_mers = openPickledFile('141110_1628_ortho_8mers_2303.txt')
handle_color = {56 : [0, 5], 88 : [1, 6], 120 : [2, 7], 152 : [3, 8], 184 : [4, 9]}

print "before adding handles"
# print_sequences()
helix_to_handle_seq_idx = {16: 0, 18: 1}
for end_num, seqidx in handle_color.items():
    for row in sorted_stap_output_ra:
        # row[3] is cadnano_start coordinate (helix, basenum), e.g. [2, 200]
        if row[3][1] == end_num and row[3][0] in (16, 18):
            row[-1] += 'TTTT' + comp_seq_FN(eight_mers[seqidx[0 if row[3][0] == 16 else 1]])

print "after adding handles"
print_sequences()
