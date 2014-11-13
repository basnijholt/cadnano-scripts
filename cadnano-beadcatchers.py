
# coding: utf-8

                # Functions:
#     - load_json(file_name, period=32)
#     - save_json(file_name)
#     - removeCrossover(helix_num, start, step, num)
#     - insertBreak(helix_num, start, step, num)
#     - insertDeletions(start, step)
#     - findStaples()
#     - colorCycle(i)
#     - resetColor()
#     - colorBased_on_helix()
#     - stapleLength()
#     - removeAllStaples()
#     - joinStaple()
#     - removeStaples(helix_num, start, step, num)
#     - insertScaffCrossover(up_helix, bot_helix, base_num)
#     - forcePath(helix_num, start, stop)
                
# In[37]:

import matplotlib.pyplot as plt
import numpy as np
import json
    
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
    
    vstrands = data['vstrands']
    num_helices = len(vstrands)
    num_bases = len(vstrands[0]['scaf'])
    idx = {}
    polarity = {}
    for helix_num in range(num_helices):
        idx[vstrands[helix_num]['num']] = helix_num
        polarity[helix_num] = vstrands[idx[helix_num]]['col'] + vstrands[idx[helix_num]]['row']
    return data, vstrands, num_helices, num_bases, idx, polarity, period

def save_json(file_name):
    """
    Saves the .json file in the correct way, so it can still be read by cadnano2.

    Parameters
    ----------
    file_name : str
        The name and location of the *.json file
    """
    data['vstrands'] = vstrands
    with open(file_name, 'wb') as outfile:
        json.dump(data, outfile)

def removeCrossover(strand, start, step, num, side='left'):
    """
    Removes a crossover of the staple strands.
    
    Parameters
    ----------
    strand : int
        the strand number of one the two strands of which you want to remove the crossover
    start : int
        the number of the basepair
    step : int
        the step interval of repeating pattern
    num : int
        number of times the pattern repeats
    side : string
        the side at which the staple makes the crossover, ___| is right, |___ is left
    """
    if (idx[strand]%2 == 0) and (side == 'left') or (idx[strand]%2 != 0) and (side == 'right'):
        next_strand = vstrands[idx[strand]]['stap'][start][2]
        # if the first strand is (even & left) or (odd & right):
        for i in range(num):
            vstrands[idx[strand]]['stap'][start+i*step][2:]= [-1,-1]
            vstrands[idx[next_strand]]['stap'][start+i*step][:2]= [-1,-1]
    else:
        next_strand = vstrands[idx[strand]]['stap'][start][0]
        for i in range(num):
            vstrands[idx[strand]]['stap'][start+i*step][:2]= [-1,-1]
            vstrands[idx[next_strand]]['stap'][start+i*step][2:]= [-1,-1]

def insertBreak(helix_num, start, step, num):
    """
    Inserts a break in the staple strand.
    
    Parameters
    ----------
    strand : int
        the strand number of where you want to insert a break
    start : int
        the number of the basepair
    step : int
        the step interval of repeating pattern
    num : int
        number of times the pattern repeats
    """
    if polarity[idx[helix_num]]%2 == 0: 
        for i in range(num):
            vstrands[idx[helix_num]]['stap'][start+i*step][:2]= [-1,-1]
            vstrands[idx[helix_num]]['stap'][start+1+i*step][2:]= [-1,-1]
    else:
        for i in range(num):
            vstrands[idx[helix_num]]['stap'][start+i*step][2:]= [-1,-1]
            vstrands[idx[helix_num]]['stap'][start+1+i*step][:2]= [-1,-1]    
            
def insertScaffBreak(helix_num, start, step, num):
    """
    Inserts a break in the scaffold strand.
    
    Parameters
    ----------
    strand : int
        the strand number of where you want to insert a break
    start : int
        the number of the basepair
    step : int
        the step interval of repeating pattern
    num : int
        number of times the pattern repeats
    """
    if polarity[idx[helix_num]]%2!=0:        
        for i in range(num):
            vstrands[idx[helix_num]]['scaf'][start+i*step][:2]= [-1,-1]
            vstrands[idx[helix_num]]['scaf'][start+1+i*step][2:]= [-1,-1]
    else:
        for i in range(num):
            vstrands[idx[helix_num]]['scaf'][start+i*step][2:]= [-1,-1]
            vstrands[idx[helix_num]]['scaf'][start+1+i*step][:2]= [-1,-1]        
        
def insertDeletions(start, step, num):
    """
    Inserts a deletion at every 'step' bases and starts at 'start'.
    
    Parameters
    ----------
    start : int
        base at which the first deletion is made
    step : int
        the distance between the deletions, in a square lattice this is 48 to compensate the undertwist.
    """
    for helix_num in idx:
        if vstrands[idx[helix_num]]['scaf'][num_bases/2] != [-1, -1, -1, -1]:
            for i in range(num):
                vstrands[idx[helix_num]]['skip'][start+i*step] = -1  

def findStaples():
    """ 
    Finds the beginning points of all staples
    
    Returns
    -------
    staples : list
        A list of tuples with the num_helix and the num_base where 
        the staple starts [(num_helix, num_base), (0, 154), (1, 14), ..., (14, 158)]
    num_staples : int
        Number of staples in the design
        
    Example
    -------
    staples, num_staples = findStaples()
        returns the list and number of staples
    """
    staples = []
    for helix_num in idx:
        for base_num in range(num_bases):
            staple = vstrands[idx[helix_num]]['stap'][base_num]
            if staple[:2] == [-1, -1] and staple[2:] != [-1, -1]:
                staples.append((helix_num, base_num))
    num_staples = len(staples)
    return staples, num_staples

def colorCycle(i):
    """ 
    This funtion returns a color, and loops back to the first one when it reaches the last one. 
    """
    colors = [13369344, 243362,  1507550, 16204552, 8947848, 12060012, 29184, 5749504, 7536862,  3355443, 11184640, 16225054]
    return colors[i%len(colors)]

def resetColor():
    """ 
    Resets all staples to the color grey 
    """
    staples, num_staples = findStaples()
    # empty all lists, so no color data remains
    for helix_num in idx:
        vstrands[idx[helix_num]]['stap_colors'] = []
    # create nested lists with color data, all in grey.
    for staple_num in range(num_staples):
        vstrands[idx[staples[staple_num][0]]]['stap_colors'].append([staples[staple_num][1], 8947848])

def colorBased_on_helix():
    """ 
    This funtions gives each staple that starts on the same helix, the same color. 
        
    Dependends
    ----------
    resetColor()
        findStaples()
    """
    resetColor()
    for helix_num in idx:
        num_staples = len(vstrands[idx[helix_num]]['stap_colors'])
        for staple_num in range(num_staples):
            vstrands[idx[helix_num]]['stap_colors'][staple_num][1] = colorCycle(helix_num)
            
def colorBased_on_length():
    """ 
    This funtions gives each staple that starts on the same helix, the same color. 
        
    Dependends
    ----------
    resetColor()
        findStaples()
    stapleLength()
        findStaples()
    """
    resetColor()
    staple_info = stapleLength()
    bins = np.array([0, 32, 40, 41, 47, 49, 1000])
    inds = np.digitize(staple_info[:,0], bins)
    i = 0
    for helix_num in idx:
        num = len(vstrands[idx[helix_num]]['stap_colors'])
        for staple_num in range(num):
            vstrands[idx[helix_num]]['stap_colors'][staple_num][1] = colorCycle(inds[i])
            i += 1
            
def stapleLength(plot=0):
    """
    Returns an array with the lengths of the staples in the structure and a histogram with the lengths of the staples.
    
    Parameters
    ----------
    plot : int
        Use stapleLength(plot=1) to plot a histogram or stapleLength() for no plot.
    
    Dependends
    ----------
    findStaples()
    
    Returns
    -------
    staple_info : array_like
        An array with the lengths of all staples [staple_length helix_num base_num]
    histogram : plot
        A plot with the staple lengths vs. the frequency
    """
    staples, num_staples = findStaples()
    staples_length = []
    staple_info = np.zeros((num_staples,6), dtype=int)
    for staple_num in range(num_staples):
        helix_num = staples[staple_num][0]
        base_num = staples[staple_num][1]
        staple = vstrands[idx[helix_num]]['stap'][base_num]  
        staple_info[staple_num, 1] = helix_num
        staple_info[staple_num, 2] = base_num
        i = 1
        while staple[2:] != [-1, -1] and i < 1000:  # the < 1000 is just for debugging
            last_helix = vstrands[idx[helix_num]]['stap'][base_num][0]
            next_helix = vstrands[idx[helix_num]]['stap'][base_num][2]
            prev_base = vstrands[idx[helix_num]]['stap'][base_num][1]
            next_base = vstrands[idx[helix_num]]['stap'][base_num][3]
            helix_num = next_helix
            base_num = next_base
            staple = vstrands[idx[next_helix]]['stap'][next_base]
            i += 1
        staple_info[staple_num,0] = i
    if plot == 1:
        plt.hist(staple_info[:,0], bins= 10)
        plt.title("Staple length")
        plt.xlabel("Length")
        plt.ylabel("Frequency")
        plt.show()
    return staple_info

def removeAllStaples():
    """ 
    This removes all Staples from the structure. 
    """
    for helix_num in idx:
        for base_num in range(num_bases):
            vstrands[idx[helix_num]]['stap'][base_num] = [-1, -1, -1, -1]

def joinStaple(helix_num, start, step, num):
    """
    Joins a break between two staples.
    
    Parameters
    ----------
    helix_num : int
        the strand number of where you want to join the staples
    start : int
        the number of the basepair
    step : int
        the step interval of repeating pattern
    num : int
        number of times the pattern repeats
    """
    for i in range(num):
        if polarity[idx[helix_num]]%2 == 0:
            vstrands[idx[helix_num]]['stap'][start+i*step][:2]= [idx[helix_num],start+1+i*step]
            vstrands[idx[helix_num]]['stap'][start+1+i*step][2:]= [idx[helix_num],start+i*step]
        else:
            vstrands[idx[helix_num]]['stap'][start+i*step][2:]= [idx[helix_num],start+1+i*step]
            vstrands[idx[helix_num]]['stap'][start+1+i*step][:2]= [idx[helix_num],start+i*step]
        
def removeStaples(helix_num, start, step, num):
    """
    Removes staples.
    
    Parameters
    ----------
    strand : int
        the strand number of where you want the staples removed
    start : int
        the number of the basepair where the staple starts. Start at the leftmost base.
    step : int
        the step interval of repeating pattern
    num : int
        number of times the pattern repeats
    """ 
    for i in range(num):
        staple = vstrands[idx[helix_num]]['stap']
        base_num = start+i*step
        if polarity[idx[helix_num]]%2!=0:
            while staple[base_num][2:] != [-1, -1]:
                staple[base_num] = [-1, -1, -1, -1]
                base_num += 1
            staple[base_num] = [-1, -1, -1, -1]
        else:
            while staple[base_num][:2] != [-1, -1]:
                staple[base_num] = [-1, -1, -1, -1]
                base_num += 1
            staple[base_num] = [-1, -1, -1, -1]

def insertScaffCrossover(up_helix, bot_helix, base_num):
    """
    Inserts a scaffold crossover between up_helix and bot_helix.
    
    Parameters
    ----------
    up_helix : int
        The number of the top helix.
    bot_helix : int
        The number of the bottom helix.
    base_num : int
        The number of the left base.
    """
    if polarity[idx[up_helix]]%2 != 0:
        vstrands[idx[up_helix]]['scaf'][base_num][:2]= [bot_helix, base_num]
        vstrands[idx[up_helix]]['scaf'][base_num+1][2:]= [bot_helix, base_num+1]
        
        vstrands[idx[bot_helix]]['scaf'][base_num][2:]= [up_helix, base_num]
        vstrands[idx[bot_helix]]['scaf'][base_num+1][:2]= [up_helix, base_num+1]
    else:
        vstrands[idx[up_helix]]['scaf'][base_num][2:]= [bot_helix, base_num]
        vstrands[idx[up_helix]]['scaf'][base_num+1][:2]= [bot_helix, base_num+1]
        
        vstrands[idx[bot_helix]]['scaf'][base_num][:2]= [up_helix, base_num]
        vstrands[idx[bot_helix]]['scaf'][base_num+1][2:]= [up_helix, base_num+1]
        
def forcePath(helix_num, start, stop):
    """
    Forces a path between staples on the same helix
    
    Parameters
    ----------
    helix_num : int
        The number of the helix
    start : int
        number of the base of the 3' (or 5') end
    stop : int
        number of the base of the 5' (or 3') end
        
    Example
    -------
    forcePath(6, 19, 210)
    >>>> forces path between staple on base 19 and 210 on helix number 6
    """
    if polarity[idx[helix_num]]%2 != 0: 
        vstrands[idx[helix_num]]['stap'][start][:2] = [helix_num, stop]
        vstrands[idx[helix_num]]['stap'][stop][2:] = [helix_num, start]
    else:
        vstrands[idx[helix_num]]['stap'][start][2:] = [helix_num, stop]
        vstrands[idx[helix_num]]['stap'][stop][:2] = [helix_num, start]


## Plate design (Incl. William's edits 4 nov)

# In[61]:

data, vstrands, num_helices, num_bases, idx, polarity, per = load_json('7x3_input.json')
##############################################

for i in [9, 11, 13]:
    removeCrossover(i, 39, per, 6, 'right')

for i in [7, 9, 11]:
    removeCrossover(i, 55, per, 6, 'right')
    
for i in [0, 2, 4, 6]:
    insertBreak(i, 39, per, 6)

for i in [14, 16, 18, 20]:
    insertBreak(i, 23, per, 6)

for i in [2, 4]:
    removeCrossover(i, 23, 1, 1, 'right')
    
for i in [9, 11, 13]:
    insertScaffBreak(i, 64, 48, 3)    

for i in [0, 5, 9, 11, 13, 14]:
    forcePath(i, 16, 207)
    
for i in [1, 3, 5, 15, 17, 19]:
    insertScaffCrossover(i, i+1, 127)

insertScaffCrossover(7, 20, 119)
insertBreak(5, 23, 1, 1)
insertBreak(7, 23, 1, 1)

# script doesn't work with forcePath between different helices.
vstrands[6]['stap'][19] = [6, 20, 6, 210]
vstrands[7]['stap'][19] = [7, 210, 7, 20]
vstrands[7]['stap'][210] = [7, 209, 7, 19]
vstrands[6]['stap'][210] = [6, 19, 6, 209]

resetColor()
insertDeletions(20, 48, 4)

##############################################
# color the staples purple that need to be made longer
for i in range(1, len(vstrands[20]['stap_colors'])):
    vstrands[20]['stap_colors'][i][1] = colorCycle(8)
for i in range(1, len(vstrands[18]['stap_colors'])):
    if i%2 == 0:
        vstrands[18]['stap_colors'][i][1] = colorCycle(8)

# color the staples for the DNA-bead attachment sites in orange
for i in range(len(vstrands[0]['stap_colors'])-1):
    vstrands[0]['stap_colors'][i][1] = colorCycle(11)
for i in range(len(vstrands[2]['stap_colors'])):
    if i%2 == 1:
        vstrands[2]['stap_colors'][i][1] = colorCycle(11)

# color the staples that need to be omitted in the first step green
for i in [2, 4, 6, 14, 16, 18, 20]:
    vstrands[idx[i]]['stap_colors'][0][1] = colorCycle(7)
for i in [2, 4, 14, 16, 18]:
    vstrands[idx[i]]['stap_colors'][-1][1] = colorCycle(7)
##############################################

save_json('7x3_output.json')


## Ruler design

# In[59]:

data, vstrands, num_helices, num_bases, idx, polarity, per = load_json('ruler_input.json')
insertBreak(0,47,32,107)
insertBreak(1,47,32,107)
insertDeletions(60, 48, 72)
resetColor()
for i in range(len(vstrands[0]['stap_colors'])):
    vstrands[0]['stap_colors'][i-1][1] = colorCycle(i%6)
    vstrands[1]['stap_colors'][i][1] = colorCycle(i%6)

save_json('ruler_output.json')

