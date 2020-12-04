#!/usr/bin/env python

# ***************************************************************************** #
# Python script for plotting and analysing data from the output(s) of 2DES
# Spectron files.
# ***************************************************************************** #


# LAST UPDATE 05/12/2019

# ***************************************************************************** #
# Improvements / Questions:
#
#   - NOTE: the axes in the spectron output are --> 
#       *) w3index | w1index | w3 | w1 | Im | Re 
#
#   - Check consistency between data in different folders
#   - Warning: !! I'm assuming the same subfolder structure for all the logdir(s) !!
#
#   - Save max and min for every snaps, and allow to plot them in this way
#      ---> not necessary: it is enough to remove cbmin/cbmax
#   
#   - Optimize the reading of multiple logdirs!!!!
#
#   - multiple sig2plot are allowed in the parser, but this is not properly
#     implemented in the code --> cycle over the signals once the data are properly reelaborated
#
#   - Use also complex values of the signal if needed!!! -> implement PS more properly?
#
#   - !!!!!!!!!!!!!!!!CHANGE OutDir setting if multiple signals are required!!!!!!!!!!!!!!!!!!!!!!!!!!!
#      UPDATE around line 400 ..  append may be solution
#
#   - MAGIC CUT 
#       * sometimes the x axis is W3 sometimes is W1.... How to make it general??!!
#
#   - CONVOLUTION --> remove first points and last points.... they are not meaningfull!
#
#   - WHEN SETTING A WINDOW, IMPOSE THAT THE STARTING AND ENDING POINTS ARE in the W1axis 
#     and W3axis vectors! This make it possible to plot a point in time even if you don't have
#     its exact coordinates (now you have to exactly select this point)!
#
#	- Add option to give arbitrary pulse shape in input
#
#   - FT of the data: 
#       * save cosines and sines in a vector, so you don't need to compute them all the times and frequencies
#       * implement with np. routines
#       * Remove exponential decay  --> DIFFICULT!
#       * In both FT2D and FTPP remove the data (very early and very late,
#         for which the time convolution is not good) BEFORE the average
#
# ***************************************************************************** #

# Inport Needed Python Modules
import sys,os,glob
import numpy       as np
import argparse    as arg
import os
#from scipy.optimize import curve_fit


# ***************************************************************************** #
# Functions:
# ***************************************************************************** #


# Check if file(s) exists
def checkfile(file):
    if (os.path.isfile(file) == False):
        print("\n File %s not found! I will stop.\n" % f)
        sys.exit()

# When reading a file, jump useless lines
def jump_lines(f,n_lines):
    for i in range(n_lines):
        line = f.readline()

def backline():        
    print('\r', end='') 

def unit_conversion(a,starting_units,ending_units):
    if starting_units == ending_units:
        return a
    elif (starting_units == 'nm' and ending_units == 'cmm1') or (starting_units == 'cmm1' and ending_units == 'nm'):
        return 1e7/a
    elif (starting_units == 'nm' and ending_units == 'eV') or (starting_units == 'eV' and ending_units == 'cmm1'):
        return 1239.84/a
    elif starting_units == 'cmm1' and ending_units == 'eV':
        return a/8065.54
    elif starting_units == 'eV' and ending_units == 'cmm1':
        return a*8065.54
    elif starting_units == 'fs' and ending_units == 'ps':
        return a/1000
    elif starting_units == 'ps' and ending_units == 'fs':
        return a*1000


def convolute_data_in_time(data,sigmat):

    global LW1,LW3, t2list

    minimum = 1e7
    maximum = -1e7

    auxiliary_data = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))]
    
    if True:
        for t in range(len(t2list)):
            print("... Convoluting time %d/%d" % (t2list[t],t2list[len(t2list)-1]), end = '')    
            backline()
            renorm = 0
            for t1 in range(len(t2list)):
                if abs(t2list[t]-t2list[t1]) < 2*sigmat:
                    renorm += 1
                    expo = np.exp(-(float(t2list[t]-t2list[t1])**2)/(2*sigmat*sigmat))
                    for i in range(LW1):
                        for j in range(LW3):
                            auxiliary_data[t][i][j]+=data[t1][i][j]*expo
            for i in range(LW1):
                for j in range(LW3):
                    auxiliary_data[t][i][j]/=renorm
                    if auxiliary_data[t][i][j] < minimum: minimum = auxiliary_data[t][i][j]
                    if auxiliary_data[t][i][j] > maximum: maximum = auxiliary_data[t][i][j]

    # Using np. routines makes it much faster but it is difficult to set up...
    else:
        gaussian=[0 for x in range(2*len(t2list))]
        gaussian=np.asarray(gaussian)

        for i in range(len(t2list)):
            gaussian[i]=np.exp(-(float(t2list[i]-t2list[len(t2list)-1])**2)/(2*sigmat*sigmat))
            gaussian[i+len(t2list)]=np.exp(-(float(t2list[i])**2)/(2*sigmat*sigmat))

        for i in range(LW1):
            print("... Convoluting time %d/%d" % (i,LW1), end = '')
            backline()
            for j in range(LW3):
                data_ij = []
                for t in reversed(range(len(t2list))):
                    data_ij.append(data[t][i][j])
                for t in range(len(t2list)):
                    data_ij.append(data[t][i][j])
                temp=np.convolve(np.asarray(data_ij),gaussian,'full')
                for t in range(len(t2list)):
                    auxiliary_data[t][i][j]=temp[t+2*len(t2list)]
                    if auxiliary_data[t][i][j] < minimum: minimum = auxiliary_data[t][i][j]
                    if auxiliary_data[t][i][j] > maximum: maximum = auxiliary_data[t][i][j]

    return auxiliary_data, minimum,maximum


def data_frequency_filter(data,CW1,CW3,sigmaW1,sigmaW3):

    global LW1,LW3, t2list, W1axis, W3axis
    minimum=1e7
    maximum=-1e7

    if sigmaW1 > 0:
        auxiliary_data = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))]
        if sigmaW3 > 0:
            # Do gaussian freqs. filtering along both axes
            for t in range(len(t2list)):
                for i in range(LW1):
                    for j in range(LW3):
                        auxiliary_data[t][i][j]=data[t][i][j]*np.exp(-((W1axis[i]-CW1)**2)/(2*sigmaW1**2))*np.exp(-((W3axis[j]-CW3)**2)/(2*sigmaW3**2))    
                        if auxiliary_data[t][i][j] < minimum: minimum = auxiliary_data[t][i][j]
                        elif auxiliary_data[t][i][j] > maximum: maximum = auxiliary_data[t][i][j]

            return auxiliary_data, minimum, maximum
        else:
            # filter only along W1
            for t in range(len(t2list)):
                for i in range(LW1):
                    for j in range(LW3):
                        auxiliary_data[t][i][j]=data[t][i][j]*np.exp(-((W1axis[i]-CW1)**2)/(2*sigmaW1**2)) 
                        if auxiliary_data[t][i][j] < minimum: minimum = auxiliary_data[t][i][j]
                        elif auxiliary_data[t][i][j] > maximum: maximum = auxiliary_data[t][i][j]
                    
            return auxiliary_data, minimum, maximum
    else:
        if sigmaW3 > 0:
            auxiliary_data = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))]
            # filter only along W3
            for t in range(len(t2list)):
                for i in range(LW1):
                    for j in range(LW3):
                        auxiliary_data[t][i][j]=data[t][i][j]*np.exp(-((W3axis[j]-CW3)**2)/(2*sigmaW3**2))
                        if auxiliary_data[t][i][j] < minimum: minimum = auxiliary_data[t][i][j]
                        elif auxiliary_data[t][i][j] > maximum: maximum = auxiliary_data[t][i][j]

            return auxiliary_data, minimum, maximum


def magiccut(data, t, y1, y2, x1, x2):

    global LW1,LW3, W1axis, W3axis
    Npoints = 1000
    sigma = 2*max(W1axis[1]-W1axis[0],W3axis[1]-W3axis[0])
    cut_vector = [0 for x in range(Npoints)]

    minimum = 1e7
    maximum = -1e7

    # There are two possibilities:
    #   - 1) m finite
    #   - 2) vertical line

    # 1)
    if x1 != x2:
        #  y = m*x + q
        m=(y1-y2)/(x1-x2)
        q=((y1+y2)-m*(x1+x2))/2.

        if abs(m) > 1:
            # move along y
            Dy = (y2-y1) / Npoints
            for i in range(Npoints):
                y = y1 + i*Dy
                x = (y - q)/m
                # convolute around (x,y)
                for j in range(LW1):
                    if abs(W1axis[j]-x) < 2*sigma:
                        expo1=np.exp(-((W1axis[j]-x)**2)/(2*sigma**2))
                        for k in range(LW3):
                            if abs(W3axis[k]-y) < 2*sigma:
                                cut_vector[i] += data[t][j][k]*expo1*np.exp(-((W3axis[k]-y)**2)/(2*sigma**2))
                if cut_vector[i] < minimum: minimum = cut_vector[i]
                elif cut_vector[i] > maximum: maximum = cut_vector[i]
        else:
            # move along x
            Dx = (x2-x1) / Npoints
            for i in range(Npoints):
                x = x1 + i*Dx
                y = m*x+q
                # convolute around (x,y)
                for j in range(LW1):
                    if abs(W1axis[j]-x) < 2*sigma:
                        expo1=np.exp(-((W1axis[j]-x)**2)/(2*sigma**2))
                        for k in range(LW3):
                            if abs(W3axis[k]-y) < 2*sigma:
                                cut_vector[i] += data[t][j][k]*expo1*np.exp(-((W3axis[k]-y)**2)/(2*sigma**2))
                if cut_vector[i] < minimum: minimum = cut_vector[i]
                elif cut_vector[i] > maximum: maximum = cut_vector[i]
    #2)           
    else:
        # x = x1 // move along y
        for i in range(1000):
            Dy = (y2-y1) / Npoints
            for i in range(Npoints):
                y = y1 + i*Dy
                x = x1
                # convolute around (x,y)
                for j in range(LW1):
                    if abs(W1axis[j]-x) < 2*sigma:
                        expo1=np.exp(-((W1axis[j]-x)**2)/(2*sigma**2))
                        for k in range(LW3):
                            if abs(W3axis[k]-y) < 2*sigma:
                                cut_vector[i] += data[t][j][k]*expo1*np.exp(-((W3axis[k]-y)**2)/(2*sigma**2))
                if cut_vector[i] < minimum: minimum = cut_vector[i]
                elif cut_vector[i] > maximum: maximum = cut_vector[i]

    return cut_vector, minimum, maximum

# ***************************************************************************** #
# MAIN function to read data: it detects the folders containing 2D maps at 
# different times, the maps dimension 
# ***************************************************************************** #


def read_spectron_out_old(folder,complex):

    global LW1, LW3, infilename, W1i_index, W1f_index, W3i_index, W3f_index
    data = [[0 for x in range(LW3)] for y in range(LW1)]

    if complex == "Re":
        index=5
    elif complex == "Im":
        index=4

    logfile = "%s/%s" % (folder,infilename)
    checkfile(logfile) # Check if the file exists
   
    print("... Reading folder: %s" % folder, end = '')
    backline()
      
    with open(logfile,'r') as f:
        #jump_lines(f,2)
        while True:
            line = f.readline()
            if not line: break
            elif "#" in line: continue
            else:
                # NOTE: the axes in the spectron output are
                # w3index w1index w3 w1 Im Re
                if int(line.split()[1]) < W1f_index and int(line.split()[1]) >= W1i_index:
                    if int(line.split()[0]) < W3f_index and int(line.split()[0]) >= W3i_index:
                        data[int(line.split()[1])-W1i_index][int(line.split()[0])-W3i_index] = -float(line.split()[index])

    return data 


def read_spectron_out(folder,complex):

    # We can take advantage of knowing the structure of the Spectron output file
    # to make the reading routine go faster.

    # NOTE: the axes in the spectron output are
    # w3index w1index w3 w1 Im Re

    global LW1, LW3, LW1original, LW3original, infilename, W1i_index, W1f_index, W3i_index, W3f_index
    data = [[0 for x in range(LW3)] for y in range(LW1)]

    if complex == "Re":
        index=5
    elif complex == "Im":
        index=4

    logfile = "%s/%s" % (folder,infilename)
    checkfile(logfile) # Check if the file exists
   
    print("... Reading folder: %s" % folder, end = '')
    backline()

    with open(logfile,'r') as f:
        # Two commented lines + a number of useless lines
        jump_lines(f,2+W3i_index*LW1original)
        for j in range(W3i_index,W3f_index):
            jump_lines(f,W1i_index)
            for i in range(W1i_index,W1f_index):
                line = f.readline()
                data[int(line.split()[1])-W1i_index][int(line.split()[0])-W3i_index] = -float(line.split()[index])
            jump_lines(f,LW1original-W1f_index)

    return data 

# ***************************************************************************** #
# Fill OPT dictionary with default values. If you want to change the default
# options you can modify here.
# ***************************************************************************** #

# cm^-1 to eV // cm^-1 to nm


PhyCon =   {
            'au2cmm1'   : 219474.63 ,       # wavenumber per au
            'eV2cmm1'   : 8065.5443 ,       # wavenumber per eV
            }

Units =   {
            'cmm1'   : 219474.63 ,       # wavenumber per au
            'nm'   : 8065.5443 ,       # wavenumber per eV
            }



# ***************************************************************************** #
#                            MAIN BODY OF THE INTERFACE
# ***************************************************************************** #


# ArgumentParser is used here to generate a user-friendly command-line interface.
parser = arg.ArgumentParser(description="Read Spectron 2DES files, plot and analyze them with a variety of tools.",formatter_class=arg.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v',help='Increase the verbosity of the output',action="count",default=0)
parser.add_argument('logdir',nargs='+',help='Path to the directory in which the data are stored. The directory should either have the Spectron output, \
    or a sequence of folders named t2_* containing 2DES maps at different t2 time (in fs). \
    If more than one directory is given, the corresponing data (same t2 time) are summed.')
parser.add_argument('-units',nargs='+',help='Set frequency and time units. Default units are: nm (for frequency axes) and fs (for time axis).', \
    default=['cmm1','fs'], choices=['nm','cmm1','eV','fs','ps'])
parser.add_argument('-sig2plot',nargs='+',help='Signal(s)/property(ies) that one wants to plot. choices=[\'2Dmap\',\'PP\',\'PPheatmap\', \
    \'Point2\',\'magiccut\',\'FT2Dmap\',\'FTPPheatmap\',\'SpecDiff\',\'SpecDiff2\']',default=['2Dmap'])
parser.add_argument('-w1i',help='W1 Initial frequency (in the units you choose)',type=float)
parser.add_argument('-w1f',help='W1 Final frequency (in the units you choose)',type=float)
parser.add_argument('-w3i',help='W3 Initial frequency (in the units you choose)',type=float)
parser.add_argument('-w3f',help='W3 Final frequency (in the units you choose)',type=float)
parser.add_argument('-t2i',help='Set t2 initial time (in the units you choose).',type=float,default=0)
parser.add_argument('-t2f',help='Set t2 final time (in the units you choose).',type=float,default=-1)
parser.add_argument('-sigmat2',help='Finite width pulse shape in time: sigma is the standard deviation of a gaussian for time \
    convolution of the spectra (in the units you choose).',type=float,default=0)
parser.add_argument('-CW1',help='Finite pulse width along W1: Central frequency (in the units you choose).',type=float,default=-1)
parser.add_argument('-sigmaW1',help='Finite pulse width along W1: standard deviation (in the units you choose).',type=float,default=-1)
parser.add_argument('-CW3',help='Finite pulse width along W3: Central frequency (in the units you choose).',type=float,default=-1)
parser.add_argument('-sigmaW3',help='Finite pulse width along W3: standard deviation (in the units you choose).',type=float,default=-1)

parser.add_argument('-info',help='Print to screen ONLY the info about maps edges and t2 time interval.',action='store_true')

group = parser.add_argument_group('Developer keywords')
# This can be usefull when fitting e.g. ESA dipoles without doing new computations
group.add_argument('-weights',nargs='+',help='If multiple data directories are supplied, one may weight differently the data of the different directories. \
    Number of weights should be equal to number of input directories.',default=[],type=float)
group.add_argument('-t2dname',help='Declare the name of the t2 dir(s), if it is not the default one. Then data will be read from #name_$fs folders',default='t2')
group.add_argument('-outfname',help='Declare the name of the output file containing the 2DES data, if it is not the default one.',default='sig-PP.dat')
group.add_argument('-cbmin',help='Set the cbrange minimum.',type=float)
group.add_argument('-cbmax',help='Set the cbrange maximum.',type=float)
group.add_argument('-w2',nargs='+',help='For ft2dmap at w2 frequency (in cmm1!!)',default=[])

args = parser.parse_args()

global LW1, LW3, W1axis, W3axis, t2list, W1i_index, W1f_index, W3i_index, W3f_index
global infilename, LW1original, LW3original


infilename = args.outfname
t2dirname = args.t2dname

W1axis = []
W3axis = []
t2list = []

# Read the units 
for u in args.units:
    if u == 'nm':
        freq_unit = 'nm'
        if len(args.units) == 1:
            time_unit = 'fs'
    elif u == 'cmm1':
        freq_unit = 'cmm1'
        if len(args.units) == 1:
            time_unit = 'fs'
    elif u == 'eV':
        freq_unit = 'eV'
        if len(args.units) == 1:
            time_unit = 'fs'
    elif u == 'fs':
        time_unit = 'fs'
        if len(args.units) == 1:
            freq_unit = 'cmm1'
    elif u == 'ps':
        time_unit = 'ps'
        if len(args.units) == 1:
            freq_unit = 'cmm1'
 
for i in range(len(args.logdir)):
    if args.logdir[i][-1] == "/":
        args.logdir[i]=args.logdir[i][:-1]


# ************************************************* #
# Read the data from the input dirs
# ************************************************* #

print("\n*********** READ THE DATA FILES ************\n")


# When multiple directories are given you can check the consistency 
# between the data in each of them
print("I will read the data in the following directories:")
print("(Assuming the sub-folders are labeled with t2 in fs)\n")
for dirname in args.logdir:
    print("\t*) "+dirname)
    t2min=10000
    t2max=-10
    nt2=0
    for i in range(10000):
        t2dirs = "%s/%s_%d" % (dirname,t2dirname,i)
        infile = "%s/%s" % (t2dirs,infilename)
        if os.path.isdir(t2dirs) == True:
            if os.path.isfile(infile) == True:
                t2list.append(i)
                nt2+=1
                if i < t2min: t2min=i
                if i > t2max: t2max=i
    print("\t   Number of t2 times: %d" % nt2)
    if nt2 > 0:
        print("\t   t2 minimal value: %d [%s]" % (unit_conversion(t2min,'fs',time_unit),time_unit))
        print("\t   t2 maximal value: %d [%s]" % (unit_conversion(t2max,'fs',time_unit),time_unit))
    else:
        t2min=0
        t2max=0

    print("Warning: I'm assuming the same structure is there for all folders\n")
    break
print()

# I'm assuming files of different folders have same structure
# If there are no t2_* dirs, check for the 2DES 
# file directly in the logdir
if nt2 == 0:
    infile =  "%s/%s" % (args.logdir[0],infilename)
    if (os.path.isfile(infile) == True):
        print("Only one file found - No t2 info")
    else:
        print("No output file found. I will stop.")
        exit(0)
else:
    infile =  "%s/t2_%d/%s" % (args.logdir[0],t2min,infilename)

if (os.path.isfile(infile) == True):
    print("Reading info from %s file, and providing them in the user specified units." % infile)
    print("I'm assuming the Spectron output is in cm^-1")
    with open(infile,'r') as f:
        lines = f.readlines()
        line = lines[-1] 
        
        # NOTE: the axes in the spectron output are
        # w3index w1index w3 w1 Im Re 
        LW3=int(line.split()[0])+1
        LW1=int(line.split()[1])+1
        
        W3max=float(line.split()[2])
        W1max=float(line.split()[3])
        line = lines[-2]
        #print(line)
        DW1=W1max-float(line.split()[3])
        line = lines[-LW1-2]
        #print(line)
        DW3=W3max-float(line.split()[2])

        W1min=W1max-LW1*DW1
        W3min=W3max-LW3*DW3
        
    print("\n\t   Number of freq.s along W1: %d" % LW1)
    print("\t   Number of freq.s along W3: %d" % LW3)
    if freq_unit != 'nm':
        print("\t   W1 freq. window: %.2lf -> %.2lf (DW1 = %.3lf) [%s]" % (unit_conversion(W1min,'cmm1',freq_unit),unit_conversion(W1max,'cmm1',freq_unit),unit_conversion(DW1,'cmm1',freq_unit),freq_unit))
        print("\t   W3 freq. window: %.2lf -> %.2lf (DW3 = %.3lf) [%s]" % (unit_conversion(W3min,'cmm1',freq_unit),unit_conversion(W3max,'cmm1',freq_unit),unit_conversion(DW3,'cmm1',freq_unit),freq_unit))
    else: 
        print("\t   W1 freq. window: %.2lf -> %.2lf [%s] (DW1 = %.3lf [%s]) " % (unit_conversion(W1max,'cmm1',freq_unit),unit_conversion(W1min,'cmm1',freq_unit),freq_unit,DW1,'cmm1'))
        print("\t   W3 freq. window: %.2lf -> %.2lf [%s] (DW3 = %.3lf [%s]) " % (unit_conversion(W3max,'cmm1',freq_unit),unit_conversion(W3min,'cmm1',freq_unit),freq_unit,DW3,'cmm1'))    
    print("\n")

LW1original=LW1
LW3original=LW3

# Print the data info to screen. 
if(args.info == True):
    exit(0)


# Convert all the inputs in the standard units (cmm1,fs)
if freq_unit == 'nm':
    if args.w1i != None:
        W1f=unit_conversion(args.w1i,freq_unit,'cmm1')
        if W1f > W1max:
            W1f=W1max
            print("w1i is too small: I will reset it to: %.2lf [%s]" % (unit_conversion(W1f,'cmm1',freq_unit),freq_unit))
    else:
        W1f=W1max
    if args.w1f != None:
        W1i=unit_conversion(args.w1f,freq_unit,'cmm1')
        if W1i < W1min:
            W1i=W1min
            print("w1f is too large: I will reset it to: %.2lf [%s]" % (unit_conversion(W1i,'cmm1',freq_unit),freq_unit))
    else: 
        W1i=W1min
    if args.w3i != None:
        W3f=unit_conversion(args.w3i,freq_unit,'cmm1')
        if W3f > W3max:
            W3f = W3max
            print("w3i is too small: I will reset it to: %.2lf [%s]" % (unit_conversion(W3f,'cmm1',freq_unit),freq_unit))
    else:
        if args.w1i != None:
            W3f=W1f
        else:
            W3f=W3max
    if args.w3f != None:
        W3i=unit_conversion(args.w3f,freq_unit,'cmm1')
        if W3i < W3min:
            W3i=W3min
            print("w3f is too large: I will reset it to: %.2lf [%s]" % (unit_conversion(W3i,'cmm1',freq_unit),freq_unit))
    else:
        if args.w1f != None:
            W3i=W1i
        else:
            W3i=W3min
else:
    if args.w1i != None:
        W1i=unit_conversion(args.w1i,freq_unit,'cmm1')
        if W1i < W1min:
            W1i=W1min
            print("w1i is too small: I will reset it to: %.2lf [%s]" % (unit_conversion(W1i,'cmm1',freq_unit),freq_unit))
    else:
        W1i=W1min
    if args.w1f != None:
        W1f=unit_conversion(args.w1f,freq_unit,'cmm1')
        if W1f > W1max:
            W1f=W1max
            print("w1f is too large: I will reset it to: %.2lf [%s]" % (unit_conversion(W1f,'cmm1',freq_unit),freq_unit))
    else: 
        W1f=W1max
    if args.w3i != None:
        W3i=unit_conversion(args.w3i,freq_unit,'cmm1')
        if W3i < W3min:
            W3i = W3min
            print("w3i is too small: I will reset it to: %.2lf [%s]" % (unit_conversion(W3i,'cmm1',freq_unit),freq_unit))
    else:
        if args.w1i != None:
            W3i=W1i
        else:
            W3i=W3min
    if args.w3f != None:
        W3f=unit_conversion(args.w3f,freq_unit,'cmm1')
        if W3f > W3max:
            W3f=W3max
            print("w3f is too large: I will reset it to: %.2lf [%s]" % (unit_conversion(W3f,'cmm1',freq_unit),freq_unit))
    else: 
        if args.w1f != None:
            W3f=W1f
        else:
            W3f=W3max

# I could insert args.wi > args.wf
if W1i > W1f:
    temp = W1f
    W1f = W1i
    W1i = temp
if W3i > W3f:
    temp = W3f
    W3f = W3i
    W3i = temp


if args.CW1 > 0:
    CW1 = unit_conversion(args.CW1,freq_unit,'cmm1')
    if CW1 < W1i or CW1 > W1f:
        CW1 = (W1f-W1i)/2.
        print("CW1 outside the selected window: I will reset it to: %.2lf [%s]" % (unit_conversion(CW1,'cmm1',freq_unit),freq_unit))
else:
    CW1 = (W1f-W1i)/2.
if args.CW3 > 0:
    CW3 = unit_conversion(args.CW3,freq_unit,'cmm1')
    if CW3 < W3i or CW3 > W3f:
        CW3 = (W3f-W3i)/2.
        print("CW3 outside the selected window: I will reset it to: %.2lf [%s]" % (unit_conversion(CW3,'cmm1',freq_unit),freq_unit))
else:
    CW3 = (W3f-W3i)/2.

if args.t2i != 0:
    t2i=unit_conversion(args.t2i,time_unit,'fs')
    if t2i < t2min:
       t2i = t2min
       print("t2i is too small: I will reset it to: %.2lf [%s]" % (unit_conversion(t2i,'fs',time_unit),time_unit))
else:
    t2i = t2min
if args.t2f > -1:
    t2f=unit_conversion(args.t2f,time_unit,'fs')
    if t2f > t2max:
        t2f = t2max
        print("t2f is too large: I will reset it to: %.2lf [%s]" % (unit_conversion(t2f,'fs',time_unit),time_unit))
else:
    t2f = t2max
if args.sigmat2 != 0:
    sigmat2=unit_conversion(args.sigmat2,time_unit,'fs')
else:
    sigmat2=0

if args.sigmaW1 > 0:
    if freq_unit == 'nm':
        sigmaW1 = unit_conversion(unit_conversion(CW1,'cmm1','nm')-args.sigmaW1,'nm','cmm1')-CW1
    else:
        sigmaW1 = unit_conversion(args.sigmaW1,freq_unit,'cmm1')
else: sigmaW1 = args.sigmaW1
if args.sigmaW3 > 0:
    if freq_unit == 'nm':
        sigmaW3 = unit_conversion(unit_conversion(CW3,'cmm1','nm')-args.sigmaW3,'nm','cmm1')-CW3
    else:
        sigmaW3 = unit_conversion(args.sigmaW3,freq_unit,'cmm1')
else: sigmaW3 = args.sigmaW3



# Weight 
if len(args.weights) == 0:
    weights = [ 1 for x in range(len(args.logdir))]
elif len(args.weights) > len(args.logdir):
    print("Warning: too many weights! I will stop.")
    exit(0)
elif len(args.weights) > 0 and len(args.weights) < len(args.logdir):
    print("Warning: missing weights! I will stop.")
    exit(0)
else:
    tot_weight = 0
    weights = [ 0 for x in range(len(args.logdir))]
    for i in range(len(args.logdir)):
        tot_weight+=abs(args.weights[i])
    if tot_weight !=0:
        for i in range(len(args.logdir)):
            weights[i]=args.weights[i]/tot_weight
            print("Folder: %s --> Weight: %lf" % (args.logdir[i],weights[i]))
    else:
        print("Error in weights assignment.. They are all zero. I will stop.")
        exit(0)

print("\n")

# Create the axis vectors
for i in range(LW1):
    W1axis.append(W1min+i*DW1)
for i in range(LW3):
    W3axis.append(W3min+i*DW3)

# Reshape the axis according to the map limits given in input!

#print(0, LW1)
#print(W1axis)
#print()

#print(0, LW3)
#print(W3axis)
#print()

removed = 0
W1i_index = 0
W1f_index = LW1
for i in range(len(W1axis)):
    if W1axis[i-removed] < W1i:
        W1axis.pop(i-removed)
        removed+=1
        W1i_index+=1
    elif W1axis[i-removed] > W1f:
        W1axis.pop(i-removed)
        removed+=1
        W1f_index-=1

LW1 = LW1 - removed

W3i_index = 0
W3f_index = LW3

removed = 0
for i in range(len(W3axis)):
    if W3axis[i-removed] < W3i:
        W3axis.pop(i-removed)
        removed+=1
        W3i_index+=1
    elif W3axis[i-removed] > W3f:
        W3axis.pop(i-removed)
        removed+=1
        W3f_index-=1

LW3 = LW3 -removed

#print(W1i_index, W1f_index)
#print(W1axis)
#print()

#print(W3i_index, W3f_index)
#print(W3axis)
#print()


# Create output dir(s) with proper name
map2D_flag = 'NO'
map2Dsd_flag = 'NO'
PP_flag = 'NO'
PPheatmap_flag = 'NO'
Point_t2_flag = 'NO'
SpecDiff_flag = 'NO'
MagicCut_flag = 'NO'
SpecDiff2_flag = 'NO'
FT2Dmap_flag = 'NO'
FTPPheatmap_flag = 'NO'

# Modify if more signals are required.
for i in args.sig2plot:
    if i.lower() == '2dmap':
        map2D_flag = 'YES'
        Outdir_root = "results/2Dmap"
    elif i.lower() == '2dmapsd':
        map2Dsd_flag = 'YES'
        Outdir_root = "results/2DmapSD"
    elif i.lower() == 'pp':
        PP_flag = 'YES'
        Outdir_root = "results/PP" 
    elif i.lower() == 'ppheatmap':
        PPheatmap_flag = 'YES'
        Outdir_root = "results/PPheat" 
    elif i.lower() == 'point2':
        Point_t2_flag = 'YES' 
        Outdir_root = "results/Point2" 
    elif i.lower() == 'specdiff':
        SpecDiff_flag = 'YES' 
        Outdir_root = "results/SpecDiff" 
    elif i.lower() == 'specdiff2':
        SpecDiff2_flag = 'YES' 
        Outdir_root = "results/SpecDiff_method3" 
    elif i.lower() == 'magiccut':
        MagicCut_flag = 'YES' 
        Outdir_root = "results/MagicCut" 
    elif i.lower() == 'ft2dmap':
        FT2Dmap_flag = 'YES' 
        Outdir_root = "results/FT2Dmap" 
    elif i.lower() == 'ftppheatmap':
        FTPPheatmap_flag = 'YES' 
        Outdir_root = "results/FTPPheat" 


# Attach the name of the dir(s) in the folder name
for j in args.logdir:
    Outdir_root = "%s_%s" % (Outdir_root,j)

if not os.path.exists("results"):
    os.makedirs("results")

if not os.path.exists(Outdir_root):
    Outdir = Outdir_root
    os.makedirs(Outdir)
else:
    i=1
    while True:
        Outdir = "%s-%d" % (Outdir_root,i) #(args.sig2plot,i)        
        if not os.path.exists(Outdir):
            os.makedirs(Outdir)
            break
        else:
            i=i+1

OutFile = "%s/INFO.txt" % Outdir
with open(OutFile,'w') as OutF:
    OutF.write(' '.join(sys.argv))
    OutF.write("\n")
# ***************************************** #


# Pre-screen some of the time according to user input
removed=0

for t in range(len(t2list)):
    if t2list[t-removed]<t2i-2*sigmat2:
        t2list.pop(t-removed)
        removed+=1
    if t2list[t-removed]>t2f+2*sigmat2:
        t2list.pop(t-removed)
        removed+=1

nt2 = nt2-removed


# Load the data from the different files
# 2Ddata_re[time][w1][w3]
# 2Ddata_im[time][w1][w3]

#data2DRe = np.array([[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))])
data2DRe = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))]
#auxiliary_data = np.array([[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))])
auxiliary_data = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(t2list))]

time=0
for t in t2list:
    dircounter = 0
    for dirname in args.logdir:
        folder = "%s/%s_%d" % (dirname,t2dirname,t)
        #auxiliary_data[time] = read_spectron_out(folder,"Re")
        #data2DRe[time]+=auxiliary_data[time]
        auxiliary_data[time] = read_spectron_out(folder,"Re")
        for i in range(LW1):
            for j in range(LW3):
                data2DRe[time][i][j] += auxiliary_data[time][i][j]*weights[dircounter]
        dircounter+=1
    time+=1

#data2DRe = data2DRe.tolist()

maximum = -1e7
minimum = 1e7
for t in range(len(t2list)):
    for i in range(LW1):
        for j in range(LW3):
            if data2DRe[t][i][j] < minimum: minimum = data2DRe[t][i][j]
            elif data2DRe[t][i][j] > maximum: maximum = data2DRe[t][i][j]

print("\n")


print("\n********* READ THE DATA FILES - DONE **********\n")

print("\n************** RE-ELABORATE DATA **************\n")

#for i in range(LW1):
#    for j in range(LW3):
#        print("%d\t%d\t%.2lf\t%.2lf\t%.2e\t" % (i,j,W1axis[i],W3axis[j],data2DRe[0][i][j]))


# Re-elaborate the data - convolute along time
# Re-compute max and min
if sigmat2 > 0:   
    data2DRe,minimum,maximum = convolute_data_in_time(data2DRe,sigmat2)

# At this point if the user have specified a restricted time range I can shorten the data list
removed=0
print("Warning: I will NOT remove from the list of times those that are not completely under the pulse shape..!")
print("Initial t2 list: ")
print(t2list)
for t in range(len(t2list)):
    if t2list[t-removed]<t2i:#+sigmat2:
        t2list.pop(t-removed)
        data2DRe.pop(t-removed)
        removed+=1
    if t2list[t-removed]>t2f:#-sigmat2:
        t2list.pop(t-removed)
        data2DRe.pop(t-removed)
        removed+=1
        
print("Final t2 list: ")
print(t2list)

# Re-elaborate the data - frequency windowing
# Re-compute max and min
if sigmaW1>0 or sigmaW3>0:
    data2DRe,minimum,maximum = data_frequency_filter(data2DRe,CW1,CW3,sigmaW1,sigmaW3)



print("Maximum and minimum values in the (Re) map: [%lf, %lf] " % (minimum, maximum))

if args.cbmin != None:
    minimum = args.cbmin
    print("Minimum of the range set by the user: %lf " % (minimum))
if args.cbmax != None:
    maximum = args.cbmax
    print("Maximum of the range set by the user: %lf " % (maximum))

print("\n********** RE-ELABORATE DATA - DONE ***********\n")

print("\n*************** PLOT/ANALYSIS *****************\n")

# Here the various signals possibilities are explored, and the 
# plotting files printed.
if map2D_flag == "YES" or map2Dsd_flag == 'YES':

    print("... Plotting 2D maps ...")
    print("... I will use the same units specified by the user in the input ...")

    fileplotname = "plot2Dmap.gp"


    for t in range(len(t2list)):
        filedataname = "data_2Dmap_t2_%d.txt" % t2list[t]

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:
            OutF.write("#Units are %s\n" % freq_unit)
            if freq_unit == 'nm':
                for i in reversed(range(LW1)):
                    for j in reversed(range(LW3)):
                        OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(W1axis[i],'cmm1','nm'),unit_conversion(W3axis[j],'cmm1','nm'),data2DRe[t][i][j]))
                    OutF.write("\n")
            else:
                for i in range(LW1):
                    for j in range(LW3):
                        OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(W1axis[i],'cmm1',freq_unit),unit_conversion(W3axis[j],'cmm1',freq_unit),data2DRe[t][i][j]))       
                    OutF.write("\n")

    # I will add the line that connects the max at every fix detection wavelength on the top of the map
    if map2Dsd_flag == 'YES':
        print("Adding the max-line (for spectral diffusion) in the map")
       

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the 2Dmap in time \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        OutF.write("tsleep = \"0.5\"\n")


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("ncontour = %d\n" % 11)
        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
  

        # Window
        if freq_unit == 'nm':
            OutF.write("mapxmin = %lf\nmapxmax = %lf\n" % (unit_conversion(W3f,'cmm1','nm'),unit_conversion(W3i,'cmm1','nm')))
            OutF.write("diagmin = %lf\n" % max(unit_conversion(W3f,'cmm1','nm'),unit_conversion(W1f,'cmm1','nm')))
            OutF.write("diagmax = %lf\n" % min(unit_conversion(W3i,'cmm1','nm'),unit_conversion(W1i,'cmm1','nm')))
            OutF.write("mapymin = %lf\nmapymax= %lf \n\n" % (unit_conversion(W1f,'cmm1','nm'),unit_conversion(W1i,'cmm1','nm')))
        else: 
            OutF.write("mapxmin = %lf\nmapxmax= %lf\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))
            OutF.write("diagmin = %lf\n" % max(unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit)))
            OutF.write("diagmax = %lf\n" % min(unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit)))
            OutF.write("mapymin = %lf\nmapymax = %lf\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit)))
        #

        if W1i == W3i and W1f == W3f:
            OutF.write("set size square\n\n")

        OutF.write("#... parameters done ...#")
        OutF.write("############# STARTING CONTOURS EVALUATION ###########")

        for t in t2list:
            filedataname = "data_2Dmap_t2_%d.txt" % t

            OutF.write("""
#... doing countours ...#    
print \"Doing coutours...\"

#Temporary options for creating contour table
unset surface
set contour base 
set cntrparam level incremental cbmin,(cbmax-cbmin)/ncontour,cbmax
set view map
unset clabel

""")
            OutF.write("set table \"tmp.contours_%d.dat\"\n" % t)
            OutF.write("splot '%s' using 2:1:3  notitle" % filedataname)

            OutF.write("""

#... contours done ...#
""")
        OutF.write("############# CONTOURS EVALUATION DONE ###########")

        OutF.write("""
#... doing map ...#

# Go back to normal options
unset table
unset contour
set surface

print \"Doing map...\"
set pm3d map interpolate 2,2
# OutF.write("set palette defined (0 \"#00008B\", 1 \"#4169E1\",2 \"#228B22\", 3 \"yellow\", 4 \"red\" )\\n")
set palette model RGB defined (cbmin "#4169E1", cbmin/2 "#00ffff",0 "#ffffff", cbmax/2 "yellow", cbmax "red")

set cbrange [cbmin:cbmax] # DO NOT SET AUTO CB

set xtics out offset 0.,0.8 scale 0.8
set ytics out offset 0.,0.1 scale 0.8
set cbtics scale 0.2 format '%2.0tx10^%T'
""")
        if freq_unit == 'nm':
            OutF.write("set xlabel 'Emission wavelength (nm)' enhanced font 'Helvetica,26' offset 0,0\n")
            OutF.write("set ylabel 'Excitation wavelength (nm)' enhanced font 'Helvetica,26' offset -1.5,0\n")
        else:
            OutF.write("set xlabel 'Emission energy (%s)' enhanced font 'Helvetica,26' offset 0,0\n" % freq_unit)
            OutF.write("set ylabel 'Excitation energy (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n" % freq_unit)

        OutF.write("""
# Here set ranges
set xrange [mapxmin:mapxmax]
set yrange [mapymin:mapymax]


set arrow front from first diagmin,diagmin to first diagmax,diagmax nohead lt 1 lc rgb '#000000'
""")

        time = 0
        for t in t2list:
            filedataname = "data_2Dmap_t2_%d.txt" % t

            ####################
            if map2Dsd_flag == 'YES':
                P = [0 for x in range(LW3)]
                for j in range(LW3):
                    if W3axis[j]>=W3i and W3axis[j]<=W3f:
                        maximum_pp = -1e7
                        for i in range(LW1):
                            if data2DRe[time][i][j] > maximum_pp: 
                                maximum_pp=data2DRe[time][i][j]
                                P[j]=unit_conversion(W1axis[i],'cmm1',freq_unit)
                                
                    if P[j] != 0:
                        OutF.write("! echo %f %f 3 >> tmp.contours_%d.dat\n" % (unit_conversion(W3axis[j],'cmm1',freq_unit),P[j],t))
            time += 1          
            ####################


            OutF.write("if(SAVE_TO_FILE==1) outfile  = '2Dmap_t2_%d'\n" % t)
            OutF.write("if(SAVE_TO_FILE == 1) set term postscript color size 10,8 enhanced font 'Helvetica,26' ; else set term x11 size 1000,800 enhanced font 'Helvetica,26'\n")
            OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

            OutF.write("set title 'Time t2 = %.2lf (fs)'\n" % float(t)) 
            OutF.write("splot '%s' using ($2):($1):($3) notitle, \"tmp.contours_%d.dat\" w l lc rgb '#000000' notitle\n" % (filedataname,t))
            OutF.write("system(\"sleep \".tsleep)\n")
            OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n")
            
        OutF.write("\n\n#... map done ...#\n!rm tmp.contours*\n")



    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)

elif PP_flag == "YES":

    print("... Plotting PP curves ...")
    print("... I will use the same units specified by the user in the input ...")

    fileplotname = "plotPP.gp"
    PP = [[0 for x in range(LW3)] for y in range(len(t2list))]

    # At this point the PP vector is already convoluted in time and windowed in frequency!
    # No further rielaboration is needed
    minimum = 1e7
    maximum = -1e7
    for t in range(len(t2list)):
        for j in range(LW3):
            for i in range(LW1):
                PP[t][j]+=data2DRe[t][i][j]
                if PP[t][j] < minimum: minimum = PP[t][j]
                elif PP[t][j] > maximum: maximum = PP[t][j]

    for t in range(len(t2list)):
        filedataname = "data_PP_t2_%d.txt" % t2list[t]

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:
            OutF.write("#Units are %s\n" % freq_unit)
            if freq_unit == 'nm':
                for j in reversed(range(LW3)):
                    OutF.write("%lf\t%.4e\n" % (unit_conversion(W3axis[j],'cmm1','nm'),PP[t][j]))
            else:
                for j in range(LW3):
                    OutF.write("%lf\t%.4e\n" % (unit_conversion(W3axis[j],'cmm1',freq_unit),PP[t][j]))       
                

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the PP curve at time %d\n\n" % t2list[t])
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        OutF.write("tsleep = \"0.5\"\n")


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
  

        # Window
        if freq_unit == 'nm':
            OutF.write("min = %lf\nmax = %lf\n" % (unit_conversion(W3f,'cmm1','nm'),unit_conversion(W3i,'cmm1','nm')))
            OutF.write("set xtics 25\n");
            
        else: 
            OutF.write("min = %lf\nmax= %lf\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))
        #

        OutF.write("set arrow from min,0 to max,0 nohead lt 1 lc rgb 'black'\n" )

        OutF.write("set sample 5000\n\n");
        OutF.write("#... parameters done ...#\n\n")

        
        OutF.write("# Here set ranges\n")
        OutF.write("set xrange [min:max]\n")
        OutF.write("set yrange [cbmin:cbmax]\n")

        if freq_unit == 'nm':
            OutF.write("set xlabel 'Emission wavelength (nm)' enhanced font 'Helvetica,26' offset 0,0\n")
        else:
            OutF.write("set xlabel 'Emission energy (%s)' enhanced font 'Helvetica,26' offset 0,0\n" % freq_unit)
        OutF.write("set ylabel 'Intensity (a.u.)' enhanced font 'Helvetica,26' offset -1.5,0\n\n")

        for t in t2list:
            filedataname = "data_PP_t2_%d.txt" % t
             

            OutF.write("if(SAVE_TO_FILE==1) outfile  = 'PP_t2_%d'\n" % t)
            OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
            OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

            OutF.write("set title 'Time t2 = %.2lf (fs)'\n" % float(t)) 
            OutF.write("plot '%s' w lp lw 4 lc rgb 'red' title ''\n" % filedataname)
            OutF.write("system(\"sleep \".tsleep)\n")
            OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
            
        OutF.write("\n\n#... PP plot done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif PPheatmap_flag == "YES":
    print("... Plotting PP heat map ...")
    print("... I will use the same units specified by the user in the input ...")

    fileplotname = "plotPPheatmap.gp"
    PP = [[0 for x in range(LW3)] for y in range(len(t2list))]

    # At this point the PP vector is already convoluted in time and windowed in frequency!
    # No further rielaboration is needed
    minimum = 1e7
    maximum = -1e7
    for t in range(len(t2list)):
        for j in range(LW3):
            for i in range(LW1):
                PP[t][j]+=data2DRe[t][i][j]
                if PP[t][j] < minimum: minimum = PP[t][j]
                elif PP[t][j] > maximum: maximum = PP[t][j]


    filedataname = "data_PPheatmap_t2.txt" 

    # print to file the data
    OutDataFile = "%s/%s" % (Outdir,filedataname)
    with open(OutDataFile,'w') as OutF:
        OutF.write("#Units are %s %s\n" % (time_unit,freq_unit))

        for t in range(len(t2list)):
            if freq_unit == 'nm':
                for j in reversed(range(LW3)):
                    OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(t2list[t],'fs',time_unit),unit_conversion(W3axis[j],'cmm1','nm'),PP[t][j]))
            else:
                for j in range(LW3):
                    OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(t2list[t],'fs',time_unit),unit_conversion(W3axis[j],'cmm1',freq_unit),PP[t][j]))       
            OutF.write("\n")

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the PP heat map \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
    
        # Window
        if freq_unit == 'nm':
            OutF.write("min = %lf\nmax = %lf\n" % (unit_conversion(W3f,'cmm1','nm'),unit_conversion(W3i,'cmm1','nm')))
            OutF.write("set ytics 25\n");
            
        else: 
            OutF.write("min = %lf\nmax= %lf\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))

        if freq_unit == 'nm':
            OutF.write("set ylabel 'Emission wavelength (nm)' enhanced font 'Helvetica,26' offset 0,0\n")
        else:
            OutF.write("set ylabel 'Emission energy (%s)' enhanced font 'Helvetica,26' offset 0,0\n" % freq_unit)
        OutF.write("set xlabel 'Time (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % time_unit)


        OutF.write("# Here set ranges\n")
        OutF.write("set yrange [min:max]\n")
        OutF.write("set xrange [0:]\n")        
        OutF.write("set cbrange [cbmin:cbmax]\n")

        OutF.write("set palette defined (-1 '#00008B', -0.5 'blue', 0.5 'cyan', 1. 'yellow', 1.5 'orange', 2 'red', 2.5 'brown')\n\n")
        OutF.write("set pm3d map\n")

        OutF.write("if(SAVE_TO_FILE==1) outfile  = 'PPheatmap'\n")
        OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
        OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

        OutF.write("set title 'Sigma t2 = %.2lf (%s)'\n" % (unit_conversion(sigmat2,'fs',time_unit),time_unit)) 
        OutF.write("splot '%s' u 1:2:3 title ''\n" % filedataname)
        OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
                
        OutF.write("\n\n#... PP heatmap plot done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


    ### ONE MAY ADD THE PLOT OF A SPECIFIC CUT OF THE PP map ALONG TIME ###


elif Point_t2_flag == "YES":

    print("... Plotting Point in time ...")
    print("... I will use the same units specified by the user in the input ...")    

    # It uses W1i W1f W3i W3f to determine the region / point

    fileplotname = "plotPoint2.gp"
    P = [0 for x in range(len(t2list))]

    # At this point the Point is already convoluted in time and windowed in frequency!
    # No further rielaboration is needed
    minimum = 1e7
    maximum = 1e-7

    for t in range(len(t2list)):
        for j in range(LW3):
            if W3axis[j]>=W3i and W3axis[j]<=W3f:
                for i in range(LW1):
                    if W1axis[i]>=W1i and W1axis[i]<=W1f:
                       P[t]+=data2DRe[t][i][j]
        if P[t] < minimum: minimum = P[t]
        if P[t] > maximum: maximum = P[t]


    filedataname = "data_Point2.txt" 

    # print to file the data
    OutDataFile = "%s/%s" % (Outdir,filedataname)
    with open(OutDataFile,'w') as OutF:
        OutF.write("#Units are %s\n" % time_unit)

        for t in range(len(t2list)):
            OutF.write("%lf\t%.2e\n" % (unit_conversion(t2list[t],'fs',time_unit),P[t]))

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the Point in time \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
    
        OutF.write("set ylabel 'Intensity (au)' enhanced font 'Helvetica,26' offset 0,0\n")
        OutF.write("set xlabel 'Time (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % time_unit)

        OutF.write("# Here set ranges\n")
        OutF.write("set yrange [cbmin:cbmax]\n")
        OutF.write("set xrange [0:]\n")        

        OutF.write("set samples 5000\n\n")        

        OutF.write("if(SAVE_TO_FILE==1) outfile  = 'Point2'\n")
        OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
        OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

        if freq_unit == 'nm':
            OutF.write("set title 'Window = [%.1lf,%.1lf] pump , [%.1lf,%.1lf] probe (%s)'\n" % (unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),freq_unit)) 
        else:
            OutF.write("set title 'Window = [%.1lf,%.1lf] pump , [%.1lf,%.1lf] probe (%s)'\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),freq_unit))

        OutF.write("plot '%s' w l lw 4 lc rgb 'red' title ''\n" % filedataname)
        OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")

        OutF.write("\n\n#... Point along time done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif SpecDiff_flag == "YES":

    # At every detection wavelength find maximum along excitation wavelength 
    # Save the points, and fit them with a linear relation
    # The angle (arctg of the angular coefficient of the line) should change in time due to spectral diffusion
    print("... Spectral diffusion analysis ...")
    print("... I will use the same units specified by the user in the input ...")    


    fileplotname = "SpecDiff.gp"
    
    # At this point the maps are already convoluted in time and windowed in frequency!
    # No further rielaboration is needed
    jmin=LW1
    jmax=0

    for t in range(len(t2list)):
        P = [0 for x in range(LW3)]
        filedataname = "data_sd_t2_%d.txt" % t2list[t]

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:
            OutF.write("#Units are %s\n" % freq_unit)
            for j in range(LW3):
                if W3axis[j]>=W3i and W3axis[j]<=W3f:
                    maximum = -1e7
                    for i in range(LW1):
                        if data2DRe[t][i][j] > maximum: 
                            maximum=data2DRe[t][i][j]
                            P[j]=W1axis[i]
                            jstep=j

                    #print(j,istep)
                    if jstep < jmin: jmin = jstep
                    if jstep > jmax: jmax = jstep

                if P[j] != 0:
                    OutF.write("%lf\t%lf\n" % (unit_conversion(W3axis[j],'cmm1',freq_unit),unit_conversion(P[j],'cmm1',freq_unit)))

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#FIT the SpecDiff angle in time t2 \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        #OutF.write("tsleep = \"0.5\"\n")

        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("set ylabel 'Intensity (au)' enhanced font 'Helvetica,26' offset 0,0\n")
        OutF.write("set xlabel 'Time (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % time_unit)

        OutF.write("# Here set ranges\n")

        if jmin==0: jmin=1
        if jmax==LW3-1: jmax=LW3-2
        if False:
            if freq_unit == 'nm':
                #OutF.write("set xrange [%lf:%lf]\n" % (unit_conversion(W3axis[jmax+1],'cmm1',freq_unit),unit_conversion(W3axis[jmin-1],'cmm1',freq_unit)))        
                #OutF.write("set yrange [%lf:%lf]\n" % (unit_conversion(W1axis[-1],'cmm1',freq_unit),unit_conversion(W1axis[0],'cmm1',freq_unit)))
                OutF.write("set xrange [%lf:%lf]\n" % (unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit)))        
                OutF.write("set yrange [%lf:%lf]\n" % (unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit)))
            else:
                #OutF.write("set xrange [%lf:%lf]\n" % (W3axis[jmin-1],W3axis[jmax+1]) )        
                #OutF.write("set yrange [%lf:%lf]\n" % (W1axis[0],W1axis[-1]))        
                OutF.write("set xrange [%lf:%lf]\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))        
                OutF.write("set yrange [%lf:%lf]\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit)))


        OutF.write("set samples 5000\n\n")        

        OutF.write("if(SAVE_TO_FILE==1) outfile  = 'SpecDiff_angle_t2'\n")
        OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
        OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

        OutF.write("set print \"angle.dat\"\n")

        time=0
        for t in t2list:
            filedataname = "data_sd_t2_%d.txt" % t

            OutF.write("m00 = 1; q00 = 50\n")

            if time == 0:
                OutF.write("m%d = m00; q%d = q00\n" % (time,time))
            else:
                OutF.write("m%d = m%d; q%d = q%d\n" % (time,time-1,time,time-1))
            OutF.write("f%d(x)=m%d*x+q%d\n" % (time,time,time) )

            OutF.write("fit f%d(x) '%s' using 1:2 via m%d,q%d \n" % (time,filedataname,time,time))
            #OutF.write("set title 'Time t2 = %.2lf (fs)'\n" % float(t)) 
            #OutF.write("plot '%s' w lp lw 4 lc rgb 'red' title '', f%d(x) \n" % (filedataname,time))
            #OutF.write("system(\"sleep \".tsleep)\n")

            OutF.write("print %d, atan(m%d)*180/3.1415\n" % (t,time))
            time+=1

        OutF.write("plot '%s' w lp lw 4 lc rgb 'red' title ''\n" % "angle.dat")
        OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
        OutF.write("\n\n#... Point along time done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif MagicCut_flag == "YES":

    print("... Plotting 2DMap Cut in time ...")
    print("... I will use the same units specified by the user in the input ...")    

    if W3i == W3f and W1i == W1f:
        print("Error, you have selected a point of the map! I will stop.")
        exit(0)

    if args.w1i == None:
        if args.w3i == None:
            args.w1i = args.w3i = W1i
        else:             
            args.w1i = args.w3i
    else:
        if args.w3i == None:
            args.w3i = args.w1i

    if args.w1f == None:
        if args.w3f == None:
            args.w1f = args.w3f = W1f
        else:             
            args.w1f = args.w3f
    else:
        if args.w3f == None:
            args.w3f = args.w1f

    fileplotname = "plot_cuts.gp"

    # Given two points in the map [ (W1i,W3i) and (W1f,W3f) ] 
    # obtain the cut along the line that connects the two points
    maximum = -1e7
    minimum = 1e7

    for t in range(len(t2list)):

        filedataname = "cut_t2_%d.txt" % t2list[t]
        cut_vector = [0 for x in range(1000)]

        print("... Processing t2 time: %d/%d (fs)" % (t2list[t],t2list[len(t2list)-1]), end = '')
        backline()

        if args.w1i > args.w1f:
            if args.w3i > args.w3f:
                cut_vector, temp_min, temp_max=magiccut(data2DRe,t,W1f,W1i,W3f,W3i)
            else:
                cut_vector, temp_min, temp_max=magiccut(data2DRe,t,W1f,W1i,W3i,W3f)
        else:
            if args.w3i > args.w3f:
                cut_vector, temp_min, temp_max=magiccut(data2DRe,t,W1i,W1f,W3f,W3i)
            else:
                cut_vector, temp_min, temp_max=magiccut(data2DRe,t,W1i,W1f,W3i,W3f)

        if temp_max > maximum: maximum = temp_max
        if temp_min < minimum: minimum = temp_min

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:

            if W3i != W3f:
                if abs((W1i-W1f)/(W3i-W3f)) <= 1:
                    Dx = (W3f-W3i) / 1000
                    x1 = W3i
                else:
                    Dx = (W1f-W1i) / 1000
                    x1 = W1i
            else: 
                Dx = (W1f-W1i) / 1000
                x1 = W1i
            if freq_unit == 'nm':
                for i in reversed(range(1000)):
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x1+i*Dx,'cmm1',freq_unit),cut_vector[i]))
            else:
                for i in range(1000):
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x1+i*Dx,'cmm1',freq_unit),cut_vector[i]))

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        if freq_unit == 'nm':
            OutF.write("#CUT between (pump,probe) points: (%lf,%lf) and (%lf,%lf) in (%s) \n\n" % (unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),freq_unit))
        else:
            OutF.write("#CUT between (pump,probe) points: (%lf,%lf) and (%lf,%lf) in (%s) \n\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),freq_unit))
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        OutF.write("tsleep = \"0.5\"\n")
    

        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("set ylabel 'Intensity (au)' enhanced font 'Helvetica,26' offset 0,0\n")
        OutF.write("set xlabel 'Energy (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % freq_unit)
        OutF.write("set yrange [%lf:%lf]\n" % (minimum,maximum))
        OutF.write("#... parameters done ...#\n\n")

        for t in t2list:
            filedataname = "cut_t2_%d.txt" % t
             
            OutF.write("if(SAVE_TO_FILE==1) outfile  = 'cut_t2_%d'\n" % t)
            OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
            OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

            OutF.write("set title 'Time t2 = %.2lf (fs)'\n" % float(t)) 
            OutF.write("plot '%s' w l lw 4 lc rgb 'red' title ''\n" % filedataname)
            OutF.write("system(\"sleep \".tsleep)\n")
    
            OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
            
        OutF.write("\n\n#... cut plot done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif SpecDiff2_flag == "YES":

    # Given y = m*x+q, which goes through (x0,y0), the perpendicular line (which goes through (x0,y0)) is y = -1/m*x + (y1+1/m*x1)

    print("... Plotting Cuts for Spectral Diffusion method 3 (in time) ...")
    print("... I will use the same units specified by the user in the input ...")    

    if W3i == W3f and W1i == W1f:
        print("Error, you have selected a point of the map! I will stop.")
        exit(0)
    elif W3i != W1i and W3f != W3f:
        print("Warning, I was expecting the points to lye on the main diagonal!\nI will force them to stay there.")
        if W3i < W1i:
            W1i = W3i
        else:
            W3i = W1i
        if W3f > W1f:
            W1f = W3f
        else:
            W3f = W1f

    fileplotname = "plot_cuts_t2.gp"

    # Given two points in the map [ (W1i,W3i) and (W1f,W3f) ] 
    # obtain the cut along the line that connects the two points
    
    maximum = -1e7
    minimum = 1e7

    FWHM_diag = [0 for x in range(len(t2list))]
    FWHM_antidiag = [0 for x in range(len(t2list))]

    for t in range(len(t2list)):

        filedataname = "diagcut_t2_%d.txt" % t2list[t]
        filedataname_anti = "antidiagcut_t2_%d.txt" % t2list[t]
        
        cut_vector = [0 for x in range(1000)]
        anticut_vector = [0 for x in range(1000)]

        print("... Processing t2 time: %d/%d (fs)" % (t2list[t],t2list[len(t2list)-1]), end = '')
        backline()


        cut_vector, temp_min, temp_max=magiccut(data2DRe,t,W1i,W1f,W3i,W3f)
        if temp_max > maximum: maximum = temp_max
        if temp_min < minimum: minimum = temp_min


        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:

            Dx = (W3f-W3i) / 1000
            x1 = W3i

            if freq_unit == 'nm':
                for i in reversed(range(1000)):
                    x = x1+i*Dx
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x,'cmm1',freq_unit),cut_vector[i]))
                    # locate x_max
                    if cut_vector[i] == temp_max:
                        xmax = x
                        ymax = x
            else:
                for i in range(1000):
                    x = x1+i*Dx
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x,'cmm1',freq_unit),cut_vector[i]))
                    # locate x_max
                    if cut_vector[i] == temp_max:
                        xmax = x
                        ymax = x


        # Find FWHM for cut and anticut
        for i in range(1000-1):
            if cut_vector[i] - temp_max/2 < 0 and cut_vector[i+1] - temp_max/2 > 0:
                FWHM_diag[t] = -(x1+i*Dx)
            elif cut_vector[i] - temp_max/2 > 0 and cut_vector[i+1] - temp_max/2 < 0:
                FWHM_diag[t] = FWHM_diag[t]+(x1+i*Dx)

        # Now build the magic cut on the anti-diagonal 
        # Given y = m*x+q, which goes through (x0,y0), the perpendicular line (which goes through (x0,y0)) is y' = -1/m*x + (y1+1/m*x1)
        # Here in particular y = x --> y' = -x + (ymax+xmax) 
        # Given the boundaries of the first line ( the window is defined by the points (W3i,W3i) (W3f,W3f), 
        # the boundaries of y' are defined by: (W3i,-W3i+(ymax+xmax) and (W3f,-W3f+(ymax+xmax)))
        anticut_vector, temp_min, temp_max=magiccut(data2DRe,t,-W3i+(ymax+xmax),-W3f+(ymax+xmax),W3i,W3f)
        if temp_max > maximum: maximum = temp_max
        if temp_min < minimum: minimum = temp_min

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname_anti)
        with open(OutDataFile,'w') as OutF:

            Dx = (W3f-W3i) / 1000
            x1 = W3i

            if freq_unit == 'nm':
                for i in reversed(range(1000)):
                    x = x1+i*Dx
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x,'cmm1',freq_unit),anticut_vector[i]))
            else:
                for i in range(1000):
                    x = x1+i*Dx
                    OutF.write("%lf\t%lf\n" % (unit_conversion(x,'cmm1',freq_unit),anticut_vector[i]))

        # Find FWHM for cut and anticut
        for i in range(1000-1):
            if anticut_vector[i] - temp_max/2 < 0 and anticut_vector[i+1] - temp_max/2 > 0:
                FWHM_antidiag[t] = -(x1+i*Dx)
            elif anticut_vector[i] - temp_max/2 > 0 and anticut_vector[i+1] - temp_max/2 < 0:
                FWHM_antidiag[t] = FWHM_antidiag[t]+(x1+i*Dx)

    # compute flattening and flattening max and min
    flattening_min = 1e7
    flattening_max = -1e7
    OutDataFile = "%s/flattening.txt" % (Outdir)
    with open(OutDataFile,'w') as OutF:
        for t in range(len(t2list)):
            if FWHM_diag[t] != 0:
                flattening = (FWHM_diag[t]-FWHM_antidiag[t])/FWHM_diag[t]
                OutF.write("%lf\t%lf\n" % (unit_conversion(t2list[t],'fs',time_unit),flattening))
                if flattening > flattening_max: flattening_max = flattening
                if flattening < flattening_min: flattening_min = flattening
            else:
                print("Warning: there should be some problem as FWHM_diag[%d (fs)] is null." % t2list[t])


    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        if freq_unit == 'nm':
            OutF.write("#CUT between (pump,probe) points: (%lf,%lf) and (%lf,%lf) in (%s) \n\n" % (unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),freq_unit))
        else:
            OutF.write("#CUT between (pump,probe) points: (%lf,%lf) and (%lf,%lf) in (%s) \n\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit),freq_unit))
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        OutF.write("tsleep = \"0.5\"\n")

        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("set ylabel 'Intensity (au)' enhanced font 'Helvetica,26' offset 0,0\n")
        OutF.write("set xlabel 'Energy (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % freq_unit)
        OutF.write("set yrange [%lf:%lf]\n" % (minimum,maximum))

        OutF.write("#... parameters done ...#\n\n")

        for t in t2list:
            filedataname = "diagcut_t2_%d.txt" % t
             
            OutF.write("if(SAVE_TO_FILE==1) outfile  = 'cut_t2_%d'\n" % t)
            OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
            OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

            OutF.write("set title 'Time t2 = %.2lf (fs)'\n" % float(t)) 
            OutF.write("plot '%s' w l lw 4 lc rgb 'red' title 'diag', 'anti%s' w l lw 4 lc rgb 'blue' title 'antidiag'\n" % (filedataname,filedataname))
            OutF.write("system(\"sleep \".tsleep)\n")
            OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
            
        OutF.write("\n\n#... cut plot done ...#\n")


    #os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))


    fileplotname = "plot_flattening_t2.gp"

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#Flattening ( (FWHM_diag-FWHM_antidiag) / FWHM_diag ) in t2 time (%s) \n\n" % (time_unit))
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("set ylabel 'Intensity (au)' enhanced font 'Helvetica,26' offset 0,0\n")
        OutF.write("set xlabel 'Time (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" % time_unit)
        OutF.write("set yrange [%lf:%lf]\n" % (flattening_min,flattening_max))

        OutF.write("#... parameters done ...#\n\n")

        filedataname = "flattening.txt" 
             
        OutF.write("if(SAVE_TO_FILE==1) outfile  = 'flattening'\n" )
        OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
        OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

        OutF.write("plot '%s' w lp lw 4 lc rgb 'red' title ''\n" % (filedataname))
        OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
            
        OutF.write("\n\n#... flattening plot done ...#\n")


    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif FT2Dmap_flag == "YES":

    # Two possibilities: 
    #   * I have the real part of the data (data2DRe) -> FT them -> take module squared (or module) of dataFT : dataFTRe^2 + dataFTIm^2 
    #   * I have both real and imaginary part of the data -> I FT them -> I take the module of the result.
    # Here I have implemented the first one

    print("... Plotting FT of 2D maps ...")
    print("... The maps are at different frequencies ...")
    print("... I will use the same units specified by the user in the input ...")

    fileplotname = "plotFT2Dmap.gp"


    # Fourier transform of the map in time

    Dt = (t2list[1]-t2list[0])

    # Qui la salvo in fs^-1
    w2list = []

    # Fourier transform the PP vector in time
    Dt = (t2list[1]-t2list[0]) #in fs
    w2list = []

    if len(args.w2) == 0:
        wmin = 2*np.pi*max(400,33356.4/(t2list[len(t2list)-1]*2*np.pi)) / 33356.4
        wmax = 2*np.pi*2000 / 33356.4
        Dw = (wmax-wmin)/(len(t2list)-1) 

        for k in range(len(t2list)):
            w = (wmin + k*Dw) #* 33356.4 / (2*np.pi)
            w2list.append(w)
    else:
        #wmin = 1e7
        #wmax = -1e7
        for i in range(len(args.w2)):
        #    if args.w2[i] < wmin: wmin = args.w2[i]
        #    if args.w2[i] > wmax: wmax = args.w2[i]
            w2list.append(2*np.pi*float(args.w2[i])/33356.4)       
        #wmin = 2*np.pi*wmin / 33356.4
        #wmax = 2*np.pi*wmax / 33356.4
        #Dw = 0

 
    dataFT2DRe = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(w2list))]
    dataFT2DIm = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(w2list))]
    module2 = [[[0 for x in range(LW3)] for y in range(LW1)] for z in range(len(w2list))]


    average = [[0 for x in range(LW3)] for y in range(LW1)]

    for j in range(LW3):
        for i in range(LW1):
            for t in range(len(t2list)):
                average[i][j]+=data2DRe[t][i][j]
            for t in range(len(t2list)):
                data2DRe[t][i][j]-=average[i][j]/len(t2list)


    indext2min = 0
    indext2max = len(t2list)

    if sigmat2 > 0:
        for t in range(len(t2list)):
            if abs(t2list[t] - t2list[0]) < 1.5*sigmat2:
                indext2min+=1
            if t2list[-1] - t2list[t] < 1.5*sigmat2:
                indext2max-=1


    maximum = -1e7
    minimum = 1e7

    for i in range(LW1):
        print("... Processing frequencies: %d/%d" % (i,LW1-1), end = '')
        backline()
        for j in range(LW3):
            for k in range(len(w2list)):
                w = w2list[k]
                dataFT2DRe[k][i][j] = (data2DRe[indext2min][i][j]*np.cos(w*t2list[indext2min]) + data2DRe[indext2max-1][i][j]*np.cos(w*t2list[indext2max-1]))/2
                dataFT2DIm[k][i][j] = (data2DRe[indext2min][i][j]*np.sin(w*t2list[indext2min]) + data2DRe[indext2max-1][i][j]*np.sin(w*t2list[indext2max-1]))/2

                for t in range(indext2min+1,indext2max-1):
                    dataFT2DRe[k][i][j]+=data2DRe[t][i][j]*np.cos(w*t2list[t])        
                    dataFT2DIm[k][i][j]+=data2DRe[t][i][j]*np.sin(w*t2list[t])

                #dataFT2DRe[k][i][j]*=2.*np.pi*Dt / 33356.4
                #dataFT2DIm[k][i][j]*=2.*np.pi*Dt / 33356.4

                module2[k][i][j] = (dataFT2DRe[k][i][j]**2+dataFT2DIm[k][i][j]**2)

                if module2[k][i][j] > maximum: maximum = module2[k][i][j]
                if module2[k][i][j] < minimum: minimum = module2[k][i][j]


    # Fourier transform of the map done

    for k in range(len(w2list)):
        filedataname = "data_FT2Dmap_w2_%.0lf.txt" % (33356.4*w2list[k]/(2.*np.pi))

        # print to file the data
        OutDataFile = "%s/%s" % (Outdir,filedataname)
        with open(OutDataFile,'w') as OutF:
            OutF.write("#Units are %s\n" % freq_unit)
            if freq_unit == 'nm':
                for i in reversed(range(LW1)):
                    for j in reversed(range(LW3)):
                        OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(W1axis[i],'cmm1','nm'),unit_conversion(W3axis[j],'cmm1','nm'),module2[k][i][j]))
                    OutF.write("\n")
            else:
                for i in range(LW1):
                    for j in range(LW3):
                        OutF.write("%lf\t%lf\t%.2e\n" % (unit_conversion(W1axis[i],'cmm1',freq_unit),unit_conversion(W3axis[j],'cmm1',freq_unit),module2[k][i][j]))       
                    OutF.write("\n")


    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the FT2Dmap in freq w2 \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)
        OutF.write("tsleep = \"0.5\"\n")

        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
 
        # Window
        if freq_unit == 'nm':
            OutF.write("mapxmin = %lf\nmapxmax = %lf\n" % (unit_conversion(W3f,'cmm1','nm'),unit_conversion(W3i,'cmm1','nm')))
            OutF.write("diagmin = %lf\n" % max(unit_conversion(W3f,'cmm1','nm'),unit_conversion(W1f,'cmm1','nm')))
            OutF.write("diagmax = %lf\n" % min(unit_conversion(W3i,'cmm1','nm'),unit_conversion(W1i,'cmm1','nm')))
            OutF.write("mapymin = %lf\nmapymax= %lf \n\n" % (unit_conversion(W1f,'cmm1','nm'),unit_conversion(W1i,'cmm1','nm')))
        else: 
            OutF.write("mapxmin = %lf\nmapxmax= %lf\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))
            OutF.write("diagmin = %lf\n" % max(unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W1i,'cmm1',freq_unit)))
            OutF.write("diagmax = %lf\n" % min(unit_conversion(W3f,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit)))
            OutF.write("mapymin = %lf\nmapymax = %lf\n" % (unit_conversion(W1i,'cmm1',freq_unit),unit_conversion(W1f,'cmm1',freq_unit)))
        #

        if W1i == W3i and W1f == W3f:
            OutF.write("set size square\n\n")

        OutF.write("#... parameters done ...#")

        OutF.write("""
#... doing map ...#

print \"Doing map...\"
set pm3d map interpolate 2,2
set palette defined (-1 '#00008B', -0.5 'blue', 0.5 'cyan', 1. 'yellow', 1.5 'orange', 2 'red', 2.5 'brown')
#set palette model RGB defined (cbmin "#4169E1", cbmin/2 "#00ffff",0 "#ffffff", cbmax/2 "yellow", cbmax "red")


set xtics out offset 0.,0.8 scale 0.8
set ytics out offset 0.,0.1 scale 0.8
set cbtics scale 0.2 format '%2.0tx10^%T'
""")

        if len(args.w2) == 0:
            OutF.write("set cbrange [cbmin:cbmax] # DO NOT SET AUTO CB\n")
        else:
            OutF.write("#set cbrange [cbmin:cbmax] # DO NOT SET AUTO CB\n")    
        if freq_unit == 'nm':
            OutF.write("set xlabel 'Emission wavelength (nm)' enhanced font 'Helvetica,26' offset 0,0\n")
            OutF.write("set ylabel 'Excitation wavelength (nm)' enhanced font 'Helvetica,26' offset -1.5,0\n")
        else:
            OutF.write("set xlabel 'Emission energy (%s)' enhanced font 'Helvetica,26' offset 0,0\n" % freq_unit)
            OutF.write("set ylabel 'Excitation energy (%s)' enhanced font 'Helvetica,26' offset -1.5,0\n" % freq_unit)

        OutF.write("""
# Here set ranges
set xrange [mapxmin:mapxmax]
set yrange [mapymin:mapymax]

set arrow front from first diagmin,diagmin to first diagmax,diagmax nohead lt 1 lc rgb '#000000'
""")

        for k in range(len(w2list)):
            filedataname = "data_FT2Dmap_w2_%.0lf.txt" % (33356.4*w2list[k]/2./np.pi)

            OutF.write("if(SAVE_TO_FILE==1) outfile  = 'FT2Dmap_w2_%.lf'\n" % (33356.4*w2list[k]/2./np.pi))
            if len(args.w2) != 0:
                OutF.write("if(SAVE_TO_FILE == 1) set term postscript color size 10,8 enhanced font 'Helvetica,26' ; else set term x11 size 1000,800 enhanced font 'Helvetica,26'\n")
            else:
                OutF.write("if(SAVE_TO_FILE == 1) set term postscript color size 10,8 enhanced font 'Helvetica,26' ; else set term x11 %d size 1000,800 enhanced font 'Helvetica,26'\n" % k+1)
            OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

            OutF.write("set title 'Freq w2 = %.2lf (cmm1)'\n" % (33356.4*w2list[k]/2./np.pi)) 
            OutF.write("splot '%s' using ($2):($1):($3) notitle\n" % (filedataname))
            #if len(args.w2) != 0:
            OutF.write("system(\"sleep \".tsleep)\n")
            OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n")
            
        OutF.write("\n\n#... map done ...#\n")



    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


elif FTPPheatmap_flag == "YES":

    if len(t2list) <= 1:
        print("Error: more then 1 point in time is required... I will stop.")
        exit(0)

    print("... Plotting FT of PP heat map ...")
    print("... I will use the same units specified by the user in the input ...")

    fileplotname = "plotFTPPheatmap.gp"

    # Create the PP vector 
    # At this point the PP vector is already convoluted in time and windowed in frequency!
    # No further rielaboration is needed

    PP = [[0 for x in range(LW3)] for y in range(len(t2list))]
    average = [0 for x in range(LW3)]

    for t in range(len(t2list)):
        for j in range(LW3):
            for i in range(LW1):
                PP[t][j]+=data2DRe[t][i][j]
            average[j]+=PP[t][j]

    for t in range(len(t2list)):
        for j in range(LW3):
            PP[t][j]-=average[j]/len(t2list)
     

    if False:

        # fit with a decreasing exp f(x) = A + B*exp(t/tau)
        # Remove the exponential from the data
        #`FT the residuals (you may need a much denser grid of points along w2, because now you will have superthin curves)
        
        #def func(x, a, b, c):
        #    return a+b*np.exp(-x/c)

        #ydata = [[0 for x in range(len(t2list))]]
        #for j in range(LW3):
        #    for t in range(len(t2list)):
        #        ydata[t]=PP[t][j]
        #    popt, pcov = curve_fit(func, t2list, ydata[t])

        #    for t in range(len(t2list)):
        #        PP[t][j]-=func(t2list[t],popt[0],popt[1],popt[2])


        # Save PP data in time
        for j in range(LW3):

            filedataname = "data_PP_freqIndex_%d.txt" % j 

            # print to file the data
            OutDataFile = "%s/%s" % (Outdir,filedataname)
            with open(OutDataFile,'w') as OutF:
                for t in range(len(t2list)):
                    OutF.write("%lf\t%.4e\n" % (t2list[t],PP[t][j]))       

        # Gnuplot fitting script

        if False:
            filefitname = "fitPP_j_along_t2.gp"

            # print to file the data
            OutDataFile = "%s/%s" % (Outdir,filefitname)
            with open(OutDataFile,'w') as OutF:
                OutF.write("set print \"fitting.dat\"\n")
                OutF.write("a = 1; tau = 100\n" )
                for j in range(LW3):
                    filedataname = "data_PP_freqIndex_%d.txt" % j 
                    OutF.write("f(x)=a*exp(-x/tau)\n" )
                    OutF.write("fit f(x) '%s' using 1:2 via a,tau \n" % (filedataname))
                    OutF.write("print a,tau\n")

            os.system("cd %s; gnuplot -persist %s " % (Outdir,filefitname))

            exit(0)

    # PP vector done 

    # Fourier transform the PP vector in time
    Dt = (t2list[1]-t2list[0]) #in fs
    # Qui la salvo in cmm1
    w2list = []

    if len(args.w2) == 0:
        wmin = 2*np.pi*max(400,33356.4/(t2list[len(t2list)-1]*2*np.pi)) / 33356.4
        wmax = 2*np.pi*2000 / 33356.4
        Dw = (wmax-wmin)/(len(t2list)-1) 

        for k in range(len(t2list)):
            w = (wmin + k*Dw) * 33356.4 / (2*np.pi)
            w2list.append(w)

    else:
        print("Warning: single frequency in PP heat map along frequencIES.. Probably it is not what you wanted to do. I will stop.")
        exit(0)
        #wmin = 2*np.pi*args.w2 / 33356.4
        #wmax = 2*np.pi*args.w2 / 33356.4
        #Dw = 0
        #w2list.append(args.w2)        


    FTPPRe = [[0 for x in range(LW3)] for y in range(len(w2list))]
    FTPPIm = [[0 for x in range(LW3)] for y in range(len(w2list))]
    module2 = [[0 for x in range(LW3)] for y in range(len(w2list))]
 
    maximum = -1e7
    minimum = 1e7

    indext2min = 0
    indext2max = len(t2list)

    if sigmat2 > 0:
        for t in range(len(t2list)):
            if abs(t2list[t] - t2list[0]) < 1.5*sigmat2:
                indext2min+=1
            if t2list[-1] - t2list[t] < 1.5*sigmat2:
                indext2max-=1


    for i in range(LW3):
        print("... Processing frequencies: %d/%d" % (i,LW3-1), end = '')
        backline()
        for k in range(len(w2list)):
            w = wmin + k*Dw
            FTPPRe[k][i] = (PP[indext2min][i]*np.cos(w*t2list[indext2min]) + PP[indext2max-1][i]*np.cos(w*t2list[indext2max-1]))/2.
            FTPPIm[k][i] = (PP[indext2min][i]*np.sin(w*t2list[indext2min]) + PP[indext2max-1][i]*np.sin(w*t2list[indext2max-1]))/2.

            for t in range(indext2min+1,indext2max-1):
                FTPPRe[k][i]+= (PP[t][i]*np.cos(w*t2list[t]))
                FTPPIm[k][i]+= (PP[t][i]*np.sin(w*t2list[t]))

            #FTPPRe[k][i]*=2.*np.pi*Dt / 33356.4
            #FTPPIm[k][i]*=2.*np.pi*Dt / 33356.4

            module2[k][i] = (FTPPRe[k][i]**2+FTPPIm[k][i]**2)

            if module2[k][i] > maximum: maximum = module2[k][i]
            if module2[k][i] < minimum: minimum = module2[k][i]

    # Fourier transform of the map done


    filedataname = "data_FTPPheatmap_w2.txt" 

    # print to file the data
    OutDataFile = "%s/%s" % (Outdir,filedataname)
    with open(OutDataFile,'w') as OutF:
        OutF.write("#Units are %s (w2 freq. - x axis) %s (det. freq./wavelength - y axis)\n" % ('cmm1',freq_unit))

        for k in range(len(w2list)):
            if freq_unit == 'nm':
                for j in reversed(range(LW3)):
                    OutF.write("%lf\t%lf\t%.2e\n" % (w2list[k],unit_conversion(W3axis[j],'cmm1','nm'),module2[k][j]))
            else:
                for j in range(LW3):
                    OutF.write("%lf\t%lf\t%.2e\n" % (w2list[k],unit_conversion(W3axis[j],'cmm1',freq_unit),module2[k][j]))       
            OutF.write("\n")

    # print to file gnuplot plotting instructions
    OutFile = "%s/%s" % (Outdir,fileplotname)
    with open(OutFile,'w') as OutF:
        # General and good for all the snaps
        OutF.write("#PLOT the FTPP heat map \n\n" )
        OutF.write("SAVE_TO_FILE = %d\n#SAVE_TO_FILE = 1 if you want to save a figure\n" % 0)


        OutF.write("\n#... doing parameters ...#\n")

        OutF.write("cbmin = %lf\ncbmax = %lf\n" % (minimum,maximum))
    
        # Window
        if freq_unit == 'nm':
            OutF.write("min = %lf\nmax = %lf\n" % (unit_conversion(W3f,'cmm1','nm'),unit_conversion(W3i,'cmm1','nm')))
            OutF.write("set ytics 25\n");
            
        else: 
            OutF.write("min = %lf\nmax= %lf\n" % (unit_conversion(W3i,'cmm1',freq_unit),unit_conversion(W3f,'cmm1',freq_unit)))

        if freq_unit == 'nm':
            OutF.write("set ylabel 'Emission wavelength (nm)' enhanced font 'Helvetica,26' offset 0,0\n")
        else:
            OutF.write("set ylabel 'Emission energy (%s)' enhanced font 'Helvetica,26' offset 0,0\n" % freq_unit)
        OutF.write("set xlabel 'W2 freq (cmm1)' enhanced font 'Helvetica,26' offset -1.5,0\n\n" )


        OutF.write("# Here set ranges\n")
        OutF.write("set yrange [min:max]\n")
        OutF.write("set xrange [%lf:%lf]\n" % (w2list[0],w2list[-1]))
        
        OutF.write("set cbrange [cbmin:cbmax]\n")

        OutF.write("set palette defined (-1 '#00008B', -0.5 'blue', 0.5 'cyan', 1. 'yellow', 1.5 'orange', 2 'red', 2.5 'brown')\n\n")
        OutF.write("set pm3d map\n")

        OutF.write("if(SAVE_TO_FILE==1) outfile  = 'FTPPheatmap'\n")
        OutF.write("if(SAVE_TO_FILE == 1) set term postscript color enhanced font 'Helvetica,26' ; else set term x11 enhanced font 'Helvetica,26'\n")
        OutF.write("if(SAVE_TO_FILE==1) print \"Save to \".outfile; set output outfile\n")

        OutF.write("set title 'Sigma t2 = %.2lf (%s)'\n" % (unit_conversion(sigmat2,'fs',time_unit),time_unit)) 
        OutF.write("splot '%s' u 1:2:3 title ''\n" % filedataname)
        OutF.write("if(SAVE_TO_FILE == 1) command=\"epstopdf -o=\".outfile.\".pdf \".outfile.\".eps\";  system(command)\n\n")
                
        OutF.write("\n\n#... FTPP heatmap plot done ...#\n")

    os.system("cd %s; gnuplot -persist %s " % (Outdir,fileplotname))
    print("\n\nFolder ----> %s" % Outdir)


print("\n\n************ PLOT/ANALYSIS - DONE *********** \n")


exit(0)








#! /usr/bin/env python

import numpy as np
import sys
import os
#from scipy.interpolate import interp1d
#from scipy.optimize import bisect

#filtertime = 200.  # fs

# list of colors
#Colors = [ "blue", "green", 'red', 'magenta', 'cyan', 'black' ]

# define the calculation directories
# when no command argument is given, use the working dir
#if sys.argv[1:]:
#   DirList = sys.argv[1:]
#else:
#   DirList = [ os.getcwd() ]
   
# define graphs
#import matplotlib.pyplot as plt
#graphA = plt.subplot("211")
#graphB = plt.subplot("212")

# define function to compute the first zero of Re(Auto) which is used to shift the final spectrum
#def FirstRootReAuto( TList, ReAutoList):
     # identify a suitable interval to look for the root
#     firstNegative = 3+next(i for i in range(len(ReAutoList)) if ReAutoList[i] < 0.0)
     # interpolate with spline in the interval
#     splineinterp = interp1d( TList[0:firstNegative], ReAutoList[0:firstNegative] )
     # find zero by bisection
#     return bisect(splineinterp, 0.0, TList[firstNegative-1])

# cycle over directories
#n = 0
#for d in DirList:

   # check the presence of a file with the autocorrelation function:
   # "auto" MCTDH-style function with 1st:t,2nd:real(Auto),3rd:img(Auto)
   # "wf_prop" vMCG-style function with 2nd:t,5th:real(Auto),6th:img(Auto)
   # "expectations.dat" new vMCG-style function with 1st:t,4th:real(Auto),5th:img(Auto)

   #cat wf_prop | awk ' {print $2,$5,$6,sqrt($5**2+$6**2)}' > auto; PyPlot auto
   #if os.path.isfile(os.path.join(d,"wf_prop")):
      # vMCG file is present
   #   f = os.path.join(d,"wf_prop")
   #   columns  = (1,4,5)
   #elif os.path.isfile(os.path.join(d,"auto")):
   #   # MCTDH file is present
   #   f = os.path.join(d,"auto")
   #   columns  = (0,1,2)
   #elif os.path.isfile(os.path.join(d,"expectations.dat")):
   #   # MCTDH file is present
   #   f = os.path.join(d,"expectations.dat")
   #   columns  = (0,3,4)
   #else:
   #   print( "error reading directory "+str(d)+": autocorrelation file is missing!" )
   #   sys.exit()

   # read cross correlation 
   #t,r,i = np.loadtxt(f, unpack=True, usecols=columns)
   #if len(t)>10001:
   #   t = t[0:10001]
   #   r = r[0:10001]
   #   i = i[0:10001]
#   if len(t)<10001:
#      r = np.array( list(r) + [0.0]*(10001-len(r)) )
#      i = np.array( list(i) + [0.0]*(10001-len(i)) )
#      t = np.array(range(len(r)))*(t[1]-t[0])

   # Compute shift in energy
#   Eshift = np.pi/2.0/FirstRootReAuto( t/0.0241887, r ) 
 
   # scale auto correlation function by the norm
   #Auto = list(r + 1j*i)

   # apply frequency filter
   #for i in range(len(Auto)):
   #    Auto[i] = Auto[i]*np.exp(-t[i]/filtertime)

tminus = [ -z for z in t[-1:0:-1]]
Autominus = [ z.conjugate() for z in Auto[-1:0:-1]]
t    = list(t[:])
Auto = Auto[:]

# compute fourier transform
Spectra = np.fft.ifft(Auto[:]+Autominus[:])*t[1]/0.0241887
#Spectra *= len(Spectra)
FreqSpacing = 2.0*np.pi/(t[1]/0.0241887*(len(t+tminus)))*219474.
FreqGrid = np.array(range(len(Spectra)))*FreqSpacing 

   #graphA.plot( FreqGrid, Spectra.real, label = d, ls="-", c=Colors[n] )
#   graphB.plot( tminus+t, np.imag(np.array(Autominus+Auto)), label = d, ls="-", c=Colors[n] )
   #graphB.plot( tminus+t, abs(np.array(Autominus+Auto)), label = d, ls="-", c=Colors[n] )

   #n += 1

#graphA.set_yscale( "log" )
# graphA.set_xlim( (3000.,17000.) )
#graphB.set_xlim( (-10. ,1000.  ) )

#graphA.legend()
#plt.show()
