#!/usr/bin/env python

# ***************************************************************************** #
# Python script for extracting quantities (energies, dipoles...) from the
# output of different Quantum Chemistry codes and prepare Spectron input files.
# ***************************************************************************** #


# ***************************************************************************** #
# References:
#
#   * Reference paper for this code ->
#
#   * Spectron -> DOI 10.1021/cr800268n
#   * NWChem ->
#   * OpenMolcas ->
#
#   * Frequency scaling factors ->
#
# ***************************************************************************** #

# LAST UPDATE 22/05/2020 (dd/mm/yyyy)

# ***************************************************************************** #
# Recently implemented:
#   - PP technique:
#       * When asking only for GSB -diagrams or SE do not compute and write the Sff and Sef spectral
#         densities (at least allow for this option, that will speed up the spectra comput!)
#         The ESAcheck == 1 condition is used to implement this.
#
# ***************************************************************************** #


# ***************************************************************************** #
# Improvements / Questions:
#
#       * NWCHEM
#   - DIFFERENT OUTPUT IN PRESENCE OF SYMMETRY!!!!!!
#   - Are "tasks" always and automatically written in NWChem? I use
#     this keyword to decide what computation has been done in a given file.
#       * OPEN-MOLCAS
#   - Set warning if computation with symmetry -> e.g. in gradients the output
#     is different, with just the symmetry inequivalent atoms coords printed.
#       * GENERAL
#   - Double check the "-free" implementation (consistency between different
#     softwares and files (of different QC outputs) order).
#   - trasform arrays in lists to use .pop() for removing list elements. NECESSARY??
#   - For PP technique:
#           * if many f states are present, divide the Spectron input in a number
#             of (ESA) computations.
#           * when I delete transitions outside the window, I should consider also
#             the reorganization energy, that may shift transitions in the window
#             even if from the (adiabatic) energetic point of view they are outside.
#             One may force the code to include more states employing a larger args.dep
#           * Target states should not be removed.. If these are outside the window
#             a warning should be generated. If these are dark (and a threshold for
#             brightness has been set) leave them anyway. Set an Hierarchy of conflicting instructions!
#           * DISORDER_INTRA_DIAG_GAUSS does not work in spectron.. DISORDER_INTRA_DIAG_GAUSS_F with
#             file provided works (at least in my version of the Spectron code)
#   - When building the SD -> set right sign according to the software output?
#     (i.e. gradient or force)
#   - Set scaling factor for frequency according to level of theory!
#     May be different for every functional. (and write scaling factor in VIB_INFO)
#   - in -free implementation: GS sometimes is named (root) "1" (when providing gradients)
#     sometimes is "0" (when providing transition dipoles)..... unify!
#   - Warning file that collects all the warnings (e.g. we are assuming our
#     molecule to be not linear) -could be the INFO file-
#   - gradient Root xx from min Root yy [additional degree of freedom for gradients
#     computed from different minima].
#     ES: G_S0[minS3], G_S1[minS0], G_S2[minS1], G_S3[minS3]~0, G_S4[minS1], G_S5[minS4]
#         and you want all of them to be centered at the minimum of S0 (note that
#         G_Sn[minSn] should be close to zero, but may be slightly different.)
#
#   - ...
#
# ***************************************************************************** #

# Inport Needed Python Modules
import sys,os,glob
import numpy       as np
import argparse    as arg
import os
import shutil

# ***************************************************************************** #
# Functions:
# ***************************************************************************** #

# Check if file(s) exists
def checkfile(files):
    for f in files:
        if (os.path.isfile(f) == False):
            print("\n File %s not found! I will stop.\n" % f)
            sys.exit()

# When reading a file, jump useless lines
def jump_lines(f,n_lines):
    for i in range(n_lines):
        line = f.readline()

# Read Nroots - OpenMolcas
def read_Nroots(f,line):
    # It assumes that "ciroot" has already been found in line
    # Open Molcas
    if "=" in line:
        Nroots = int(line.split()[1])
    # Molcas
    else:
        line = f.readline()
        Nroots = int(line.split()[1])
    return Nroots

# Read Natoms - OpenMolcas
def read_Natoms_atoms(f,line):
    # It assumes that "**** Cartesian Coordinates"
    # has already been found in line
    atoms = []

    jump_lines(f,3)
    
    while True:
        line = f.readline()
        if len(line.strip()) == 0: break
        Natoms=int(line.split()[0])
        atoms.append(line.split()[1].replace(str(Natoms),''))

    return Natoms,atoms


# Read energies RASSCF or MS-PT2 - OpenMolcas
def read_energies_RAS_MSPT2(f,line,Nroots,PT2,PT2check):
    ener = []

    if Nroots == 0:
        print("Error: number of Roots is zero")
        exit(0)

    # Extract state energies - MS/CASPT2
    elif (PT2check == 'YES') and ("Total MS-CASPT2 energies:" in line) and (PT2 == 'MS'):
        print (" ... Reading MS-CASPT2 energies ...")
        print("Warning: if energies are at MS-PT2 level and dipoles/gradients at RASSCF (or SS-PT2) level \n there may be a wrong  assignment of states properties (reordering issue)...")

        for i in range(Nroots):
            line = f.readline()
            ener.append(float(line.split()[6]))
            ener[i]*=PhyCon['au2cmm1']

    # Extract state energies - CASSCF
    elif (PT2check == 'NO') and ("Final state energy(ies)" in line):
        print  (" ... Reading CASSCF energies ...")
        jump_lines(f,2)
        for i in range(Nroots):
            line = f.readline()
            ener.append(float(line.split()[7]))
            ener[i]*=PhyCon['au2cmm1']

    return ener


# Read energies SS-PT2 - OpenMolcas
def read_energies_SSPT2(f,line,Nroots,PT2,PT2check,ener,total_pt2_en):

    if Nroots == 0:
        print("Error: number of Roots is zero")
        exit(0)

    if total_pt2_en == 0:
        print (" ... Reading SS-CASPT2 energies ...")
        ener = [0 for x in range(Nroots)]

    root_num=int(line.split()[-1])

    if root_num > Nroots:
        print("Warning: #PT2root > #Roots. Inconsistency detected, I will stop!")
        exit(0)
        # Weaker condition
        # print("Warning: #PT2root > #Roots. Skip roots %d" % root_num)
    else:
        #first occurrence of this SS-PT2 root
        if ener[root_num-1] == 0:
            while True:
                line = f.readline()
                if "Total energy:" in line:
                    ener[root_num-1]=float(line.split()[-1])*PhyCon['au2cmm1']
                elif ("Reference weight" in line) or (not line):
                    break
            #print(ener)
        # This SS-PT2 root has already been found
        # Check if the present value is compatible with the one already found
        else:
            while True:
                line = f.readline()
                if "Total energy:" in line:
                    if np.absolute(ener[root_num-1]+float(line.split()[-1])*PhyCon['au2cmm1']) != 0:
                        if np.absolute((ener[root_num-1]-float(line.split()[-1])*PhyCon['au2cmm1'])/(ener[root_num-1]+float(line.split()[-1])*PhyCon['au2cmm1']))>0.0005:
                            print("ERROR: inconsistency between energy (SS-PT2) values for the same root!!!")
                            exit(0)
                    break
                elif ("Total CASPT2 energies" in line) or (not line):
                    print("Something went wrong, energy not found")
                    exit(0)

    return ener

# Read dipoles - OpenMolcas
def read_dips(f,line,Nroots):
    
    # It assumes that "FOR THE RASSCF" has been found in line
    print(" ... Reading dipole matrix ...")
    
    # initialize dipo list
    if Nroots != 0:
        dipo = [[[0 for x in range(Nroots)] for y in range(Nroots)] for z in range(3)]
    else:
        print("Error: number of Roots is zero")
        exit(0)
        
    # FORMATTING #
    ndip = Nroots%4
    if ndip != 0:
        nblocks = Nroots//4 + 1 # // is floor division
    else:
        nblocks = Nroots//4
        ndip = 4
    ##############
            
    jump_lines(f,2)
    line = f.readline()
                
    for comp in range(3):
        # This if statement is required to handle computations with
        # symmetry in which null blocks of the dipole matrix may be removed
                    
        if ("PROPERTY:" in line and int(line.split()[4])==comp+1):
            i=comp
            jump_lines(f,3)
                        
            for j in range(nblocks):
                for k in range(Nroots):
                    line = f.readline()
                    if j == nblocks - 1:
                        for l in range(ndip):
                            dipo[i][k][(nblocks-1)*4+l] = float(line.split()[l+1])
                    else:
                        for l in range(4):
                            dipo[i][k][j*4+l] = float(line.split()[l+1])
                if j < nblocks-1:
                    jump_lines(f,6)
            
            jump_lines(f,2)
            line = f.readline()

    return dipo


# Read modes and frequencies - OpenMolcas
def read_modes_freqs(f,line,Natoms,atoms):
    if Natoms == 0:
        print("Error: number of Atoms is zero")
        exit(0)
    
    # It assumes that "Numerical differentiation is finished!" has been found in line
    print  (" ... Reading Normal Modes ...")
    print  (" ... assuming non-linear molecule ...")
    
    freqs = []
    
    # NB this is correct only for non-linear molecules
    Nvibmodes=3*Natoms-6
    
    modes = [[[0 for x in range(Nvibmodes)] for y in range(Natoms)] for z in range(3)]
    
    jump_lines(f,12)
    
    # FORMATTING #
    nmo = Nvibmodes%6
    if nmo != 0:
        nblocks = Nvibmodes//6 + 1 # // is floor division
    else:
        nblocks = Nvibmodes//6
        nmo = 6
    ##############

    counter=0

    for j in range(nblocks):
        # Read and save frequencies
        line = f.readline()
        if j == nblocks - 1:
            for l in range(nmo):
                counter+=1
                freqs.append(float(line.split()[l+1]))
            jump_lines(f,4)
    
            for k in range(Natoms):
                for m in range(3):
                    line = f.readline()
                    for l in range(nmo):
                        modes[m][k][l+j*6] = float(line.split()[l+2])

        else:
            for l in range(6):
                counter+=1
                if line.split()[l+1][0] != 'i':
                    freqs.append(float(line.split()[l+1]))
                else:
                    print("Warning: imaginary frequency detected. I'll ignore mode nr %d." % counter)
                    freqs.append(-float(line.split()[l+1][1:]))

            jump_lines(f,4)

            for k in range(Natoms):
                for m in range(3):
                    line = f.readline()
                    for l in range(6):
                        modes[m][k][l+j*6] = float(line.split()[l+2])
        
        if j < nblocks-1:
            jump_lines(f,4)

    # In molcas the modes are not normalized.. to normalize one
    # should multiply modes[m][k][l+j*6] by np.sqrt(red_mass[l+j*6])
    # , where red_mass has to be saved before

    # Mass weighted normal modes + normalization
    # Obtaining unitless normal mode coordinates
    for m in range(Nvibmodes):
        norm = 0
        for k in range(3):
            for l in range(Natoms):
                modes[k][l][m] *= np.sqrt(AtomMass[atoms[l]])
                norm += modes[k][l][m]**2
        for k in range(3):
            for l in range(Natoms):
                modes[k][l][m] /= np.sqrt(norm)

     # Check mode normalization
     #for m in range(Nvibmodes):
     #   tot=0
     #   for k in range(3):
     #       for l in range(Natoms):
     #           tot+=modes[k][l][m]*modes[k][l][m]
     #   print(tot)

    return Nvibmodes,modes,freqs

# Read gradients - OpenMolcas
def read_grads(f,line,Nroots,Natoms,total_grads,grads,gradPT2check,grad_pt2_list):
    
    global grad_list
    
    if (Nroots == 0) or (Natoms == 0):
        print("Error: number of Roots or number of Atoms is zero")
        exit(0)
    if total_grads == 0:
        if "Lagrangian multipliers are calculated for state no" in line:
            print  (" ... Reading CASSCF-ANALYTICAL gradients ...")
        elif "Numerical gradient, root" in line:
            if gradPT2check == 'NO': print  (" ... Reading CASSCF-NUMERICAL gradients ...")
            else: print  (" ... Reading CASPT2-NUMERICAL gradients ...")
        grads = [[[0 for x in range(Natoms)] for y in range(3)] for z in range(Nroots)]
        grad_list = [0 for x in range(Nroots)]
    
    if gradPT2check == 'NO': grad_num = int(line.split()[-1])
    else:
        grad_num = grad_pt2_list[0]
        grad_pt2_list.pop(0)

    if grad_num > Nroots:
        print("Warning: #Grads > #Roots. Skip grad number %d" % grad_num)
    else:
        #first occurrence of this gradient
        if grad_list[grad_num-1] == 0:
            grad_list[grad_num-1]=1
            # Analytical gradients
            if "Lagrangian multipliers are calculated for state no" in line:
                while True:
                    line = f.readline()
                    if "Molecular gradients" in line:
                        jump_lines(f,7)
                        for l in range(Natoms):
                            line = f.readline()
                            for k in range(3):
                                grads[grad_num-1][k][l] = float(line.split()[k+1])
                        break
                    elif not line: break
            # Numerical gradients
            elif "Numerical gradient, root":
                jump_lines(f,3)
                for l in range(Natoms):
                    line = f.readline()
                    for k in range(3):
                        grads[grad_num-1][k][l] = float(line.split()[k+1])

        # This gradient number has already been found
        # Check if the present value is compatible with the one already found
        else:
            # Analytical gradients
            if "Lagrangian multipliers are calculated for state no" in line:
                while True:
                    line = f.readline()
                    if "Molecular gradients" in line:
                        jump_lines(f,7)
                        for l in range(Natoms):
                            line = f.readline()
                            for k in range(3):
                                if np.absolute(grads[grad_num-1][k][l]+float(line.split()[k+1]))!=0:
                                    if np.absolute(2*(grads[grad_num-1][k][l]-float(line.split()[k+1]))/(grads[grad_num-1][k][l]+float(line.split()[k+1])))>0.05:
                                        print("ERROR: inconsistency between gradients for the same root!!!")
                                        exit(0)
                        break
                    elif not line: break
            # Numerical gradients
            elif "Numerical gradient, root":
                jump_lines(f,3)
                for l in range(Natoms):
                    line = f.readline()
                    for k in range(3):
                        if np.absolute(grads[grad_num-1][k][l]+float(line.split()[k+1]))>0:
                            if np.absolute(2*(grads[grad_num-1][k][l]-float(line.split()[k+1]))/(grads[grad_num-1][k][l]+float(line.split()[k+1])))>0.05:
                                print("ERROR: inconsistency between gradients for the same state!!!")
                                exit(0)


    return grads


# ***************************************************************************** #
# MAIN function to read data: it detects what kind of QC software is used;
# Energies should be transform in cm^-1 and dipole moments in eA;
# Frequencies are in [cm^-1], gradients in [Hartree/Bohr] and normal
# modes coordinates [Unitless].
# The Nvibmodes=3*Natoms-6 or 3*Natoms-5 (in linear molecules)
# The dipoles moments are organized in a matrix d[i][j][k]
#       - i runs over the three dipole components 0,1,2 (for x,y and z)
#       - j and k run over all possibile state indices;
#           * if j==k one is considering the permanent dipole of state j(=k)
#           * if j!=k one is considering trans. dipoles from st. j to st. k
# The modes are organized in a matrix M[i][j][k]
#       - i runs over the three cartesian components 0,1,2 (for x,y and z)
#       - j run over the atom number
#       - k run over the mode number
# The gradients are organized in a matrix G[i][j][k]
#       - i run over the state number
#       - j runs over the three cartesian components 0,1,2 (for x,y and z)
#       - k run over the atom number
# ***************************************************************************** #

def readlog(logfile,free):

    global Natoms, Nroots
    Natoms=0
    Nroots=0

    # Check if the file exists
    checkfile(logfile)

    Ener = []; Dipoles = [[[]]] ; Gradients = [[[]]] ; Modes = [[[]]]; Frequencies = []
    
    # Check only the first file to find which QC Software has been used
    with open(logfile[0],'r') as f:
        while True:
            line = f.readline()
            if not line: break
            
            elif "MOLCAS" in line and free==0:
                OPT['read'] = 'molcas'
                OPT['logfile'] = logfile
                print  (" ... Processing MOLCAS Input ...")
                print  (" ... Reading in %s log file(s) ..." % OPT['logfile'])
                Ener,Dipoles,Gradients,Modes,Frequencies = readmollog(OPT['logfile'])
                break
            elif "Entering Gaussian System, Link 0=g16" in line or "Entering Gaussian System, Link 0=g09" in line and free==0:
                OPT['read'] = 'gaussian'
                OPT['logfile'] = logfile
                print  (" ... Processing GAUSSIAN (v.09/v.16) Input(s) ...")
                print  (" ... Reading in %s log file(s) ..." % OPT['logfile'])
                Ener,Dipoles,Gradients,Modes,Frequencies = readgaulog(OPT['logfile'])
                break
            elif "nwchem" in line and free==0:
                OPT['read'] = 'NWChem'
                OPT['logfile'] = logfile
                print  (" ... Processing NWCHEM Input(s) ...")
                print  (" ... Reading in %s log file(s) ..." % OPT['logfile'])
                Ener,Dipoles,Gradients,Modes,Frequencies = readnwclog(OPT['logfile'])
                break
            elif free!=0:
                OPT['read'] = 'FreeReading'
                OPT['logfile'] = logfile
                print  (" ... Processing More General Input(s) ...")
                print  (" ... Reading in %s log file(s) ..." % OPT['logfile'])
                Ener,Dipoles,Gradients,Modes,Frequencies = readfreelog(OPT['logfile'],free)
                break

    print  ("\n")
    return Ener,Dipoles,Gradients,Modes,Frequencies


# ***************************************************************************** #
# Read into freeformat .log files to extract states energy and dipoles ...
# ***************************************************************************** #

def readfreeformat(logfile):

    global Nroots
    global Natoms
    global GRADcheck
    
    ener = []
    dipo = [[[]]]
    grads = [[[]]]
    GRADcheck='NO' # It may have read gradients from other QC files

    for fs in logfile:
        with open(fs,'r') as f:
            line = f.readline()
            # Search for energy input
            #if "FreeFormat - energies" in line:
            if "free" in line.lower() and "format" in line.lower() and "ener" in line.lower():
                print("\nReading energies from external file. These are assumed to be in cm^-1.\n")
                while True :
                    line = f.readline()
                    if not line: break
                    else: ener.append(float(line.split()[0]))
            # Search for dipo input
            #if "FreeFormat - dipoles" in line:
            if "free" in line.lower() and "format" in line.lower() and "dip" in line.lower():
                print("\nReading dipoles from external file.\n")
                dipo = [[[0 for x in range(Nroots)] for y in range(Nroots)] for z in range(3)]
                while True :
                    line = f.readline()
                    if not line: break
                    else:
                        for m in range(3):
                            dipo[m][int(line.split()[0])][int(line.split()[1])]=float(line.split()[2+m])
                            dipo[m][int(line.split()[1])][int(line.split()[0])]=dipo[m][int(line.split()[0])][int(line.split()[1])]
            if "free" in line.lower() and "format" in line.lower() and "grad" in line.lower():
                
                # If this is the first file in which you detected a gradient -> initialize grads to zero
                if GRADcheck == 'NO':
                    grads = [[[0 for x in range(Natoms)] for y in range(3)] for z in range(Nroots)]
                    print("\nReading gradients from external file. GS should be 'root 1.'\n")
            
                GRADcheck='YES'

                if Natoms == 0:
                    print("Warning: the Natoms was not specified in any other file. This is unexpected, I will stop.")
                    exit(0) #Weaker condition -> remove exit(0), and provide Natoms
                line = f.readline()
                
                # Expected format ".... gradient, root  n"
                # GS root # is 1
                if "gradient" in line.lower() and "root" in line.lower():
                    grad_num = int(line.split()[-1])
                else:
                    print("Warning: I was expecting a line specifying the root number. No such line detected. I will stop.\n")
                    exit(0)

                for l in range(Natoms):
                    line = f.readline()
                    for k in range(3):
                        grads[grad_num-1][k][l] = float(line.split()[k+1])

                ##### TO BE IMPLEMENTED --> GRADIENTS FROM DIFFERENT MINIMA #####
                #if "gradient" in line.lower() and "root" in line.lower():
                #    grad_num = int(line.split()[2])
                #    if "from" in line.lower():
                #       grad_reference[grad_num-1] = int(line.split()[-1])
                #       should declare grad_reference = [] full of Nroots zeros.
                #    else:
                #       grad_reference[grad_num-1] = 0
                #    #the reference numbering should be consistent --> GS is 'root 1'
                #    #syntax : gradient root 1 from root 3
                #else: ...
                #
                #for l in range(Natoms): --> this is the same but grads is grads[state][referece][k][l]
                #                        -> not necessary. Just introduce grads[state][state] if necessary
                #
                ##We should rewrite the gradients here according to the grad_reference
                # for reference
                #   for states =! reference
                #       grads[state][referece][k][l]-=grads[reference][referece][k][l]
                #   grads[reference][referece][k][l] = 0
                #
                #
                ##Now, fix the reference to the GSmin
                # for states
                #      if grad_reference[state] == grad_reference[0]
                #            grads[state][k][l]-=grads[0][k][l]
                #            grad_reference[state] = 0
                #-> now all the states that were from refrence grad_reference[0] are from reference state 0
                #
                #
                ##Now cycle over the states "i" with reference 0 and change the states "j" with reference "i" in this way:
                # finished = 0
                # while finished != Nroots
                # finished = 0
                # for states i
                #   if grad_referece[i] != 0
                #       if grad_reference[grad_reference[i]] == 0
                #           grads[i]=grads[i]+grads[grad_reference[i]]
                #           grad_reference[i]=0
                #           finished +=1
                #   else
                #       finished +=1

    return ener,dipo,grads


# ***************************************************************************** #
# Read into MOLCAS .log files to extract states energy, dipoles ...
# ***************************************************************************** #

def readmollog(logfile):
    
    # ******************************************** #
    # Improvements :
    #   * ...
    #   * ...
    # ******************************************** #

    # Variables
    global Nroots, Nvibmodes, Natoms, atoms, GRADcheck, FREQcheck, grad_list
    ener = []; freqs = []; dipo = [[[]]]; grads = [[[]]]; modes = [[[]]];
    grad_pt2_list = []

    # Initialization
    total_grads = 0
    total_pt2_en = 0
    PT2='MS'
    PT2check='NO'
    gradPT2check='NO'
    GRADcheck='NO'
    FREQcheck='NO'
    
    if Natoms == 0:
        atoms = []
    
    Natoms_check=0
    Nroots_check=0
    
    file_RASSCF = []
    file_SSRASPT2 = []
    file_MSRASPT2 = []
    file_TDMs = []
    file_GRADs = []
    file_FREQs = []
    
    # Open the file(s one by one) and check the input
    # Recognize sections and keywords
    for fs in logfile:
        
        # Search for RASSCF input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "-- -----------------" in line: break
                elif "&rasscf" in line.lower():
                    file_RASSCF.append(fs)

        # Search for RASPT2 input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "-- -----------------" in line: break
                elif "&caspt2" in line.lower():
                    PT2check='YES'
                    while True:
                        line = f.readline()
                        if "nomu" in line.lower():
                            PT2='SS'
                            file_SSRASPT2.append(fs)
                            break
                        elif "-- -----------------" in line:
                            file_MSRASPT2.append(fs)
                            break
                    break

        # Search for TDMs input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "-- -----------------" in line: break
                elif "&rassi" in line.lower():
                    while True:
                        line = f.readline()
                        if "Mltpl 1" in line:
                            file_TDMs.append(fs)
                            break
                        elif "-- -----------------" in line:
                            break
                    break

        # Search for GRADs input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "-- -----------------" in line: break
                elif "&alaska" in line.lower() or "&slapaf" in line.lower():
                    file_GRADs.append(fs)
                    GRADcheck='YES'
                        
        # Search for FREQs input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "-- -----------------" in line: break
                elif "&mckinley" in line.lower():
                    file_FREQs.append(fs)
                    FREQcheck='YES'

        ###### Check atoms number and order consistency ######
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if not line: break
                elif "**** Cartesian Coordinates" in line and Natoms==0:
                    atoms = []
                    Natoms,atoms = read_Natoms_atoms(f,line)
                    break
                elif "**** Cartesian Coordinates" in line and Natoms!=0:
                    atoms_check = []
                    Natoms_check,atoms_check = read_Natoms_atoms(f,line)
                    break
            if Natoms != 0 and Natoms_check != 0:
                if Natoms != Natoms_check:
                    print ("Error: inconsistent number of atoms in different input files")
                    exit(0)
                else:
                    if atoms != atoms_check:
                        print ("Error: inconsistent order of atoms in different input files")
                        exit(0)
        ######################################################

    if PT2check == 'NO':
        print("Warning: NO PT2 file found...")


    # Now read info from various files
    # Info from RASSCF file - if multiple files: check consistency
    if len(file_RASSCF) != 0:
        for fs in file_RASSCF:
            with open(fs,'r') as f:
                while True :
                    line = f.readline()
                    if not line: break
                    # Extract number of roots - CIROOTS
                    elif "ciro" in line.lower():
                        if Nroots == 0:
                            Nroots=read_Nroots(f,line)
                        else:
                            Nroots_check=read_Nroots(f,line)
                            if Nroots != Nroots_check:
                                print("Error: inconsistent number of roots in different input files")
                                exit(0)
                    # Extract energies if there is no PT2 file
                    elif (PT2check == 'NO') and ("Final state energy(ies)" in line):
                        ener = read_energies_RAS_MSPT2(f,line,Nroots,PT2,PT2check)


    # Info from RASPT2 file
    if PT2check == 'YES':
        if PT2 == 'MS':
            if len(file_MSRASPT2) > 1:
                print("\nWarning, multiple MS-PT2 files detected. I will stop.\n (Btw, did you insert the same file multiple times?\n")
                exit(0)
            else:
                with open(file_MSRASPT2[0],'r') as f:
                    while True :
                        line = f.readline()
                        if not line: break
                        # Extract MS-PT2 energies
                        elif "Total MS-CASPT2 energies:" in line:
                            ener = read_energies_RAS_MSPT2(f,line,Nroots,PT2,PT2check)
        else:
            if len(file_MSRASPT2) != 0 and len(file_SSRASPT2) != 0:
                print("\nWarning, MS-PT2 and SS-PT2 files detected. Incompatible inputs. I will stop.\n")
                exit(0)
            for fs in file_SSRASPT2:
                with open(fs,'r') as f:
                    while True :
                        line = f.readline()
                        if not line: break
                        # Extract SS-PT2 energies
                        # It handles situations in which the gradients were not computed for all the states
                        elif "Compute H0 matrices for state" in line:
                            ener = read_energies_SSPT2(f,line,Nroots,PT2,PT2check,ener,total_pt2_en)
                            total_pt2_en+=1

    # Info from TDMs file
    if len(file_TDMs) > 1:
        print("Warning: more than one dipole file provided. One is expected. I will stop.")
        exit(0)
    elif len(file_TDMs) == 1:
        with open(file_TDMs[0],'r') as f:
            while True :
                line = f.readline()
                if not line: break
                # Extracts the dipole moments (length)
                elif "FOR THE RASSCF" in line:
                    dipo=read_dips(f,line,Nroots)

    # Info from GRADs file
    if GRADcheck == 'YES':
        for fs in file_GRADs:
            with open(fs,'r') as f:
                while True :
                    line = f.readline()
                    if not line: break
                    # Check if gradients are ar the caspt2 level
                    elif "&caspt2" in line.lower():
                        gradPT2check='YES'
                    # If PT2 gradients: check the roots
                    elif gradPT2check=='YES' and "Compute H0 matrices for state" in line:
                        grad_pt2_list.append(int(line.split()[-1]))
                    # Extracts analytical and numerical gradients
                    # It handles situations in which the gradients were not computed for all the states
                    elif "Lagrangian multipliers are calculated for state no" in line or "Numerical gradient, root" in line:
                        grads = read_grads(f,line,Nroots,Natoms,total_grads,grads,gradPT2check,grad_pt2_list)
                        total_grads+=1
                if total_grads != 0 and len(grad_pt2_list) != 0:
                    print("Warning: some gradients are at PT2 level, some at RASSCF. I will stop.")
                    exit(0)

    # Info from FREQs and MODEs file
    if FREQcheck == 'YES':
        if len(file_FREQs) > 1:
            print("Warning: more than one freq. file provided. One is expected. I will stop.")
            exit(0)
        with open(file_FREQs[-1],'r') as f:
            while True :
                line = f.readline()
                if not line: break
                # Extracts the normal modes (frequencies and coordinates)
                elif "Numerical differentiation is finished!" in line:
                    Nvibmodes,modes,freqs=read_modes_freqs(f,line,Natoms,atoms)


    # CHECKS
    if total_grads < Nroots and GRADcheck=='YES':
        print("\nWarning: less gradients than roots.....\n")
    elif total_grads > Nroots and GRADcheck=='YES':
        print("Error: more gradients than roots.....")

    # File is closed
    return ener,dipo,grads,modes,freqs


# ***************************************************************************** #
# Read into NWChem .log files to extract states energies, dipoles ...
# ***************************************************************************** #

def readnwclog(logfile):

    # ******************************************** #
    # Improvements :
    #   * Update for imag. or neg. frequencies
    #   * Read multiple "=== echo of input deck ==="
    #   * ...
    # ******************************************** #


    # initialize lists
    global Nroots, Nvibmodes, Natoms, atoms, GRADcheck, FREQcheck
    ener = []; dipo = [[[]]]; grads = [[[]]]; modes = [[[]]]; freqs = []; atoms = []
    
    GRADcheck='NO'
    FREQcheck='NO'

    inside = False  #true when we are inside output block
    start_tag = "Convergence criterion met"
    end_tag = "Excited state energy"
    max_osc_search = 4  #max number of lines after root energy to look for transition dipole moments

    file_EN_TDMs = []
    file_GRADs = []
    file_FREQs = []

    for fs in logfile:
        # Search for energy and TDMs
        with open(fs,'r') as f:
            already_flag = False
            while True :
                line = f.readline()
                if "========" in line and "echo of input deck" not in line: break
                elif already_flag == False and "tddft" in line and "gradient" not in line and "#" not in line:
                    file_EN_TDMs.append(fs)
                    already_flag = True

        # Search for grads
        with open(fs,'r') as f:
            already_flag = False
            while True :
                line = f.readline()
                if "========" in line and "echo of input deck" not in line: break
                elif already_flag == False and "grad" in line and "#" not in line:
                    GRADcheck='YES'
                    file_GRADs.append(fs)
                    already_flag = True


        # Search for FREQs input
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if "========" in line and "echo of input deck" not in line: break
                elif "freq" in line and "#" not in line:
                    FREQcheck='YES'
                    file_FREQs.append(fs)

    # Now read info from various files
    # Info from energy and TDM file
    if len(file_EN_TDMs) != 0:
        with open(file_EN_TDMs[0],'r') as f:
            while True :
                line = f.readline()
                if not line: break

                # Extract number (and type) of atoms
                elif "geometry" in line and "print xyz" in line:
                    while True:
                        line = f.readline()
                        if "end" in line:
                            break
                        else:
                            atoms.append(line.split()[0])
                            Natoms+=1

                # Extract number of roots
                elif "nroots" in line and "#" not in line:
                    Nroots_temp = int(line.split()[1])
                    if Nroots_temp+1 > Nroots:
                        Nroots = Nroots_temp
                        Nroots += 1         # Add the GS
                        dipo = [[[0 for x in range(Nroots)] for y in range(Nroots)] for z in range(3)]
                        ener = [0 for x in range(Nroots)]
                
                elif start_tag in line:
                    inside = True
                elif end_tag in line:
                    inside = False

                # Extract energies & transition dipoles
                elif inside and "Root" in line and "singlet" in line and "eV" in line:
                    line_strip = line.strip()
                    line_split = line_strip.split()
                    try:
                        line_start = line_split[0]     # contains "Root"
                        line_n = line_split[1]         # contains root number (int)
                        line_ev_tag = line_split[-1]   # contains "eV"
                        line_e = line_split[-2]        # contains excitation energy in eV
                    except:
                        raise Exception ("Failed to parse data line for root: {0}".format(line_strip))
                    
                    if line_start == "Root" and line_ev_tag == "eV":
                            try:
                                state_numb = int(line_n)
                                #ener.append(float(line_e)*PhyCon['eV2cmm1'])
                                ener[state_numb]=float(line_e)*PhyCon['eV2cmm1']
                                if float(line_e) < 0.0:
                                    raise Exception ("Invalid negative energy: {0}".format(float(line_e)))
                            except:
                                raise Exception ("Failed to convert root values: {0}".format(line_strip))
                    else:
                        raise Exception ("Unexpected format for root: {0}".format(line_strip))
                    
                    # Now look at the transition dipole moments, which will be a few
                    # lines down (though the exact position may vary it seems).
                    ioscline = -1
                    while True:
                        ioscline = ioscline + 1
                        if ioscline >= max_osc_search:
                            raise Exception ("Failed to find oscillator strength after looking {0} lines.".format(ioscline))

                        line = f.readline()
                        line_strip = line.strip()
                        line_split = line_strip.split()
                        if "Transition Moments" in line and "X" in line and not "XX" in line and not "YY" in line:
                            try:
                                tdmx = float(line_split[-5])     # contains TDMx (float)
                                tdmy = float(line_split[-3])     # contains TDMy (float)
                                tdmz = float(line_split[-1])     # contains TDMz (float)
                            except:
                                raise Exception ("Failed to convert transition moments: {0}".format(line_strip))
                            dipo[0][0][state_numb]=tdmx
                            dipo[0][state_numb][0]=dipo[0][0][state_numb]
                            dipo[1][0][state_numb]=tdmy
                            dipo[1][state_numb][0]=dipo[1][0][state_numb]
                            dipo[2][0][state_numb]=tdmz
                            dipo[2][state_numb][0]=dipo[2][0][state_numb]
                            break
                # END extract energies & transition dipoles

    # Extract frequencies and normal modes
    # Frequencies are in [cm^-1], modes are unitless and orthonormal
    if FREQcheck == 'YES':
        with open(file_FREQs[0],'r') as f:
            # Read properties
            while True :
                line = f.readline()
                if not line: break

                # Extract number (and type) of atoms (if not previously done)
                elif Natoms == 0 and "geometry" in line and "print xyz" in line:
                    while True:
                        line = f.readline()
                        if "end" in line:
                            break
                        else:
                            atoms.append(line.split()[0])
                            Natoms+=1

                elif "MASS-WEIGHTED PROJECTED HESSIAN" in line:
                    if Natoms == 0:
                        print("Error: number of Atoms is zero")
                        exit(0)
                    print  (" ... Reading Normal Modes ...")

                    # Note that global rotation and translations are still written in
                    # the NWChem output, as the first 6 zero frequency modes
                    Nvibmodes=3*Natoms-6 # valid for non-linear molecules
                    modes = [[[0 for x in range(Nvibmodes)] for y in range(Natoms)] for z in range(3)]

                    while True:
                        line = f.readline()
                        if not line: break
                        elif "NORMAL MODE EIGENVECTORS IN CARTESIAN COORDINATES" in line:
                            
                            jump_lines(f,(7+(Nvibmodes+6)+3))

                            # FORMATTING #
                            nmo = (Nvibmodes+6)%6 #skip the first 6 columns, that contains global rot and transl
                            if nmo != 0:
                                nblocks = Nvibmodes//6 + 1 # // is floor division
                            else:
                                nblocks = Nvibmodes//6
                                nmo = 6
                            ##############

                            for j in range(nblocks):
                                # Read and save frequencies
                                line = f.readline()
                                if j == nblocks - 1:
                                    for l in range(nmo):
                                        freqs.append(float(line.split()[l+1]))
                                    jump_lines(f,1)
                            
                                    for k in range(Natoms):
                                        for m in range(3):
                                            line = f.readline()
                                            for l in range(nmo):
                                                modes[m][k][l+j*6] = float(line.split()[l+1])

                                else:
                                    for l in range(6):
                                        freqs.append(float(line.split()[l+1]))
                                    jump_lines(f,1)

                                    for k in range(Natoms):
                                        for m in range(3):
                                            line = f.readline()
                                            for l in range(6):
                                                modes[m][k][l+j*6] = float(line.split()[l+1])
                                
                                if j < nblocks-1:
                                    jump_lines(f,3)


                            # Mass weighted normal modes + normalization
                            # Obtaining unitless normal mode coordinates
                            for m in range(Nvibmodes):
                                norm = 0
                                for k in range(3):
                                    for l in range(Natoms):
                                        modes[k][l][m] *= np.sqrt(AtomMass[atoms[l]])
                                        norm += modes[k][l][m]**2
                                for k in range(3):
                                    for l in range(Natoms):
                                        modes[k][l][m] /= np.sqrt(norm)

                            # Check mode normalization
                            #for m in range(Nvibmodes):
                            #    tot=0
                            #    for k in range(3):
                            #        for l in range(Natoms):
                            #           tot+=modes[k][l][m]*modes[k][l][m]
                            #   print(tot)

                # END extract frequencies and normal modes


    # Extract ES gradients (units are [Hartree/Bohr])
    if GRADcheck == 'YES':
        with open(file_GRADs[0],'r') as f:
            grads = [[[0 for x in range(Natoms)] for y in range(3)] for z in range(Nroots)]
            # Read properties
            while True :
                line = f.readline()
                if not line: break

                elif "NWChem TDDFT Gradient Module" in line:
                    print  (" ... Reading TDDFT gradients ...")
                    
                    while True:
                        line = f.readline()
                        if not line: break
                        
                        # The root number is written also in other places. Choose the best one.
                        elif "Root" in line:
                            root_num = int(line.split()[-1])
                    
                        # Extract the gradient of Root = root_num
                        # grads[root_num][0][i] and not grads[root_num-1][0][i] because should
                        # include the GS, whose gradient is assumed to be zero.
                        elif "TDDFT ENERGY GRADIENTS" in line:
                            jump_lines(f,3)
                            for i in range(Natoms):
                                line = f.readline()
                                line_strip = line.strip()
                                line_split = line_strip.split()
                                grads[root_num][0][i]=float(line_split[-3])
                                grads[root_num][1][i]=float(line_split[-2])
                                grads[root_num][2][i]=float(line_split[-1])
                            break
                # END extract ES gradients



        #end while

    # File is closed
    return ener,dipo,grads,modes,freqs


# ***************************************************************************** #
# Read into gdvh36 .log files to extract states energies, dipoles ...
# ***************************************************************************** #

def readgaulog(logfile):
    
    # ******************************************** #
    # Improvements :
    #   * Now it only reads gradients.. generalize!
    #   * Scale the freq. according to the functional
    #   * ...
    # ******************************************** #
    
    # initialize lists
    global Nroots, Nvibmodes, Natoms, atoms, GRADcheck, FREQcheck
    ener = []; dipo = [[[]]]; modes = [[[]]]; grads = [[[]]]; freqs = [];
    
    try:
        GRADcheck
    except:
        GRADcheck = 'NO'
    try:
        FREQcheck
    except:
        FREQcheck = 'NO'

    if Natoms == 0:
        atoms = []
    
    file_EN_TDMs = []
    file_GRADs = []
    file_FREQs = []

    HP_modes_check = False

    # Frequency scaling factor -> wB97XD
    #freq_scaling = 0.989 #0.975 ?

    # If the file exists then open it and read line by line
    with open(logfile[0],'r') as f:
        # Read properties
        while True:
            line = f.readline()
            if not line: break

            # Extract number of roots
            #elif "nstate=" in line:

            # Extracts the dipole moments (length)
            #elif "electric dipole" in line:

            # Extracts the site energy values
            #elif "nm " in line:


            # Extract number (and type) of atoms
            elif "Input orientation" in line and Natoms==0:
                jump_lines(f,4)
                while True:
                    line = f.readline()
                    if "-----------------------------------" in line:
                        break
                    else:
                        atoms.append(int(line.split()[1]))
                        Natoms+=1

            elif "hpmodes" in line.lower() and "freq" in line.lower():
                HP_modes_check = True

            # Extract frequencies and normal modes
            # Frequencies are in [cm^-1]
            elif "Harmonic frequencies (cm**-1)" in line:
                FREQcheck = 'YES'
                if Natoms == 0:
                    print("Error: number of Atoms is zero")
                    exit(0)

                if len(freqs) == 0:
                    Nvibmodes=3*Natoms-6 # valid for non-linear molecules
                    modes = [[[0 for x in range(Nvibmodes)] for y in range(Natoms)] for z in range(3)]

                jump_lines(f,3)
                line = f.readline()
                if HP_modes_check == True and Nvibmodes>3 and int(line.split()[-1]) == 3:
                    # This happens when you have both high and low precision modes outputs
                    # I will skip the low precision one.
                    break
                # If this does not happend  --> I have found the high precision output.
                # I will use this.

                jump_lines(f,1)

                # HIGH PRECISION OUTPUT FORMATTING
                if HP_modes_check == True:
                    print  (" ... Reading Normal Modes (High Precision) ...")

                    # FORMATTING #
                    nmo = (Nvibmodes)%5
                    if nmo != 0:
                        nblocks = Nvibmodes//5 + 1 # // is floor division
                    else:
                        nblocks = Nvibmodes//5
                        nmo = 5
                    ##############

                    counter=0
                    
                    for j in range(nblocks):
                        # Read and save frequencies
                        line = f.readline()
                        if j == nblocks - 1:
                            for l in range(nmo):
                                counter+=1
                                freqs.append(float(line.split()[l+2]))
                            jump_lines(f,4)
                    
                            for k in range(Natoms):
                                for m in range(3):
                                    line = f.readline()
                                    for l in range(nmo):
                                        modes[m][k][l+j*5] = float(line.split()[l+3])
                        else:
                            for l in range(5):
                                counter+=1
                                if float(line.split()[l+2]) < 0:
                                    print("Warning: negative frequency detected. I will ignore mode nr. %d" % counter)
                                freqs.append(float(line.split()[l+2]))
                            jump_lines(f,4)

                            for k in range(Natoms):
                                for m in range(3):
                                    line = f.readline()
                                    for l in range(5):
                                        modes[m][k][l+j*5] = float(line.split()[l+3])
                        
                        if j < nblocks-1:
                            jump_lines(f,2)

                # LOW PRECISION OUTPUT FORMATTING
                else:
                    print  (" ... Reading Normal Modes (Low Precision) ...")

                    # FORMATTING #
                    nmo = (Nvibmodes)%3
                    if nmo != 0:
                        nblocks = Nvibmodes//3 + 1 # // is floor division
                    else:
                        nblocks = Nvibmodes//3
                        nmo = 3
                    ##############

                    counter=0
                    for j in range(nblocks):
                        # Read and save frequencies
                        line = f.readline()
                        if j == nblocks - 1:
                            for l in range(nmo):
                                counter+=1
                                #append freqs and scale value
                                freqs.append(float(line.split()[l+2]))
                            jump_lines(f,4)
                    
                            for k in range(Natoms):
                                line = f.readline()
                                for l in range(nmo):
                                    for m in range(3):
                                        modes[m][k][l+j*3] = float(line.split()[l*3+2+m])
                        else:
                            for l in range(3):
                                counter+=1
                                if float(line.split()[l+2]) < 0:
                                    print("Warning: negative frequency detected. I will ignore mode nr. %d" % counter)
                                freqs.append(float(line.split()[l+2]))
                            jump_lines(f,4)

                            for k in range(Natoms):
                                line = f.readline()
                                for l in range(3):
                                    for m in range(3):
                                        modes[m][k][l+j*3] = float(line.split()[l*3+2+m])
                        
                        if j < nblocks-1:
                            jump_lines(f,2)


                # Mass weighted normal modes + normalization
                # Obtaining unitless normal mode coordinates
                for m in range(Nvibmodes):
                    norm = 0
                    for k in range(3):
                        for l in range(Natoms):
                            modes[k][l][m] *= np.sqrt(AtomMass[atoms[l]])
                            norm += modes[k][l][m]**2
                    for k in range(3):
                        for l in range(Natoms):
                            modes[k][l][m] /= np.sqrt(norm)

                # Check mode normalization
                #for m in range(Nvibmodes):
                #    tot=0
                #    for k in range(3):
                #        for l in range(Natoms):
                #           tot+=modes[k][l][m]*modes[k][l][m]
                #   print(tot)
                

            # Extracts gradients
            #

        #end while
        
    # File is closed
    return ener,dipo,grads,modes,freqs



# ***************************************************************************** #
# Read into generic .log files to extract states energies, dipoles, gradients ...
# ***************************************************************************** #

def readfreelog(logfile,free):

    #****************** GENERAL IDEA ********************#
    # Once it has recognized which type of QC software   #
    # produced which file, it uses the proper readlog    #
    # function (molcas / gaussian / nwchem) to extract   #
    # the quantities of interest.                        #
    #****************************************************#

    print("\nThe user should check the consistency (Nroots, Natoms, atom order...) of the files coming from different QC softwares\n")
    # The root number in the free format gradients may change if one provides full input for energies or free format energies..

    # initialize lists
    global Nroots, Nvibmodes, Natoms, atoms, GRADcheck, FREQcheck
    ener = []; dipo = [[[]]]; grads = [[[]]]; modes = [[[]]]; freqs = []; atoms = []
    ener1 = []; dipo1 = [[[]]]; grads1 = [[[]]]; modes1 = [[[]]]; freqs1 = []
    ener2 = []; dipo2 = [[[]]]; grads2 = [[[]]]; modes2 = [[[]]]; freqs2 = []
    ener3 = []; dipo3 = [[[]]]; grads3 = [[[]]]; modes3 = [[[]]]; freqs3 = []
    
    GRADcheck='NO'
    FREQcheck='NO'
    Natoms=0
    Nroots=free

    file_OpenMolcas = []
    file_Gaussian = []
    file_NWChem = []
    file_freeformat = []

    for fs in logfile:
        # Identify the used softwares
        with open(fs,'r') as f:
            while True :
                line = f.readline()
                if not line: break
                elif "molcas" in line.lower():
                    file_OpenMolcas.append(fs)
                    break
                elif "Entering Gaussian System, Link 0=g16" in line or "Entering Gaussian System, Link 0=g09" in line:
                    file_Gaussian.append(fs)
                    break
                elif "nwchem" in line:
                    file_NWChem.append(fs)
                    break
                # The free format file is a file similar to the Spectron input file,
                # that can be used to fill the missing info or to force value of
                # desired quantities.
                # Headers (no case sensitive):
                #  * for energy in free format:  Free Format - energies [ener]
                #  * for dipoles in free format: Free Format - dipoles [dipo/dips/dipl]
                #  * for gradients in free format: Free Format - gradients [grad]
                elif "free" in line.lower() and "format" in line.lower():
                    file_freeformat.append(fs)
                    break

    # Since the number of roots has already been given, whatever order with which we read the files is ok.

    if len(file_OpenMolcas) != 0:
        ener1,dipo1,grads1,modes1,freqs1 = readmollog(file_OpenMolcas)
    if len(file_Gaussian) != 0:
        ener2,dipo2,grads2,modes2,freqs2 = readgaulog(file_Gaussian)
    if len(file_NWChem) != 0:
        ener3,dipo3,grads3,modes3,freqs3 = readnwclog(file_NWChem)
    if len(file_freeformat) != 0:
        ener,dipo,grads = readfreeformat(file_freeformat)

    # If redundant info are given (e.g. energies both as an external file and with Gaussian)
    # there is an hierarchy of the method that is followed: on the top there are the
    # externally provided files in free format, then OpenMolcas' files, then Gaussian's
    # and then NWChem's.

    # check energies
    if len(ener) != 0:
        if len(ener1) !=0 or len(ener2) !=0 or len(ener3) !=0:
            print("\nInconsistency detected: same info from multiple files - ener\n")
        print("I will use external energy inputs.")
    elif len(ener1) !=0:
        if len(ener2) !=0:
            print("\nInconsistency detected: same info from multiple files - ener\n")
        elif len(ener3) !=0:
            print("\nInconsistency detected: same info from multiple files - ener\n")
        print("I will use OpenMolcas energy inputs.")
        ener = ener1
    elif len(ener2) != 0:
        if len(ener3) !=0:
            print("\nInconsistency detected: same info from multiple files - ener\n")
        print("I will use Gaussian energy inputs.")
        ener = ener2
    elif len(ener3) != 0:
        print("I will use NWChem energy inputs.")
        ener = ener3
    else:
        print("No energy provided. I will stop.")
        exit(0)

    # check dipoles
    if len(dipo) == 3:
        if len(dipo1) == 3 or len(dipo2) == 3 or len(dipo3) == 3:
            print("\nInconsistency detected: same info from multiple files - dipo\n")
        print("I will use external dipole inputs.")
    elif len(dipo1) == 3:
        if len(dipo2) == 3:
            print("\nInconsistency detected: same info from multiple files - dipo\n")
        elif len(dipo3) == 3:
            print("\nInconsistency detected: same info from multiple files - dipo\n")
        print("I will use OpenMolcas dipole inputs.")
        dipo = dipo1
    elif len(dipo2) == 3:
        if len(dipo3) == 3:
            print("\nInconsistency detected: same info from multiple files - dipo\n")
        print("I will use Gaussian dipole inputs.")
        dipo = dipo2
    elif len(dipo3) == 3:
        print("I will use NWChem dipole inputs.")
        dipo = dipo3
    else:
        print("No dipoles provided. I will stop.")
        exit(0)

    # check gradients
    if len(grads) == Nroots:
        if len(grads1) == Nroots or len(grads2) == Nroots or len(grads3) == Nroots:
            print("\nInconsistency detected: same info from multiple files - grads\n")
        print("I will use external gradient inputs.")
    elif len(grads1) == Nroots:
        if len(grads2) == Nroots:
            print("\nInconsistency detected: same info from multiple files - grads\n")
        elif len(grads3) == Nroots:
            print("\nInconsistency detected: same info from multiple files - grads\n")
        print("I will use OpenMolcas gradient inputs.")
        grads = grads1
    elif len(grads2) == Nroots:
        if len(grads3) == Nroots:
            print("\nInconsistency detected: same info from multiple files - grads\n")
        print("I will use Gaussian gradient inputs.")
        grads = grads2
    elif len(grads3) == Nroots:
        print("I will use NWChem gradient inputs.")
        grads = grads3

    # check modes
    if len(modes1) == 3:
        if len(modes2) == 3 or len(modes3) == 3:
            print("\nInconsistency detected: same info from multiple files - modes\n")
        print("I will use OpenMolcas normal modes inputs.")
        modes = modes1
    elif len(modes2) ==3:
        if len(modes3) ==3:
            print("\nInconsistency detected: same info from multiple files - modes\n")
        print("I will use Gaussian normal modes inputs.")
        modes = modes2
    else:
        print("I will use NWChem normal modes inputs.")
        modes = modes3

    # check frequencies
    if len(freqs1) != 0:
        if len(freqs2) !=0 or len(freqs3) !=0:
            print("\nInconsistency detected: same info from multiple files - freqs\n")
        print("I will use OpenMolcas freq.s inputs.")
        freqs = freqs1
    elif len(freqs2) !=0:
        if len(freqs3) !=0:
            print("\nInconsistency detected: same info from multiple files - freqs\n")
        print("I will use Gaussian freq.s inputs.")
        freqs = freqs2
    else:
        print("I will use NWChem freq.s inputs.")
        freqs = freqs3

    return ener,dipo,grads,modes,freqs


# ***************************************************************************** #
# Fill OPT dictionary with default values. If you want to change the default
# options you can modify here.
# ***************************************************************************** #

OPT =    {
            'read'      : 'molcas',         # Select MOLCAS/GAUSSIAN/NWChem
            'verbosity' : 0 ,               # Verbosity of the output on screen
            'logfile'   : None,             # Molcas/NWChem/Gaussian LogFile
        }

PhyCon =   {
            'au2cmm1'   : 219474.63 ,       # wavenumber per au
            'eV2cmm1'   : 8065.5443 ,       # wavenumber per eV
            }

atommass=[1.007825,4.002603,7.016004,9.012182,11.009305,12.000000,14.003074,15.994914,18.998403,19.992440,22.989769,23.985041,26.981538,27.976926,30.973761,31.972071,34.968852,39.962383,38.963706,39.962590]

AtomMass = {
    'H' : atommass[0], 1 : atommass[0], 'He' : atommass[1], 2 : atommass[1], 'Li' : atommass[2], 3 : atommass[2], 'Be' : atommass[3], 4 : atommass[3], 'B' : atommass[4], 5 : atommass[4], 'C' : atommass[5], 6 : atommass[5], 'N' : atommass[6], 7 : atommass[6], 'O' : atommass[7], 8 : atommass[7], 'F' : atommass[8], 9 : atommass[8], 'Ne' : atommass[9], 10 : atommass[9], 'Na' : atommass[10], 11 : atommass[10], 'Mg' : atommass[11], 12 : atommass[11], 'Al' : atommass[12], 13 : atommass[12], 'Si' : atommass[13], 14 : atommass[13], 'P' : atommass[14], 15 : atommass[14], 'S' : atommass[15], 16 : atommass[15], 'Cl' : atommass[16], 'CL' : atommass[16], 17 : atommass[16], 'Ar' : atommass[17], 18 : atommass[17], 'K' : atommass[18], 19 : atommass[18], 'Ca' : atommass[19], 20 : atommass[19],
            }



# ***************************************************************************** #
#                           MAIN BODY OF THE INTERFACE
# ***************************************************************************** #


# ArgumentParser is used here to generate a user-friendly command-line interface.
parser = arg.ArgumentParser(description="Read NWChem/OpenMolcas/Cobramm/Gaussian .log files, extract energies and dipole moments (and normal modes, gradients, etc.) and build Spectron Inputs for linear and non-linear spectroscopy.",formatter_class=arg.ArgumentDefaultsHelpFormatter)
parser.add_argument('-v',help='Increase the verbosity of the output',action="count",default=0)
parser.add_argument('logfile',nargs='+',help='OpenMolcas/NWChem/Cobramm/Gaussian output log file(s) to be processed')
parser.add_argument('-free',help='Read data from outputs of different QC codes: provide number of roots here (including the GS). External files with energies ([cm^-1], including the GS), dipoles (#i #f x y z) (GS is state 0), gradients ([Hartree/Bohr], Atm x y z) (GS is \'root 1\') may also be provided. These shoud have the headers <<Free Format - energies/dipoles/gradients>>',type=int,default=0)
parser.add_argument('-sig',help='Signal type. Supported signals are LA/PP(+KI+KII)',default='PP', choices=['LA','PP'])
parser.add_argument('-w1i',help='W1 Initial frequency (cm^-1)',type=float)
parser.add_argument('-w1f',help='W1 Final frequency (cm^-1)',type=float)
parser.add_argument('-w3i',help='W3 Initial frequency (cm^-1)',type=float)
parser.add_argument('-w3f',help='W3 Final frequency (cm^-1)',type=float)
parser.add_argument('-tgst',nargs='+',help='Target specific states: list them in any order (GS is state 0)',default=[])
parser.add_argument('-thr',help='Set threshold for bright transitions (only |d|>thr*|dmax| will be considered)',default=0,type=float)
parser.add_argument('-dep',help='Dephasing (cm^-1)',default=150,type=float)
parser.add_argument('-stdis',help='Set (same) static disorder for all states (cm^-1)',default=0,type=float)
parser.add_argument('-nshots',help='Set number of shots',default=1,type=float)
parser.add_argument('-t2',help='Set time t2 in non-linear spectra [-sig PP]',default=0,type=float)
parser.add_argument('-tmp',help='Set bath temperature (K)',default=300,type=float)
parser.add_argument('-diagrams',nargs='+',help='Decide to produce only ESA and/or GSB and/or SE diagrams [-sig PP]',default=['ESA','GSB','SE'])

group = parser.add_argument_group('Developer keywords')
group.add_argument('-nw1',help='Number of samples along W1',type=int,default=200)
group.add_argument('-nw3',help='Number of samples along W3',type=int)
group.add_argument('-fscale',help='Scaling factor for normal modes frequency: w -> w * fscale',default=1,type=float)
group.add_argument('-SDdamp',help='Intra-Molecular Spectral Density damping parameter (cm^-1)',default=5,type=float)
group.add_argument('-lOBO',help='Coupling to the low frequency spectral density (Overdamped Brownian Oscillator) (cm^-1)',default=150,type=float)
group.add_argument('-tOBO',help='Dephasing time of the low frequency spectral density (fs)',default=400,type=float)
# Convert tOBO (fs) in wOBO (cm^-1) !
group.add_argument('-SDfmin',help='Discard normal modes below this value (in cm^-1)',default=300,type=float)
group.add_argument('-SDfmax',help='Discard normal modes above this value (in cm^-1)',default=3000,type=float)
group.add_argument('-tag',help='Add this tag at the end of the folder\'s name',default='')
group.add_argument('--overwrite',help='Use for debugging purposes',action='store_true')
## Add POLARIZATIONS -pol
## Add PULSES SHAPE  -ps

args = parser.parse_args()
OPT['verbosity']=args.v



# List of input files
# It does work also with /*
InLogFile = args.logfile

if args.sig == 'LA':
    directory = 'LA%s/' % args.tag
elif args.sig == 'PP':
    directory = 'PP%s/' % args.tag

# Create proper folder name #
if args.overwrite is True:
    directory = 'aaaTEST_FOLDER_FOR_DEBUG/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        shutil.rmtree(directory)
        os.makedirs(directory)
elif not os.path.exists(directory):
    os.makedirs(directory)
else:
    i=1
    while True:
        if args.sig == 'LA':
            directory = 'LA%s-%d/' % (args.tag,i)
        elif args.sig == 'PP':
            directory = 'PP%s-%d/' % (args.tag,i)
        
        if not os.path.exists(directory):
            os.makedirs(directory)
            break
        else:
            i=i+1

# Info file. Move down to account for variable values
# changes during execution
OutFile = "%sINFO.txt" % directory
dict_args = parser.parse_known_args()[0].__dict__
with open(OutFile,'w') as OutF:
    OutF.write(' '.join(sys.argv))
    OutF.write("\n\n")
    OutF.write("Value for 'key' is set to: XX\n")
    # Print all the args values (default/non-default)
    for key in dict_args:
        OutF.write("  %s  --->\t%s\n" % (key, dict_args[key]))

print("\n********** READ INPUT FILE(S) **********\n")
# Energies are here assumed to be in cm^-1, dipoles in electron charges times Angstrom (eA)
energies,dipo,grads,modes,freqs = readlog(InLogFile,args.free)
# Use normal modes frequency scaling factor if specified by the user.
for i in range(len(freqs)):
    freqs[i] *= args.fscale

if len(energies) == 0 or len(dipo) == 0:
    print("No energies and/or no dipoles detected. I will stop.")
    exit(0)

SDcheck = 'NO'
if GRADcheck == 'YES' and FREQcheck == 'YES':
    SDcheck = 'YES'

# Reorder transitions in increasing energy
# May be necessary in RASSCF and SS-PT2 computations
state_indexes = [x for x in range(Nroots)]
list1, list2 = (list(t) for t in zip(*sorted(zip(energies, state_indexes))))

#Create trasnformation matrix
transf = np.array([[0 for x in range(Nroots)] for y in range(Nroots)])
for i in range(Nroots):
    for j in range(Nroots):
        if j == list2[i]: transf[i][j]=1

# Reorder target states
sorted_tgst_list = [] # Necessary for the module that look at tgst
if len(args.tgst) != 0:
    tgst_list = [0 for x in range(Nroots)]
    sorted_tgst_list = [0 for x in range(Nroots)]
    for i in range(Nroots):
        for j in args.tgst:
            if i == int(j): tgst_list[i] = 1
    sorted_tgst_list = transf.dot(tgst_list)

# Reorder energies
sorted_energies = [0 for x in range(Nroots)]
sorted_energies = transf.dot(energies)
# Subtract the GS energy
sorted_energies=sorted_energies-sorted_energies[0]

new_dips = [[[0 for x in range(Nroots)] for y in range(Nroots)] for z in range(3)]
#reorder dipole moments according to new energy order
new_dips[0]=(transf.dot(dipo[0]).dot(transf.transpose()))
new_dips[1]=(transf.dot(dipo[1]).dot(transf.transpose()))
new_dips[2]=(transf.dot(dipo[2]).dot(transf.transpose()))

print("******* READ INPUT FILE(S) - DONE ******* \n")

#### SQUARED WINDOW IS ASSUMED if args.w3i and args.w3f are not set  ####
if args.w3i == None:
    args.w3i=args.w1i
if args.w3f == None:
    args.w3f=args.w1f
if args.nw3 == None:
    args.nw3 = args.nw1


if SDcheck == 'YES':
    print("\n********* COMPUTE DHO DISPLACEMENTS **********\n")

    displ = [[0 for x in range(Nvibmodes)] for y in range(Nroots)]
    
    if args.sig != 'LA':
        print("\nThe present implementation assumes that all the gradients\nare computed from the same reference geometry\n (tipically the optimized GS geometry).")
        # One may have gradients of different states computed from different positions (e.g. minimum of an e-manifold state).
        # If necessary the interface may be updated to handle that "mixed" situation.
    
    # Remove very low and very high frequecies
    modes_removed=0
    counts_modes_removed_below_threshold=0
    print("\nModes with frequency below %.1lf cm^-1 and above %.1lf cm^-1 are removed.\n" % (args.SDfmin,args.SDfmax))
    for m in range(Nvibmodes):
        if (freqs[m-modes_removed] < args.SDfmin) or (freqs[m-modes_removed] > args.SDfmax):
            if (freqs[m-modes_removed] < args.SDfmin): counts_modes_removed_below_threshold+=1
            freqs.pop(m-modes_removed)
            for j in range(Natoms):
                modes[0][j].pop(m-modes_removed)
                modes[1][j].pop(m-modes_removed)
                modes[2][j].pop(m-modes_removed)
            modes_removed+=1

    Nvibmodes=Nvibmodes-modes_removed

    # eq. 20 in SI of Azo paper: d_i = -W^{-2} P M^{-1/2} f_i with
    #    d_i: array of the displacements on state i [Bohr*sqrt(amu)]
    #    W: normal mode frequencies [cm-1]
    #    P: normal mode matrix [unitless]
    #    M: atomic masses [amu]
    #    f: gradient of i-th state [Hartree/Bohr]
    # input units: Hartree / (Bohr * sqrt(amu) * (cm-1)**2)
    # output units: Bohr*sqrt(amu)
    # gradient [Hartree/Bohr] -> [4.36*1e-18 J / 5.3*1e-11 m]
    # normal modes: unitless
    # mass [sqrt(amu)] -> [sqrt(0.166)*1e-13 kg]
    # frequency [cm-1] -> [2.997*1e+10 1/s]
    # conversion from m*sqrt(kg) to Bohr*sqrt(amu) -> 1.889*1e+10*sqrt(6.022)*1e+13
    # SIGN depends on whether one uses the gradient (sign is -) OR the force (sign is +) as it is the force that is projected on the normal modes
    # this implementation expects gradients

    # Reorder grads to be consistent with energy sorting
    # --> all the other quantities (HR, reorg_en) will be reordered consistently
    # Implemented with grads[st][k][l] --> grads[list2[st]][k][l]
    displ_const = (1.889*1e+10*np.sqrt(6.022)*1e+13)*(4.36*1e-18/(5.29*1e-11))/(np.sqrt(0.166)*1e-13*((2.0*np.pi*2.997*1e+10)**2))

    for st in range(Nroots):
        for m in range(Nvibmodes):
            for k in range(3):
                for l in range(Natoms):
                    displ[st][m]-=displ_const*grads[list2[st]][k][l]*modes[k][l][m]/(np.sqrt(AtomMass[atoms[l]])*(freqs[m]**2))
            # rescale displacement with respect to GS displacement. This should be very small or zero if we are at the GS minimum, and the GS opt was done at the same level of theory of the gradients.
            if st != 0:
                displ[st][m]-=displ[0][m]


    print("\n******* COMPUTE DHO DISPLACEMENTS - DONE *******\n")


print("\n******** PREPARE SPECTRON-INPUT *********\n")

OutFile_evals_name = "SOS_evals.txt"
OutFile_edips_name = "SOS_edips.txt"
if SDcheck == 'YES':
    OutFile_sd_name = "spectral_densities.txt"
    OutFile_sbcoup_name = "sb_couplings.txt"
    OutFile_lifet_name = "transport_rates.txt"
    #OutFile_stdis_name = "static_disorder.txt"

OutFile_evals = "%s%s" % (directory,OutFile_evals_name)
OutFile_edips = "%s%s" % (directory,OutFile_edips_name)
if SDcheck == 'YES':
    OutFile_sd = "%s%s" % (directory,OutFile_sd_name)
    OutFile_sbcoup = "%s%s" % (directory,OutFile_sbcoup_name)
    OutFile_lifet = "%s%s" % (directory,OutFile_lifet_name)
    #OutFile_stdis = "%s%s" % (directory,OutFile_stdis_name)


### WHICH DIAGRAMS WERE SELECTED ###
ESAcheck=0
GSBcheck=0
SEcheck=0
for diag in args.diagrams:
    if diag == "ESA":
        ESAcheck=1
    elif diag == "SE":
        SEcheck=1
    elif diag == "GSB":
        GSBcheck=1

### IF NSHOTS = 1 and STATIC DIS != 0: WARNING ###
if args.stdis != 0 and args.nshots == 1:
    print("WARNING:")
    print("You must increase the number of shots to be performed (-nshots ).")
    print("I will set them to 100.\n")
    args.nshots = 100


### CREATE FIRST EXC. MANIFOLD ###

#trasform arrays in lists to use .pop() for removing list elements
sorted_energies=sorted_energies.tolist()
new_dips[0]=new_dips[0].tolist()
new_dips[1]=new_dips[1].tolist()
new_dips[2]=new_dips[2].tolist()


# If a list of target states has been specified, remove all the states
# that do not belong to that list
if len(sorted_tgst_list) != 0:
    if OPT['read']=='molcas':
        print("Warning: if energies are at MS-PT2 level and dipoles/gradients at RASSCF (or SS-PT2) level \n there may be a wrong  assignment of states properties (reordering issue)...")
    removed=0
    for i in range(1,Nroots):
        check_ij=0
        if(int(sorted_tgst_list[i])==1):
            check_ij=1
        if check_ij == 0:
            sorted_energies.pop(i-removed)
            new_dips[0].pop(i-removed)
            new_dips[1].pop(i-removed)
            new_dips[2].pop(i-removed)
            if SDcheck == 'YES':
                displ.pop(i-removed)
            for j in range(Nroots-removed-1): #-1 because I have already removed 1 row
                new_dips[0][j].pop(i-removed)
                new_dips[1][j].pop(i-removed)
                new_dips[2][j].pop(i-removed)
            removed+=1
    Nroots_reduced=Nroots-removed
else:
    Nroots_reduced=Nroots

first_exc_energies = []

if args.w1i-args.dep/2. <= 0:
    nummodes=-1 #remove the GS if it has been included in the 1st exc. manifold (probably the window has not been set properly).
else:
    nummodes=0

for i in range(Nroots_reduced):
    #NB could be necessary to include also lower lying states that are not in the window, but that can be populated via population transfer and from which there can be bright ESA signals in the window
    if args.w1i-args.dep/2. <= sorted_energies[i] <= args.w1f+args.dep/2.:
        first_exc_energies.append(float(sorted_energies[i]))
        nummodes+=1

esnumst=0
if args.sig != 'LA' and ESAcheck == 1:
    # Usually 1st exc. manifold and 2nd exc. manifold are well separated; if they are not (i.e. some higher lying states of the 1st exc. manifold can be reached from lower lying states of the same manifold, in the considered window) one should include some states in both manifolds (or add these states in both manifolds, if Spectron complains about having one single state in both manifolds). This is a very rare case, so is not accounted for here. One can for example split the computation in two: one containing all GSB and SE signals, the other only the ESA.
    for j in range(nummodes+1,Nroots_reduced): #nummodes+1 to accelerate (NB since I have not yet removed the states lying outside the window, by starting from nummodes+1 I'm probably accounting for some of them, but still speeding up the computation).
        for i in range(nummodes): #considering only the transitions from the previously selected states of the first exc. manif.
            if args.w3i-args.dep/2. <= sorted_energies[j]-first_exc_energies[i] <= args.w3f+args.dep/2.:
                esnumst+=1
                break

# First exc. manifold + Second exc. manifold + GS
esnumst+=nummodes+1

# Delete the transitions lying below the window
removed=0
for i in range(1,Nroots_reduced):
    if sorted_energies[i-removed] < args.w1i-args.dep/2.:
        sorted_energies.pop(i-removed)
        new_dips[0].pop(i-removed)
        new_dips[1].pop(i-removed)
        new_dips[2].pop(i-removed)
        if SDcheck == 'YES':
            displ.pop(i-removed)
        for j in range(Nroots_reduced-removed-1): #-1 because I have already removed 1 row
            new_dips[0][j].pop(i-removed)
            new_dips[1][j].pop(i-removed)
            new_dips[2][j].pop(i-removed)
        removed+=1

Nroots_reduced-=removed

print(" ... %2d Roots removed from first exc. manifold\n\t\t which now has %2d states ..." % (removed, nummodes))

if args.sig != 'LA' and ESAcheck == 1:
    #check if some higher lying states are outside the probed window when reached from whatever initial state
    removed=0
    for j in range(nummodes+1,Nroots_reduced):
        flag=0
        for i in range(nummodes):
            if (sorted_energies[j-removed]-first_exc_energies[i] > args.w3f+args.dep/2.) or (sorted_energies[j-removed]-first_exc_energies[i] < args.w3i-args.dep/2.):
                flag+=1
        # If the signal is always outside the window (i.e. from whatever starting state): remove it from the list
        if flag == nummodes:
            sorted_energies.pop(j-removed)
            new_dips[0].pop(j-removed)
            new_dips[1].pop(j-removed)
            new_dips[2].pop(j-removed)
            if SDcheck == 'YES':
                displ.pop(j-removed)
            for k in range(Nroots_reduced-removed-1): #-1 because I have already removed 1 row
                new_dips[0][k].pop(j-removed)
                new_dips[1][k].pop(j-removed)
                new_dips[2][k].pop(j-removed)
            removed+=1

    Nroots_reduced-=removed
    print(" ... %2d Roots removed from second exc. manifold\n\t\t which now has %2d states ..." % (removed,esnumst-1-nummodes))
elif args.sig == 'LA':
    Nroots_reduced = nummodes+1
    
print(" ... Original n. of states: %2d, Actual n. of states: %2d ..." % (Nroots,Nroots_reduced))


# ************************************************* #
# Compute SPECTRAL DENSITIES from the displacements
# Part 1...
# ************************************************* #

if SDcheck == 'YES':

    # Here I use only displacements from the GS (taken as the reference) to the e states.
    # With these I will compute only HR_ee.
    
    # Huang-Ryes Factors per mode -> S_ii=omega_i*d_i**2 / 2 hbar [unitless]
    # eq. 14 from Azo SI
    # input units: Bohr**2 * sqrt(amu)**2 * cm-1
    # output units: unitless
    # d [Bohr*sqrt(amu)] -> [5.3*1e-11 m * sqrt(0.166)*1e-13 kg]
    # frequency [cm-1] -> [2*pi*2.997*1e+10 1/s]
    # hbar = 1.05*1e-34

    HR_ee = [[0 for x in range(Nvibmodes)] for y in range(nummodes)]
   
    # Reorg. energy --> lambda_i = S_i * omega_i [cm^-1]
    # eq. 14 from Azo SI
    reorg_en_ee = [[0 for x in range(Nvibmodes)] for y in range(nummodes)]
    total_reorg_en_ee = [0 for x in range(nummodes)]

    # Build and Write Undamped SD --> C(omega) = Sum_i S_i * omega_i **2 delta(omega-omega_i)
    undamped_SDh = [[0 for x in range(Nvibmodes)] for y in range(nummodes)]
    hr_constant = 2*np.pi*2.997*1e+10*((5.29*1e-11*np.sqrt(0.166)*1e-13)**2)/(2*1.05*1e-34)

    for st in range(nummodes):
        for m in range(Nvibmodes):
            #HR_ee[st][m]= (2*np.pi*freqs[m]*2.997*1e+10)*((displ[st+1][m]*5.29*1e-11*np.sqrt(0.166)*1e-13)**2)/(2*1.05*1e-34);
            HR_ee[st][m] = freqs[m]*(displ[st+1][m]**2)*hr_constant
            reorg_en_ee[st][m] = HR_ee[st][m]*freqs[m]
            total_reorg_en_ee[st]+=reorg_en_ee[st][m]
            undamped_SDh[st][m]=HR_ee[st][m]*freqs[m]**2

    if args.sig == 'PP' and ESAcheck == 1:
    
        # Here I have:
        #    * displacements from the GS (taken as the reference) to the e states.
        #    * displacements from the GS (taken as the reference) to the f states.
        # With these I will compute:
        #    * HR_ff
        #    * HR_ef_p (positive)
        #    * HR_ef_n (negative)

        # Huang-Ryes Factors per mode -> S[ef]_i=omega_i*d[e]_i*d[f]_i / 2 hbar [unitless]
        # eq. 14 from Azo SI
        # input units: Bohr**2 * sqrt(amu)**2 * cm-1
        # output units: unitless
        # d [Bohr*sqrt(amu)] -> [5.3*1e-11 m * sqrt(0.166)*1e-13 kg]
        # frequency [cm-1] -> [2*pi*2.997*1e+10 1/s]
        # hbar = 1.05*1e-34

        HR_ff = [[0 for x in range(Nvibmodes)] for y in range(esnumst-1-nummodes)]
        HR_ef_p = [[[0 for x in range(Nvibmodes)] for y in range(esnumst-1-nummodes)] for z in range(nummodes)]
        HR_ef_n = [[[0 for x in range(Nvibmodes)] for y in range(esnumst-1-nummodes)] for z in range(nummodes)]
        use_HR_ef_p = [[0 for x in range(esnumst-1-nummodes)] for y in range(nummodes)] # 0 -> not use , 1 -> use
        use_HR_ef_n = [[0 for x in range(esnumst-1-nummodes)] for y in range(nummodes)] # 0 -> not use , 1 -> use
        
        # Reorg. energy --> lambda_ij = S[ef]_i * omega_i [cm^-1]
        # eq. 14 from Azo SI
        reorg_en_ff = [[0 for x in range(Nvibmodes)] for y in range(esnumst-1-nummodes)]
        total_reorg_en_ff = [0 for x in range(esnumst-1-nummodes)]
        reorg_en_ef = [[[0 for x in range(Nvibmodes)] for y in range(esnumst-1-nummodes)] for z in range(nummodes)]
        total_reorg_en_ef = [[0 for x in range(esnumst-1-nummodes)] for y in range(nummodes)]
        
        # Build and Write Undamped SD --> C(omega) = Sum_i S_i * omega_i **2 delta(omega-omega_i)
        # NOT NECESSARY

        for st in range(nummodes):
            for st1 in range(esnumst-1-nummodes):
                for m in range(Nvibmodes):
                    displs_product = displ[st+1][m]*displ[1+nummodes+st1][m]
                    if st ==0:
                        #HR_ff[st1][m]= (2*np.pi*freqs[m]*2.997*1e+10)*((displ[1+nummodes+st1][m]*5.29*1e-11*np.sqrt(0.166)*1e-13)**2)/(2*1.05*1e-34);
                        HR_ff[st1][m]= freqs[m]*(displ[1+nummodes+st1][m]**2)*hr_constant
                        reorg_en_ff[st1][m] = HR_ff[st1][m]*freqs[m]
                        total_reorg_en_ff[st1] += reorg_en_ff[st1][m]
                    if displs_product > 0:
                        #HR_ef_p[st][st1][m] = (2*np.pi*freqs[m]*2.997*1e+10)*displs_product*((5.29*1e-11*np.sqrt(0.166)*1e-13)**2)/(2*1.05*1e-34);
                        HR_ef_p[st][st1][m] = freqs[m]*displs_product*hr_constant
                        use_HR_ef_p[st][st1]=1
                    elif displs_product < 0:
                        # I save it as positive (revert the sign) but in the HR_ef_n factor
                        HR_ef_n[st][st1][m] = -freqs[m]*displs_product*hr_constant
                        use_HR_ef_n[st][st1]=1

                    # HR_ef can be either + or -, this justifies the following expression
                    reorg_en_ef[st][st1][m] = (HR_ef_p[st][st1][m]+HR_ef_n[st][st1][m])*freqs[m]
                    total_reorg_en_ef[st][st1]+=reorg_en_ef[st][st1][m]


    #*********** Write to file the SDs if verbose output ************#
    if OPT['verbosity'] > 0:
        print("\nSave to file SDs: not implemented yet.\n")
        # COPY IT AND ADAPT HERE FROM PREVIOUS VERSION.
    #*********** DONE - Write to file the SDs if verbose output ************#



# ************************************************* #
# Apply THRESHOLD for transition brightness
# Store correctly formatted Energies and trans. dips.
# ************************************************* #
# We don't delete first excitation manifold states which are dark,
# because they can be populated through energy trasfer and could have
# bright ESAs in the window of interest


# find transition dipoles module max, and discard transitions accordingly
dip_module2 = []
maxx = -1000.
for i in range(nummodes):
    #i+1 to account for the GS
    dip_module2.append(float(new_dips[0][0][i+1]**2+new_dips[1][0][i+1]**2+new_dips[2][0][i+1]**2))
    if dip_module2[i] > maxx:
        maxx = dip_module2[i]

threshold = np.sqrt(maxx)*float(args.thr)

# Energies and transition from GS to first exc. manifold
with open(OutFile_evals,'w') as OutF_ene:
    OutF_ene.write("0\n")
    with open(OutFile_edips,'w') as OutF_dip:
        for i in range(nummodes):
            OutF_ene.write("%.2f\n" % first_exc_energies[i] )
            if np.sqrt(dip_module2[i]) > threshold:
                #i+1 to account for the GS
                OutF_dip.write(" 0 %2d %8.4f %8.4f %8.4f\n" % (i+1,new_dips[0][0][i+1],new_dips[1][0][i+1],new_dips[2][0][i+1]))

count_sd_ef = 0

# Only for multidimensional techniques
if (args.sig != 'LA') and (esnumst > nummodes+1):
    with open(OutFile_evals,'a') as OutF_ene:
        for i in range(nummodes):
            dip_module2 = []
            maxx = -1000.
            k=0
            for j in range(nummodes+1,Nroots_reduced):
                dip_module2.append(float(new_dips[0][i+1][j]**2+new_dips[1][i+1][j]**2+new_dips[2][i+1][j]**2))
                if dip_module2[k] > maxx:
                    maxx = dip_module2[k]
                k+=1
            threshold = np.sqrt(maxx)*float(args.thr)
            with  open(OutFile_edips,'a') as OutF_dip:
                k=0
                for j in range(nummodes+1,Nroots_reduced):
                    if i==0: OutF_ene.write("%.2f\n" % sorted_energies[j] )
                    if np.sqrt(dip_module2[k]) > threshold:
                        OutF_dip.write("%2d %2d %8.4f %8.4f %8.4f\n" % (i+1,j,new_dips[0][i+1][j],new_dips[1][i+1][j],new_dips[2][i+1][j]) )
                        count_sd_ef+=2
                    # If the transition e->f is dark, I can avoid to consider that ef spectral density.
                    else:
                        use_HR_ef_p[i][j-nummodes-1]=0
                        use_HR_ef_n[i][j-nummodes-1]=0
                    k+=1


# ************************************************* #
# Compute SPECTRAL DENSITIES from the displacements
# Part 2...
# ************************************************* #

if SDcheck == 'YES':
    
    #*********** Write spectral density ************#
    # Write spectral density file with correct format ( Nfreq*Fstep must be equal to Maxfreq)
    # Number of frequencies
    # Freq. step
    # column of values ...

    # Build Damped SD --> C(omega)=Sum_i S_i*omega_i**3*omega*gamma/((omega**2-omega_i**2)**2+2gamma**2*omega**2)
    # GAMMA is the damping parameter

    # I need to compute SD_ee, SD_ff, SD_ef_p, SD_ef_n
    # SD_ef_p(/n) is computed only if use_HR_ef_p(/n)[st][st1] != 0 (which is the flag that tells that the transition between st -> st1 is dark

    Dw = 1.0 # [cm^-1]
    wmax=int(max(args.w1f/Dw,args.w3f/Dw))
    #wmax=2000
    gamma = args.SDdamp #cm^-1

    SDh_ee = [[0 for x in range(wmax)] for y in range(nummodes)]

    print(" ")
    for st in range(nummodes):
        print("Processing SDh_ee number %d/%d" % (st+1,nummodes))
        for w in range(0,wmax):
            # consider only the range in which SD is effectively different from zero.
            if w < args.SDfmax: #should be ~ the same value of the freq cutoff when extracting the normal modes
                for m in range(Nvibmodes):
                    SDh_ee[st][w]+=(2*np.sqrt(2)*HR_ee[st][m]*(freqs[m]**3)*gamma*w*Dw)/(((w*Dw)**2-freqs[m]**2)**2+2*(gamma**2)*((w*Dw)**2))

    if args.sig == 'PP' and ESAcheck == 1:
        
        SDh_ff = [[0 for x in range(wmax)] for y in range(esnumst-1-nummodes)]
        SDh_ef_p = [[[0 for x in range(wmax)] for y in range(esnumst-1-nummodes)] for z in range(nummodes)] #plus sign
        SDh_ef_n = [[[0 for x in range(wmax)] for y in range(esnumst-1-nummodes)] for z in range(nummodes)] #minus sign

        l=0
        k=0
        for st in range(nummodes):
            for st1 in range(esnumst-1-nummodes):
                print("Processing SDh_ff %d/%d" % (st1+1,esnumst-1-nummodes))
                for w in range(0,wmax):
                    # consider only the range in which SD is effectively different from zero.
                    if w < args.SDfmax: #should be ~ the same value of the freq cutoff when extracting the normal modes
                        for m in range(Nvibmodes):
                            if st == 0:
                                SDh_ff[st1][w]+=(2*np.sqrt(2)*HR_ff[st1][m]*(freqs[m]**3)*gamma*w*Dw)/(((w*Dw)**2-freqs[m]**2)**2+2*(gamma**2)*((w*Dw)**2))
                            if use_HR_ef_p[st][st1] != 0:
                                if(w==0 and m==0):
                                    print("Processing SDh_ef_p %de->%df (%d/%d)" % (st+1,st1+1,k+1,count_sd_ef//2))
                                    k+=1
                                SDh_ef_p[st][st1][w]+=(2*np.sqrt(2)*HR_ef_p[st][st1][m]*(freqs[m]**3)*gamma*w*Dw)/(((w*Dw)**2-freqs[m]**2)**2+2*(gamma**2)*((w*Dw)**2))
                            if use_HR_ef_n[st][st1] != 0:
                                if(w==0 and m==0):
                                    print("Processing SDh_ef_n %de->%df (%d/%d)" % (st+1,st1+1,l+1,count_sd_ef//2))
                                    l+=1
                                SDh_ef_n[st][st1][w]+=(2*np.sqrt(2)*HR_ef_n[st][st1][m]*(freqs[m]**3)*gamma*w*Dw)/(((w*Dw)**2-freqs[m]**2)**2+2*(gamma**2)*((w*Dw)**2))
            

    # Add a phenomenological OBO inter-molecular SD
    # The same for all the transitions
    print("\nAdding a (common) OBO inter-molecular (low-frequency) SD")
    SDl = [0 for x in range(wmax)]

    lOBO = args.lOBO #[cm^-1]
    wOBO = 33356.4/args.tOBO #[cm^-1]

    for w in range(0,wmax):
        SDl[w]=2*lOBO*(w*Dw)*wOBO/((w*Dw)**2+wOBO**2)


    # Write to file the SDs
    # NSDs counts the number of spectral densities
    NSDs=0
    with open(OutFile_sd,'w') as OutF:
        OutF.write("%d\n%.2lf\n" % (wmax,Dw))
        NSDs+=1
        # Write SDl
        for w in range(0,wmax):
            OutF.write("%.2lf\n" % (SDl[w]))
        # Write SDh_ee
        for st in range(nummodes):
            OutF.write("%d\n%.2lf\n" % (wmax,Dw))
            NSDs+=1
            for w in range(0,wmax):
                OutF.write("%.2lf\n" % (SDh_ee[st][w]))
        # If PP write the other
        if args.sig == 'PP' and ESAcheck == 1:
            # Write SDh_ff
            for st1 in range(esnumst-1-nummodes):
                OutF.write("%d\n%.2lf\n" % (wmax,Dw))
                NSDs+=1
                for w in range(0,wmax):
                    OutF.write("%.2lf\n" % (SDh_ff[st1][w]))
            # Write SDh_ef when the transition e->f is not dark
            for st in range(nummodes):
                for st1 in range(esnumst-1-nummodes):
                    if use_HR_ef_p[st][st1]!=0:
                        OutF.write("%d\n%.2lf\n" % (wmax,Dw))
                        NSDs+=1
                        for w in range(0,wmax):
                            OutF.write("%.2lf\n" % (SDh_ef_p[st][st1][w]))
                    if use_HR_ef_n[st][st1]!=0:
                        OutF.write("%d\n%.2lf\n" % (wmax,Dw))
                        NSDs+=1
                        for w in range(0,wmax):
                            OutF.write("%.2lf\n" % (SDh_ef_n[st][st1][w]))
            

    # Write to file "coupling-triangles"
    # For every SD you have a triangle of couplings which is esnumst * esnumst
    # Every position is identified by two indexes, which are the indexes
    # of the given transition that is coupled with that SD.

    # Example: SDl (blank lines are not printed, but make here easyer to distinguish the sectors)
    # de = 1, df = 2 --> de*de = 1 , df*df = 4 , de*df = 2

    # Here it is implemented the correlated+anticorrelated option for ESA (i.e. round shaped peaks)
    # If you want to modify that manually, just multiply everything by the same scaling factor.

    # GS 0
    #
    # e  0   1
    #    0   0 1
    #    0   0 0 1
    #
    # f  0   0 0 0   4
    #    0   0 2 0   0 4
    #
    #    GS  e       f
    
    # SDh_ee come first;
    # IF PP we have SDh_ff first, and then the various SDh_ef_p, SDh_ef_n, SDh_ef_p, SDh_ef_n ..
    with open(OutFile_sbcoup,'w') as OutF_sbcoup:
        # Low frequency SDl
        # Implemented to have circular ESA (sum of correlated - de=1/sqrt(2), df=sqrt(2) -
        # and anticorrelated - de=1/sqrt(2) df=0 - ESA peaks)..
        # May not be the correct choice in specific cases, but without additional information this is probably the best.
        for i in range(esnumst):
            for j in range(i+1):
                # ee region
                if i == j and i != 0 and i<nummodes+1:
                    OutF_sbcoup.write("1 ")
                # ff region
                elif args.sig == 'PP' and i == j and i != 0 and i>=nummodes+1:
                    OutF_sbcoup.write("2 ") # ("0 ") for peak anticorrelated #("4 ") for peak correlated
                # ef region
                elif args.sig == 'PP' and i >= nummodes+1 and j!= 0 and j < nummodes+1:
                    if use_HR_ef_p[j-1][i-1-nummodes]!=0 and use_HR_ef_n[j-1][i-1-nummodes] != 0:
                        # they are both zero when the transition was removed because it was dark
                        # May be necessary to change this if I remove selectively some SDs_ef_p/n
                        OutF_sbcoup.write("1 ") # ("0 ") for peak anticorrelated #("2 ") for peak correlated
                    else:
                        OutF_sbcoup.write("0 ")
                else:
                    OutF_sbcoup.write("0 ")
            OutF_sbcoup.write("\n")
        # e-manifold SDhs
        for sd in range(nummodes):
            for i in range(esnumst):
                for j in range(i+1):
                    if i == j and i != 0 and i == sd+1:
                        OutF_sbcoup.write("1 ")
                    else:
                        OutF_sbcoup.write("0 ")
                OutF_sbcoup.write("\n")
        # IF PP write the other SDs
        if args.sig == 'PP' and ESAcheck == 1:
            # f-manifold SDhs
            for sd in range(esnumst-1-nummodes):
                for i in range(esnumst):
                    for j in range(i+1):
                        if i == j and i == sd+nummodes+1:
                            OutF_sbcoup.write("1 ")
                        else:
                            OutF_sbcoup.write("0 ")
                    OutF_sbcoup.write("\n")
            for st in range(nummodes):
                for st1 in range(esnumst-1-nummodes):
                    if use_HR_ef_p[st][st1]!=0:
                        # ed-manifold SDhs_plus
                        for i in range(esnumst):
                            for j in range(i+1):
                                if i != 0 and i >= nummodes+1 and j!= 0 and j < nummodes+1:
                                    if st+1 == j and st1+1+nummodes == i:
                                        OutF_sbcoup.write("1 ")
                                    else:
                                        OutF_sbcoup.write("0 ")
                                else:
                                    OutF_sbcoup.write("0 ")
                            OutF_sbcoup.write("\n")
                    if use_HR_ef_n[st][st1]!=0:
                        # ed-manifold SDhs_minus
                        for i in range(esnumst):
                            for j in range(i+1):
                                if i != 0 and i >= nummodes+1 and j!= 0 and j < nummodes+1:
                                    if st+1 == j and st1+1+nummodes == i:
                                        OutF_sbcoup.write("-1 ")
                                    else:
                                        OutF_sbcoup.write("0 ")
                                else:
                                    OutF_sbcoup.write("0 ")
                            OutF_sbcoup.write("\n")


    # Write to file the matrix of rates and lifetimes [set all to 0 -> infinite lifetime]
    with open(OutFile_lifet,'w') as OutF_lifet:
        for i in range(nummodes+1):
            OutF_lifet.write("0.0 "*(nummodes+1))
            OutF_lifet.write("\n")


    # Static disorder --> in case set also Nshots
    #with open(OutFile_stdis,'w') as OutF_stdis:
    #   for i in range(esnumst): --> if PP # or nummodes --> if LA
    #


    # Write to file info about intra-molecular SDs
    # The mode-index used here is the same of the mode-index
    # in the output of the frequency computation to
    # simplify the modes identification.
    OutFile = "%sVIB_INFO.txt" % directory
    with open(OutFile,'w') as OutF:
        OutF.write("ONLY displacements different from zero are printed\n")
        OutF.write("Modes with frequency smaller than % .2f and greater than % .2f [cm^-1] are deleted\n" % (args.SDfmin,args.SDfmax) )
        # ee
        OutF.write("*******************************    e-manifold    *******************************\n\n")
        for st in range(0,nummodes):
            OutF.write("State nr. %d (Total Reorg. Energy = %.2lf)\n\n" % (st+1,total_reorg_en_ee[st]))
            OutF.write(" mode nr.   Freq.[cm^-1]   Displ [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
            for m in range(Nvibmodes):
                if(displ[st+1][m]!=0): OutF.write("   %02d\t     % .2f \t        % .2e \t     %.2e\t     %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[st+1][m], HR_ee[st][m], reorg_en_ee[st][m]))
            OutF.write("\n")
            
            OutF.write("Main modes (Mode Reorg. En. > 15% of the total) \n\n")
            OutF.write(" mode nr.   Freq.[cm^-1]   Displ [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
            for m in range(Nvibmodes):
                if(displ[st+1][m]!=0 and reorg_en_ee[st][m] > 0.15*total_reorg_en_ee[st]): OutF.write("   %02d\t     % .2f \t        % .2e \t     %.2e\t     %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[st+1][m], HR_ee[st][m], reorg_en_ee[st][m]))
            OutF.write("\n\n")

        if args.sig == 'PP' and ESAcheck == 1:
            # ff
            OutF.write("*******************************    f-manifold    *******************************\n\n")
            for st in range(0,esnumst-1-nummodes):
                OutF.write("State nr. %d (Total Reorg. Energy = %.2lf)\n\n" % (st+1,total_reorg_en_ff[st]))
                OutF.write(" mode nr.   Freq.[cm^-1]   Displ [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
                for m in range(Nvibmodes):
                    if(displ[1+nummodes+st][m]!=0): OutF.write("   %02d\t     % .2f \t        % .2e \t     %.2e\t     %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[1+nummodes+st][m], HR_ff[st][m], reorg_en_ff[st][m]))
                OutF.write("\n")

                OutF.write("Main modes (Mode Reorg. En. > 15% of the total) \n\n")
                OutF.write(" mode nr.   Freq.[cm^-1]   Displ [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
                for m in range(Nvibmodes):
                    if(displ[1+nummodes+st][m]!=0 and reorg_en_ff[st][m] > 0.15*total_reorg_en_ff[st]): OutF.write("   %02d\t     % .2f \t        % .2e \t     %.2e\t     %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[1+nummodes+st][m], HR_ff[st][m], reorg_en_ff[st][m]))
                OutF.write("\n\n")


            # ef
            OutF.write("***************************************    ef-manifold    ***************************************\n\n")
            for st in range(0,nummodes):
                for st1 in range(0,esnumst-1-nummodes):
                    if use_HR_ef_p[st][st1]!=0 or use_HR_ef_n[st][st1]!=0:
                        OutF.write("Transition %d_e -> %d_f  (Total Reorg. Energy = %.2lf)\n\n" % (st+1,st1+1,total_reorg_en_ef[st][st1]))
                        OutF.write(" mode nr.   Freq.[cm^-1]    Displ_e    Displ_f [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
                        for m in range(Nvibmodes):
                            if(displ[st+1][m]*displ[st1+1+nummodes][m]!=0):
                                OutF.write("   %02d\t     % .2f\t   % .2e \t      % .2e \t    %.2e\t      %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[st+1][m],displ[st1+1+nummodes][m], HR_ef_p[st][st1][m]+HR_ef_n[st][st1][m], reorg_en_ef[st][st1][m]))
                        OutF.write("\n")

                        OutF.write("Main modes (Mode Reorg. En. > 15% of the total) \n\n")
                        OutF.write(" mode nr.   Freq.[cm^-1]    Displ_e    Displ_f [Bohr*sqrt(amu)]    HR [n.a.]\t  Reorg. En. [cm^-1]\n")
                        for m in range(Nvibmodes):
                            if(displ[st+1][m]*displ[st1+1+nummodes][m]!=0 and reorg_en_ef[st][st1][m] > 0.15*total_reorg_en_ef[st][st1] ): OutF.write("   %02d\t     % .2f\t   % .2e \t      % .2e \t    %.2e\t      %.2e\n" % (m+1+counts_modes_removed_below_threshold,freqs[m],displ[st+1][m],displ[st1+1+nummodes][m], HR_ef_p[st][st1][m]+HR_ef_n[st][st1][m], reorg_en_ef[st][st1][m]))
                        OutF.write("\n")
        OutF.write("*************************************************************************************************\n")


########################################################################################


print("\n")
print(" ... Electric energies printed to file            ----> %s" % OutFile_evals)
print(" ... Electric transition dipoles printed to file  ----> %s" % OutFile_edips)
if SDcheck == 'YES':
    print(" ... Formatted spectral density printed to file   ----> %s" % OutFile_sd)
    print(" ... System-bath coupling printed to file         ----> %s" % OutFile_sbcoup)
    print(" ... Transport rate matrix printed to file        ----> %s" % OutFile_lifet)
    #print(" ... States Static disorder printed to file       ----> %s" % OutFile_stdis)

# ************************************************* #
# Producing the Spectron main input file
# ************************************************* #

OutFile = "%sinput.com" % directory
print(" ... Creating the Spectron input file             ----> %s." % OutFile)

with open(OutFile,'w') as OutF:
    OutF.write("$REGISTRATION\n\n")
    OutF.write("%s\n" % args.sig)
    OutF.write("\n$END\n")

    if args.sig == 'LA':
        OutF.write("\n$SYSTEM\n\n")
        OutF.write("NUMMODES %d\n" % nummodes)
        OutF.write("ES_NUMST %d\n" % (nummodes+1))
        OutF.write("ES_EVALS %s\n" % OutFile_evals_name)
        OutF.write("ES_EDIPS %s\n\n" % OutFile_edips_name)
        if SDcheck == 'YES':
            OutF.write("ES_LAMBDA %s\n" % OutFile_sbcoup_name)
            OutF.write("DISORDER_INTRA_DIAG_GAUSS %d\n" % 0)
            #OutF.write("DISORDER_INTRA_DIAG_GAUSS_F %s\n" % OutFile_stdis_name)

            OutF.write("\nTRANSPORT %d\n" % 1)
            OutF.write("TRANSPORT_RATES %s\n" % OutFile_lifet_name)
            OutF.write("PRINT_POPULATION_M %d\n" % 1)
        else:
            OutF.write("CONST_DEPHASING %d\n" % args.dep)
        OutF.write("\n$END\n")

        OutF.write("\n$%s\n\n" % args.sig)
        if SDcheck == 'YES':
            OutF.write("CAL_METHOD SOS_CGF_C\n\n")
        else:
            OutF.write("CAL_METHOD SOS_CGF\n\n")
        OutF.write("INI_FREQ %.2f\n" % args.w1i)
        OutF.write("FIN_FREQ %.2f\n" % args.w1f)
        OutF.write("NUM_FREQ %d  \n" % args.nw1)
        OutF.write("NUM_SHOTS %d\n" % args.nshots)
        OutF.write("\nOUT_FILE %s \n" % "sig-LA.dat")
        OutF.write("\n$END\n")
            
        if SDcheck == 'YES':
            OutF.write("\n$BATH\n\n")
            OutF.write("OSCILLATORS_NUM %d\n" % NSDs) #Number of SDs
            OutF.write("TEMPERATURE %d\n" % args.tmp)
            OutF.write("BATH_MODEL MM_Continuous_spectral_density\n")
            OutF.write("SPECTRAL_DENSITIES %s\n" % OutFile_sd_name)
            OutF.write("\n$END\n")

    elif args.sig == 'PP':
        OutF.write("\n$SYSTEM\n\n")
        OutF.write("NUMMODES %d\n" % nummodes)
        OutF.write("ES_NUMST %d\n" % esnumst)
        OutF.write("ES_EVALS %s\n" % OutFile_evals_name)
        OutF.write("ES_EDIPS %s\n\n" % OutFile_edips_name)
        if SDcheck == 'YES':
            OutF.write("ES_LAMBDA %s\n" % OutFile_sbcoup_name)
            OutF.write("DISORDER_INTRA_DIAG_GAUSS %d\n" % args.stdis)
            #OutF.write("DISORDER_INTRA_DIAG_GAUSS_F %s\n" % OutFile_stdis_name)

            OutF.write("\nTRANSPORT %d\n" % 1)
            OutF.write("TRANSPORT_RATES %s\n" % OutFile_lifet_name)
            OutF.write("PRINT_POPULATION_M %d\n" % 1)
        else:
            OutF.write("CONST_DEPHASING %d\n" % args.dep)
        OutF.write("\n$END\n")
        
        OutF.write("\n$%s\n\n" % args.sig)
        if SDcheck == 'YES':
            OutF.write("CAL_METHOD SOS_CGF_C\n\n")
        else:
            OutF.write("CAL_METHOD SOS_CGF\n\n")
        OutF.write("INI_FREQ1 %.2f\n" % args.w1i)
        OutF.write("FIN_FREQ1 %.2f\n" % args.w1f)
        OutF.write("NUM_FREQ1 %d\n" % args.nw1)
        OutF.write("INI_FREQ3 %.2f\n" % args.w3i)
        OutF.write("FIN_FREQ3 %.2f\n" % args.w3f)
        OutF.write("NUM_FREQ3 %d\n" % args.nw3)

        OutF.write("\nNUM_SHOTS %d\n" % args.nshots)
        OutF.write("DEL_TIME2 %d\n" % args.t2)

        OutF.write("\nESA_CONTRIBUTIONS %d\n" % ESAcheck)
        OutF.write("GSB_CONTRIBUTIONS %d\n" % GSBcheck)
        OutF.write("SE_CONTRIBUTIONS %d\n" % SEcheck)
        OutF.write("\nOUT_FILE %s \n" % "sig-PP.dat")

        OutF.write("\n$END\n")

        if SDcheck == 'YES':
            OutF.write("\n$BATH\n\n")
            OutF.write("OSCILLATORS_NUM %d\n" % NSDs) #Number of SDs
            OutF.write("TEMPERATURE %d\n" % args.tmp)
            OutF.write("BATH_MODEL MM_Continuous_spectral_density\n")
            OutF.write("SPECTRAL_DENSITIES %s\n" % OutFile_sd_name)
            OutF.write("\n$END\n")


print("\n\n***** PREPARE SPECTRON-INPUT - DONE ***** \n")
