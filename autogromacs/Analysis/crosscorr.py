
import os
import correlationplus
import prody
from correlationplus.calculate import *
from correlationplus.visualize import *
from collections import Counter

import numpy as np


def atom_group_to_prody(atom_group: mda.AtomGroup) -> prody.AtomGroup:
    filename = "ATOM_CNV$!.pdb"
    # No in-memory buffer is accepted here, so we write a file x.x
    atom_group.write(filename)
    output_group = prody.parsePDB(filename)
    os.remove(filename)
    return(output_group)


def trajectory_cross_correlation(universe: mda.Universe, output_file):

    atom_group_p = atom_group_to_prody(universe.atoms)

    crosscorr_matrix = trajectory_cross_correlation_calculate(universe)
    trajectory_cross_correlation_visualize(atom_group_p, output_file)

    return crosscorr_matrix

def trajectory_cross_correlation_calculate(universe: mda.Universe):
    """ atom_group should be CA only!"""

    method = "GNM"
    cut_off = 10
    num_mod = 100

    out_file = "cross.dat"

    calcMDnDCC_(
        universe,
        normalized=True,
        saveMatrix=True,
        out_file=out_file
    )


# Ripped directly from https://github.com/tekpinar/correlationplus/blob/master/correlationplus/calculate.py
def calcMDnDCC_(
        universe, startingFrame=0, endingFrame=(-1),
               normalized=True, alignTrajectory=True,
               saveMatrix=True, out_file="DCC"):
    """
        Calculate normalized dynamical cross-correlations when a topology
        and a trajectory file is given.
    Parameters
    ----------
    topology: string
        A PDB file.
    trajectory: string
        A trajectory file in dcd, xtc or trr format.
    startingFrame: int
        You can specify this value if you want to exclude some initial frames
        from your cross-correlation calculations. Default value is 0.
    endingFrame: int
        You can specify this value if you want to calculate cross-correlation
        calculations up to a certain ending frame. Default value is -1 and it
        indicates the last frame in your trajectory.
    normalized: bool
        Default value is True and it means that the cross-correlation matrix
        will be normalized.
    alignTrajectory: bool
        Default value is True and it means that all frames in the trajectory
        will be aligned to the initial frame.
    saveMatrix: bool
        If True, cross-correlation matrix will be written to an output file.
    out_file: string
        Output file name for the cross-correlation matrix.
        Default value is DCC and the file extension is .dat.
    Returns
    -------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.
    """
    # Create the universe (That sounds really fancy :)

    # Create an atomgroup from the alpha carbon selection
    calphas = universe.select_atoms("protein and name CA")
    N = calphas.n_atoms
    print(f"@> Parsed {N} Calpha atoms.")
    # Set your frame window for your trajectory that you want to analyze
    # startingFrame = 0
    if endingFrame == -1:
        endingFrame = universe.trajectory.n_frames
    skip = 1

    # Perform Calpha alignment first
    if alignTrajectory:
        print("@> Aligning only Calpha atoms to the initial frame!")
        alignment = align.AlignTraj(universe, universe, select="protein and name CA", in_memory=True)
        alignment.run()

    Rvector = []
    # Iterate through the universe trajectory
    for timestep in universe.trajectory[startingFrame:endingFrame:skip]:
        Rvector.append(calphas.positions.flatten())
    ##############################################

    # I reassign this bc in ccMatrix calculation, we may skip some part of the trajectory!
    N_Frames = len(Rvector)

    R_average = np.mean(Rvector, axis=0)
    print("@> Calculating cross-correlation matrix:")
    ccMatrix = DCCmatrixCalculation(N, np.array(Rvector), R_average)

    # Do the averaging
    ccMatrix = ccMatrix / float(N_Frames)

    if normalized:
        cc_normalized = np.zeros((N, N), np.double)
        for i in range(0, N):
            for j in range (i, N):
                cc_normalized[i][j] = ccMatrix[i][j] / ((ccMatrix[i][i] * ccMatrix[j][j]) ** 0.5)
                cc_normalized[j][i] = cc_normalized[i][j]

        if saveMatrix:
            np.savetxt(out_file, cc_normalized, fmt='%.6f')

        return cc_normalized
    else:
        for i in range(0, N):
            for j in range(i + 1, N):
                ccMatrix[j][i] = ccMatrix[i][j]
        if saveMatrix:
            np.savetxt(out_file, ccMatrix, fmt='%.6f')
        return ccMatrix


def trajectory_cross_correlation_visualize(atom_group, output_file):

    visualizemapApp(
        inp_file="cross.dat",
        out_file=output_file,
        sel_type="ndcc",
        atom_group=atom_group,
        vmin_fltr=0.75,
        vmax_fltr=None,
        dmin_fltr=0.0,
        dmax_fltr=9999.9,
        cyl_rad=None
    )




# Ripped off directly from
# https://github.com/tekpinar/correlationplus/blob/master/correlationplus/scripts/visualize.py
# This function had to be modified in order to fit with our system.
def visualizemapApp(
        inp_file, out_file, sel_type, atom_group: prody.AtomGroup,
        vmin_fltr, vmax_fltr,
        dmin_fltr, dmax_fltr, cyl_rad):

    ##########################################################################
    selectedAtoms = atom_group.select('protein and name CA')
    ##########################################################################
    minColorBarLimit = 0.0
    maxColorBarLimit = 1.0
    # Read data file and assign to a numpy array
    if sel_type.lower() == "ndcc":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True,
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = np.loadtxt(inp_file, dtype=float)
        # Check the data range in the matrix.
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        if minCorrelationValue < 0.0:
        # Assume that it is an nDCC file
            minColorBarLimit = -1.0
        if maxCorrelationValue > 1.00001:
            print("This correlation map is not normalized!")
            # TODO: At this point, one can ask the user if s/he wants to normalize it!
            sys.exit(-1)
        else:
            maxColorBarLimit = 1.0
    elif sel_type.lower() == "dcc":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True,
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = np.loadtxt(inp_file, dtype=float)
        # Check the data range in the matrix.
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    elif sel_type.lower() == "absndcc":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = np.absolute(parseSparseCorrData(inp_file, selectedAtoms, \
                                                        Ctype=True,
                                                        symmetric=True,
                                                        writeAllOutput=False))
        else:
            ccMatrix = np.absolute(np.loadtxt(inp_file, dtype=float))
        minColorBarLimit = 0.0
        maxColorBarLimit = 1.0
    elif sel_type.lower() == "lmi":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True,
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=False)
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
        #minColorBarLimit = 0.0

    elif sel_type.lower() == "nlmi":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True,
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = convertLMIdata2Matrix(inp_file, writeAllOutput=False)
        #minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = 0.0

        # Ideally, it is supposed to be 1 but I used 1.00001 to avoid
        # rounding problems
        if maxCorrelationValue > 1.00001:
            print("This LMI map is not normalized!")
            # TODO: At this point, one can ask the user if s/he wants to normalize it!
            sys.exit(-1)
        else:
            maxColorBarLimit = 1.0

    elif sel_type.lower() == "coeviz":
        ccMatrix = np.loadtxt(inp_file, dtype=float)
        minColorBarLimit = 0.0
        maxColorBarLimit = 1.0

    elif sel_type.lower() == "evcouplings":
        ccMatrix = parseEVcouplingsScores(inp_file, selectedAtoms, False)
        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    elif sel_type.lower() == "generic":
        # Check if the data type is sparse matrix
        data_file = open(inp_file, 'r')
        allLines = data_file.readlines()
        data_file.close()

        # Read the first line to determine if the matrix is sparse format
        words = allLines[0].split()

        # Read the 1st line and check if it has three columns
        if (len(words) == 3):
            ccMatrix = parseSparseCorrData(inp_file, selectedAtoms, \
                                            Ctype=True,
                                            symmetric=True,
                                            writeAllOutput=False)
        else:
            ccMatrix = np.loadtxt(inp_file, dtype=float)

        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    elif sel_type.lower() == "eg":
        # The data type is elasticity graph
        ccMatrix = parseElasticityGraph(inp_file, selectedAtoms, \
                                            writeAllOutput=False)

        minCorrelationValue = np.min(ccMatrix)
        maxCorrelationValue = np.max(ccMatrix)
        minColorBarLimit = minCorrelationValue
        maxColorBarLimit = maxCorrelationValue
    else:
        print("@> ERROR: Unknown data type: Type can only be ndcc, absndcc, lmi,\n")
        print("@>        coeviz or evcouplings. If you have your data in full \n")
        print("@>        matrix format and your data type is none of the options\n")
        print("@>        mentionned, you can set data type 'generic'.\n")
        sys.exit(-1)

    # Set vmin_fltr and vmax_fltr
    if (vmin_fltr == None):
        vmin_fltr = minColorBarLimit
    if (vmax_fltr == None):
        vmax_fltr = maxColorBarLimit

    print(f"""@> Min. value filter: {vmin_fltr}""")
    print(f"""@> Max. value filter: {vmax_fltr}""")

    ##########################################################################
    # Call overall correlation calculation

    overallCorrelationMap(ccMatrix, minColorBarLimit, maxColorBarLimit,
                          out_file, " ", selectedAtoms)

    plotDistributions = True
    VMDcylinderRadiusScale = 0.5
    PMLcylinderRadiusScale = 0.3

    if (cyl_rad == None):
        if sel_type.lower() == "evcouplings":
            VMDcylinderRadiusScale = 0.02
            PMLcylinderRadiusScale = 0.02
        else:
            VMDcylinderRadiusScale = 0.5
            PMLcylinderRadiusScale = 0.3
        print(f"""@> VMD Cylinder radius: {VMDcylinderRadiusScale}""")
        print(f"""@> PyMol Cylinder radius: {PMLcylinderRadiusScale}""")

    else:
        VMDcylinderRadiusScale = float(cyl_rad)
        PMLcylinderRadiusScale = float(cyl_rad)
        print(f"""@> Cylinder radius: {cyl_rad}""")

    if plotDistributions:
        if sel_type.lower() == "ndcc":
            distanceDistribution(ccMatrix, out_file, "nDCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        elif sel_type.lower() == "dcc":
            distanceDistribution(ccMatrix, out_file, "DCC", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif sel_type.lower() == "absndcc":
            distanceDistribution(ccMatrix, out_file, "Abs(nDCC)",
                                 selectedAtoms, absoluteValues=True, writeAllOutput=False)

        elif sel_type.lower() == "lmi":
            distanceDistribution(ccMatrix, out_file, "LMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)
        elif sel_type.lower() == "nlmi":
            distanceDistribution(ccMatrix, out_file, "nLMI", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)
        elif sel_type.lower() == "coeviz":
            distanceDistribution(ccMatrix, out_file, "CoeViz", selectedAtoms,
                                 absoluteValues=True, writeAllOutput=True)

        elif sel_type.lower() == "evcouplings":
            distanceDistribution(ccMatrix, out_file, "EVcoupling Score", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)

        elif sel_type.lower() == "generic":
            distanceDistribution(ccMatrix, out_file, "Correlation", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        elif sel_type.lower() == "eg":
            distanceDistribution(ccMatrix, out_file, "Force Constants", selectedAtoms,
                                 absoluteValues=False, writeAllOutput=True)
        else:
            print("Warning: Unknows correlation data.\n")
            print("         Correlations can be dcc, ndcc, absndcc, lmi,\n")
            print("         nlmi, coeviz or evcouplings!\n")

    ##########################################################################
    # Check number of chains. If there are multiple chains, plot inter and
    # intra chain correlations
    chains = Counter(selectedAtoms.getChids()).keys()
    saveMatrix = False
    plotChains = True
    if len(chains) > 1 and plotChains:
        intraChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)
        interChainCorrelationMaps(ccMatrix, minColorBarLimit, maxColorBarLimit,
                                  out_file, " ", selectedAtoms, saveMatrix)

    # Here, we can filter some correlation values closer than a distance.
    # Typically, it is supposed to filter out the correlation within the
    # same secondary structure etc.
    filterByDistance = True
    if filterByDistance:
        disMinValue = float(dmin_fltr)
        disMaxValue = float(dmax_fltr)
        ccMatrix = filterCorrelationMapByDistance(ccMatrix, out_file, " ",
                                                  selectedAtoms,
                                                  disMinValue, disMaxValue,
                                                  absoluteValues=False,
                                                  writeAllOutput=False)

    # Overall projection
    projectCorrelationsOntoProteinVMD(
        out_file,
        ccMatrix,
        out_file,
        selectedAtoms,
        vminFilter=float(vmin_fltr),
        vmaxFilter=float(vmax_fltr),
        cylinderRadiusScaler=VMDcylinderRadiusScale,
        absoluteValues=True,
        writeAllOutput=True
    )

    projectCorrelationsOntoProteinPyMol(
        out_file,
        ccMatrix,
        out_file,
        selectedAtoms,
        vminFilter=float(vmin_fltr),
        vmaxFilter=float(vmax_fltr),
        cylinderRadiusScaler=PMLcylinderRadiusScale,
        absoluteValues=True, writeAllOutput=True
    )
