
import argparse

import MDAnalysis
from MDAnalysis.analysis import diffusionmap, align, rms

import numpy.linalg

import matplotlib.pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", dest='Main', nargs="*")


    return parser.parse_args()


def analyzeMD(arguments):


    print(arguments.Main)
    us = [
        MDAnalysis.Universe(FP + ".gro", FP + ".trr")
        for FP in arguments.Main
    ]
    # can access via segid (4AKE) and atom name
    # we take the first atom named N and the last atom named C
    #nterm = u.select_atoms('segid 4AKE and name N')[0]
    #cterm = u.select_atoms('segid 4AKE and name C')[-1]

    bb = [u.select_atoms('protein and backbone') for u in us] # a selection (AtomGroup)

    aligner = align.AlignTraj(us[0], us[1], select='name CA',
                              in_memory=True).run()

    matrix = diffusionmap.DistanceMatrix(us[0], select='name CA').run()

    plt.imshow(matrix.dist_matrix, cmap='viridis')
    plt.xlabel('Frame')
    plt.ylabel('Frame')
    plt.colorbar(label=r'RMSD ($\AA$)')

    plt.show()

    for i, ts in enumerate(us[0].trajectory):     # iterate through all frames
        #r = cterm.position - nterm.position # end-to-end vector from atom positions
        print(dir(ts))
        print(ts.positions)
        print(i)
    #r = bb.position
    #d = numpy.linalg.norm(r)  # end-to-end distance
    #rgyr = bb.radius_of_gyration()  # method of AtomGroup
    #print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
    #      ts.frame, d, rgyr))


def main():
    arguments = parse_arguments()
    analyzeMD(arguments)


if __name__ == "__main__":
    main()

