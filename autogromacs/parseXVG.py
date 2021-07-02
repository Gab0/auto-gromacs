#!/bin/python

"""

This module parses simple and multi-data .xvg arrays.

"""

import argparse
from typing import Optional, List
from . import mdplots


class XVGParserError(Exception):
    pass


def parse(filepath: str):

    AT = List[List[float]]
    NCOL: Optional[int] = None
    Arrays: List[AT] = []
    Array: AT = []
    with open(filepath) as f:
        for i, line in enumerate(f.readlines()):
            if not line:
                continue
            if line.startswith(('#', '@')):
                continue
            if line.startswith('&'):
                Arrays.append(Array)
                Array = []
                continue
            try:
                row = [
                    float(x)
                    for x in line.strip().split()
                ]
            except Exception as e:
                raise XVGParserError(f"Parser error {e}")

            ncol = len(row)
            if NCOL is None:
                NCOL = ncol
            if NCOL != ncol:
                raise XVGParserError(
                    f"Invalid column length of {ncol}," +
                    f"while {NCOL} was expected at line {i + 1}."
                )
            Array.append(row)

    if Arrays:
        return Arrays

    return [Array]


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", dest="filepaths", nargs='*')
    parser.add_argument("-y", "--ytransform", type=float, default=1.0)
    parser.add_argument("-l", "--label-prefix", dest="label_prefix")
    parser.add_argument("-t", "--totaltime", type=float, default=1.0)
    return parser.parse_args()


def main():

    arguments = parse_arguments()

    for filepath in arguments.filepaths:
        values = parse(filepath)

        YValues = [[v[-1] * arguments.ytransform for v in z] for z in values]
        Labels = [
            f"{arguments.label_prefix}{i + 1}"
            for i, _ in enumerate(YValues)
        ]
        Times = [round(arguments.totaltime) for _, _ in enumerate(YValues)]

        mdplots.show_rms_series_stacked(
            YValues,
            Labels,
            Times,
            filepath + ".png",
            "RMSF"
        )
