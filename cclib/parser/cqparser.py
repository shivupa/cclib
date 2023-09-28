# -*- coding: utf-8 -*-
#
# Copyright (c) 2023, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for CQ_Dev output files"""

import itertools
import math
import re
import datetime
import numpy

from packaging.version import parse as parse_version, Version

from cclib.parser import logfileparser
from cclib.parser import utils


class CQDev(logfileparser.Logfile):
    """A CQDev log file."""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="CQDev", *args, **kwargs)

    def __str__(self):
        """Return a string representation of the object."""
        return f"CQDev log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'CQDev("{self.filename}")'

    def before_parsing(self):
        pass

    def after_parsing(self):
        super(CQDev, self).after_parsing()

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        # Extract the atomic numbers and coordinates of the atoms.
        if "Geometry:" in line:
            convertor = lambda x: utils.convertor(x, "bohr", "Angstrom")
            self.skip_lines(inputfile, ["dashes", "cols", "blank"])
            atomelements = []
            atomcoords = []
            line = next(inputfile)
            while list(set(line.strip())) != [""]:
                entry = line.split()
                atomelements.append(entry[1])
                atomcoords.append([convertor(float(value)) for value in entry[3:]])
                line = next(inputfile)
            self.append_attribute("atomcoords", atomcoords)

            # We calculate and handle atomnos no matter what, since in
            # the case of fragment calculations the atoms may change,
            # along with the charge and spin multiplicity.
            self.atomnos = []
            self.atomelements = []
            for atomelement in atomelements:
                self.atomelements.append(atomelement.split("-")[0])
                if atomelement == "GH":
                    self.atomnos.append(0)
                else:
                    self.atomnos.append(self.table.number[atomelement])
            self.natom = len(self.atomnos)
            self.atommap = self.generate_atom_map()
            self.formula_histogram = self.generate_formula_histogram()
