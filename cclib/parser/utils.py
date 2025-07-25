# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Utilities often used by cclib parsers and scripts"""

import re
from itertools import accumulate
from math import sqrt
from typing import List, Sequence

import numpy
import periodictable


def find_package(package: str) -> bool:
    """Check if a package exists without importing it.

    Derived from https://stackoverflow.com/a/14050282
    """
    from importlib.util import find_spec

    module_spec = find_spec(package)
    return module_spec is not None and module_spec.loader is not None


_found_scipy = find_package("scipy")
if _found_scipy:
    import scipy.spatial


def symmetrize(m: numpy.ndarray, use_triangle: str = "lower") -> numpy.ndarray:
    """Symmetrize a square NumPy array by reflecting one triangular
    section across the diagonal to the other.
    """

    if use_triangle not in ("lower", "upper"):
        raise ValueError
    if not len(m.shape) == 2:
        raise ValueError
    if not (m.shape[0] == m.shape[1]):
        raise ValueError

    dim = m.shape[0]

    ms = m.copy()

    if use_triangle == "lower":
        lower_indices = numpy.tril_indices(dim, k=-1)
        upper_indices = (lower_indices[1], lower_indices[0])
        ms[upper_indices] = ms[lower_indices]
    if use_triangle == "upper":
        upper_indices = numpy.triu_indices(dim, k=1)
        lower_indices = (upper_indices[1], upper_indices[0])
        ms[lower_indices] = ms[upper_indices]

    return ms


def convertor(value: float, fromunits: str, tounits: str) -> float:
    """Convert from one set of units to another.

    Sources:
        NIST 2010 CODATA (http://physics.nist.gov/cuu/Constants/index.html)
        Documentation of GAMESS-US or other programs as noted
    """

    _convertor = {
        "time_au_to_fs": lambda x: x * 0.02418884,
        "fs_to_time_au": lambda x: x / 0.02418884,
        "Angstrom_to_bohr": lambda x: x * 1.8897261245,
        "bohr_to_Angstrom": lambda x: x * 0.5291772109,
        "wavenumber_to_eV": lambda x: x / 8065.54429,
        "wavenumber_to_hartree": lambda x: x / 219474.6313708,
        "wavenumber_to_kcal/mol": lambda x: x / 349.7550112,
        "wavenumber_to_kJ/mol": lambda x: x / 83.5934722814,
        "wavenumber_to_nm": lambda x: 1e7 / x,
        "wavenumber_to_Hz": lambda x: x * 29.9792458,
        "eV_to_wavenumber": lambda x: x * 8065.54429,
        "eV_to_hartree": lambda x: x / 27.21138505,
        "eV_to_kcal/mol": lambda x: x * 23.060548867,
        "eV_to_kJ/mol": lambda x: x * 96.4853364596,
        "hartree_to_wavenumber": lambda x: x * 219474.6313708,
        "hartree_to_eV": lambda x: x * 27.21138505,
        "hartree_to_kcal/mol": lambda x: x * 627.50947414,
        "hartree_to_kJ/mol": lambda x: x * 2625.4996398,
        "kcal/mol_to_wavenumber": lambda x: x * 349.7550112,
        "kcal/mol_to_eV": lambda x: x / 23.060548867,
        "kcal/mol_to_hartree": lambda x: x / 627.50947414,
        "kcal/mol_to_kJ/mol": lambda x: x * 4.184,
        "kJ/mol_to_wavenumber": lambda x: x * 83.5934722814,
        "kJ/mol_to_eV": lambda x: x / 96.4853364596,
        "kJ/mol_to_hartree": lambda x: x / 2625.49963978,
        "kJ/mol_to_kcal/mol": lambda x: x / 4.184,
        "nm_to_wavenumber": lambda x: 1e7 / x,
        # Taken from GAMESS docs, "Further information",
        # "Molecular Properties and Conversion Factors"
        "Debye^2/amu-Angstrom^2_to_km/mol": lambda x: x * 42.255,
        # Conversion for charges and multipole moments.
        "e_to_coulomb": lambda x: x * 1.602176565 * 1e-19,
        "e_to_statcoulomb": lambda x: x * 4.80320425 * 1e-10,
        "coulomb_to_e": lambda x: x * 0.6241509343 * 1e19,
        "statcoulomb_to_e": lambda x: x * 0.2081943527 * 1e10,
        "ebohr_to_Debye": lambda x: x * 2.5417462300,
        "ebohr2_to_Buckingham": lambda x: x * 1.3450341749,
        "ebohr2_to_Debye.ang": lambda x: x * 1.3450341749,
        "ebohr3_to_Debye.ang2": lambda x: x * 0.7117614302,
        "ebohr4_to_Debye.ang3": lambda x: x * 0.3766479268,
        "ebohr5_to_Debye.ang4": lambda x: x * 0.1993134985,
        # FIXME I'm not sure this is correct.
        "hartree/bohr2_to_mDyne/angstrom": lambda x: x * 8.23872350 / 0.5291772109,
    }

    return _convertor[f"{fromunits}_to_{tounits}"](value)


def _get_rmat_from_vecs(a, b):
    """Get rotation matrix from two 3D vectors, a and b
    Args:
       a (numpy.ndarray): 3d vector with shape (3,0)
       b (numpy.ndarray): 3d vector with shape (3,0)
    Returns:
       numpy.ndarray
    """
    a_ = a / numpy.linalg.norm(a, 2)
    b_ = b / numpy.linalg.norm(b, 2)
    v = numpy.cross(a_, b_)
    s = numpy.linalg.norm(v, 2)
    c = numpy.dot(a_, b_)
    # skew-symmetric cross product of v
    vx = numpy.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rmat = numpy.identity(3) + vx + numpy.matmul(vx, vx) * ((1 - c) / s**2)
    return rmat


def get_rotation(a, b):
    """Get rotation part for transforming a to b, where a and b are same positions with different orientations
    If one atom positions, i.e (1,3) shape array, are given, it returns identify transformation

    Args:
        a (numpy.ndarray): positions with shape(N,3)
        b (numpy.ndarray): positions with shape(N,3)
    Returns:
        A scipy.spatial.transform.Rotation object
    """
    if not _found_scipy:
        raise ImportError("You must install `scipy` to use this function")

    assert a.shape == b.shape
    if a.shape[0] == 1:
        return scipy.spatial.transform.Rotation.from_euler("xyz", [0, 0, 0])
    # remove translation part
    a_ = a - a[0]
    b_ = b - b[0]
    if hasattr(scipy.spatial.transform.Rotation, "align_vectors"):
        r, _ = scipy.spatial.transform.Rotation.align_vectors(b_, a_)
    else:
        if numpy.linalg.matrix_rank(a_) == 1:
            # in the case of linear molecule, e.g. O2, C2H2
            idx = numpy.argmax(numpy.linalg.norm(a_, ord=2, axis=1))
            rmat = _get_rmat_from_vecs(a_[idx], b_[idx])
            r = scipy.spatial.transform.Rotation.from_dcm(rmat)
        else:
            # scipy.spatial.transform.Rotation.match_vectors has bug
            # Kabsch Algorithm
            cov = numpy.dot(b_.T, a_)
            V, S, W = numpy.linalg.svd(cov)
            if (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0:
                S[-1] = -S[-1]
                V[:, -1] = -V[:, -1]
            rmat = numpy.dot(V, W)
            r = scipy.spatial.transform.Rotation.from_dcm(rmat)
    return r


def skip_until_no_match(inputfile, regex):
    """Skip lines that match a regex. First non-matching line is returned.

    This method allows to skip a variable number of lines, allowing for example,
    to parse sections that might have different whitespace/spurious lines for
    different versions of the software.
    """
    line = next(inputfile)
    while re.match(regex, line):
        line = next(inputfile)
    return line


def str_contains_only(string, chars):
    """Checks if string contains only the specified characters."""
    return all([c in chars for c in string])


class PeriodicTable:
    """Allows conversion between element name and atomic no."""

    def __init__(self) -> None:
        self.element = [None]
        self.number = {}

        for e in periodictable.elements:
            if e.symbol != "n":
                self.element.append(e.symbol)
                self.number[e.symbol] = e.number

        # Add common placeholder atoms.  These are not ghost atoms, which
        # still have basis functions associated with a parent element.
        self.number["-"] = 0


class WidthSplitter:
    """Split a line based not on a character, but a given number of field
    widths.
    """

    def __init__(self, widths: Sequence[int]) -> None:
        self.start_indices = [0] + list(accumulate(widths))[:-1]
        self.end_indices = list(accumulate(widths))

    def split(self, line: str, truncate: bool = True) -> List[str]:
        """Split the given line using the field widths passed in on class
        initialization.
        """
        elements = [
            line[start:end].strip() for (start, end) in zip(self.start_indices, self.end_indices)
        ]
        # Handle lines that contain fewer fields than specified in the
        # widths; they are added as empty strings, so remove them.
        if truncate:
            while len(elements) and elements[-1] == "":
                elements.pop()
        return elements


def _dim_from_tblock_size(x: int) -> int:
    """Given the number of elements in a symmetric matrix lower triangle,
    including the diagonal, compute the dimension (length of one side) of the
    full matrix.
    """
    r1 = 0.5 * (sqrt(8 * x + 1) - 1)
    r2 = 0.5 * (-sqrt(8 * x + 1) - 1)
    m = max(r1, r2)
    if _BUILTIN_FLOAT(round(m)) != m:
        raise RuntimeError(f"The number of elements ({x}) isn't possible for a matrix triangle")
    return int(m)


def block_to_matrix(block: numpy.ndarray) -> numpy.ndarray:
    """Convert a flattened symmetric matrix lower triangle to its full
    symmetrized form.
    """
    assert len(block.shape) == 1
    dim = _dim_from_tblock_size(block.shape[0])
    mat = numpy.empty(shape=(dim, dim), dtype=block.dtype)
    mat[numpy.tril_indices_from(mat)] = block
    return symmetrize(mat)


_BUILTIN_FLOAT = float


# This is at the bottom of the file so it doesn't interfere with type hints in
# any of the file's other functions.
def float(number: str) -> float:
    """Convert a string to a float.

    This method should perform certain checks that are specific to cclib,
    including avoiding the problem with Ds instead of Es in scientific notation.
    Another point is converting string signifying numerical problems (*****)
    to something we can manage (Numpy's NaN).
    """

    if list(set(number)) == ["*"]:
        return numpy.nan

    return _BUILTIN_FLOAT(number.replace("D", "E"))
