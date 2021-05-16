'''
A python implementation of ChemHammer, for use with serialisation over
clusters. This has been compiled into an executable by the numba library.

Copyright (C) 2019  Cameron Hargreaves

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

--------------------------------------------------------------------------------
'''

import numpy as np

import os

from app.network_simplex import network_simplex as ns_jit

from collections import Counter, OrderedDict
import re
import json


def gen_numpy_arrays(comp1):
    # Parse both formulae into their ratios
    if str(comp1) is "((C H3)4 N) (Cu Cd (C N)4)) (C Cl4)":
        comp1 = "((C H3)4 N) (Cu Cd (C N4)) (C Cl4)"

    comp1 = _parse_formula(comp1)
    comp1 = _normalise_composition(comp1)

    comp1_labels = []
    comp1_demands = []

    for k, v in comp1.items():
        comp1_labels.append(_get_position(k))
        comp1_demands.append(v)

    # Turn these lists into numpy arrays and sort these
    source_labels = np.array(comp1_labels, dtype=np.int64)
    source_ord = np.argsort(source_labels)
    source_labels = source_labels[source_ord]
    source_demands = np.array(comp1_demands, dtype=np.float64)[source_ord]

    return source_labels, source_demands


def min_flow_dist(comp1_labels, comp1_demands, comp2_labels, comp2_demands):
    '''
    For two given atom labels and associated ratios calculate the min flow 
    distance between these
    '''
    source_labels = np.array(comp1_labels, dtype=np.int64)
    source_demands = np.array(comp1_demands, dtype=np.float64)
    sink_labels = np.array(comp2_labels, dtype=np.int64)
    sink_demands = np.array(comp2_demands, dtype=np.float64)
    
    flow_cost = ns_jit(source_labels, source_demands, sink_labels, sink_demands)

    return flow_cost


def min_flow_dist2(comp1, comp2):
    comp1_labels, comp1_demands = gen_numpy_arrays(comp1)
    comp2_labels, comp2_demands = gen_numpy_arrays(comp2)

    flow_cost = ns_jit(comp1_labels, comp1_demands, comp2_labels, comp2_demands)
    return flow_cost


def _get_position(element):
        """
        Return either the x, y coordinate of an elements position, or the
        x-coordinate on the Pettifor numbering system as a 2-dimensional
        """
        atomic_num = _get_atomic_num(element)
        try:
            atom_info = _get_periodic_tab()['elements'][atomic_num]
        except:
            atom_info = {'mod_petti_num' : 1} #If there's an unknown element make equal to H

        return atom_info['mod_petti_num']


def _get_atomic_num(element_string):
    """ Return atomic number from element string """
    for i, element in enumerate(_get_periodic_tab()['elements']):
        if element['symbol'] == element_string:
            return i


def _get_periodic_tab():
    """
    Just load from file
    """
    with open('./ElementData.json') as json_data:
            periodic_data = json.load(json_data)

    return periodic_data


def _normalise_composition(composition):
    """ Sum up the numbers in our counter to get total atom count """
    normed_comp = {}
    atom_count =  sum(composition.values(), 0.0)

    for atom in composition.keys():
        normed_comp[atom] = composition[atom] / atom_count

    return normed_comp


def _is_balanced(formula):
    """Check if all sort of brackets come in pairs."""
    # Very naive check, just here because you always need some input checking
    c = Counter(formula)
    return c['['] == c[']'] and c['{'] == c['}'] and c['('] == c[')']


def _dictify(tuples):
    """Transform tuples of tuples to a dict of atoms."""
    res = dict()
    for atom, n in tuples:
        try:
            res[atom] += float(n or 1)
        except KeyError:
            res[atom] = float(n or 1)
    return res


def _fuse(mol1, mol2, w=1):
    """ Fuse 2 dicts representing molecules. Return a new dict. """
    return {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w for atom in set(mol1) | set(mol2)}


def _parse(formula):
    """
    Return the molecule dict and length of parsed part.
    Recurse on opening brackets to parse the subpart and
    return on closing ones because it is the end of said subpart.
    """
    ATOM_REGEX = '([A-Z][a-z]*)(\d*\.*\d*)'
    OPENERS = '({['
    CLOSERS = ')}]'

    q = []
    mol = {}
    i = 0

    while i < len(formula):
        # Using a classic loop allow for manipulating the cursor
        token = formula[i]

        if token in CLOSERS:
            # Check for an index for this part
            m = re.match('\d+\.*\d*', formula[i+1:])
            if m:
                weight = float(m.group(0))
                i += len(m.group(0))
            else:
                weight = 1

            submol = _dictify(re.findall(ATOM_REGEX, ''.join(q)))
            return _fuse(mol, submol, weight), i

        elif token in OPENERS:
            submol, l = _parse(formula[i+1:])
            mol = _fuse(mol, submol)
            # skip the already read submol
            i += l + 1
        else:
            q.append(token)

        i += 1

    # Fuse in all that's left at base level
    return _fuse(mol, _dictify(re.findall(ATOM_REGEX, ''.join(q)))), i


def _parse_formula(formula):
    """Parse the formula and return a dict with occurences of each atom."""
    if not _is_balanced(formula):
        raise ValueError("Your brackets not matching in pairs ![\{]$[&?)]\}!]")

    return _parse(formula)[0]


if __name__ == "__main__":
    import time

    time_start = time.time()

    a = "(Li1.104 Mn0.896) Mn O4"
    b = "Li1.3Al0.3Ti1.7(PO4)3"

    # source_labels, source_demands, sink_labels, sink_demands = gen_numpy_arrays(a, b)
    # This is meant to be 24.315
    print(f'network simplex {min_flow_dist(gen_numpy_arrays(a), gen_numpy_arrays(b))}')
    print(time.time() - time_start)

    time_start = time.time()

    for i in range(1, 101):
        min_flow_dist(gen_numpy_arrays(f"Li Ba{i} B3 S6"), gen_numpy_arrays("Li1.3Al0.3Ti1.7(PO4)3"))

    print(time.time() - time_start)
