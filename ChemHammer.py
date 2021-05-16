"""
A class to compute the atomic similarity between two structural compounds based
on the minimal cost flow between these.

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

Python Parser Source: https://github.com/Zapaan/python-chemical-formula-parser

Periodic table JSON data: https://github.com/Bowserinator/Periodic-Table-JSON,
updated to include the Pettifor number and modified Pettifor number from
https://iopscience.iop.org/article/10.1088/1367-2630/18/9/093011

"""

import re
import json
import os
import sys
import warnings

import urllib.request

from copy import deepcopy
from collections import Counter, OrderedDict
from math import sqrt
from functools import partial

import numpy as np

from scipy.spatial.distance import pdist, squareform
from scipy.optimize import linear_sum_assignment

from ortools.graph import pywrapgraph


# from MinFlowFinder import MinFlowFinder

def main():
    x = ChemHammer("H")
    y = ChemHammer("KF")
    z = ChemHammer("He")

    print(f"Verbose Minimum Flow between {x.formula} and {y.formula} is: {x.min_flow_dist(y)}\n")
    print(f"Verbose Minimum Flow between {y.formula} and {z.formula}  is: {y.min_flow_dist(z)}\n")
    print(f"Verbose Minimum Flow between {x.formula} and {z.formula}  is: {x.min_flow_dist(z)}\n")

class ChemHammer():
    ATOM_REGEX = '([A-Z][a-z]*)(\d*\.*\d*)'
    OPENERS = '({['
    CLOSERS = ')}]'

    # As the current optimization solver only takes in ints we must multiply
    # all floats to capture the decimal places
    FP_MULTIPLIER = 100000000

    # How much the change in element position affects the distance metric
    # Advsied to use 0.1 for pettifor numbers as they have much larger variance
    DIST_MOD = 2

    # How much the addition of a new element affects the distance metric
    LEVENSH_MOD = 0.1

    def __init__(self, formula, metric="mod_petti"):
        self.formula = formula.replace(" ", "")
        self.periodic_tab = self._get_periodic_tab()
        self.composition = self._parse_formula(self.formula)
        self.normed_composition = self._normalise_composition(self.composition)
        self.distance_metric = metric

    def min_flow_dist(self, comp2, comp1 = None, verbose=False):
        comp1 = comp1 if comp1 is not None else self.normed_composition

        if isinstance(comp2, str):
            comp2 = self._parse_formula(comp2)
            comp2 = self._normalise_composition(comp2)

        if isinstance(comp2, ChemHammer):
            comp2 = comp2.normed_composition

        start_nodes, end_nodes, labels, capacities, costs, supplies = self._generate_parameters(comp1, comp2)

        # Google OR-tools only take integer values, so we multiply our floats
        # by self.FP_MULTIPLIER and cast to int
        capacities = [int(x * self.FP_MULTIPLIER) for x in capacities]
        supplies = [int(x * self.FP_MULTIPLIER) for x in supplies]

        # Due to rounding errors, the two supplies may no longer be equal to one
        # another. We add the difference to the largest value in the smaller set
        # to allow this to process and keep the error minimised
        source_tot = sum([x for x in supplies if x > 0])
        sink_tot = -sum([x for x in supplies if x < 0])

        while sink_tot < source_tot:
            supplies[supplies.index(min(supplies))] -= 1
            sink_tot = -sum([x for x in supplies if x < 0])

        while source_tot < sink_tot:
            supplies[supplies.index(max(supplies))] += 1
            source_tot = sum([x for x in supplies if x > 0])

        while sum(capacities) < source_tot or sum(capacities) < sink_tot:
            capacities[capacities.index(max(capacities))] += 1

        # Instantiate a SimpleMinCostFlow solver
        min_cost_flow = pywrapgraph.SimpleMinCostFlow()

        # Add each arc.
        for i in range(0, len(start_nodes)):
            min_cost_flow.AddArcWithCapacityAndUnitCost(start_nodes[i], end_nodes[i],
                                                        capacities[i], costs[i])

        # Add node supplies.
        for i in range(len(supplies)):
            min_cost_flow.SetNodeSupply(i, supplies[i])

        feasibility_status = min_cost_flow.Solve()

        if feasibility_status == min_cost_flow.OPTIMAL:
            dist = min_cost_flow.OptimalCost() / self.FP_MULTIPLIER

            if verbose:
                print(f"Arc{' ' * 9}Flow{' ' * 13}Capacity{' ' * 8}Pettifor Dist.   Cost")
                for i in range(min_cost_flow.NumArcs()):
                    cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i)
                    print('%-2s -> %-2s    %-15s  %-15s %-16s %-5s' % (
                        labels[min_cost_flow.Tail(i)].split('_')[0],
                        labels[min_cost_flow.Head(i)].split('_')[0],
                        min_cost_flow.Flow(i) / self.FP_MULTIPLIER,
                        min_cost_flow.Capacity(i) /self.FP_MULTIPLIER,
                        costs[i],
                        cost / self.FP_MULTIPLIER ))
                print(f"Total Cost: {min_cost_flow.OptimalCost() / self.FP_MULTIPLIER}\n")
            return dist

        else:
            return "Infeasible solution"

    def _generate_parameters(self, source, sink):
        start_nodes = []
        start_labels = []
        end_nodes = []
        end_labels = []

        capacities = []
        costs = []
        supply_tracker = OrderedDict()

        for i, key_value_source in enumerate(source.items()):
            for j, key_value_sink in enumerate(sink.items()):
                start_nodes.append(i)
                start_labels.append(key_value_source[0])
                end_nodes.append(j + len(source))
                end_labels.append(key_value_sink[0])
                capacities.append(min(key_value_source[1], key_value_sink[1]))
                costs.append(abs(self._get_position(key_value_source[0]) - self._get_position(key_value_sink[0])))

        for lab in start_labels:
            supply_tracker[lab + "_source"] = source[lab]

        for lab in end_labels:
            supply_tracker[lab + "_sink"] = -sink[lab]

        labels = list(supply_tracker.keys())
        supplies = list(supply_tracker.values())

        return start_nodes, end_nodes, labels, capacities, costs, supplies

    def euclidean_dist(self, comp2, comp1 = None):
        """
        Simply take the euclidean distance between two vectors excluding atom
        similarity. Here we take the normalised vector and assume that the input
        is a string formula
        """
        comp1 = comp1 if comp1 is not None else self.normed_composition

        if isinstance(comp2, str):
            comp2 = self._parse_formula(comp2)
            comp2 = self._normalise_composition(comp2)

        dist = 0

        # Loop through the union of  keys
        for key in list(comp1.keys()) + list(comp2.keys()):
            # Simply take the distance between if element is present in both
            if key in comp1 and key in comp2:
                dist += (comp1[key] - comp2[key]) ** 2

            # If its' not a shared element try and take distance from each dict
            elif key in comp1:
                dist += comp1[key] ** 2

            elif key in comp2:
                dist += comp2[key] ** 2

            else:
                print("Key not in either, strange bug occurred")

        return sqrt(dist)

    def hamming_dist(self, comp2, comp1=None, flag=True):
        """
        DEPRECATED, use min_flow_dist()

        This is similar to euclidean distance, however adds a further distance
        metric depending on manhattan distance between closest neighbours
        """
        comp1 = comp1 if comp1 is not None else self.normed_composition

        if isinstance(comp2, str):
            comp2 = self._parse_formula(comp2)
            comp2 = self._normalise_composition(comp2)

        if isinstance(comp2, ChemHammer):
            comp2 = comp2.composition

        if len(comp1) != len(comp2) and flag:
            warnings.warn("Must have equal numbers of unique elements for Hamming distance, use levenshtein distance")

        # Pair up each of the atoms so that the overall distances are minimised
        comp1_pos = self._return_positions(comp1)
        comp2_pos = self._return_positions(comp2)
        pairwise_matches = self._pairwise_dist(comp1_pos, comp2_pos)

        dist = 0
        # Loop through the union of  keys
        for key, value in pairwise_matches.items():
            dist += abs(comp1[key] - comp2[value[0]]) + value[1] * self.DIST_MOD

        if not flag:
            return dist, pairwise_matches

        else:
            return dist

    def levenshtein_dist(self, comp2, comp1 = None):
        """
        DEPRECATED, use min_flow_dist()

        Similar methodology to Hamming distance except we match the values from
        the root composition to the test compostion and then add the weights of
        the leftover elements to the distance metric
        """
        self.testing_formula = comp2

        comp1 = comp1 if comp1 is not None else self.normed_composition

        if isinstance(comp2, str):
            comp2 = self._parse_formula(comp2)
            comp2 = self._normalise_composition(comp2)


        # Calculate the hamming dist first
        dist, pairwise_matches = self.hamming_dist(comp2, flag=False)

        # Now mop  up the remaining elements that weren't mapped to one another
        for element, distribution in comp1.items():
            if element not in pairwise_matches:
                dist += distribution * self.LEVENSH_MOD

        for element, distribution in comp2.items():
            if element not in [value[0] for key, value in pairwise_matches.items()]:
                dist += distribution * self.LEVENSH_MOD

        return dist

    def _get_periodic_tab(self):
        """
        Attempt to load periodic data from the same folder, else download
        it from the web
        """
        try:
            with open('./ElementData.json') as json_data:
                periodic_data = json.load(json_data)
            return periodic_data

        except FileNotFoundError as e:
            print(f"ELement lookup table failed to load due to {e}")
            print("Attempting to download from the web, please allow firewall access")
            url = 'https://raw.githubusercontent.com/SurgeArrester/ChemHammer/master/ElementData.json?token=ABW3CRPQZQDGZHESYOVZLJK5JMV5U'

            try:
                response = urllib.request.urlopen(url)
                data = response.read()      # a `bytes` object
                data = data.decode('utf-8')
                periodic_data = json.loads(data)
                return periodic_data

            except urllib.error.URLError as e:
                print(f"Failed to download due to {e}, ensure you are connected to the internet and firewall access is granted, or download the ElementData.json file from {url} and save into this directory. Exiting.")
                sys.exit()

        except Exception as e:
            print(f"Failed due to {e}")

    def _is_balanced(self, formula):
        """Check if all sort of brackets come in pairs."""
        # Very naive check, just here because you always need some input checking
        c = Counter(formula)
        return c['['] == c[']'] and c['{'] == c['}'] and c['('] == c[')']

    def _dictify(self, tuples):
        """Transform tuples of tuples to a dict of atoms."""
        res = dict()
        for atom, n in tuples:
            try:
                res[atom] += float(n or 1)
            except KeyError:
                res[atom] = float(n or 1)
        return res

    def _fuse(self, mol1, mol2, w=1):
        """ Fuse 2 dicts representing molecules. Return a new dict. """
        return {atom: (mol1.get(atom, 0) + mol2.get(atom, 0)) * w for atom in set(mol1) | set(mol2)}

    def _parse(self, formula):
        """
        Return the molecule dict and length of parsed part.
        Recurse on opening brackets to parse the subpart and
        return on closing ones because it is the end of said subpart.
        """
        q = []
        mol = {}
        i = 0

        while i < len(formula):
            # Using a classic loop allow for manipulating the cursor
            token = formula[i]

            if token in self.CLOSERS:
                # Check for an index for this part
                m = re.match('\d+\.*\d*', formula[i+1:])
                if m:
                    weight = float(m.group(0))
                    i += len(m.group(0))
                else:
                    weight = 1

                submol = self._dictify(re.findall(self.ATOM_REGEX, ''.join(q)))
                return self._fuse(mol, submol, weight), i

            elif token in self.OPENERS:
                submol, l = self._parse(formula[i+1:])
                mol = self._fuse(mol, submol)
                # skip the already read submol
                i += l + 1
            else:
                q.append(token)

            i += 1

        # Fuse in all that's left at base level
        return self._fuse(mol, self._dictify(re.findall(self.ATOM_REGEX, ''.join(q)))), i

    def _parse_formula(self, formula):
        """Parse the formula and return a dict with occurences of each atom."""
        if not self._is_balanced(formula):
            raise ValueError("Your brackets not matching in pairs ![{]$[&?)]}!]")

        return self._parse(formula)[0]

    def _normalise_composition(self, input_comp):
        """ Sum up the numbers in our counter to get total atom count """
        composition = deepcopy(input_comp)
        # check it has been processed
        if isinstance(composition, str):
            composition = self._parse_formula(composition)

        atom_count =  sum(composition.values(), 0.0)

        for atom in composition:
            composition[atom] /= atom_count

        return composition

    def _get_atomic_num(self, element_string):
        """ Return atomic number from element string """
        for i, element in enumerate(self.periodic_tab['elements']):
            if element['symbol'] == element_string:
                return i

    def _get_position(self, element, metric=None):
        """
        Return either the x, y coordinate of an elements position, or the
        x-coordinate on the Pettifor numbering system as a 2-dimensional
        """
        metric = metric if metric is not None else self.distance_metric

        atomic_num = self._get_atomic_num(element)
        atom_info = self.periodic_tab['elements'][atomic_num]

        if metric == "mod_petti":
            return atom_info['mod_petti_num']

        elif metric == "petti":
            return atom_info['petti_num']

        elif metric == "manhattan":
            return (atom_info['xpos'], atom_info['ypos'])

    def _return_positions(self, composition):
        """ Return a dictionary of associated positions for each element """
        element_pos = {}

        for element in composition:
            element_pos[element] = self._get_position(element, metric="manhattan")

        return element_pos

    def _pairwise_dist(self, comp1_orig, comp2_orig):
        """
        Return matched pairs of closest elements
        TODO: Refactor the code to use a simpler datastructure than dictionaries
        of dictionarys
        """
        # Make copies of these to avoid changing original
        comp1 = deepcopy(comp1_orig)
        comp2 = deepcopy(comp2_orig)

        pairing_dict = {}

        # Check every element in the table, add to pair and pop from both lists if theres an exact match
        for element in list(comp1.keys()):
            if element in comp2:
                # Remove from both lists if we have a shared element and update with distance 0
                pairing_dict[element] = [element, 0]
                comp1.pop(element)
                comp2.pop(element)

        # List the elements we're using for use in a square matrix, create
        # a list of the coordinates of these on the periodic table, and make
        # a square distance matrix
        comp1_elements = list(comp1.keys())
        comp2_elements = list(comp2.keys())
        elements = comp1_elements + comp2_elements

        kw_func = partial(self._get_position, metric="manhattan")
        coords = list(map(kw_func, elements))

        if len(comp1) > 0 and len(comp2) > 0:
            dist_matrix = np.array(squareform(pdist(coords, metric="cityblock")))

            # As we only want to match those from opposing compositions we will
            # look at the top right quadrant only
            dist_cost = dist_matrix[:len(comp1), len(comp1):]

            # Use the minimum weight matching algorithm for bipartite graphs to find
            # the best combination of these
            row_ind, col_ind = linear_sum_assignment(dist_cost)

            for i, _ in enumerate(row_ind):
                pairing_dict[comp1_elements[row_ind[i]]] = [comp2_elements[col_ind[i]], dist_cost[row_ind[i]][col_ind[i]]]

        return pairing_dict

if __name__ == "__main__":
    main()