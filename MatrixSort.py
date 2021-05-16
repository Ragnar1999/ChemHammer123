"""
Author: Cameron Hargreaves

This code will take in a square distance matrix and then use a hierarchical
clustering algorithm to sort these points into a dendogram
"""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from fastcluster import linkage
from sklearn import datasets
import matplotlib.pyplot as plt

def main():
    iris = datasets.load_iris()
    iris.data.shape
    dist_mat = squareform(pdist(iris.data))

    N = len(iris.data)
    X = iris.data[np.random.permutation(N),:]

    dist_mat = squareform(pdist(X))

    methods = ["ward", "single", "average", "complete"]
    for method in methods:
        print("Method:\t", method)

        sorted_mat = DistanceMatrixSorter(dist_mat, method)
        # ordered_dist_mat, res_order, res_linkage = compute_serial_matrix(dist_mat,method)

        plt.pcolormesh(sorted_mat.ordered_dist_mat)
        plt.xlim([0,N])
        plt.ylim([0,N])
        plt.show()

class DistanceMatrixSorter():
    """
    This class takes in a distance matrix and sorts this via the dendrogram
    hierarchical clustering method

    Sorting methods are "ward", "single", "average", "complete"
    """
    def __init__(self, dist_matrix, method="complete"):
        self.dist_matrix = dist_matrix
        self.method = method
        ordered_dist_mat, sorted_index, dendogram = self.compute_serial_matrix(dist_matrix, method)
        self.ordered_dist_mat = ordered_dist_mat
        self.sorted_index = sorted_index
        self.dendogram = dendogram

    def compute_serial_matrix(self, dist_mat,method="ward"):
        '''
            input:
                - dist_mat is a distance matrix
                - method = ["ward","single","average","complete"]
            output:
                - seriated_dist is the input dist_mat,
                but with re-ordered rows and columns
                according to the seriation, i.e. the
                order implied by the hierarchical tree
                - res_order is the order implied by
                the hierarhical tree
                - res_linkage is the hierarhical tree (dendrogram)

            compute_serial_matrix transforms a distance matrix into
            a sorted distance matrix according to the order implied
            by the hierarchical tree (dendrogram)
        '''
        N = len(dist_mat)
        flat_dist_mat = squareform(dist_mat)

        # Create a hierarchical dendogram
        res_linkage = linkage(flat_dist_mat, method=self.method, preserve_input=False)
        # Sort this dendogram to get the resultant 1D sorted indices
        res_order = self.seriation(res_linkage, N, N + N-2)

        # Placeholder for the sorted distance matrix
        seriated_dist = np.zeros((N, N))
        a, b = np.triu_indices(N, k=1)

        # Create the upper right triangle of the distance matrix
        seriated_dist[a, b] = dist_mat[[res_order[i] for i in a], [res_order[j] for j in b]]
        # Mirror this into the bottom left triangle
        seriated_dist[b, a] = seriated_dist[a,b]

        return seriated_dist, res_order, res_linkage

    def seriation(self, Z, N, cur_index):
        '''
            input:
                - Z is a hierarchical tree (dendrogram)
                - N is the number of points given to the clustering process
                - cur_index is the position in the tree for the recursive traversal
            output:
                - order implied by the hierarchical tree Z

            Recursively works its' way through the hierarchical tree and adds each
            of the leaves in order
        '''
        if cur_index < N:
            return [cur_index]

        else:
            left = int(Z[cur_index-N, 0])
            right = int(Z[cur_index-N, 1])
            return (self.seriation(Z, N, left) + self.seriation(Z, N, right))

if __name__ == "__main__":
    main()
