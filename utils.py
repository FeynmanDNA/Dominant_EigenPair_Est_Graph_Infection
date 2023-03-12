import networkx as nx
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt


def adj_mat_ground_truth(G, input_is_matrix=False):
    if not input_is_matrix:
        print("#nodes, #edges = ", G.number_of_nodes(), G.number_of_edges())

        # NOTE: networkx.linalg.spectrum.adjacency_spectrum
        # Returns eigenvalues of the adjacency matrix of G.
        adj_spectrum = nx.adjacency_spectrum(G, weight='weight')  # weightstring or None, optional (default=’weight’)
        print("adj_spectrum =\n", adj_spectrum)
        print("sorted adj_spectrum = ", np.sort(adj_spectrum))

        A = nx.adjacency_matrix(G)
    else:
        # when direct import matrix data
        A = G

    print("A.shape = ", A.shape)
    # print("adj A = ")
    # print(A)
    # if want to see the full rows of A:
    # for line in A.todense():
    #     print(*line)

    # eigen ground-truth, might be too memory demanding
    if not input_is_matrix:
        print(A.todense(), A.todense().shape)
        w, v = LA.eig(A.todense())
    else:
        print(A)
        w, v = LA.eig(A)
    print("w:\n", w)
    print("v:\n", v)
    """
    The normalized (unit “length”) eigenvectors,
    such that the column v[:,i] is the eigenvector
    corresponding to the eigenvalue w[i].
    """
    # index of dominant eigenvalue
    dom_index_1 = np.argmax(w)
    dom_index_2 = np.argmin(w)
    print("dom_index 1 and 2 = ", dom_index_1, dom_index_2)
    dom_eigVec_1 = v[:, dom_index_1]
    dom_eigVec_2 = v[:, dom_index_2]
    with np.printoptions(precision=5, suppress=True, threshold=100):
        print("Sorted eigenvalues w/o abs:")
        print(np.sort(w))
        print("Sorted abs(eigenvalues):")
        print(np.sort(abs(w)))
        print("Dominant Eigenvector1:\n", dom_eigVec_1)
        print("Dominant Eigenvector2:\n", dom_eigVec_2)
    # verify the dominant eigenpair works
    print("A @ dom_eigVec_1=\n", A @ dom_eigVec_1)
    print("np.max(w)*dom_eigVec_1=\n", np.max(w)*dom_eigVec_1)
    print("A @ dom_eigVec_2=\n", A @ dom_eigVec_2)
    print("np.min(w)*dom_eigVec_2=\n", np.min(w)*dom_eigVec_2)

    # alternative ways to find the max/min eigenvalue
    if not input_is_matrix:
        e = np.linalg.eigvals(A.todense())
    else:
        e = np.linalg.eigvals(A)
    print("Largest eigenvalue:", max(e))
    print("Smallest eigenvalue:", min(e))
    return A, max(e), dom_eigVec_1

def angle_dot(a, b):
    """
    from https://stackoverflow.com/questions/64501805/dot-product-and-angle-in-degrees-between-two-numpy-arrays
    and
    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    """
    # Return the real part of the complex argument for symmetric graphs
    # dom_eigVec_1 = dom_eigVec_1.real  # optional

    # NOTE: dot for 2D input is matrix multiplication, not a dot product.
    # dot two 1D arrays takes a dot product and produces a scalar result
    # Use np.ravel (for a 1D view) or np.ndarray.flatten (for a 1D copy)
    a = np.ravel(a)
    b = np.ravel(b)
    # print("a, b shape = ", a.shape, b.shape)
    dot_product = np.dot(a, b)
    prod_of_norms = np.linalg.norm(a) * np.linalg.norm(b)
    # print("dot_product, prod_of_norms = ", dot_product, prod_of_norms)
    # print("dot_product / prod_of_norms = ", dot_product / prod_of_norms)
    angle = round(
        np.degrees(
            np.arccos(
                # clip the range to prevent NAN
                np.clip(
                    (dot_product / prod_of_norms).real,
                    -1.0, 1.0
                )
            )
        ),
        1
    )
    print("dot-prod angle: ", angle)
    return round(dot_product, 1), angle


def get_abs_angle(angle):
    """
    if angle near 180, get the difference to 180
    if angle near 0, get the difference to 0
    """
    dist_to_180 = abs(angle - 180)
    dist_to_0 = abs(angle - 0)
    return min(dist_to_180, dist_to_0)


def draw_with_options(G):
    options = {
        "font_size": 15,
        "node_size": 500,
        "node_color": "white",
        "edgecolors": "black",
    #     "linewidths": 5,
    #     "width": 5,
    }
    nx.draw(G, with_labels=True, **options)
    plt.show()
    nx.draw_circular(G, with_labels=True, **options)
    plt.show()

    pos = nx.spring_layout(G)
    nodes = nx.draw_networkx_nodes(
        G,
        pos,
        # node_color="indigo"
    )
    edges = nx.draw_networkx_edges(
        G,
        pos,
        arrowstyle="->",
        arrowsize=10,
        width=2,
    )
    plt.show()
