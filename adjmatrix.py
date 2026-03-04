from itertools import permutations
import numpy as np 
import networkx as nx
from networkx.algorithms import isomorphism


def assing_nodes(edge_to_nodes, i, node_id):
    """ Assign node to a type 2 segment."""
    m1, m2 = edge_to_nodes[i]
    if m1 == -1 and m2 == -1:
        edge_to_nodes[i] = [node_id, node_id + 1]
        return node_id + 2
    elif m1 != -1 and m2 == -1:
        edge_to_nodes[i] = [m1, node_id]
        return node_id + 1
    elif m1 == -1 and m2 != -1:
        edge_to_nodes[i] = [node_id, m2]
        return node_id + 1
    return node_id



def classify_neighbors(array_connected_branches, i):
    """ Classify neighbors according to connectivity"""
    group1, group2 = [], []
    for neighbor in array_connected_branches[i]:
        segments = list(set(array_connected_branches[neighbor]) & set(array_connected_branches[i]))
        if not group1:
            segments.append(neighbor)
            group1.append(segments)
        elif neighbor not in group1[0] and not group2:
            segments.append(neighbor)
            group2.append(segments)
        elif group2 and neighbor not in group2[0]:
            segments.append(neighbor)
            group1.append(segments)
    
    return group1, group2

def upload_nodes(edge_to_nodes, group, nodo):
    """ Upload nodes of the edges in a group. """
    df =False
    if len(group) >= 1:       
        for segment in group[0]:
            n1, n2 = edge_to_nodes[segment]
            if n1 != -1 and n2 != -1:
                df = True
                continue
            elif n1 != -1 and n1 != nodo:
                edge_to_nodes[segment] = [n1, nodo]
            elif n2 != -1 and n2 != nodo:
                edge_to_nodes[segment] = [nodo, n2]
            else:
                edge_to_nodes[segment] = [nodo, -1]
    else:
        print("possible mistake in graph")
    return df

def process_segments(array_connected_branches, seg_type):
    """ Process the edges and assign nodes. """
    m = len(array_connected_branches)
    edge_to_nodes = {i: [-1, -1] for i in range(m)}
    node_id = 1

    for i in range(m):
        if seg_type[i] == 2:
            node_id = assing_nodes(edge_to_nodes, i, node_id)
            group1, group2 = classify_neighbors(array_connected_branches, i)
            df = upload_nodes(edge_to_nodes, group2, edge_to_nodes[i][1])
            if df == True:
                df = upload_nodes(edge_to_nodes, group1, edge_to_nodes[i][1])
            else:
                df = upload_nodes(edge_to_nodes, group1, edge_to_nodes[i][0])

        elif seg_type[i] == 1:
            edge_to_nodes[i] = [edge_to_nodes[i][0],node_id]
            node_id += 1
        else:
            #type 0:
            edge_to_nodes[0] =  [1,2]
    if 2 not in seg_type:
        for v in edge_to_nodes.values():
            if -1 in v:
                v[v == -1] = node_id
            
    return edge_to_nodes



def find_adj_matrix_exact_pattern(AM1, AM2, tol=1e-12, max_isomorphisms=None, refine_node_signature=True):
    """
    Fidns a permutation of AM2 such as:
      - The pattern of non-zero elements (representing an edge) matches exactly with AM1
      - Among all valid permutations, minimises the sum of |AM1 - AM2_permuted| over the non-zero indices of AM1.

    Args:
      AM1, AM2: np.ndarray (n x n), weigthed adjencency matrixes.
      tol: tolerance to consider an element "non-zero" 
      max_isomorphisms: if not None, it limits the number of isomorphs to be explored to that number 
      refine_node_signature: if True, it adds a signature (grade + neighbors) to reduce the pairing candidates.

    Returns:
      (min_cost, best_AM2_permuted, best_perm)
      - min_cost: float (infinite or None if there is no isomorphism).
      - best_AM2_permuted: AM2 permutaded to align with AM1 (or None if there is no solution).
      - best_perm: list such as best_perm[i] = j (node i of AM1 corresponds to node j of AM2).
    """
    AM1 = np.asarray(AM1)
    AM2 = np.asarray(AM2)
    if AM1.shape != AM2.shape or AM1.ndim != 2 or AM1.shape[0] != AM1.shape[1]:
        raise ValueError("AM1 and AM2 must be square matrices of the same dimension.")

    n = AM1.shape[0]
    mask1 = np.abs(AM1) > tol
    mask2 = np.abs(AM2) > tol

    # Detect directed/undirected by symmetry (if both symmetrical -> undirected)
    directed = not (np.allclose(AM1, AM1.T, atol=tol) and np.allclose(AM2, AM2.T, atol=tol))

    # Constructing binary graphs (only indicating the existence of edges)
    if directed:
        G1 = nx.DiGraph()
        G2 = nx.DiGraph()
    else:
        G1 = nx.Graph()
        G2 = nx.Graph()

    G1.add_nodes_from(range(n))
    G2.add_nodes_from(range(n))

    if directed:
        for i in range(n):
            for j in range(n):
                if mask1[i, j]:
                    G1.add_edge(i, j)
                if mask2[i, j]:
                    G2.add_edge(i, j)
    else:
        for i in range(n):
            for j in range(i, n):
                if mask1[i, j]:
                    G1.add_edge(i, j)
                if mask2[i, j]:
                    G2.add_edge(i, j)

    # If different numbers of edges => impossible
    if G1.number_of_edges() != G2.number_of_edges():
        return None, None, None

    # Add nodal signatures to speed up the process 
    if refine_node_signature:
        if directed:
            outdeg1 = np.array([G1.out_degree(i) for i in range(n)], dtype=int)
            indeg1  = np.array([G1.in_degree(i) for i in range(n)], dtype=int)
            outdeg2 = np.array([G2.out_degree(i) for i in range(n)], dtype=int)
            indeg2  = np.array([G2.in_degree(i) for i in range(n)], dtype=int)

            def node_signature_directed(G, AM_mask):
                sigs = {}
                for i in range(n):
                    out_neigh = sorted([int(G.out_degree(j)) for j in G.successors(i)])
                    in_neigh  = sorted([int(G.in_degree(j)) for j in G.predecessors(i)])
                    sigs[i] = (int(G.out_degree(i)), int(G.in_degree(i)), tuple(out_neigh), tuple(in_neigh))
                return sigs

            sigs1 = node_signature_directed(G1, mask1)
            sigs2 = node_signature_directed(G2, mask2)

            for i in range(n):
                G1.nodes[i]['sig'] = sigs1[i]
                G2.nodes[i]['sig'] = sigs2[i]

            node_match = lambda a, b: a.get('sig') == b.get('sig')

        else:
            deg1 = np.array([G1.degree(i) for i in range(n)], dtype=int)
            deg2 = np.array([G2.degree(i) for i in range(n)], dtype=int)

            sigs1 = {}
            sigs2 = {}
            for i in range(n):
                neigh_degs_1 = tuple(sorted(int(G1.degree(j)) for j in G1.neighbors(i)))
                neigh_degs_2 = tuple(sorted(int(G2.degree(j)) for j in G2.neighbors(i)))
                sigs1[i] = (int(G1.degree(i)), neigh_degs_1)
                sigs2[i] = (int(G2.degree(i)), neigh_degs_2)

            for i in range(n):
                G1.nodes[i]['sig'] = sigs1[i]
                G2.nodes[i]['sig'] = sigs2[i]

            node_match = lambda a, b: a.get('sig') == b.get('sig')
    else:
        node_match = None

    # Selecting the correct matcher class
    if directed:
        GM = isomorphism.DiGraphMatcher(G1, G2, node_match=node_match)
    else:
        GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)

    min_cost = float('inf')
    best_perm = None
    best_AM2 = None

    # Iterate isomorphisms (VF2); there may be many if the graph is very symmetrical.
    for k, mapping in enumerate(GM.isomorphisms_iter()):
        # mapping: node of G1 -> node of G2
        perm = [mapping[i] for i in range(n)]
        AM2p = AM2[np.ix_(perm, perm)]

        # redundant verification: patterns must match exactly
        if not np.array_equal(mask1, np.abs(AM2p) > tol):
            continue  # por seguridad, aunque VF2 ya debería garantizar esto

        # calculate cost only on non-zero indices of AM1
        mask = mask1
        cost = np.sum(np.abs(AM1[mask] - AM2p[mask]))

        if cost < min_cost:
            min_cost = cost
            best_perm = perm.copy()
            best_AM2 = AM2p.copy()

        if (max_isomorphisms is not None) and (k + 1 >= max_isomorphisms):
            break

    if best_perm is None:
        return None, None, None
    else:
        return float(min_cost), best_AM2, best_perm

