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
    Encuentra una permutación de AM2 tal que:
      - el patrón de índices no nulos (presencia de aristas) coincida EXACTAMENTE con AM1,
      - entre todas las permutaciones válidas minimiza la suma de |AM1 - AM2_permuted| sobre los índices no nulos de AM1.

    Args:
      AM1, AM2: np.ndarray (n x n), matrices de adyacencia (pesadas).
      tol: tolerancia para considerar un elemento "no nulo".
      max_isomorphisms: si no es None, limitamos la exploración a ese número de isomorfismos.
      refine_node_signature: si True, añadimos una firma (grado + vecinos) para reducir candidatos de emparejado.

    Returns:
      (min_cost, best_AM2_permuted, best_perm)
      - min_cost: float (infinito o None si no hay isomorfismo).
      - best_AM2_permuted: AM2 permutada para alinear con AM1 (o None si no hay solución).
      - best_perm: lista tal que best_perm[i] = j (nodo i de AM1 corresponde a nodo j de AM2).
    """
    AM1 = np.asarray(AM1)
    AM2 = np.asarray(AM2)
    if AM1.shape != AM2.shape or AM1.ndim != 2 or AM1.shape[0] != AM1.shape[1]:
        raise ValueError("AM1 y AM2 deben ser matrices cuadradas de la misma dimensión")

    n = AM1.shape[0]
    mask1 = np.abs(AM1) > tol
    mask2 = np.abs(AM2) > tol

    # Detectar dirigido/undirigido por simetría (si ambas simétricas -> no dirigido)
    directed = not (np.allclose(AM1, AM1.T, atol=tol) and np.allclose(AM2, AM2.T, atol=tol))

    # Construir grafos binarios (sólo indicando existencia de arista)
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

    # Si números de aristas distintos => imposible
    if G1.number_of_edges() != G2.number_of_edges():
        return None, None, None

    # Añadir firmas nodales para podar (opcional pero suele acelerar mucho)
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

    # Seleccionar la clase de matcher correcta
    if directed:
        GM = isomorphism.DiGraphMatcher(G1, G2, node_match=node_match)
    else:
        GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)

    min_cost = float('inf')
    best_perm = None
    best_AM2 = None

    # Iterar isomorfismos (VF2); puede ser muchos si el grafo es muy simétrico
    for k, mapping in enumerate(GM.isomorphisms_iter()):
        # mapping: nodo de G1 -> nodo de G2
        perm = [mapping[i] for i in range(n)]
        AM2p = AM2[np.ix_(perm, perm)]

        # comprobación redundante: patrones deben coincidir exactamente
        if not np.array_equal(mask1, np.abs(AM2p) > tol):
            continue  # por seguridad, aunque VF2 ya debería garantizar esto

        # calcular coste solamente sobre índices no nulos de AM1
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

# def find_adj_matrix(AM1,AM2):
#     n = len(AM1)
#     best_perm = None
#     min_cost = float('inf')
    

#     # Probar todas las permutaciones posibles de nodos
#     for perm in permutations(range(n)):  
#         # Reordenar filas y columnas de A2 según la permutación
#         AM2_permuted = AM2[np.ix_(perm, perm)]
#         cost = 0
        
#         # Calcular la diferencia de pesos
#         for i in range(n):
#             for j in range(n):
#                 if AM1[i][j] !=0 and AM2_permuted[i][j] !=0:
#                     cost = cost + np.abs(AM1[i][j] - AM2_permuted[i][j])
#                 elif AM1[i][j] ==0  and AM2_permuted[i][j] ==0:
#                     continue
#                 else:
#                     cost = cost + 1000

#         # Guardar la mejor permutación
#         if cost < min_cost:
#             min_cost = cost
#             best_perm = perm
            
#     best_AM2 = AM2[np.ix_(best_perm, best_perm)]
#     # Imprimir la mejor asignación de nodos
#     node_mapping = {i: best_perm[i] for i in range(n)}

#     #print("Asignación óptima de nodos:", node_mapping)
#     return min_cost,best_AM2