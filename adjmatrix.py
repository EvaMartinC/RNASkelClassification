from itertools import permutations
import numpy as np 
import networkx as nx

def asignar_nodos(eje_a_nodos, i, nodo_id):
    """ Asigna nodos a un segmento del tipo 2."""
    m1, m2 = eje_a_nodos[i]
    if m1 == -1 and m2 == -1:
        eje_a_nodos[i] = [nodo_id, nodo_id + 1]
        return nodo_id + 2
    elif m1 != -1 and m2 == -1:
        eje_a_nodos[i] = [m1, nodo_id]
        return nodo_id + 1
    elif m1 == -1 and m2 != -1:
        eje_a_nodos[i] = [nodo_id, m2]
        return nodo_id + 1
    return nodo_id

def clasificar_vecinos(array_connected_branches, i):
    """ Clasifica los vecinos en dos grupos según la conectividad."""
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

def actualizar_nodos(eje_a_nodos, group, nodo):
    """ Actualiza los nodos de los segmentos en un grupo."""
    if len(group) >= 1:       
        for segment in group[0]:
            n1, n2 = eje_a_nodos[segment]
            if n1 != -1 and n2 != -1:
                continue
            elif n1 != -1 and n1 != nodo:
                eje_a_nodos[segment] = [n1, nodo]
            elif n2 != -1 and n2 != nodo:
                eje_a_nodos[segment] = [nodo, n2]
            else:
                eje_a_nodos[segment] = [nodo, -1]
    else:
        print("possible mistake in graph")

def procesar_segmentos(array_connected_branches, seg_type):
    """ Procesa los segmentos y asigna nodos."""
    m = len(array_connected_branches)
    eje_a_nodos = {i: [-1, -1] for i in range(m)}
    nodo_id = 1
    for i in range(m):
        if seg_type[i] == 2:
            nodo_id = asignar_nodos(eje_a_nodos, i, nodo_id)
            group1, group2 = clasificar_vecinos(array_connected_branches, i)
            actualizar_nodos(eje_a_nodos, group1, eje_a_nodos[i][0])
            actualizar_nodos(eje_a_nodos, group2, eje_a_nodos[i][1])
        else:
            eje_a_nodos[i] = [eje_a_nodos[i][0],nodo_id]
            nodo_id += 1
    for v in eje_a_nodos.values():
        if -1 in v:
            v[v == -1] = nodo_id
    return eje_a_nodos


def find_adj_matrix(Gw,Gw2):
    AM1= nx.adjacency_matrix(Gw,nodelist=list(range(1,len(Gw)+1))).toarray()
    AM2= nx.adjacency_matrix(Gw2,nodelist=list(range(1,len(Gw2)+1))).toarray()
    n = len(AM1)
    best_perm = None
    min_cost = float('inf')

    # Probar todas las permutaciones posibles de nodos
    for perm in permutations(range(n)):  
        #print(perm)
        # Reordenar filas y columnas de A2 según la permutación
        AM2_permuted = AM2[np.ix_(perm, perm)]
        
        # Calcular la diferencia de pesos
        cost = np.sum(np.abs(AM1 - AM2_permuted))
        # Guardar la mejor permutación
        if cost < min_cost:
            min_cost = cost
            best_perm = perm
            
    best_AM2 = AM2[np.ix_(best_perm, best_perm)]
    # Imprimir la mejor asignación de nodos
    node_mapping = {i: best_perm[i] for i in range(n)}
    #print("Asignación óptima de nodos:", node_mapping)
    return AM1,best_AM2