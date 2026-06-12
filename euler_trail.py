"""
euler_trail.py
==============
Standard cell layout generation using Euler trail algorithm.

Background
----------
In CMOS standard cell layout, transistors are arranged in a single row
for both PFET (pull-up) and NFET (pull-down) networks. The key insight is:

  If adjacent transistors share a drain/source net, no extra wire is needed
  between them. This minimizes cell width and routing complexity.

Finding such an optimal ordering is equivalent to finding an Euler trail
on a graph where:
  - Vertex = drain/source net  (VDD, VSS, Z, net_0, ...)
  - Edge   = transistor        (labeled by its gate net: A1, A2, Z_neg, ...)

An Euler trail visits every edge exactly once, which maps directly to
placing every transistor exactly once in a sequence where adjacent
transistors always share a net.

Input
-----
  Nangate_15nm.sp  : SPICE netlist containing transistor-level descriptions
                     of standard cells (AND2, NOR3, DFFRNQ, etc.)

Output
------
  placement_result.txt : optimal PFET / NFET transistor ordering per cell
  routing_pin_result.txt : (x, y, layer) pin coordinates for routing

Reference
---------
  Euler trail condition:
    - Eulerian circuit  : all vertices have even degree
    - Eulerian trail    : exactly 2 vertices have odd degree
    - Otherwise         : add dummy edges to pair up odd-degree vertices
"""

import copy
from itertools import permutations
from collections import deque


# ---------------------------------------------------------------------------
# Graph data class
# ---------------------------------------------------------------------------

class CellGraph:
    """
    Represents one standard cell as two graphs (PFET and NFET).

    Attributes
    ----------
    name        : cell name, e.g. 'AND2_X1'
    p_vertex    : list of net names used as drain/source in PFET network
    n_vertex    : list of net names used as drain/source in NFET network
    edge        : list of gate net names shared by both networks
    connect_list: [pfet_connections, nfet_connections]
                  each entry is [gate_net, drain_net, source_net]
    pmatrix     : adjacency matrix for PFET graph
    nmatrix     : adjacency matrix for NFET graph
    ans_p       : candidate Euler trail edge sequences for PFET
    ans_n       : candidate Euler trail edge sequences for NFET
    visited     : vertex visit sequences corresponding to ans_p / ans_n
    """
    def __init__(self, name):
        self.name = name
        self.p_vertex    = []
        self.n_vertex    = []
        self.edge        = []
        self.width       = [[], []]
        self.connect_list = [[], []]
        self.pmatrix     = []
        self.nmatrix     = []
        self.ans_p       = []
        self.ans_n       = []
        self.visited     = [[], []]


# ---------------------------------------------------------------------------
# Core Euler trail
# ---------------------------------------------------------------------------

def euler_path(matrix, start, n_vertices, visited):
    """
    Hierholzer's algorithm for Euler trail / circuit.

    Recursively traverses all edges from `start`, removing each edge as it
    is used (to prevent revisiting). Appends vertices to `visited` in
    reverse post-order, which produces the correct Euler sequence.

    Parameters
    ----------
    matrix    : adjacency matrix (modified in-place — pass a deep copy)
    start     : index of the starting vertex
    n_vertices: total number of vertices
    visited   : list that accumulates the trail (output)
    """
    for neighbor in range(n_vertices):
        while matrix[start][neighbor]:
            # Remove the edge before recursing (undirected: remove both directions)
            matrix[start][neighbor] -= 1
            matrix[neighbor][start] -= 1
            euler_path(matrix, neighbor, n_vertices, visited)
    # All edges from `start` are exhausted — record this vertex
    visited.append(start)


def vertex_degree(n_vertices, matrix):
    """
    Compute the degree of every vertex.

    Returns
    -------
    degree_list : list[int], degree_list[i] = degree of vertex i
    """
    degree_list = []
    for v in range(n_vertices):
        degree = sum(matrix[v][u] for u in range(n_vertices))
        degree_list.append(degree)
    return degree_list


# ---------------------------------------------------------------------------
# Trail → edge/vertex sequence conversion
# ---------------------------------------------------------------------------

def find_path(visit, edge, vertex, pn_flag):
    """
    Convert a vertex visit sequence into an ordered list of edge indices.

    Given the Euler trail as a sequence of vertex indices, identify which
    edge (transistor) connects each consecutive pair of vertices.

    Parameters
    ----------
    visit    : vertex index sequence from euler_path()
    edge     : global edge (gate net) list
    vertex   : [pfet_edge_info, nfet_edge_info]
               each entry: [edge_idx, vertex_a_idx, vertex_b_idx]
    pn_flag  : 'pfet' or 'nfet'

    Returns
    -------
    ans : list of edge indices representing the ordered transistor sequence
    """
    pn = 1 if pn_flag == 'nfet' else 0
    ans = []

    # Track how many times each edge and each connection entry can be used
    edge_use     = [[0] * len(edge)        for _ in range(2)]
    edge_num_use = [[0] * len(vertex[pn])  for _ in range(2)]

    for entry_idx in range(len(vertex[pn])):
        edge_idx = vertex[pn][entry_idx][0]
        edge_use[pn][edge_idx] += 1

    for step in range(len(visit) - 1):
        src, dst = visit[step], visit[step + 1]
        for entry_idx, entry in enumerate(vertex[pn]):
            e_idx, va, vb = entry
            if edge_use[pn][e_idx] == 0 or edge_num_use[pn][entry_idx] != 0:
                continue
            # Match edge in either direction (undirected)
            if (va == src and vb == dst) or (vb == src and va == dst):
                edge_use[pn][e_idx]         -= 1
                edge_num_use[pn][entry_idx]  = 1
                ans.append(entry_idx)
                break
    return ans


def euler_to_vertex(visited, vertex_names):
    """Map vertex index sequence → vertex name sequence."""
    return [vertex_names[i] for i in visited]


def ans_to_edge(ans, edge_names, edge_info):
    """Map edge-entry index sequence → gate net name sequence."""
    return [edge_names[edge_info[i][0]] for i in ans]


# ---------------------------------------------------------------------------
# Odd-degree handling (dummy edge insertion)
# ---------------------------------------------------------------------------

def change_info(start, matrix, edge, vertex, pn_flag, i, j):
    """
    Pair two odd-degree vertices by inserting a dummy edge between them.

    An Euler trail exists only when exactly 0 or 2 vertices have odd degree.
    When more odd-degree vertices exist, we add dummy edges to make them even,
    then remove the dummy positions from the final answer.

    Parameters
    ----------
    start   : current list of odd-degree vertex indices
    matrix  : adjacency matrix to modify
    edge    : edge (gate) list — must already contain 'dummy'
    vertex  : [pfet_connections, nfet_connections]
    pn_flag : 'pfet' or 'nfet'
    i, j    : indices into `start` — the two vertices to connect

    Returns
    -------
    new_start  : remaining odd-degree vertices after pairing
    new_vertex : updated connection list with dummy entry appended
    new_matrix : updated adjacency matrix
    """
    pn = 0 if pn_flag == 'pfet' else 1

    new_start  = copy.deepcopy(start)
    new_matrix = copy.deepcopy(matrix)
    new_vertex = copy.deepcopy(vertex)

    dummy_idx = edge.index('dummy')

    # Connect the two odd-degree vertices with a dummy edge
    new_matrix[start[i]][start[j]] += 1
    new_matrix[start[j]][start[i]] += 1
    new_vertex[pn].append([dummy_idx, start[i], start[j]])

    # Remove the paired vertices from the odd list
    del new_start[new_start.index(start[i])]
    del new_start[new_start.index(start[j])]

    return new_start, new_vertex, new_matrix


# ---------------------------------------------------------------------------
# Same-connection detection (handles parallel transistors)
# ---------------------------------------------------------------------------

def same_connect(vertex_list, pn_flag):
    """
    Find groups of edges that share identical endpoint pairs.

    In cells like NAND2, two PFET transistors connect the same two nets
    (VDD and Z_neg) — they are parallel. Their positions in the Euler trail
    can be swapped without affecting correctness. This function identifies
    such interchangeable groups so we can enumerate all valid permutations.

    Returns
    -------
    ans : list of groups, each group is a list of edge-entry indices
          that are interchangeable with each other
    """
    pn = 0 if pn_flag == 'pfet' else 1
    entries = vertex_list[pn]
    ans = []
    ans_cnt = 0

    for i in range(len(entries) - 1):
        group = [i]
        already_grouped = any(i in g for g in ans)
        if already_grouped:
            continue
        for j in range(i + 1, len(entries)):
            endpoints_i = entries[i][1:3]
            endpoints_j = entries[j][1:3]
            # Same endpoints in either order → parallel transistors
            if endpoints_i == endpoints_j or endpoints_i == endpoints_j[::-1]:
                group.append(j)
        if len(group) > 1:
            ans.append(group)

    return ans


def same_location(ans, same_groups):
    """
    Find the positions of each interchangeable group within the trail.

    Returns
    -------
    locations : list of position lists, one per group in same_groups
    """
    locations = []
    for group in same_groups:
        pos_list = []
        for member in group:
            for idx, edge_entry in enumerate(ans):
                if member == edge_entry:
                    pos_list.append(idx)
                    break
        locations.append(pos_list)
    return locations


def plus_ans(same_groups, locations, ans, total_ans,
             total_vertex_num, vertex_num, visited, visit, pn_flag):
    """
    Enumerate all valid permutations of interchangeable transistor groups.

    For each group of parallel transistors, generate all orderings and
    append them as additional candidate trails.
    """
    storage = [ans]

    if pn_flag == 'pfet':
        total_vertex_num[0].append(vertex_num[0])
        visited[0].append(visit[0])
    else:
        total_vertex_num[1].append(vertex_num[1])
        visited[1].append(visit[1])

    for g_idx, group in enumerate(same_groups):
        # Process one permutation at a time — avoids loading all into memory
        for perm in permutations(group, len(group)):
            if list(perm) == group:  # skip the original ordering
                continue
            for base in storage[:]:
                candidate = copy.deepcopy(base)
                for k in range(len(locations[g_idx])):
                    candidate[locations[g_idx][0] + k] = perm[k]
                if candidate not in storage:
                    storage.append(candidate)
                    if pn_flag == 'pfet':
                        total_vertex_num[0].append(vertex_num[0])
                        visited[0].append(visit[0])
                    else:
                        total_vertex_num[1].append(vertex_num[1])
                        visited[1].append(visit[1])

    for candidate in storage:
        total_ans.append(candidate)


# ---------------------------------------------------------------------------
# Main trail search
# ---------------------------------------------------------------------------

def odd_num_ans(start, start_idx, graph, total_vertex_num,
                vertex_num, matrix, edge, same_groups, pn_flag):
    """
    Run Euler trail from one starting vertex and store all valid results.

    Parameters
    ----------
    start     : list of candidate starting vertices
    start_idx : which entry in `start` to use this call
    graph     : CellGraph instance (results stored here)
    ...
    """
    visit = [[], []]
    pn = 0 if pn_flag == 'pfet' else 1

    mat_copy = copy.deepcopy(matrix)
    euler_path(mat_copy, start[start_idx],
               len(graph.p_vertex) if pn_flag == 'pfet' else len(graph.n_vertex),
               visit[pn])

    trail_edges = find_path(visit[pn], edge, vertex_num, pn_flag)
    locations   = same_location(trail_edges, same_groups)

    total_ans = graph.ans_p if pn_flag == 'pfet' else graph.ans_n
    plus_ans(same_groups, locations, trail_edges, total_ans,
             total_vertex_num, vertex_num, graph.visited, visit, pn_flag)


def ans_oddnum(pn_flag, odd_num, start, graph, vertex_num_list, total_vertex_num):
    """
    Dispatch Euler trail search based on the number of odd-degree vertices.

    Euler trail exists when odd_num <= 2. For higher odd counts, we pair
    up odd-degree vertices with dummy edges (2 at a time) and retry.

    odd_num == 0 or 1 : Eulerian circuit — any vertex can be start
    odd_num == 2       : Eulerian trail  — must start from an odd vertex
    odd_num == 4,6,8   : Add (odd_num/2 - 1) dummy edges, reducing to 2
    """
    matrix = graph.pmatrix if pn_flag == 'pfet' else graph.nmatrix
    copy_edge = copy.deepcopy(graph.edge)

    def run(start_list, vlist, mat, edge_list):
        for s_idx in range(len(start_list)):
            sg = same_connect(vlist, pn_flag)
            odd_num_ans(start_list, s_idx, graph, total_vertex_num,
                        vlist, mat, edge_list, sg, pn_flag)

    if odd_num <= 3:
        run(start, vertex_num_list, matrix, graph.edge)

    elif odd_num == 4:
        copy_edge.append('dummy')
        for i in range(len(start) - 1):
            for j in range(i + 1, len(start)):
                s2, v2, m2 = change_info(start, matrix, copy_edge,
                                         vertex_num_list, pn_flag, i, j)
                run(s2, v2, m2, copy_edge)

    elif odd_num == 6:
        copy_edge.append('dummy')
        for i in range(len(start) - 1):
            for j in range(i + 1, len(start)):
                s2, v2, m2 = change_info(start, matrix, copy_edge,
                                         vertex_num_list, pn_flag, i, j)
                for ii in range(len(s2) - 1):
                    for jj in range(ii + 1, len(s2)):
                        s3, v3, m3 = change_info(s2, m2, copy_edge, v2, pn_flag, ii, jj)
                        run(s3, v3, m3, copy_edge)

    elif odd_num == 8:
        copy_edge.append('dummy')
        for i in range(len(start) - 1):
            for j in range(i + 1, len(start)):
                s2, v2, m2 = change_info(start, matrix, copy_edge,
                                         vertex_num_list, pn_flag, i, j)
                for ii in range(len(s2) - 1):
                    for jj in range(ii + 1, len(s2)):
                        s3, v3, m3 = change_info(s2, m2, copy_edge, v2, pn_flag, ii, jj)
                        for iii in range(len(s3) - 1):
                            for jjj in range(iii + 1, len(s3)):
                                s4, v4, m4 = change_info(s3, m3, copy_edge,
                                                         v3, pn_flag, iii, jjj)
                                run(s4, v4, m4, copy_edge)

    return copy_edge


# ---------------------------------------------------------------------------
# PFET ↔ NFET trail matching
# ---------------------------------------------------------------------------

def Is_include(name, ans, connect_info, case_num, odd_num, case_odd, max_graph):
    """
    Check whether the PFET Euler trail is compatible with the NFET network.

    After finding a PFET Euler trail (ans), verify that the same gate ordering
    can also be applied to the NFET network. If the two networks share the same
    gate sequence, the transistors can be placed in a single row with matching
    columns — which is the goal of standard cell layout.

    Parameters
    ----------
    name         : cell name (some cells need special handling)
    ans          : PFET Euler trail as a list of edge-entry indices
    connect_info : NFET connection list [edge_idx, vertex_a, vertex_b]
    case_num     : recursion case index (handles parallel transistors)
    odd_num      : number of odd-degree vertices in NFET graph
    case_odd     : odd vertex case index
    max_graph    : total number of edges (used for dummy edge index)

    Returns
    -------
    (include, ans_vertex, dummy_loc, ans_loc)
    include    : 1 if compatible, 0 if not
    ans_vertex : NFET vertex visit sequence matching the PFET trail
    dummy_loc  : positions where dummy transistors were inserted
    ans_loc    : positions of dummy edges in the original trail
    """
    odd = odd_num
    cnt = 1
    cnt_odd = 0
    loc = []
    tmp_list = []
    initial_ans_loc = []

    if len(ans) != len(connect_info):
        if len(connect_info) in ans:
            initial_ans_loc.append(ans.index(len(connect_info)))
        if len(ans) - 1 in ans:
            del ans[ans.index(len(ans) - 1)]
    initial_len_ans = len(ans)

    visit = deque()

    if case_num == 2:
        return 0, [], loc, initial_ans_loc
    if len(ans) == 1:
        return 1, connect_info[0][1:], loc, initial_ans_loc

    for i in range(1, len(ans)):
        connect_new = connect_info[ans[i]][1:]

        if (connect_new[case_num] in connect_info[ans[i-1]][1:]) and \
           (len(visit) == 0 or connect_new[case_num] != visit[len(visit)-1]):
            visit.append(connect_new[case_num])
            connect_new.remove(connect_new[case_num])
            cnt += 1
            if i == 1:
                if visit[0] == connect_info[ans[0]][1]:
                    visit.appendleft(connect_info[ans[0]][2])
                elif visit[0] == connect_info[ans[0]][2]:
                    visit.appendleft(connect_info[ans[0]][1])
            if i == (len(ans) - 1):
                if visit[-1] == connect_info[ans[i]][1]:
                    visit.append(connect_info[ans[i]][2])
                elif visit[-1] == connect_info[ans[i]][2]:
                    visit.append(connect_info[ans[i]][1])

        elif (connect_new[1-case_num] in connect_info[ans[i-1]][1:]) and \
             (len(visit) == 0 or connect_new[1-case_num] != visit[len(visit)-1]):
            visit.append(connect_new[1-case_num])
            connect_new.remove(connect_new[1-case_num])
            cnt += 1
            if i == 1:
                if visit[0] == connect_info[ans[0]][1]:
                    visit.appendleft(connect_info[ans[0]][2])
                elif visit[0] == connect_info[ans[0]][2]:
                    visit.appendleft(connect_info[ans[0]][1])
            if i == (len(ans) - 1):
                if visit[-1] == connect_info[ans[i]][1]:
                    visit.append(connect_info[ans[i]][2])
                elif visit[-1] == connect_info[ans[i]][2]:
                    visit.append(connect_info[ans[i]][1])

        elif (odd - 2 > 1) or (name == 'DFFRNQ_X1' and odd - 2 >= 0):
            odd = odd - 2
            tmp_info = connect_info[ans[i-1]][1:]
            if len(visit) == 0:
                tmp_list.append(max_graph)
                tmp_list.append(tmp_info[1-case_odd])
                tmp_list.append(connect_info[ans[i]][2])
                visit.append(tmp_info[case_odd])
                visit.append(tmp_info[1-case_odd])
                visit.append(connect_info[ans[i]][2])
            else:
                tmp_list.append(max_graph)
                tmp_list.append(tmp_info[1-tmp_info.index(visit[-1])])
                tmp_list.append(connect_info[ans[i]][2])
                visit.append(tmp_info[1-tmp_info.index(visit[-1])])
                if i == (len(ans) - 1):
                    next_v = connect_info[ans[i]][1:]
                    for last in next_v:
                        visit.append(last)
                else:
                    next_v = connect_info[ans[i+1]][1:]
                    same_value = list(set(connect_new).intersection(next_v))
                    visit.append(connect_new[1-connect_new.index(same_value[0])])
                i = i + cnt_odd
                loc.append(i)
                cnt_odd += 1
                cnt += 1
        else:
            include, ans_pvertex, dummy_loc, ans_loc = Is_include(
                name, ans, connect_info, case_num+1, odd_num, case_odd, max_graph)
            if include == 1:
                return 1, ans_pvertex, loc, initial_ans_loc
            else:
                return 0, [], loc, initial_ans_loc

    if cnt == len(ans):
        if len(tmp_list) != 0:
            connect_info.append(tmp_list)
        for i, dummy in enumerate(loc):
            ans.insert(dummy, initial_len_ans)
        return 1, list(visit), loc, initial_ans_loc
    else:
        return 0, [], loc, initial_ans_loc


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_placement(ans_edge, ans_pvertex, ans_nvertex,
                    edge_info, graph, result, dummy_loc, ans_loc):
    """
    Write one valid placement solution to the result file.

    Format
    ------
    Name : <cell_name>
    PMOS : <drain/source net sequence>
    NMOS : <drain/source net sequence>
    Gate : <gate net sequence>
    """
    result.write(f'Name : {graph.name}\n')

    def write_row(label, vertex_seq, vertex_names, dummy_positions, insert_loc):
        result.write(f'{label} :')
        for idx, v in enumerate(vertex_seq):
            result.write(f' {vertex_names[v]}')
            # Insert duplicate at dummy location to mark filler transistor
            if idx in dummy_positions:
                result.write(f' {vertex_names[v]}')
        result.write('\n')

    if len(ans_pvertex) == len(ans_nvertex):
        # Balanced — write directly
        result.write('PMOS :')
        for v in ans_pvertex:
            result.write(f' {graph.p_vertex[v]}')
        result.write('\n')
        result.write('NMOS :')
        for v in ans_nvertex:
            result.write(f' {graph.n_vertex[v]}')
        result.write('\n')
    elif len(ans_pvertex) < len(ans_nvertex):
        # PFET side needs a filler (dummy) transistor
        write_row('PMOS', ans_pvertex, graph.p_vertex, dummy_loc, ans_loc)
        result.write('NMOS :')
        for v in ans_nvertex:
            result.write(f' {graph.n_vertex[v]}')
        result.write('\n')
    else:
        # NFET side needs a filler (dummy) transistor
        result.write('PMOS :')
        for v in ans_pvertex:
            result.write(f' {graph.p_vertex[v]}')
        result.write('\n')
        write_row('NMOS', ans_nvertex, graph.n_vertex, dummy_loc, ans_loc)

    result.write('Gate :')
    for ans_cnt in ans_edge:
        if ans_cnt >= len(edge_info) or len(edge_info[ans_cnt]) == 0:
            result.write(' dummy')
        elif int(edge_info[ans_cnt][0]) >= len(graph.edge):
            result.write(' dummy')
        else:
            result.write(' %s' % graph.edge[edge_info[ans_cnt][0]])
    result.write('\n')
    result.write('\n')


def routing_all(pmos, nmos, gate):
    """
    Identify shared nets between PMOS, NMOS, and Gate pin sequences.

    For each net that appears in more than one location, record the
    (pmos_indices, nmos_indices, gate_indices) that share it.
    This determines which pins need to be connected by routing wires.

    Returns
    -------
    ans_list : list of [pmos_locs, nmos_locs, gate_locs] per shared net
    """
    ans_list     = []
    used_p       = [False] * len(pmos)
    used_n       = [False] * len(nmos)
    used_g       = [False] * len(gate)

    for pi, pnet in enumerate(pmos):
        if used_p[pi]:
            continue
        used_p[pi] = True
        grp_p, grp_n, grp_g = [pi], [], []

        for pi2, pnet2 in enumerate(pmos):
            if not used_p[pi2] and pnet == pnet2 and pnet != 'VDD':
                grp_p.append(pi2)
                used_p[pi2] = True

        for ni, nnet in enumerate(nmos):
            if not used_n[ni] and pnet == nnet and pnet != 'VSS':
                grp_n.append(ni)
                used_n[ni] = True

        for gi, gnet in enumerate(gate):
            if not used_g[gi] and pnet == gnet:
                grp_g.append(gi)
                used_g[gi] = True

        # Only record if the net appears in at least 2 different locations
        empty = sum(1 for g in [grp_p, grp_n, grp_g] if not g)
        multi = sum(1 for g in [grp_p, grp_n, grp_g] if len(g) > 1)
        if empty < 2 or multi > 0:
            ans_list.append([grp_p, grp_n, grp_g])

    # Also check nets that only appear in NMOS
    used_p = [False] * len(pmos)
    used_n = [False] * len(nmos)
    used_g = [False] * len(gate)

    for ni, nnet in enumerate(nmos):
        if used_n[ni]:
            continue
        used_n[ni] = True
        grp_p, grp_n, grp_g = [], [ni], []

        for pi, pnet in enumerate(pmos):
            if not used_p[pi] and nnet == pnet and pnet != 'VDD':
                grp_p.append(pi)
                used_p[pi] = True

        for ni2, nnet2 in enumerate(nmos):
            if not used_n[ni2] and nnet == nnet2 and nnet != 'VSS':
                grp_n.append(ni2)
                used_n[ni2] = True

        for gi, gnet in enumerate(gate):
            if not used_g[gi] and nnet == gnet:
                grp_g.append(gi)
                used_g[gi] = True

        empty = sum(1 for g in [grp_p, grp_n, grp_g] if not g)
        multi = sum(1 for g in [grp_p, grp_n, grp_g] if len(g) > 1)
        if (empty < 2 or multi > 0) and [grp_p, grp_n, grp_g] not in ans_list:
            ans_list.append([grp_p, grp_n, grp_g])

    return ans_list


def write_routing(routing_list, result):
    """
    Write pin coordinates for each shared net.

    Coordinate system
    -----------------
    PMOS pins : y = 6..9  (top rail)
    NMOS pins : y = 0..3  (bottom rail)
    Gate pins : x = odd column, y = 4..5 (middle)

    Each transistor occupies 2 x-columns (even = diffusion, odd = gate).
    """
    for grp_p, grp_n, grp_g in routing_list:
        for p in grp_p:
            for layer in range(4):          # 4 metal layers per diffusion pin
                result.write(f'{2*p} {layer+6} 0 ')

        for n in grp_n:
            for layer in range(4):
                result.write(f'{2*n} {layer} 0 ')

        for g in grp_g:
            for layer in range(2):          # 2 layers per gate pin
                result.write(f'{2*g+1} {layer+4} 0 ')

        result.write('\n')
