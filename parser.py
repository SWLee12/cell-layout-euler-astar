"""
parser.py
=========
Parses a SPICE netlist (.sp) and builds CellGraph objects.

SPICE format (Nangate_15nm.sp)
-------------------------------
Each standard cell block looks like:

  .SUBCKT AND2_X1 ...   ← cell starts
  * X1.Cellname: AND2_X1.
  * PININFO ...          ← transistor list begins after this
  M_1 Z_neg A1 VDD VDD pfet W=0.028   ← transistor line
  ...
  .ENDS                  ← cell ends

Transistor line fields (space-separated)
  [0] instance  M_1
  [1] drain     Z_neg
  [2] gate      A1
  [3] source    VDD
  [4] bulk      VDD
  [5] type      pfet / nfet
  [6] W=...     width
"""

from euler_trail import CellGraph


def load_cells(spice_file):
    """
    Parse a SPICE netlist and return a list of CellGraph objects.

    Parameters
    ----------
    spice_file : path to .sp file

    Returns
    -------
    cells : list[CellGraph]
    """
    cells = []
    in_cell = False
    in_transistors = False
    current = None

    with open(spice_file, 'r') as f:
        for line in f:
            # ── Cell start ──────────────────────────────────────────────
            # Line format: "* Cellname:   AND2_X1.   *"
            if 'Cellname' in line and not in_cell:
                if any(x in line for x in ['DFFSNQ_X1', 'SDFFSNQ_X1', '_X2', '_X4', '_X6', '_X8', '_X12', '_X16']):
                    continue
                # Extract cell name: last token ending with '.'
                parts = line.split()
                name_token = next((p for p in parts if p.endswith('.')), None)
                if name_token is None:
                    continue
                in_cell = True
                cell_name = name_token[:-1]   # strip trailing '.'
                current = CellGraph(cell_name)
                cells.append(current)

            # ── Transistor block begins ──────────────────────────────────
            elif 'PININFO' in line and in_cell:
                in_transistors = True

            # ── Transistor line ──────────────────────────────────────────
            elif line.startswith('M_') and in_cell and in_transistors:
                parts = line.split()
                # parts: [instance, drain, gate, source, bulk, type, W=...]
                drain  = parts[1]
                gate   = parts[2]
                source = parts[3]
                tr_type = parts[5]   # 'pfet' or 'nfet'

                # Register gate net as a graph edge (once per unique net)
                if gate not in current.edge:
                    current.edge.append(gate)

                if tr_type == 'pfet':
                    current.connect_list[0].append([gate, drain, source])
                    if drain  not in current.p_vertex:
                        current.p_vertex.append(drain)
                    if source not in current.p_vertex:
                        current.p_vertex.append(source)

                elif tr_type == 'nfet':
                    current.connect_list[1].append([gate, drain, source])
                    if drain  not in current.n_vertex:
                        current.n_vertex.append(drain)
                    if source not in current.n_vertex:
                        current.n_vertex.append(source)

            # ── Cell end ─────────────────────────────────────────────────
            elif '.ENDS' in line and in_cell:
                in_cell = False
                in_transistors = False

    return cells


def build_adjacency_matrices(cells):
    """
    Build PFET and NFET adjacency matrices for every cell.

    Matrix entry [i][j] = number of transistors connecting vertex i and j.
    (Parallel transistors → value > 1)

    Modifies cells in-place.
    """
    for cell in cells:
        n_p = len(cell.p_vertex)
        n_n = len(cell.n_vertex)
        cell.pmatrix = [[0] * n_p for _ in range(n_p)]
        cell.nmatrix = [[0] * n_n for _ in range(n_n)]

        # PFET matrix
        for gate, drain, source in cell.connect_list[0]:
            i = cell.p_vertex.index(drain)
            j = cell.p_vertex.index(source)
            cell.pmatrix[i][j] += 1
            cell.pmatrix[j][i] += 1

        # NFET matrix
        for gate, drain, source in cell.connect_list[1]:
            i = cell.n_vertex.index(drain)
            j = cell.n_vertex.index(source)
            cell.nmatrix[i][j] += 1
            cell.nmatrix[j][i] += 1
