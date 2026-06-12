# Standard Cell Layout Generation using Euler Trail and A* Routing

This project implements a standard cell layout generation pipeline for CMOS digital circuits.
Given a SPICE netlist, it automatically determines the optimal transistor placement order using the Euler trail algorithm, then connects pins using A* global routing.

This work is based on UROP research conducted at Kookmin University under Prof. Heechun Park,
and was presented at the IEIE SoC Conference 2023 (Outstanding Presentation Paper Award).

---

## Background

In CMOS standard cell layout, transistors are arranged in a single row for both PFET (pull-up)
and NFET (pull-down) networks. The key insight is:

> If adjacent transistors share a drain/source net, no extra wire is needed between them.
> This minimizes cell width and routing complexity.

Finding such an optimal ordering is equivalent to finding an **Euler trail** on a graph where:
- **Vertex** = drain/source net (VDD, VSS, Z, net_0, ...)
- **Edge** = transistor (labeled by its gate net: A1, A2, Z_neg, ...)

---

## Pipeline

```
Nangate_15nm.sp (SPICE netlist)
        │
        ▼ parser.py
Parse transistor connections → build CellGraph (vertex, edge, adjacency matrix)
        │
        ▼ euler_trail.py
Find Euler trail → determine optimal PFET/NFET transistor ordering
        │
        ▼ main2.py
Write placement_result.txt
        │
        ▼ euler_trail.py (routing_all)
Calculate pin coordinates
        │
        ▼ main2.py
Write routing_pin_result.txt
        │
        ▼ astar.py
Connect pins via A* global routing (obstacle-aware)
```

---

## File Structure

```
.
├── main2.py               # Entry point — runs full pipeline
├── parser.py              # SPICE netlist parser → CellGraph objects
├── euler_trail.py         # Euler trail algorithm + pin coordinate generation
├── astar.py               # A* global routing
├── Nangate_15nm.sp        # Input: SPICE netlist (NanGate 15nm OCL)
├── placement_result.txt   # Output: transistor placement order per cell
└── routing_pin_result.txt # Output: pin coordinates for routing
```

---

## How to Run

**Requirements:** Python 3.x (no external libraries needed)

```bash
# Step 1: Euler trail placement
python main2.py
# → generates placement_result.txt, routing_pin_result.txt

# Step 2: A* routing
python astar.py
# → connects pins in routing_pin_result.txt
```

---

## Algorithm Details

### Euler Trail (transistor placement)

A CMOS cell is modeled as a graph:
- Each transistor becomes an **edge** (labeled by its gate net)
- Each drain/source net becomes a **vertex**

An Euler trail visits every edge exactly once — this maps directly to placing every transistor
exactly once in a sequence where adjacent transistors always share a net.

**Odd-degree handling:**
An Euler trail exists only when 0 or 2 vertices have odd degree.
When more odd-degree vertices exist, dummy edges are inserted to pair them up,
then removed from the final result.

**PFET ↔ NFET matching:**
After finding a PFET Euler trail, the algorithm checks whether the same gate ordering
is compatible with the NFET network. Only matching combinations are kept as valid placements.

### A* Routing (pin connection)

After placement, pins sharing the same net must be physically connected.
A* finds the shortest path between pins on a routing grid, avoiding already-routed nets as obstacles.

- **g(n):** actual cost from start pin to current cell
- **h(n):** Manhattan distance from current cell to target pin
- **f(n) = g(n) + h(n):** total estimated cost

---

## Output Example (AND2_X1)

**placement_result.txt:**
```
Name : AND2_X1
PMOS : Z_neg VDD Z_neg VDD
NMOS : net_0 VSS Z_neg Z
Gate : A2 A1 Z_neg
```

**routing_pin_result.txt:**
```
Name :AND2_X1
7 10 2
2
0 6 0 0 7 0 0 8 0 0 9 0 1 0 0 2 0 0 3 0 0
4 6 0 4 7 0 4 8 0 4 9 0 6 0 0 6 1 0 6 2 0 6 3 0 1 4 0 1 5 0
```

---

## Known Limitations & Future Work

### X2 and above cells are excluded
**Why**: Parallel transistors in X2+ cells cause combinatorial explosion in permutation search.
X2 has 2× transistors in parallel, so permutations grow as N! — e.g. BUF_X12 would require 12! = 479,001,600 iterations, exceeding available memory.

**Fix**: Merge parallel transistors into a single representative before Euler trail search, then duplicate the result. This reduces the problem back to X1 complexity regardless of drive strength.

### A* routing may fail in congested grids
**Why**: When many nets are already routed, a pin may become completely surrounded by obstacles and unreachable.

**Fix**: Implement rip-up and reroute — temporarily remove blocking nets, reroute the failed net, then reroute the removed nets with updated constraints.

### LEF/DEF physical constraints not enforced
**Why**: Pin coordinates are computed from placement order only, without referencing actual metal layer rules (width, spacing, via enclosure).

**Fix**: Integrate tech.lef design rules into the routing cost function so the A* path respects real DRC constraints.

---

## Environment

- Language: Python 3.13
- PDK: NanGate 15nm OCL
- Input format: SPICE netlist (.sp)
- No external dependencies required
