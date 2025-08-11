# Optimizing Traffic Locality using Hilbert Curve Mapping

An analysis of using Hilbert space-filling curves to optimize network communication costs on high-performance computing (HPC) torus networks.

**Authors:** Nikolaos Lazaridis, Polyvios Pratikakis, Fabien Chaix  
*Department of Computer Science, University of Crete*

---

## 1. Overview

In high-performance computing (HPC) systems, communication efficiency is critically dependent on **traffic locality**—placing frequently communicating nodes physically close to one another. Poor locality leads to increased network latency and energy consumption, bottlenecking overall performance.

This project investigates the effectiveness of using **Hilbert space-filling curve mapping** to improve traffic locality on 2D and 3D torus networks. We developed a C++ analysis tool to model a torus architecture, process application communication patterns, and evaluate the quality of different node mappings. The evaluation compares the communication cost before and after reordering nodes with a Hilbert curve, using traffic-weighted Euclidean, Manhattan, and Chebyshev distance metrics.

Our analysis demonstrates that while Hilbert curves can significantly improve communication locality for certain traffic patterns, their effectiveness is not universal and depends heavily on the intrinsic communication structure of the application being mapped.

## 2. Key Concepts

### Torus Networks
A torus network is a popular interconnect topology in HPC systems. It arranges nodes in a multi-dimensional grid (e.g., 2D or 3D) with **wrap-around connections**. These links connect nodes on the periphery of the grid to nodes on the opposite side, creating a continuous network without physical edges. This is analogous to the surface of a donut.


*Figure 1: Illustration of a 3D Mesh Torus Network from the paper.*

### The Problem: Traffic Locality
- **Good Locality:** Nodes that communicate frequently are placed close together in the network. This minimizes the distance data has to travel.
- **Poor Locality:** Nodes that communicate frequently are far apart, leading to high-cost communication.

The goal of a locality-aware mapping is to arrange application tasks onto the physical nodes to minimize the total communication cost.

### The Proposed Solution: Hilbert Curve Mapping
A Hilbert curve is a continuous, fractal, space-filling curve that can "unroll" a multi-dimensional space (like a 3D grid) into a 1D line. Its key property is **locality preservation**: points that are close in the 1D sequence tend to be close in the original multi-dimensional space.

By traversing the torus grid along a Hilbert curve, we can generate a new, linear ordering of the nodes that aims to keep neighboring nodes together.

## 3. Methodology

Our C++ analysis tool implements the following pipeline:

1.  **Model a Torus Grid:** Define a torus of arbitrary dimensions (e.g., 2x2x4 or 2x8x16).
2.  **Parse Traffic Data:** Ingest a traffic matrix that specifies the communication volume between pairs of logical nodes for a given application (e.g., "CG", "SP").
3.  **Calculate Baseline Cost:** For the default "canonical" mapping (where logical node `i` is on physical node `i`), calculate the total communication cost. The cost is the sum of traffic-weighted distances for all communicating pairs:
    ```
    Total Communication Cost = Σ (Traffic(u,v) * TorusDistance(Pu, Pv))
    ```
    We calculate this cost using three torus-aware distance metrics:
    - **Manhattan (L1):** Hop count. `Δx + Δy + Δz`
    - **Euclidean (L2):** Straight-line distance. `sqrt(Δx² + Δy² + Δz²)`
    - **Chebyshev (L∞):** Maximum distance in any single dimension. `max(Δx, Δy, Δz)`

4.  **Apply Hilbert Mapping:**
    - Generate a Hilbert index for each physical coordinate in the torus grid.
    - Reorder the grid nodes according to their Hilbert index. This creates a new mapping of logical nodes to physical coordinates.

5.  **Calculate Hilbert-Mapped Cost:** Recalculate the `Total Communication Cost` using the new Hilbert-based node mapping.

6.  **Compare Results:** Analyze the percentage change in communication cost to determine if the Hilbert mapping improved or degraded traffic locality.

## 4. Key Findings

The effectiveness of Hilbert curve mapping is **highly dependent on the application's communication pattern**. It is not a universal solution.

### Positive Impact (e.g., "CG" Workload)
For workloads with strong inherent spatial locality (e.g., communication along diagonals or between nearby neighbors), **Hilbert mapping provides significant improvements.**

- The "CG" (Conjugate Gradient) workload saw a **reduction in total communication cost of up to 50%**.
- **Why?** The Hilbert curve's linearization aligns well with these localized communication patterns, effectively shortening the average path length for high-traffic pairs. This trend was also observed for `MiniAero`, `XHPCG`, and `XHPL` workloads.

### Negative Impact (e.g., "SP" Workload)
For workloads with more dispersed, irregular communication, **Hilbert mapping is detrimental, dramatically increasing communication costs.**

- The "SP" (Scalar Pentadiagonal) workload saw an **increase in total communication cost of over 140%** on the larger torus.
- **Why?** These patterns often rely on the torus's **wrap-around links** to create "shortcuts" between distant nodes. Hilbert reordering preserves general locality but disrupts these specific, long-range connections, forcing communication to travel a much longer path along the curve. This trend was also observed for `DT`, `EP`, `LU`, and `MG` workloads.

### Results Summary

The charts below show the percentage change in total and maximum traffic-weighted distance after applying Hilbert mapping. Negative values indicate improvement (cost reduction).


*(a) Percentage change in total communication cost.*


*(b) Percentage change in maximum single-pair communication cost.*

#### Example Comparison (Total Cost, 2x8x16 Torus):

| Benchmark  | Communication Pattern | Euclidean | Manhattan | Chebyshev | Outcome    |
| :--------- | :-------------------- | :-------: | :-------: | :-------: | :--------- |
| **CG_D_256** | Localized, Diagonal   | `-44.6%`  | `-32.2%`  | `-50.7%`  | **Improved** |
| **SP_D_256** | Dispersed, Irregular  | `+143.0%` | `+146.1%` | `+102.9%` | **Degraded** |
| **XHPL_512** | Localized             | `-49.5%`  | `-49.5%`  | `-52.6%`  | **Improved** |
| **LU_D_256** | Dispersed             | `+107.3%` | `+164.5%` | `+85.0%`  | **Degraded** |

---

## Steps to Compile and Run

1. mkdir build
2. cd build
3. cmake ..
4. make
5. make run_all (run the executable and generate the csv files)
6. python -m http.server "[port]"
7. open html index.html

---
