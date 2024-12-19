"""
Implementation of Prim's Minimum Spanning Tree algorithm.

This module provides functions to find the minimum spanning tree of a weighted
undirected graph using Prim's algorithm with a priority queue implementation.
"""

import heapq
from typing import List, Tuple, Dict


def create_adjacency_list(
    vertices: int, edges: List[List[float]]
) -> List[List[Tuple[int, float]]]:
    """
    Create an adjacency list representation of the graph.

    Args:
        vertices: Number of vertices in the graph
        edges: List of edges, where each edge is [u, v, weight]

    Returns:
        Adjacency list where each entry contains (vertex, weight) tuples
    """
    adj = [[] for _ in range(vertices)]
    for u, v, weight in edges:
        adj[int(u)].append((int(v), weight))
        adj[int(v)].append((int(u), weight))  # Add both directions for undirected graph
    return adj


def prim_mst(
    vertices: int, edges: List[List[float]]
) -> Tuple[float, List[List[float]]]:
    """
    Implements Prim's algorithm to find the Minimum Spanning Tree.

    The algorithm works by:
    1. Starting from vertex 0
    2. Using a priority queue to select minimum weight edges
    3. Adding vertices to the MST one at a time

    Args:
        vertices: Number of vertices in the graph
        edges: List of edges, where each edge is [u, v, weight]

    Returns:
        Tuple containing:
        - Total weight of the MST
        - List of edges in the MST, where each edge is [u, v, weight]
    """
    # Create adjacency list representation
    adj = create_adjacency_list(vertices, edges)

    # Initialize data structures
    visited = [False] * vertices
    mst_edges = []
    total_weight = 0

    # Priority queue entries: (weight, current_vertex, parent_vertex)
    pq = [(0, 0, -1)]  # Start with vertex 0

    while pq:
        weight, current, parent = heapq.heappop(pq)

        # Skip if vertex already visited
        if visited[current]:
            continue

        # Add vertex to MST
        visited[current] = True
        total_weight += weight

        # Add edge to MST if it has a valid parent
        if parent != -1:
            mst_edges.append([parent, current, weight])

        # Add edges to unvisited neighbors
        for neighbor, edge_weight in adj[current]:
            if not visited[neighbor]:
                heapq.heappush(pq, (edge_weight, neighbor, current))  # type: ignore[arg-type]

    print(f"Total weight of MST (Prim's): {total_weight:.6f}")
    return total_weight, mst_edges
