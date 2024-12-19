"""
Implementation of Kruskal's Minimum Spanning Tree algorithm.

This module provides a Graph class that implements Kruskal's algorithm
for finding the minimum spanning tree of a weighted undirected graph.
"""

from typing import List, Tuple


class Graph:
    """
    A class representing an undirected weighted graph and implementing Kruskal's algorithm.

    Attributes:
        vertices (int): Number of vertices in the graph
        graph (List[List[int, int, float]]): List of edges with weights [u, v, weight]
    """

    def __init__(self, vertices: int):
        """
        Initialize the graph with a given number of vertices.

        Args:
            vertices: Number of vertices in the graph
        """
        self.vertices = vertices
        self.graph: List[List[float]] = []

    def add_edge(self, u: int, v: int, weight: float) -> None:
        """
        Add an edge to the graph.

        Args:
            u: First vertex
            v: Second vertex
            weight: Weight of the edge
        """
        self.graph.append([u, v, weight])

    def find_parent(self, parent: List[int], i: int) -> int:
        """
        Find the parent of a vertex using path compression.

        Args:
            parent: List storing parent information
            i: Vertex to find parent for

        Returns:
            Parent vertex
        """
        if parent[i] != i:
            parent[i] = self.find_parent(parent, parent[i])
        return parent[i]

    def union(self, parent: List[int], rank: List[int], x: int, y: int) -> None:
        """
        Union of two subsets using rank to keep tree height small.

        Args:
            parent: List storing parent information
            rank: List storing rank information
            x: Root of first set
            y: Root of second set
        """
        if rank[x] < rank[y]:
            parent[x] = y
        elif rank[x] > rank[y]:
            parent[y] = x
        else:
            parent[y] = x
            rank[x] += 1

    def kruskal_mst(self) -> List[List[float]]:
        """
        Implements Kruskal's algorithm to find the Minimum Spanning Tree.

        The algorithm works by:
        1. Sorting edges by weight
        2. Taking edges in ascending order of weight
        3. Adding edge if it doesn't create a cycle

        Returns:
            List of edges in the MST, where each edge is [u, v, weight]
        """
        result = []  # Stores the MST edges
        i = 0  # Index for sorted edges
        e = 0  # Number of edges in MST

        # Sort edges by weight
        self.graph = sorted(self.graph, key=lambda item: item[2])

        # Initialize parent and rank arrays
        parent = list(range(self.vertices))
        rank = [0] * self.vertices

        # Build MST by adding edges that don't create cycles
        while e < self.vertices - 1 and i < len(self.graph):
            u, v, weight = self.graph[i]
            i += 1

            # Find parents (sets) of vertices
            x = self.find_parent(parent, int(u))
            y = self.find_parent(parent, int(v))

            # Add edge if it doesn't create cycle
            if x != y:
                e += 1
                result.append([u, v, weight])
                self.union(parent, rank, x, y)

        total_weight = sum(weight for _, _, weight in result)
        print(f"Total weight of MST (Kruskal's): {total_weight:.6f}")
        return result
