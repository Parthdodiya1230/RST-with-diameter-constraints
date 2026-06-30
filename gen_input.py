#!/usr/bin/env python3
import random
import sys
import os
import networkx as nx
from collections import defaultdict, deque

def gen_instance(n, seed=42, density_factor=3):
    random.seed(seed)
    edge_set = set()
    
    # Backbone for connectivity
    for i in range(1, n):
        u, v = i, i + 1
        edge_set.add((min(u, v), max(u, v)))

    # Extra edges
    target = min(n * (n - 1) // 2, density_factor * n)
    attempts = 0
    while len(edge_set) < target and attempts < 200000:
        u = random.randint(1, n)
        v = random.randint(1, n)
        if u != v:
            edge_set.add((min(u, v), max(u, v)))
        attempts += 1

    edges = sorted(edge_set)
    m = len(edges)
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    def bfs_tree(root):
        visited = {root}
        queue = deque([root])
        tree_edges = []
        while queue:
            cur = queue.popleft()
            for nb in sorted(adj[cur]):
                if nb not in visited:
                    visited.add(nb)
                    tree_edges.append((cur, nb))
                    queue.append(nb)
        return tree_edges

    ts = bfs_tree(1)
    tt = bfs_tree(n)

    # THE FIX: Set the limit to the max possible nodes to guarantee a YES instance.
    d = n 

    lines = [f"{n} {m}"]
    for u, v in edges: lines.append(f"{u} {v}")
    for u, v in ts: lines.append(f"{u} {v}")
    for u, v in tt: lines.append(f"{u} {v}")
    lines.append(str(d))
    return "\n".join(lines)

def main():
    if len(sys.argv) > 1:
        sizes = [int(x) for x in sys.argv[1:]]
    else:
        sizes = list(range(1, 16, 1))   

    os.makedirs("inputs", exist_ok=True)
    print(f"Generating GUARANTEED YES inputs for n = {sizes}")
    
    for n in sizes:
        fname = f"inputs/input_n{n:04d}.txt"
        content = gen_instance(n, seed=42+n)
        with open(fname, "w") as f:
            f.write(content)
        m = int(content.split("\n")[0].split()[1])
        print(f"  n={n:4d}  m={m:5d}  d={n:3d}  -> {fname}")
    print("Done.")

if __name__ == "__main__":
    main()