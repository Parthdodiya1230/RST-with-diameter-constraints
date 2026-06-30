#!/usr/bin/env python3
import subprocess
import sys
import os
import csv
import time
import random
from collections import defaultdict, deque

# 1. Custom sizes: 2 through 10, then 20, 30, 40, 50
SIZES = list(range(4, 11)) + [15, 20, 30, 40, 50]  
N_SEEDS = 2
TIMEOUT_SECONDS = 100.0  # 2. Hard 100-second timeout

SOURCES = {
    "optimized":   ("rst_cons.cpp",   "rst_opt"),
    "unoptimized": ("rst_unopt.cpp",  "rst_unopt"),
}
OUTPUT_CSV  = "results.csv"
OUTPUT_PLOT = "scalability.png"

def compile_binary(src, binary):
    if os.path.exists(binary): return True
    print(f"Compiling {src} -> {binary} ...", end=" ", flush=True)
    r = subprocess.run(["g++", "-O3", "-std=c++17", "-o", binary, src], capture_output=True, text=True)
    if r.returncode != 0:
        print("FAILED")
        return False
    print("OK")
    return True

def generate_guaranteed_yes_instance(n, seed):
    random.seed(seed)
    edge_set = set()
    
    for i in range(1, n): edge_set.add((min(i, i + 1), max(i, i + 1)))

    target = min(n * (n - 1) // 2, 3 * n)
    attempts = 0
    while len(edge_set) < target and attempts < 10000:
        u, v = random.randint(1, n), random.randint(1, n)
        if u != v: edge_set.add((min(u, v), max(u, v)))
        attempts += 1

    edges = sorted(edge_set)
    m = len(edges)
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    def bfs_tree(root):
        visited, queue, tree_edges = {root}, deque([root]), []
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
    d = n 

    lines = [f"{n} {m}"]
    for u, v in edges: lines.append(f"{u} {v}")
    for u, v in ts: lines.append(f"{u} {v}")
    for u, v in tt: lines.append(f"{u} {v}")
    lines.append(str(d))
    return "\n".join(lines)

def run_binary(binary, input_file):
    t0 = time.perf_counter()
    try:
        # Pass the timeout limit directly to the OS process
        r = subprocess.run([f"./{binary}", input_file], capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        output = r.stdout.upper()
        if "YES" in output: answer = "YES"
        elif "NO" in output: answer = "NO"
        else: answer = "ERR"
        return answer, elapsed_ms
    except subprocess.TimeoutExpired:
        # If it hits 50s, return the specific TIMEOUT flag
        return "TIMEOUT", TIMEOUT_SECONDS * 1000.0

def median(lst):
    s = sorted(lst)
    n = len(s)
    return s[n // 2] if n % 2 == 1 else (s[n//2 - 1] + s[n//2]) / 2.0

def main():
    print("=" * 60)
    print(f"RST Benchmark: N=2..10, 20..50 | Limit: {TIMEOUT_SECONDS}s")
    print("=" * 60)

    for label, (src, binary) in SOURCES.items():
        if not compile_binary(src, binary): sys.exit(1)

    os.makedirs("inputs", exist_ok=True)
    print(f"\n[3] Running benchmark ({N_SEEDS} seeds per n)...")
    print(f"    {'n':>5}  {'opt_ms':>12}  {'unopt_ms':>12}  {'speedup':>10}  {'answer_opt':>11}  {'answer_unopt':>13}")
    print("    " + "-" * 70)

    rows = []
    
    # Track if a binary has timed out so we can skip larger graphs
    opt_timed_out = False
    unopt_timed_out = False

    for n in SIZES:
        opt_times, unopt_times, opt_answers, unopt_answers = [], [], [], []

        for seed in range(1, N_SEEDS + 1):
            inp_str = generate_guaranteed_yes_instance(n, seed + n)
            tmp = f"inputs/tmp_n{n}_s{seed}.txt"
            with open(tmp, "w") as f: f.write(inp_str)

            # OPTIMIZED RUN
            if opt_timed_out:
                opt_answers.append("TIMEOUT")
                opt_times.append(TIMEOUT_SECONDS * 1000.0)
            else:
                ans, ms = run_binary(SOURCES["optimized"][1], tmp)
                opt_answers.append(ans)
                opt_times.append(ms)
                if ans == "TIMEOUT": opt_timed_out = True

            # UNOPTIMIZED RUN
            if unopt_timed_out:
                unopt_answers.append("TIMEOUT")
                unopt_times.append(TIMEOUT_SECONDS * 1000.0)
            else:
                ans, ms = run_binary(SOURCES["unoptimized"][1], tmp)
                unopt_answers.append(ans)
                unopt_times.append(ms)
                if ans == "TIMEOUT": unopt_timed_out = True

            os.remove(tmp)

        opt_med = median(opt_times)
        unopt_med = median(unopt_times)
        
        # Calculate speedup (handle division by zero or timeouts)
        if opt_answers[0] == "TIMEOUT" and unopt_answers[0] == "TIMEOUT": speedup = 1.0
        elif opt_med > 0: speedup = unopt_med / opt_med
        else: speedup = 0.0

        opt_ans_str = opt_answers[0] if opt_answers else "?"
        unopt_ans_str = unopt_answers[0] if unopt_answers else "?"

        print(f"    {n:>5}  {opt_med:>12.2f}  {unopt_med:>12.2f}  {speedup:>9.2f}x  {opt_ans_str:>11}  {unopt_ans_str:>13}")

        rows.append({
            "n": n, "opt_ms": round(opt_med, 3), "unopt_ms": round(unopt_med, 3),
            "speedup": round(speedup, 4), "answer_opt": opt_ans_str, "answer_unopt": unopt_ans_str,
        })

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ns = [r["n"] for r in rows]
        opt_ms = [r["opt_ms"] for r in rows]
        unopt_ms = [r["unopt_ms"] for r in rows]

        plt.figure(figsize=(10, 6))
        plt.plot(ns, unopt_ms, "s--", color="#d62728", linewidth=3, markersize=8, label="Without OPT (Evaluates Full Graph)")
        plt.plot(ns, opt_ms, "o-", color="#1f77b4", linewidth=3, markersize=8, label="With OPT (Early-Exit BFS)")
        
        # Draw a clear red line representing the timeout ceiling
        plt.axhline(TIMEOUT_SECONDS * 1000.0, color='red', linestyle='--', linewidth=2, label="100 Second Timeout Ceiling")
        
        plt.yscale('log')
        plt.xlabel("Number of vertices (N)", fontsize=12, fontweight='bold')
        plt.ylabel("Execution Time (ms)", fontsize=12, fontweight='bold')
        plt.title(f"Execution Time vs Network Size (Max {TIMEOUT_SECONDS}s)", fontsize=14, fontweight='bold')
        plt.legend(fontsize=12, loc='lower right')
        plt.grid(True, which="both", alpha=0.3)

        # Better X-axis ticks to show the gap
        plt.xticks([2, 5, 10, 20, 30, 40, 50])

        plt.tight_layout()
        plt.savefig(OUTPUT_PLOT, dpi=200)
        print(f"\nSaved plot to {OUTPUT_PLOT}")

    except Exception as e:
        print(f"Plot error: {e}")

if __name__ == "__main__":
    main()