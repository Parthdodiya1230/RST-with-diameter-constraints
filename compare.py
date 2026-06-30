#!/usr/bin/env python3
"""
compare.py — Side-by-side comparison of naive vs. constrained RST reconfiguration.

Usage:
    python3 compare.py                                        # default file names
    python3 compare.py output_constrained.txt output_naive.txt
    python3 compare.py --help

Reads the structured output files produced by:
    rst_constrained input.txt output_constrained.txt
    rst_naive       input.txt output_naive.txt

Prints:
    1. Summary statistics (steps, violations, max/avg diameter)
    2. Side-by-side step table with diameter and colour-coded violations
    3. ASCII bar-chart of diameter over steps
    4. Verdict for dynamic-network use-case
"""

import sys
import os
from collections import defaultdict

# ──────────────────────────────────────────────────────────────
# ANSI colours
# ──────────────────────────────────────────────────────────────
try:
    import shutil
    _cols = shutil.get_terminal_size().columns
    USE_COLOUR = sys.stdout.isatty()
except Exception:
    _cols = 100
    USE_COLOUR = False

def _c(code, s):
    return f"\033[{code}m{s}\033[0m" if USE_COLOUR else s

def red(s):    return _c("31;1", s)
def green(s):  return _c("32;1", s)
def yellow(s): return _c("33;1", s)
def cyan(s):   return _c("36;1", s)
def bold(s):   return _c("1", s)
def dim(s):    return _c("2", s)

# ──────────────────────────────────────────────────────────────
# PARSER
# ──────────────────────────────────────────────────────────────
def parse_output(filepath):
    """Parse an output_*.txt file.  Returns a dict."""
    if not os.path.exists(filepath):
        print(f"[ERROR] File not found: {filepath}")
        return None

    data = {
        "method": "?",
        "result": "?",
        "d_bound": None,
        "n": None,
        "m": None,
        "graph_edges": [],
        "steps": [],      # list of dicts
    }

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("METHOD"):
                data["method"] = line.split()[1]
            elif line.startswith("RESULT"):
                data["result"] = line.split()[1]
            elif line.startswith("D_BOUND"):
                data["d_bound"] = int(line.split()[1])
            elif line.startswith("N_VERTICES"):
                data["n"] = int(line.split()[1])
            elif line.startswith("N_EDGES"):
                data["m"] = int(line.split()[1])
            elif line.startswith("GRAPH_EDGES"):
                data["graph_edges"] = line.split()[1:]
            elif line.startswith("STEP"):
                # STEP i LABEL lbl DIAM d VIOLATES v [ADDED a REMOVED r] EDGES e1 e2 ...
                tokens = line.split()
                step = {}
                i = 0
                while i < len(tokens):
                    t = tokens[i]
                    if t == "STEP":   step["idx"]      = int(tokens[i+1]);  i += 2
                    elif t == "LABEL":
                        step["label"] = tokens[i+1];  i += 2
                    elif t == "DIAM":   step["diam"]   = int(tokens[i+1]);  i += 2
                    elif t == "VIOLATES": step["violates"] = int(tokens[i+1]); i += 2
                    elif t == "ADDED":  step["added"]  = tokens[i+1];       i += 2
                    elif t == "REMOVED":step["removed"]= tokens[i+1];       i += 2
                    elif t == "EDGES":
                        step["edges"] = tokens[i+1:]; break
                    else:
                        i += 1
                data["steps"].append(step)

    return data

# ──────────────────────────────────────────────────────────────
# FORMATTING HELPERS
# ──────────────────────────────────────────────────────────────
def fmt_diam(diam, d_bound, width=4):
    s = str(diam).rjust(width)
    if diam > d_bound:
        return red(f"{s} !")
    return green(s + "  ")

def fmt_label(label, max_w=18):
    return label[:max_w].ljust(max_w)

def fmt_edges(edges, max_e=5):
    shown = edges[:max_e]
    rest  = len(edges) - max_e
    s = " ".join(shown)
    if rest > 0:
        s += f" … +{rest}"
    return s

# ──────────────────────────────────────────────────────────────
# STATS
# ──────────────────────────────────────────────────────────────
def stats(data):
    d = data["d_bound"]
    steps  = data["steps"]
    diams  = [s["diam"] for s in steps]
    viols  = [s for s in steps if s.get("violates",0)]
    return {
        "n_trees":   len(steps),
        "n_steps":   len(steps) - 1,
        "violations":len(viols),
        "max_diam":  max(diams) if diams else 0,
        "avg_diam":  sum(diams)/len(diams) if diams else 0,
        "d_bound":   d,
    }

# ──────────────────────────────────────────────────────────────
# BAR CHART
# ──────────────────────────────────────────────────────────────
def bar_chart(steps_c, steps_n, d_bound, width=60):
    max_diam = max(
        max(s["diam"] for s in steps_c) if steps_c else 0,
        max(s["diam"] for s in steps_n) if steps_n else 0,
        d_bound
    )

    def bar(diam, w=width):
        filled = round(diam / max_diam * w)
        b = "█" * filled
        if diam > d_bound:
            return red(b)
        return green(b)

    print(bold("\n  Diameter per step — bar chart"))
    print(dim(f"  (max = {max_diam}, bound d = {d_bound}, █ = proportional to diameter)"))
    print()

    n_c = len(steps_c)
    n_n = len(steps_n)
    rows = max(n_c, n_n)

    # Header
    C_W = 35; N_W = 35
    print(f"  {'Step':<5}  {cyan('Algorithm 3 (constrained)'):<{C_W+10}}  {yellow('Naive'):<{N_W+10}}")
    print(f"  {'─'*5}  {'─'*C_W}  {'─'*N_W}")

    for i in range(rows):
        c_part = n_part = ""
        if i < n_c:
            s = steps_c[i]
            d = s["diam"]
            c_part = f"d={d:2d} {bar(d, 20)}"
        if i < n_n:
            s = steps_n[i]
            d = s["diam"]
            n_part = f"d={d:2d} {bar(d, 20)}"
        print(f"  {i:<5}  {c_part:<50}  {n_part}")

    print()
    # Legend
    print(f"  {green('█')} diam ≤ d = {d_bound}   {red('█')} diam > d  (constraint violated)")
    print()

# ──────────────────────────────────────────────────────────────
# SIDE-BY-SIDE TABLE
# ──────────────────────────────────────────────────────────────
def side_by_side_table(c_data, n_data):
    d = c_data["d_bound"]
    c_steps = c_data["steps"]
    n_steps = n_data["steps"]
    rows = max(len(c_steps), len(n_steps))

    COL_W = 52

    # Header
    sep = "─" * COL_W
    print(bold("  Step-by-step comparison"))
    print()
    header = f"  {'#':<4}  {cyan('Algorithm 3 (constrained)'):<{COL_W}}  {yellow('Naive')}"
    print(header)
    print(f"  {'─'*4}  {sep}  {sep}")

    for i in range(rows):
        c_str = n_str = dim("(no step)")

        if i < len(c_steps):
            s = c_steps[i]
            diam_str = fmt_diam(s["diam"], d)
            op = ""
            if s.get("added","none") != "none":
                op += f" +{s['added']}"
            if s.get("removed","none") != "none":
                op += f" -{s['removed']}"
            c_str = f"[{diam_str}] {s['label'][:20]:<20}{op}"

        if i < len(n_steps):
            s = n_steps[i]
            diam_str = fmt_diam(s["diam"], d)
            op = ""
            if s.get("added","none") != "none":
                op += f" +{s['added']}"
            if s.get("removed","none") != "none":
                op += f" -{s['removed']}"
            n_str = f"[{diam_str}] {s['label'][:20]:<20}{op}"

        print(f"  {i:<4}  {c_str:<{COL_W+20}}  {n_str}")

    print(f"\n  Legend: [diam] {green('= ok ≤ d')}   [diam {red('!')}] {red('= VIOLATION > d')}")
    print()

# ──────────────────────────────────────────────────────────────
# VERDICT
# ──────────────────────────────────────────────────────────────
def verdict(c_st, n_st):
    d     = c_st["d_bound"]
    c_vio = c_st["violations"]
    n_vio = n_st["violations"]
    c_steps = c_st["n_steps"]
    n_steps = n_st["n_steps"]

    print(bold("━" * 68))
    print(bold("  COMPARISON VERDICT"))
    print(bold("━" * 68))
    print()

    # Algorithm 3 summary
    if c_vio == 0:
        print(f"  {green('✓')} Algorithm 3   — {c_steps} steps, {green('0 violations')}, "
              f"max diam = {c_st['max_diam']} ≤ d={d}")
        print(f"      Every intermediate spanning tree satisfies diam ≤ {d}.")
        print(f"      {green('SAFE')} for dynamic networks with a hop-count SLA of {d}.")
    else:
        print(f"  {yellow('~')} Algorithm 3   — {c_steps} steps, {yellow(str(c_vio)+' violations')}, "
              f"max diam = {c_st['max_diam']}")
        print(f"      (Unexpected violations — check input/d_bound)")

    print()

    # Naive summary
    if n_vio == 0:
        print(f"  {yellow('~')} Naive approach — {n_steps} steps, 0 violations (by luck, no guarantee)")
        print(f"      Naive reached Tt without violating d={d} on THIS input.")
        print(f"      {yellow('WARNING')}: No mathematical guarantee — a different input or edge")
        print(f"      ordering could easily produce violations.")
    else:
        print(f"  {red('✗')} Naive approach — {n_steps} steps, {red(str(n_vio)+' violation(s)')}, "
              f"max diam = {n_st['max_diam']} (peak +{n_st['max_diam']-d} over bound)")
        print(f"      {n_vio} intermediate tree(s) have diameter > {d}.")
        print(f"      {red('UNSAFE')}: during these steps, network paths exceed the SLA.")
        print(f"      Packets in a live network would traverse {n_st['max_diam']}-hop paths")
        print(f"      instead of the promised ≤ {d} hops.")

    print()
    print(bold("  Why this matters for dynamic networks:"))
    print(f"  Every step in the sequence is a real network state — the topology")
    print(f"  during a live reconfiguration. If diameter > {d} at any step,")
    print(f"  routing tables see longer paths, QoS guarantees break, and")
    print(f"  latency-sensitive traffic silently degrades.")
    print()

    if c_vio == 0 and n_vio > 0:
        saved = n_steps - c_steps
        print(f"  Algorithm 3 costs {abs(saved)} {'more' if saved<0 else 'fewer'} steps but "
              f"{green('guarantees safety')} at every step.")
    elif c_vio == 0 and n_vio == 0:
        print(f"  On this instance naive was also safe — but only by coincidence.")
        print(f"  Algorithm 3 provides a provable guarantee for ALL valid inputs.")

    print()
    print(bold("━" * 68))

# ──────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────
def main():
    args = sys.argv[1:]
    if "--help" in args or "-h" in args:
        print(__doc__)
        return

    c_file = args[0] if len(args) >= 1 else "output_constrained.txt"
    n_file = args[1] if len(args) >= 2 else "output_naive.txt"

    print()
    print(bold("━" * 68))
    print(bold("  RST RECONFIGURATION — CONSTRAINED vs. NAIVE COMPARISON"))
    print(bold("━" * 68))
    print(f"  Constrained: {c_file}")
    print(f"  Naive:       {n_file}")
    print()

    c_data = parse_output(c_file)
    n_data = parse_output(n_file)

    if c_data is None or n_data is None:
        print("Run rst_constrained and rst_naive first to generate the output files.")
        return

    # Consistency check
    if c_data["d_bound"] != n_data["d_bound"]:
        print(f"  {yellow('WARNING')}: d_bound mismatch: "
              f"constrained={c_data['d_bound']}, naive={n_data['d_bound']}")

    d = c_data["d_bound"]
    print(f"  Graph      : {c_data['n']} vertices, {c_data['m']} edges")
    print(f"  Diameter d : {d}")
    print(f"  Result     : Algorithm 3 → {green(c_data['result']) if c_data['result']=='YES' else red(c_data['result'])}")
    print()

    c_st = stats(c_data)
    n_st = stats(n_data)

    # Summary row
    print(bold("  ┌─────────────────────┬─────────────────┬─────────────────┐"))
    print(bold("  │ Metric              │ Algorithm 3     │ Naive           │"))
    print(bold("  ├─────────────────────┼─────────────────┼─────────────────┤"))
    def row(label, cv, nv, good_fn=None):
        cs = str(cv).ljust(15)
        ns = str(nv).ljust(15)
        if good_fn:
            cs = good_fn(cv, "c").ljust(15 + 10)
            ns = good_fn(nv, "n").ljust(15 + 10)
        print(f"  │ {label:<19} │ {cs} │ {ns} │")

    def viol_fmt(v, side):
        s = str(v)
        if v == 0: return green(s + " ✓")
        return red(s + " ✗")

    row("Steps",            c_st["n_steps"],  n_st["n_steps"])
    row("Trees in sequence",c_st["n_trees"],  n_st["n_trees"])
    row("Violations",       c_st["violations"],n_st["violations"], viol_fmt)
    row("Max diameter",     c_st["max_diam"], n_st["max_diam"])
    row(f"Avg diameter",    f"{c_st['avg_diam']:.1f}", f"{n_st['avg_diam']:.1f}")
    row("d_bound",          c_st["d_bound"],  n_st["d_bound"])
    print(bold("  └─────────────────────┴─────────────────┴─────────────────┘"))
    print()

    side_by_side_table(c_data, n_data)
    bar_chart(c_data["steps"], n_data["steps"], d)
    verdict(c_st, n_st)


if __name__ == "__main__":
    main()