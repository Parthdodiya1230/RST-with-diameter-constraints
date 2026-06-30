import sys
import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def parse_input(filepath):
    """Parses input.txt to extract the base graph, Ts, and Tt."""
    tokens = []
    with open(filepath, 'r') as f:
        for line in f:
            if '#' in line:
                line = line.split('#')[0]
            tokens.extend(line.split())
    
    if not tokens:
        return None
        
    n = int(tokens[0])
    m = int(tokens[1])
    ptr = 2
    
    graph_edges = []
    for _ in range(m):
        graph_edges.append((int(tokens[ptr]), int(tokens[ptr+1])))
        ptr += 2
        
    Ts_edges = []
    for _ in range(n-1):
        Ts_edges.append((int(tokens[ptr]), int(tokens[ptr+1])))
        ptr += 2
        
    Tt_edges = []
    for _ in range(n-1):
        Tt_edges.append((int(tokens[ptr]), int(tokens[ptr+1])))
        ptr += 2
        
    d = float(tokens[ptr])
    
    return n, m, graph_edges, Ts_edges, Tt_edges, d

def parse_output(filepath):
    """Parses output_constrained.txt to extract the steps."""
    data = {"d_bound": 0, "steps": []}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts: continue
            
            if parts[0] == "D_BOUND":
                data["d_bound"] = int(parts[1])
            elif parts[0] == "STEP":
                try:
                    step = {}
                    step["idx"] = int(parts[1])
                    step["label"] = parts[parts.index("LABEL")+1]
                    step["diam"] = int(parts[parts.index("DIAM")+1])
                    step["violates"] = int(parts[parts.index("VIOLATES")+1])
                    
                    add_idx = parts.index("ADDED") + 1
                    step["added"] = tuple(map(int, parts[add_idx].split('-'))) if parts[add_idx] != "none" else None
                    
                    rem_idx = parts.index("REMOVED") + 1
                    step["removed"] = tuple(map(int, parts[rem_idx].split('-'))) if parts[rem_idx] != "none" else None
                    
                    edges_idx = parts.index("EDGES") + 1
                    step["edges"] = [tuple(map(int, e.split('-'))) for e in parts[edges_idx:]]
                    data["steps"].append(step)
                except ValueError:
                    continue
    return data

def draw_graph_state(ax, G_base, pos, tree_edges, added_edge, removed_edge, title_text, title_color="black", title_size=15):
    """Helper function to draw a single frame/state."""
    ax.clear()
    
    # Draw background full graph
    nx.draw_networkx_edges(G_base, pos, ax=ax, edge_color='lightgray', alpha=0.6, width=2.0)
    
    # Draw current spanning tree
    current_tree = [e for e in tree_edges if e != removed_edge]
    nx.draw_networkx_edges(G_base, pos, ax=ax, edgelist=current_tree, edge_color='dodgerblue', width=4)
    
    # Draw nodes
    nx.draw_networkx_nodes(G_base, pos, ax=ax, node_color='white', edgecolors='black', node_size=600, linewidths=1.5)
    nx.draw_networkx_labels(G_base, pos, ax=ax, font_size=11, font_weight="bold", font_family='sans-serif')
    
    # Highlight added edge (Green)
    if added_edge:
        nx.draw_networkx_edges(G_base, pos, ax=ax, edgelist=[added_edge], edge_color='limegreen', width=5)
        
    # Highlight removed edge (Red, Dashed)
    if removed_edge:
        nx.draw_networkx_edges(G_base, pos, ax=ax, edgelist=[removed_edge], edge_color='crimson', width=5, style='dashed')
        
    ax.set_title(title_text, fontsize=title_size, fontweight='bold', color=title_color, pad=15)
    ax.axis('off')

def main():
    input_file = "input.txt" if len(sys.argv) < 2 else sys.argv[1]
    output_file = "output_constrained.txt" if len(sys.argv) < 3 else sys.argv[2]
    
    print(f"Parsing {input_file}...")
    in_data = parse_input(input_file)
    if not in_data:
        print("Failed to read input file.")
        return
        
    n, m, graph_edges, Ts_edges, Tt_edges, d = in_data
    
    print(f"Parsing {output_file}...")
    out_data = parse_output(output_file)
    if not out_data["steps"]:
        print("Failed to read steps from output file.")
        return

    # Build the base graph
    G_base = nx.Graph()
    G_base.add_edges_from(graph_edges)
    
    # Generate Layout
    try:
        import scipy
        pos = nx.kamada_kawai_layout(G_base)
    except ImportError:
        pos = nx.spring_layout(G_base, seed=42) 

    # =================================================================
    # ANTI-OVERLAP NUDGE
    # =================================================================
    if 11 in pos and 8 in pos and 2 in pos:
        pos[11][0] -= 0.15  
        pos[11][1] -= 0.15  
    # =================================================================
    
    fig, ax = plt.subplots(figsize=(9, 7))
    fig.subplots_adjust(top=0.82) 
    
    # ---------------------------------------------------------
    # 1. Generate Static Images
    # ---------------------------------------------------------
    print("Generating Graph.png...")
    ax.clear()
    nx.draw_networkx_edges(G_base, pos, ax=ax, edge_color='gray', width=1.5)
    nx.draw_networkx_nodes(G_base, pos, ax=ax, node_color='white', edgecolors='black', node_size=600)
    nx.draw_networkx_labels(G_base, pos, ax=ax, font_size=11, font_weight="bold")
    ax.set_title("Full Network Graph", fontsize=24, fontweight='bold') 
    ax.axis('off')
    plt.savefig("Graph.png", dpi=150)
    
    print("Generating Ts.png...")
    draw_graph_state(ax, G_base, pos, Ts_edges, None, None, f"Initial Spanning Tree (Ts)\nDiameter Constraint: d={int(d)}", title_size=24)
    plt.savefig("Ts.png", dpi=150)

    print("Generating Tt.png...")
    draw_graph_state(ax, G_base, pos, Tt_edges, None, None, f"Target Spanning Tree (Tt)\nDiameter Constraint: d={int(d)}", title_size=24)
    plt.savefig("Tt.png", dpi=150)
    
    # ---------------------------------------------------------
    # 2. Generate Animated Sequence (MP4 VIDEO)
    # ---------------------------------------------------------
    print("Generating sequence_animation.mp4 (Encoding video)...")
    
    anim_frames = []
    
    first_step = out_data["steps"][0]
    anim_frames.append({
        "edges": first_step["edges"],
        "added": None,
        "removed": None,
        "title": f"Step 0: Initial Position\nDiameter: {first_step['diam']} (Limit: {out_data['d_bound']})",
        "color": "black"
    })
    
    for i in range(1, len(out_data["steps"])):
        prev_step = out_data["steps"][i-1]
        curr_step = out_data["steps"][i]
        
        status = "VIOLATION (Exceeds d)" if curr_step["violates"] else "SAFE (≤ d)"
        base_color = "crimson" if curr_step["violates"] else "black"
        
        if curr_step["added"]:
            anim_frames.append({
                "edges": prev_step["edges"], 
                "added": curr_step["added"], 
                "removed": None,
                "title": f"Step {curr_step['idx']}a: Forming Cycle\nAdding: {curr_step['added']} (Green)",
                "color": "black"
            })
        
        if curr_step["added"] and curr_step["removed"]:
            anim_frames.append({
                "edges": prev_step["edges"], 
                "added": curr_step["added"], 
                "removed": curr_step["removed"], 
                "title": f"Step {curr_step['idx']}b: Breaking Cycle\nRemoving: {curr_step['removed']} (Red)",
                "color": "black"
            })
            
        anim_frames.append({
            "edges": curr_step["edges"], 
            "added": None, 
            "removed": None, 
            "title": f"Step {curr_step['idx']}c: Tree Stabilized | {status}\nDiameter: {curr_step['diam']} (Limit: {out_data['d_bound']})",
            "color": base_color
        })
        
    anim_frames.extend([anim_frames[-1]] * 3) 
    
    def update(frame_idx):
        frame_data = anim_frames[frame_idx]
        draw_graph_state(ax, G_base, pos, 
                         frame_data["edges"], 
                         frame_data["added"], 
                         frame_data["removed"], 
                         frame_data["title"], 
                         frame_data["color"],
                         title_size=16)

    # Note: Interval is mostly ignored by ffmpeg writer, fps controls speed.
    anim = animation.FuncAnimation(fig, update, frames=len(anim_frames))
    
    # Save as MP4 using FFmpeg (0.8 fps = 1.25 seconds per frame)
    anim.save("sequence_animation.mp4", writer='ffmpeg', fps=1, dpi=150)
    print(f"Done! Generated sequence_animation.mp4 successfully.")
    
    plt.close()

if __name__ == "__main__":
    main()