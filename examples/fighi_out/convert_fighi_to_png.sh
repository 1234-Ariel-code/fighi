#!/usr/bin/env bash
# ==========================================================
# convert_fighi_to_png.sh
# Converts all FIGHI graph exports (.gml, .cyjs, .hyper)
# in your output directory into PNG visualizations
# ==========================================================

set -euo pipefail

# --- Target directory (your FIGHI output folder) ---
GRAPH_DIR="/home/arielghislain.kemogn/fighi_extended_package/fighi_ext/fighi_out"
OUT_DIR="${GRAPH_DIR}/figures"
mkdir -p "$OUT_DIR"

echo "[INFO] Converting FIGHI files in: $GRAPH_DIR"
echo "[INFO] Output PNGs will be saved to: $OUT_DIR"

# --- Python inline script ---
python3 - <<'PYCODE'
import json, os, sys
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path

graph_dir = Path("/example/fighi_out")
out_dir = graph_dir / "figures"
out_dir.mkdir(exist_ok=True)

def draw_and_save(G, title, out_path):
    print(f"[INFO] Rendering {title} -> {out_path}")
    pos = nx.spring_layout(G, seed=42)
    w = [G[e[0]][e[1]].get('weight', 1.0) for e in G.edges()]
    nx.draw(
        G, pos,
        with_labels=True, node_size=400, width=w,
        font_size=7, edge_color="gray", node_color="skyblue"
    )
    plt.title(title, fontsize=10)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[SAVED] {out_path}")

# ---- GML ----
for gml_file in graph_dir.glob("*.gml"):
    G = nx.read_gml(gml_file)
    out_path = out_dir / f"{gml_file.stem}.png"
    draw_and_save(G, f"GML: {gml_file.name}", out_path)

# ---- CYJS ----
for cyjs_file in graph_dir.glob("*.cyjs"):
    with open(cyjs_file) as f:
        data = json.load(f)
    G = nx.Graph()
    for n in data["elements"]["nodes"]:
        nid = n["data"]["id"]
        G.add_node(nid, **{k:v for k,v in n["data"].items() if k!="id"})
    for e in data["elements"]["edges"]:
        src, tgt = e["data"]["source"], e["data"]["target"]
        G.add_edge(src, tgt, **{k:v for k,v in e["data"].items() if k not in ("id","source","target")})
    out_path = out_dir / f"{cyjs_file.stem}.png"
    draw_and_save(G, f"CYJS: {cyjs_file.name}", out_path)

# ---- HYPER ----
for hyper_file in graph_dir.glob("*.hyper"):
    with open(hyper_file) as f:
        H = json.load(f)
    nodes = H.get("nodes") or H.get("node_names") or []
    hedges = H.get("hyperedges") or H.get("edges") or []

    B = nx.Graph()
    for n in nodes:
        B.add_node(f"n:{n}", kind="node")
    for i, he in enumerate(hedges):
        if isinstance(he, (list, tuple)) and len(he) == 3 and isinstance(he[1], list):
            _, members, w = he
        else:
            members, w = he, 1.0
        hname = f"h:{i}"
        B.add_node(hname, kind="hyperedge", weight=w)
        for m in members:
            name = nodes[m] if isinstance(m, int) else m
            B.add_edge(hname, f"n:{name}", weight=w)
    out_path = out_dir / f"{hyper_file.stem}.png"
    draw_and_save(B, f"HYPER: {hyper_file.name}", out_path)

print("[DONE] All FIGHI visualizations saved to:", out_dir)
PYCODE
