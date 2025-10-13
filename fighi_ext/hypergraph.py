import csv, json

def export_csv(path, rows, header):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

def export_gml(path, nodes, hyperedges):
    with open(path, "w", encoding="utf-8") as f:
        f.write("graph [\n  directed 0\n")
        for i, name in enumerate(nodes):
            f.write(f"  node [ id {i} label \"{name}\" type \"node\" ]\n")
        hid_offset = len(nodes)
        for j,(hid, members, w) in enumerate(hyperedges):
            nid = hid_offset + j
            f.write(f"  node [ id {nid} label \"h{hid}\" type \"hyperedge\" weight {w} ]\n")
            for m in members:
                f.write(f"  edge [ source {nid} target {m} ]\n")
        f.write("]\n")

def export_hyper_json(path, nodes, hyperedges):
    data = {
        "nodes": [{"id": i, "name": n} for i,n in enumerate(nodes)],
        "hyperedges": [{"id": hid, "members": members, "weight": w} for hid,(members,w) in enumerate([(h[1],h[2]) for h in hyperedges])]
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

def export_cytoscape_cyjs(path, nodes, hyperedges):
    # represent hyperedge as an extra node connected to its members
    elements = {"nodes": [], "edges": []}
    for i, n in enumerate(nodes):
        elements["nodes"].append({"data": {"id": f"n{i}", "label": n, "type":"node"}})
    offset = len(nodes)
    for j,(hid, members, w) in enumerate(hyperedges):
        hid_str = f"h{hid}"
        jid = offset + j
        elements["nodes"].append({"data": {"id": f"n{jid}", "label": hid_str, "type":"hyperedge", "weight": w}})
        for m in members:
            elements["edges"].append({"data": {"id": f"e{jid}_{m}", "source": f"n{jid}", "target": f"n{m}"}})
    cyjs = {"format_version":"1.0", "generated_by":"FIGHI", "elements": elements}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cyjs, f, indent=2)

def export_graphml(path, nodes, hyperedges):
    """
    Simple GraphML (nodes + hyperedge-nodes connected to member SNPs).
    Nodes 'n{i}' for SNPs, 'h{j}' for hyperedges; edge attrs include weight.
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
        f.write('  <graph edgedefault="undirected">\n')
        # SNP nodes
        for i, n in enumerate(nodes):
            f.write(f'    <node id="n{i}"><data key="label">{n}</data></node>\n')
        # hyperedge nodes + edges to members
        offset = len(nodes)
        for j,(hid, members, w) in enumerate(hyperedges):
            f.write(f'    <node id="h{j}"><data key="type">hyperedge</data><data key="weight">{w}</data></node>\n')
            for m in members:
                f.write(f'    <edge source="h{j}" target="n{m}"/>\n')
        f.write('  </graph>\n</graphml>\n')

