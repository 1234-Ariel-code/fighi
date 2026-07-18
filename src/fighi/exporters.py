from __future__ import annotations

import csv
import json
from pathlib import Path
from xml.etree import ElementTree as ET

from .search import InteractionResult


def _graph_elements(
    interactions: list[InteractionResult], graph_top: int
) -> tuple[list[str], list[InteractionResult]]:
    chosen = sorted(
        (item for item in interactions if item.significant),
        key=lambda item: (item.q_value, item.p_value),
    )[:graph_top]
    features = sorted({feature for item in chosen for feature in item.features})
    return features, chosen


def export_hypergraph_json(
    path: Path, interactions: list[InteractionResult], graph_top: int
) -> None:
    features, chosen = _graph_elements(interactions, graph_top)
    payload = {
        "schema": "fighi-hypergraph-1.0",
        "nodes": [{"id": f"s{index}", "label": name, "type": "feature"} for index, name in enumerate(features)],
        "hyperedges": [
            {
                "id": f"h{index}",
                "members": list(item.features),
                "order": item.order,
                "fi_gain": item.fi_gain,
                "p_value": item.p_value,
                "q_value": item.q_value,
                "stability": item.stability,
            }
            for index, item in enumerate(chosen)
        ],
    }
    path.write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")


def export_cytoscape(
    path: Path, interactions: list[InteractionResult], graph_top: int
) -> None:
    features, chosen = _graph_elements(interactions, graph_top)
    nodes = [
        {"data": {"id": f"s{index}", "label": name, "type": "feature"}}
        for index, name in enumerate(features)
    ]
    feature_ids = {name: f"s{index}" for index, name in enumerate(features)}
    edges = []
    for index, item in enumerate(chosen):
        hyperedge_id = f"h{index}"
        nodes.append(
            {
                "data": {
                    "id": hyperedge_id,
                    "label": " × ".join(item.features),
                    "type": "hyperedge",
                    "order": item.order,
                    "fi_gain": item.fi_gain,
                    "p_value": item.p_value,
                    "q_value": item.q_value,
                }
            }
        )
        for member in item.features:
            edges.append(
                {
                    "data": {
                        "id": f"e{index}_{feature_ids[member]}",
                        "source": hyperedge_id,
                        "target": feature_ids[member],
                    }
                }
            )
    payload = {
        "format_version": "1.0",
        "generated_by": "FIGHI 1.0.0",
        "elements": {"nodes": nodes, "edges": edges},
    }
    path.write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")


def export_gml(path: Path, interactions: list[InteractionResult], graph_top: int) -> None:
    features, chosen = _graph_elements(interactions, graph_top)

    def clean(value: str) -> str:
        return value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")

    lines = ["graph [", "  directed 0"]
    feature_ids = {name: index for index, name in enumerate(features)}
    for name, node_id in feature_ids.items():
        lines.extend(
            ["  node [", f"    id {node_id}", f'    label "{clean(name)}"', '    type "feature"', "  ]"]
        )
    offset = len(features)
    for index, item in enumerate(chosen):
        node_id = offset + index
        lines.extend(
            [
                "  node [",
                f"    id {node_id}",
                f'    label "{clean(" × ".join(item.features))}"',
                '    type "hyperedge"',
                f"    fi_gain {item.fi_gain:.17g}",
                f"    q_value {item.q_value:.17g}",
                "  ]",
            ]
        )
        for member in item.features:
            lines.extend(
                ["  edge [", f"    source {node_id}", f"    target {feature_ids[member]}", "  ]"]
            )
    lines.append("]")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_graphml(path: Path, interactions: list[InteractionResult], graph_top: int) -> None:
    features, chosen = _graph_elements(interactions, graph_top)
    namespace = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", namespace)
    root = ET.Element(f"{{{namespace}}}graphml")
    for key_id, target, name, value_type in [
        ("label", "node", "label", "string"),
        ("type", "node", "type", "string"),
        ("fi_gain", "node", "fi_gain", "double"),
        ("q_value", "node", "q_value", "double"),
    ]:
        ET.SubElement(
            root,
            f"{{{namespace}}}key",
            {"id": key_id, "for": target, "attr.name": name, "attr.type": value_type},
        )
    graph = ET.SubElement(root, f"{{{namespace}}}graph", {"id": "FIGHI", "edgedefault": "undirected"})
    feature_ids = {name: f"s{index}" for index, name in enumerate(features)}
    for name, node_id in feature_ids.items():
        node = ET.SubElement(graph, f"{{{namespace}}}node", {"id": node_id})
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "label"}).text = name
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "type"}).text = "feature"
    for index, item in enumerate(chosen):
        node_id = f"h{index}"
        node = ET.SubElement(graph, f"{{{namespace}}}node", {"id": node_id})
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "label"}).text = " × ".join(item.features)
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "type"}).text = "hyperedge"
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "fi_gain"}).text = str(item.fi_gain)
        ET.SubElement(node, f"{{{namespace}}}data", {"key": "q_value"}).text = str(item.q_value)
        for member_index, member in enumerate(item.features):
            ET.SubElement(
                graph,
                f"{{{namespace}}}edge",
                {"id": f"e{index}_{member_index}", "source": node_id, "target": feature_ids[member]},
            )
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def export_edge_list(path: Path, interactions: list[InteractionResult], graph_top: int) -> None:
    _, chosen = _graph_elements(interactions, graph_top)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["hyperedge_id", "feature", "order", "fi_gain", "q_value"])
        for index, item in enumerate(chosen):
            for member in item.features:
                writer.writerow([f"h{index}", member, item.order, item.fi_gain, item.q_value])


def export_all(outdir: Path, interactions: list[InteractionResult], graph_top: int) -> dict[str, str]:
    paths = {
        "hypergraph_json": outdir / "fighi_hypergraph.json",
        "cytoscape": outdir / "fighi_cytoscape.cyjs",
        "gml": outdir / "fighi_hypergraph.gml",
        "graphml": outdir / "fighi_hypergraph.graphml",
        "edge_list": outdir / "fighi_hypergraph_edges.tsv",
    }
    export_hypergraph_json(paths["hypergraph_json"], interactions, graph_top)
    export_cytoscape(paths["cytoscape"], interactions, graph_top)
    export_gml(paths["gml"], interactions, graph_top)
    export_graphml(paths["graphml"], interactions, graph_top)
    export_edge_list(paths["edge_list"], interactions, graph_top)
    return {name: str(path) for name, path in paths.items()}

