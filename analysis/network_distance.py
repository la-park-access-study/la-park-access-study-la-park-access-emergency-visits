#!/usr/bin/env python3
"""network_distance.py: shortest-path step of analysis/02_network_distance.R (Python 3 standard library only).

Called by 02_network_distance.R with two arguments: a directory holding the files it wrote, and the name of a park set
(all, open, merged or pna2016):
  nodes.csv                 node_id, component                (nodes of the used road network)
  edges.csv                 edge_id, u, v, len_m              (undirected road edges with their length in meters)
  park_nodes_<set>.csv      park_id, node_id, gap_m           (network nodes inside or within 60 m of a park; gap = node-to-polygon distance)
  park_edges_<set>.csv      park_id, edge_id, gap_m, pos_m    (nearest point on the nearest road edge to each park, and every road crossing of the park boundary)
  origin_edges.csv          GEOID, edge_id, gap_m, pos_m      (nearest point on the nearest road edge to each tract's population-weighted center)
  origin_parks_<set>.csv    GEOID, park_id, gap_m             (parks that a tract center lies inside or within 60 m of; reached directly)
Writes network_result_<set>.csv: GEOID, net3_m, net1_m, origin_gap_m, component, parks.

Every park access point and every tract center becomes a node inserted into its road edge at position pos_m.
A multi-source run of Dijkstra's algorithm, started from all park access points at once with the park-to-road
gap as the starting distance, gives every node the distances to its three nearest distinct parks (each node
is settled at most three times, once per distinct park). A tract's distance to a park is then the center-to-road
gap plus the distance of its inserted node.
"""
import csv, heapq, os, sys
from collections import defaultdict

D = sys.argv[1]; SET = sys.argv[2]
K = 3                                  # number of nearest distinct parks kept per node

def read(name):
    with open(os.path.join(D, name), newline="") as f:
        return list(csv.DictReader(f))

nodes = read("nodes.csv"); edges = read("edges.csv")
N = max(int(r["node_id"]) for r in nodes) + 1
comp = [0] * N
for r in nodes:
    comp[int(r["node_id"])] = int(r["component"])
elen = {}; eu = {}; ev = {}
for r in edges:
    i = int(r["edge_id"]); elen[i] = float(r["len_m"]); eu[i] = int(r["u"]); ev[i] = int(r["v"])
print(f"network: {len(nodes)} nodes, {len(edges)} edges")

# --- nodes inserted into edges: park access points and tract centers ---------------------------
inserted = defaultdict(list)           # edge_id -> [(pos_m, new_node_id)]
next_id = N
access = []                            # (node_id, park_id, gap_m)
for r in read(f"park_nodes_{SET}.csv"):
    access.append((int(r["node_id"]), int(r["park_id"]), float(r["gap_m"])))
n_band = len(access)
for r in read(f"park_edges_{SET}.csv"):
    inserted[int(r["edge_id"])].append((float(r["pos_m"]), next_id))
    access.append((next_id, int(r["park_id"]), float(r["gap_m"]))); next_id += 1
origins = []                           # (GEOID, node_id, gap_m, component)
for r in read("origin_edges.csv"):
    e = int(r["edge_id"])
    inserted[e].append((float(r["pos_m"]), next_id))
    origins.append((r["GEOID"], next_id, float(r["gap_m"]), comp[eu[e]])); next_id += 1
direct = defaultdict(dict)             # GEOID -> {park_id: straight-line gap} for parks the center lies inside or within the band of
for r in read(f"origin_parks_{SET}.csv"):
    direct[r["GEOID"]][int(r["park_id"])] = float(r["gap_m"])
print(f"park access points: {n_band} network nodes near parks + {len(access) - n_band} points on the nearest edge; "
      f"tract centers: {len(origins)}; nodes inserted into {len(inserted)} edges")

# --- adjacency list, splitting edges at the inserted nodes --------------------------------------
adj = [[] for _ in range(next_id)]
def link(a, b, w):
    adj[a].append((b, w)); adj[b].append((a, w))
for i, L in elen.items():
    if i not in inserted:
        link(eu[i], ev[i], L)
    else:
        chain = [(0.0, eu[i])] + sorted(inserted[i]) + [(L, ev[i])]
        for (p0, a), (p1, b) in zip(chain, chain[1:]):
            link(a, b, max(p1 - p0, 0.0))
print(f"graph built: {next_id} nodes")

# --- multi-source Dijkstra with the K nearest distinct parks per node ---------------------------
lab_park = [[] for _ in range(next_id)]     # parks settled at each node, nearest first
lab_dist = [[] for _ in range(next_id)]
pq = [(gap, node, park) for node, park, gap in access]
heapq.heapify(pq)
settled = 0
while pq:
    d, u, p = heapq.heappop(pq)
    lp = lab_park[u]
    if len(lp) >= K or p in lp:
        continue
    lp.append(p); lab_dist[u].append(d); settled += 1
    for v, w in adj[u]:
        lv = lab_park[v]
        if len(lv) < K and p not in lv:
            heapq.heappush(pq, (d + w, v, p))
print(f"shortest paths: {settled} node-park labels settled")

# --- results per tract ----------------------------------------------------------------------------
rows = []
for geoid, node, gap, c in origins:
    ds = lab_dist[node]; ps = lab_park[node]
    if len(ds) < K:
        raise SystemExit(f"tract {geoid}: only {len(ds)} parks reachable from its network part")
    best = {p: gap + d for p, d in zip(ps, ds)}                     # via the road network
    for p, sl in direct[geoid].items():                              # directly, for parks the center lies in or beside
        if sl < best.get(p, float("inf")): best[p] = sl
    top = sorted(best.items(), key=lambda kv: kv[1])[:K]
    tot = [d for _, d in top]
    rows.append(dict(GEOID=geoid, net3_m=f"{sum(tot) / K:.3f}", net1_m=f"{tot[0]:.3f}", origin_gap_m=f"{gap:.3f}",
                     component=c, parks=" ".join(str(p) for p, _ in top)))
with open(os.path.join(D, f"network_result_{SET}.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print(f"wrote network_result_{SET}.csv for {len(rows)} tracts")
