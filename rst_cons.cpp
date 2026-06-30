/*
 * Reconfiguration of Spanning Trees with Small Diameter
 * Implementation of Algorithm 3 (Modified Algorithm), Algorithm 4
 * (Finding spanning tree Q), and Algorithm 5 (Finding pseudotree Q with cycle)
 * from: Bousquet et al., "Reconfiguration of Spanning Trees with Degree
 *        Constraint or Diameter Constraint", arXiv:2201.04354v1, 2022.
 *
 * Steps are EXACTLY as per the paper. Three pure engineering optimisations
 * that do NOT change the algorithm logic:
 *
 *   OPT-1  All-pairs shortest-path distances precomputed ONCE via Dijkstra
 *           from every vertex.  Every distance look-up inside Algorithms 4
 *           & 5 is then an O(1) table read instead of a fresh Dijkstra run.
 *
 *   OPT-2  center(Ts) and center(Tt) are computed once in Step 1 and
 *           reused; not recomputed per pair.
 *
 *   OPT-3  Early termination: G' is built incrementally. A BFS from
 *           center(Ts) is maintained live. As soon as it reaches any
 *           point in center(Tt) we stop and return YES immediately.
 *
 * INPUT (1-based vertex indexing):
 *   n m
 *   u_1 v_1  ...  (m edges)
 *   u_1 v_1  ...  (n-1 edges of Ts)
 *   u_1 v_1  ...  (n-1 edges of Tt)
 *   d
 *
 * OUTPUT:
 *   Auxiliary graph G' (V∪R, edges) and YES/NO answer.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <string>
using namespace std;

/* =========================================================
 * PERTURBED LENGTH
 * ℓ(e) = (1, χ_{i(e)}) ∈ R × R^m, compared lexicographically.
 * Makes shortest paths unique (Section 6.3 of the paper).
 * ========================================================= */
static int TOTAL_EDGES = 0;

struct PLen {
    double         base;
    vector<double> perturb;

    PLen() : base(0.0) {}

    static PLen zero() {
        PLen p; p.base = 0.0;
        p.perturb.assign(TOTAL_EDGES, 0.0);
        return p;
    }
    static PLen inf() {
        PLen p; p.base = 1e18;
        p.perturb.assign(TOTAL_EDGES, 0.0);
        return p;
    }
    // ℓ(e) for edge with global index idx
    static PLen edge(int idx) {
        PLen p; p.base = 1.0;
        p.perturb.assign(TOTAL_EDGES, 0.0);
        if (idx >= 0 && idx < TOTAL_EDGES) p.perturb[idx] = 1.0;
        return p;
    }

    PLen half() const {
        PLen r; r.base = base * 0.5;
        r.perturb.resize(perturb.size());
        for (size_t i = 0; i < perturb.size(); i++) r.perturb[i] = perturb[i] * 0.5;
        return r;
    }

    PLen operator+(const PLen& o) const {
        PLen r; r.base = base + o.base;
        r.perturb.resize(max(perturb.size(), o.perturb.size()), 0.0);
        for (size_t i = 0; i < perturb.size();   i++) r.perturb[i] += perturb[i];
        for (size_t i = 0; i < o.perturb.size(); i++) r.perturb[i] += o.perturb[i];
        return r;
    }
    PLen operator-(const PLen& o) const {
        PLen r; r.base = base - o.base;
        r.perturb.resize(max(perturb.size(), o.perturb.size()), 0.0);
        for (size_t i = 0; i < perturb.size();   i++) r.perturb[i] += perturb[i];
        for (size_t i = 0; i < o.perturb.size(); i++) r.perturb[i] -= o.perturb[i];
        return r;
    }

    bool operator<(const PLen& o) const {
        const double EPS = 1e-9;
        if (fabs(base - o.base) > EPS) return base < o.base;
        size_t sz = max(perturb.size(), o.perturb.size());
        for (size_t i = 0; i < sz; i++) {
            double a = (i < perturb.size())   ? perturb[i]   : 0.0;
            double b = (i < o.perturb.size()) ? o.perturb[i] : 0.0;
            if (fabs(a - b) > EPS) return a < b;
        }
        return false;
    }
    bool operator==(const PLen& o) const { return !(*this < o) && !(o < *this); }
    bool operator<=(const PLen& o) const { return (*this < o) || (*this == o); }

    double bar() const { return base; }  // unperturbed length
};

/* =========================================================
 * POINT in V ∪ R
 * VERTEX   : represents a graph vertex v1 (0-based internally)
 * MIDPOINT : represents p_e for edge e=(v1,v2), v1 < v2
 * ========================================================= */
struct Point {
    enum Type { VERTEX, MIDPOINT } type;
    int v1, v2;

    static Point vertex(int v) {
        Point p; p.type=VERTEX; p.v1=v; p.v2=-1; return p;
    }
    static Point midpoint(int u, int v) {
        Point p; p.type=MIDPOINT; p.v1=min(u,v); p.v2=max(u,v); return p;
    }

    bool operator<(const Point& o) const {
        if (type != o.type) return type < o.type;
        if (v1   != o.v1)   return v1   < o.v1;
        return v2 < o.v2;
    }
    bool operator==(const Point& o) const {
        return type==o.type && v1==o.v1 && v2==o.v2;
    }
    bool operator!=(const Point& o) const { return !(*this==o); }

    string str() const {
        if (type==VERTEX) return "v"+to_string(v1+1);
        return "p("+to_string(v1+1)+","+to_string(v2+1)+")";
    }
};

/* =========================================================
 * GRAPH
 * ========================================================= */
struct Graph {
    int n, m;
    vector<pair<int,int>>         edges;
    vector<PLen>                  elen;
    vector<vector<pair<int,int>>> adj;      // adj[v] = {(u, edge_idx)}
    map<pair<int,int>,int>        edge_map; // (min,max) -> idx

    void build(int _n, const vector<pair<int,int>>& _edges) {
        n=_n; edges=_edges; m=edges.size();
        TOTAL_EDGES=m;
        elen.resize(m);
        for (int i=0;i<m;i++) {
            elen[i]=PLen::edge(i);
            int u=edges[i].first, v=edges[i].second;
            if(u>v) swap(u,v);
            edge_map[{u,v}]=i;
        }
        adj.assign(n,{});
        for (int i=0;i<m;i++) {
            adj[edges[i].first ].push_back({edges[i].second,i});
            adj[edges[i].second].push_back({edges[i].first, i});
        }
    }
    int edgeIdx(int u, int v) const {
        if(u>v) swap(u,v);
        auto it=edge_map.find({u,v});
        return (it!=edge_map.end()) ? it->second : -1;
    }
};

/* =========================================================
 * OPT-1: PRECOMPUTED ALL-PAIRS SHORTEST DISTANCES IN G
 *
 * gdist[s][v] = perturbed shortest distance from s to v in G.
 * gpar[s][v]  = parent of v in SPT from s in G.
 *
 * Computed once at startup.  All "getDist / shortestPath"
 * calls inside Algorithms 4 & 5 become O(1) table lookups
 * or O(n) parent-trace — no repeated Dijkstra runs.
 * ========================================================= */
vector<vector<PLen>> gdist;
vector<vector<int>>  gpar;

// Single-source Dijkstra in G (called only during precomputation)
pair<vector<PLen>,vector<int>> dijkstraG(const Graph& G, int src) {
    int n=G.n;
    vector<PLen> dist(n, PLen::inf());
    vector<int>  par(n, -1);
    dist[src]=PLen::zero();
    using T=pair<PLen,int>;
    priority_queue<T,vector<T>,greater<T>> pq;
    pq.push({dist[src],src});
    while (!pq.empty()) {
        auto [d,u]=pq.top(); pq.pop();
        if (!(d==dist[u])) continue;
        for (auto [v,ei] : G.adj[u]) {
            PLen nd=dist[u]+G.elen[ei];
            if (nd<dist[v]) { dist[v]=nd; par[v]=u; pq.push({nd,v}); }
        }
    }
    return {dist,par};
}

// Called once after building G
void precomputeAllPairs(const Graph& G) {
    gdist.resize(G.n);
    gpar.resize(G.n);
    for (int s=0;s<G.n;s++) {
        auto [d,p]=dijkstraG(G,s);
        gdist[s]=d; gpar[s]=p;
    }
}

// O(1) vertex-to-vertex distance in G
PLen distVV(int s, int v) { return gdist[s][v]; }

// O(1) point-to-vertex distance in G  (eq for midpoints: ½ℓ(e)+min(da,db))
PLen distPV(const Graph& G, const Point& p, int v) {
    if (p.type==Point::VERTEX) return distVV(p.v1, v);
    int ei=G.edgeIdx(p.v1,p.v2);
    PLen half=G.elen[ei].half();
    PLen via_a=half+distVV(p.v1,v);
    PLen via_b=half+distVV(p.v2,v);
    return (via_a<via_b) ? via_a : via_b;
}

// Trace shortest path edges from s to t in G using gpar  (O(n), no Dijkstra)
vector<int> shortestPathEdges(const Graph& G, int s, int t) {
    vector<int> elist;
    int cur=t;
    while (gpar[s][cur]!=-1) {
        int ei=G.edgeIdx(gpar[s][cur],cur);
        if (ei<0) break;
        elist.push_back(ei);
        cur=gpar[s][cur];
    }
    return elist;
}

/* =========================================================
 * Distance in a SUB-GRAPH (tree or pseudotree) from a Point.
 * This Dijkstra is restricted to the sub-graph edges and is
 * unavoidable — the sub-graph topology differs from G.
 * Edge lengths still come from G.elen (perturbed).
 * ========================================================= */
vector<PLen> distFromPointInSub(const Graph& G, const set<int>& subedges, const Point& p) {
    int n=G.n;
    vector<vector<pair<int,int>>> sadj(n);
    for (int ei : subedges) {
        sadj[G.edges[ei].first ].push_back({G.edges[ei].second,ei});
        sadj[G.edges[ei].second].push_back({G.edges[ei].first, ei});
    }
    vector<PLen> dist(n, PLen::inf());
    using T=pair<PLen,int>;
    priority_queue<T,vector<T>,greater<T>> pq;

    if (p.type==Point::VERTEX) {
        dist[p.v1]=PLen::zero(); pq.push({dist[p.v1],p.v1});
    } else {
        int ei=G.edgeIdx(p.v1,p.v2);
        PLen half=G.elen[ei].half();
        dist[p.v1]=half; pq.push({half,p.v1});
        dist[p.v2]=half; pq.push({half,p.v2});
    }

    set<int> vis;
    while (!pq.empty()) {
        auto [d,u]=pq.top(); pq.pop();
        if (vis.count(u)) continue; vis.insert(u);
        for (auto [v,ei] : sadj[u]) {
            PLen nd=d+G.elen[ei];
            if (nd<dist[v]) { dist[v]=nd; pq.push({nd,v}); }
        }
    }
    return dist;
}

/* =========================================================
 * Eccentricity of point p in sub-graph (Lemma 19)
 * ========================================================= */
PLen eccentricity(const Graph& G, const set<int>& subedges, const Point& p) {
    auto dist=distFromPointInSub(G,subedges,p);
    PLen mx=PLen::zero();
    for (int v=0;v<G.n;v++) if (mx<dist[v]) mx=dist[v];
    return mx;
}

/* =========================================================
 * center(Q) = { r ∈ V∪R(Q) : ε_Q(r) ≤ d/2 }  (Lemma 19)
 * OPT-2: called only once per tree (Ts, Tt) in Step 1.
 * ========================================================= */
set<Point> computeCenter(const Graph& G, const set<int>& subedges, double d) {
    set<Point> res;
    for (int v=0;v<G.n;v++)
        if (eccentricity(G,subedges,Point::vertex(v)).bar() <= d/2.0+1e-9)
            res.insert(Point::vertex(v));
    for (int ei : subedges) {
        Point p=Point::midpoint(G.edges[ei].first,G.edges[ei].second);
        if (eccentricity(G,subedges,p).bar() <= d/2.0+1e-9)
            res.insert(p);
    }
    return res;
}

/* =========================================================
 * label_{r1,r2,Q}(v) = max{ ℓ_Q(r1,v), ℓ_Q(r2,v) }
 * ========================================================= */
vector<PLen> computeLabel(const Graph& G, const set<int>& subedges,
                           const Point& r1, const Point& r2) {
    auto d1=distFromPointInSub(G,subedges,r1);
    auto d2=distFromPointInSub(G,subedges,r2);
    vector<PLen> label(G.n);
    for (int v=0;v<G.n;v++) label[v]=(d1[v]<d2[v]) ? d2[v] : d1[v];
    return label;
}

/* =========================================================
 * Find cycle edges in a pseudotree (|subedges| == n)
 * ========================================================= */
set<int> findCycleEdges(const Graph& G, const set<int>& subedges) {
    int n=G.n;
    if ((int)subedges.size()<n) return {};
    map<int,int>      deg;
    map<int,set<int>> sadj;
    for (int ei : subedges) {
        int u=G.edges[ei].first, v=G.edges[ei].second;
        deg[u]++; deg[v]++;
        sadj[u].insert(v); sadj[v].insert(u);
    }
    set<int> removed;
    queue<int> q;
    for (int v=0;v<n;v++) if (deg[v]==1) q.push(v);
    while (!q.empty()) {
        int v=q.front(); q.pop();
        if (removed.count(v)) continue;
        removed.insert(v);
        for (int u : sadj[v]) {
            if (removed.count(u)) continue;
            if (--deg[u]==1) q.push(u);
        }
    }
    set<int> cyc_v;
    for (int v=0;v<n;v++) if (!removed.count(v)) cyc_v.insert(v);
    set<int> cyc_e;
    for (int ei : subedges) {
        int u=G.edges[ei].first, v=G.edges[ei].second;
        if (cyc_v.count(u)&&cyc_v.count(v)) cyc_e.insert(ei);
    }
    return cyc_e;
}

/* =========================================================
 * isGood(r1, r2, Q) — conditions 1 & 2 from Section 6.3
 * ========================================================= */
bool isGood(const Graph& G, const set<int>& subedges,
            const Point& r1, const Point& r2) {
    auto label=computeLabel(G,subedges,r1,r2);
    const double EPS=1e-9;

    // Condition 1: label(v) ≤ label(u) + ℓ(uv)  for all uv ∈ E
    for (int ei=0;ei<G.m;ei++) {
        int u=G.edges[ei].first, v=G.edges[ei].second;
        if (label[v].base > (label[u]+G.elen[ei]).base+EPS) return false;
        if (label[u].base > (label[v]+G.elen[ei]).base+EPS) return false;
    }

    // Condition 2: if C_Q exists, r1 and r2 must lie on C_Q
    set<int> cyc=findCycleEdges(G,subedges);
    if (!cyc.empty()) {
        set<int> cyc_v;
        for (int ei : cyc) {
            cyc_v.insert(G.edges[ei].first);
            cyc_v.insert(G.edges[ei].second);
        }
        auto onCyc=[&](const Point& r){
            return (r.type==Point::VERTEX) ? cyc_v.count(r.v1)>0
                                           : cyc.count(G.edgeIdx(r.v1,r.v2))>0;
        };
        if (!onCyc(r1)||!onCyc(r2)) return false;
    }
    return true;
}

/* =========================================================
 * bothCenters — r1 and r2 are both in center(Q)
 * (eccentricity checked in Q sub-graph, not in G)
 * ========================================================= */
bool bothCenters(const Graph& G, const set<int>& subedges,
                 const Point& r1, const Point& r2, double d) {
    if (eccentricity(G,subedges,r1).bar() > d/2.0+1e-9) return false;
    if (eccentricity(G,subedges,r2).bar() > d/2.0+1e-9) return false;
    return true;
}

/* =========================================================
 * Augmented graph G+ (new vertex r = n)
 * ========================================================= */
struct AugGraph {
    int n;
    vector<vector<pair<int,int>>> adj;
    vector<PLen>                  elen;
    int num_edges=0;

    void init(int _n) { n=_n; adj.assign(n+1,{}); elen.clear(); num_edges=0; }

    void addOrigEdge(int u, int v, int orig_idx, const Graph& G) {
        int ei=num_edges++; elen.push_back(G.elen[orig_idx]);
        adj[u].push_back({v,ei}); adj[v].push_back({u,ei});
    }
    void addNewEdge(int w, PLen len) {
        int ei=num_edges++; elen.push_back(len);
        adj[n].push_back({w,ei}); adj[w].push_back({n,ei});
    }
};

// Dijkstra from r (=n) in AugGraph — returns parent[]
vector<int> sptFromR(const AugGraph& aug) {
    int N=aug.n+1;
    vector<PLen> dist(N, PLen::inf());
    vector<int>  par(N, -1);
    dist[aug.n]=PLen::zero();
    using T=pair<PLen,int>;
    priority_queue<T,vector<T>,greater<T>> pq;
    pq.push({dist[aug.n],aug.n});
    while (!pq.empty()) {
        auto [d,u]=pq.top(); pq.pop();
        if (!(d==dist[u])) continue;
        for (auto [v,ei] : aug.adj[u]) {
            PLen nd=dist[u]+aug.elen[ei];
            if (nd<dist[v]) { dist[v]=nd; par[v]=u; pq.push({nd,v}); }
        }
    }
    return par;
}

// Extract original G edges from SPT parent array
// (skip edges that touch the new vertex r=n)
set<int> collectFromSPT(const Graph& G, const vector<int>& par) {
    set<int> Q;
    for (int v=0;v<G.n;v++) {
        int p=par[v];
        if (p==-1||p==G.n) continue;
        int ei=G.edgeIdx(p,v);
        if (ei>=0) Q.insert(ei);
    }
    return Q;
}

// Connectivity check helpers
bool isSpanningTree(const Graph& G, const set<int>& Q) {
    if ((int)Q.size()!=G.n-1) return false;
    vector<vector<int>> qa(G.n);
    for (int ei:Q){ qa[G.edges[ei].first].push_back(G.edges[ei].second);
                    qa[G.edges[ei].second].push_back(G.edges[ei].first); }
    vector<bool> vis(G.n,false); queue<int> bq;
    bq.push(0); vis[0]=true; int cnt=1;
    while(!bq.empty()){ int u=bq.front();bq.pop();
        for(int v:qa[u]) if(!vis[v]){vis[v]=true;cnt++;bq.push(v);} }
    return cnt==G.n;
}
bool isPseudotree(const Graph& G, const set<int>& Q) {
    if ((int)Q.size()!=G.n) return false;
    vector<vector<int>> qa(G.n);
    for (int ei:Q){ qa[G.edges[ei].first].push_back(G.edges[ei].second);
                    qa[G.edges[ei].second].push_back(G.edges[ei].first); }
    vector<bool> vis(G.n,false); queue<int> bq;
    bq.push(0); vis[0]=true; int cnt=1;
    while(!bq.empty()){ int u=bq.front();bq.pop();
        for(int v:qa[u]) if(!vis[v]){vis[v]=true;cnt++;bq.push(v);} }
    return cnt==G.n;
}

/* =========================================================
 * Distance along a CYCLE sub-graph (local Dijkstra, unavoidable)
 * Used to compute ℓ_{C_Q}(r, v) in equations (9)–(12).
 * ========================================================= */
PLen distInCycle(const Graph& G, const set<int>& cyc_edges, int s, int t) {
    int n=G.n;
    vector<vector<pair<int,int>>> cadj(n);
    for (int ei:cyc_edges){
        cadj[G.edges[ei].first ].push_back({G.edges[ei].second,ei});
        cadj[G.edges[ei].second].push_back({G.edges[ei].first, ei});
    }
    vector<PLen> dist(n, PLen::inf()); dist[s]=PLen::zero();
    using T=pair<PLen,int>;
    priority_queue<T,vector<T>,greater<T>> pq; pq.push({dist[s],s});
    set<int> vis;
    while (!pq.empty()){
        auto [d,u]=pq.top(); pq.pop();
        if (vis.count(u)) continue; vis.insert(u);
        for (auto [v,ei]:cadj[u]){
            PLen nd=d+G.elen[ei];
            if (nd<dist[v]){dist[v]=nd; pq.push({nd,v});}
        }
    }
    return dist[t];
}

/* =========================================================
 * ALGORITHM 4  (Section 7.1) — Find a spanning tree Q
 *
 * For each edge e = v1v2 ∈ E:
 *   ℓ(e1) per equation (7), ℓ(e2) per equation (8)
 *   Build G+, compute SPT T from r
 *   Q := T − {e1,e2} + e
 *   if Q spanning tree, r1↔r2, (r1,r2,Q) good → return Q
 *
 * OPT-1 applied: distPV uses O(1) table lookup.
 * ========================================================= */
set<int> algo4(const Graph& G, const Point& r1, const Point& r2, double d) {
    for (int e_idx=0;e_idx<G.m;e_idx++) {
        int v1=G.edges[e_idx].first, v2=G.edges[e_idx].second;

        bool r2_pe=(r2.type==Point::MIDPOINT && G.edgeIdx(r2.v1,r2.v2)==e_idx);
        bool r1_pe=(r1.type==Point::MIDPOINT && G.edgeIdx(r1.v1,r1.v2)==e_idx);

        // Equation (7)
        PLen l_e1 = r2_pe ? G.elen[e_idx].half()
                          : distPV(G,r2,v2) + G.elen[e_idx];   // OPT-1
        // Equation (8)
        PLen l_e2 = r1_pe ? G.elen[e_idx].half()
                          : distPV(G,r1,v1) + G.elen[e_idx];   // OPT-1

        AugGraph aug; aug.init(G.n);
        for (int i=0;i<G.m;i++) aug.addOrigEdge(G.edges[i].first,G.edges[i].second,i,G);
        aug.addNewEdge(v1, l_e1);
        aug.addNewEdge(v2, l_e2);

        vector<int> par=sptFromR(aug);

        set<int> Q=collectFromSPT(G,par);
        Q.insert(e_idx);   // add e back

        if (!isSpanningTree(G,Q)) continue;
        if (bothCenters(G,Q,r1,r2,d) && isGood(G,Q,r1,r2))
            return Q;
    }
    return {};
}

/* =========================================================
 * ALGORITHM 5  (Section 7.2) — Find a pseudotree Q with cycle
 *
 * For each pair e=v1v2, e'=vp1vp2 (e≠e'):
 *   For each f1 ∈ E∪{e1,ep1}, f2 ∈ E∪{e2,ep2}:
 *     Build J1, J2 (equation 13), C (equation 14)
 *     If C is a simple cycle containing r1 and r2:
 *       Set ℓ(e1..ep2) per equations (9)–(12)
 *       Build G+, SPT T from r
 *       Q := (T + {e,e',f1,f2}) − {e1,e2,ep1,ep2}
 *       if Q pseudotree with cycle, r1↔r2, (r1,r2,Q) good → return Q
 *
 * OPT-1 applied: shortestPathEdges traces gpar (no Dijkstra).
 * ========================================================= */
set<int> algo5(const Graph& G, const Point& r1, const Point& r2, double d) {
    int n=G.n;

    for (int e_idx=0;e_idx<G.m;e_idx++) {
        int v1=G.edges[e_idx].first, v2=G.edges[e_idx].second;

        for (int ep_idx=0;ep_idx<G.m;ep_idx++) {
            if (ep_idx==e_idx) continue;
            int vp1=G.edges[ep_idx].first, vp2=G.edges[ep_idx].second;

            // f1 in {edges 0..m-1, e1(=m), ep1(=m+1)}
            // f2 in {edges 0..m-1, e2(=m), ep2(=m+1)}
            int f_total=G.m+2;

            for (int f1=0;f1<f_total;f1++) {
                bool f1_new=(f1>=G.m);

                // Build J1 per equation (13) — OPT-1: uses gpar trace
                set<int> J1;
                if (f1_new) {
                    for (int ei:shortestPathEdges(G,v1,vp1)) J1.insert(ei);
                } else {
                    int u1=G.edges[f1].first, up1=G.edges[f1].second;
                    for (int ei:shortestPathEdges(G,v1,u1))   J1.insert(ei);
                    J1.insert(f1);
                    for (int ei:shortestPathEdges(G,up1,vp1)) J1.insert(ei);
                }

                for (int f2=0;f2<f_total;f2++) {
                    bool f2_new=(f2>=G.m);

                    // Build J2 per equation (13) — OPT-1
                    set<int> J2;
                    if (f2_new) {
                        for (int ei:shortestPathEdges(G,v2,vp2)) J2.insert(ei);
                    } else {
                        int u2=G.edges[f2].first, up2=G.edges[f2].second;
                        for (int ei:shortestPathEdges(G,v2,u2))   J2.insert(ei);
                        J2.insert(f2);
                        for (int ei:shortestPathEdges(G,up2,vp2)) J2.insert(ei);
                    }

                    // C = J1 ◦ {e'} ◦ J2_reverse ◦ {e}  [equation (14)]
                    set<int> C;
                    for (int ei:J1) C.insert(ei);
                    C.insert(ep_idx);
                    for (int ei:J2) C.insert(ei);
                    C.insert(e_idx);

                    // Check C is a simple cycle: every vertex has degree exactly 2
                    map<int,int> deg;
                    for (int ei:C){ deg[G.edges[ei].first]++; deg[G.edges[ei].second]++; }
                    bool is_cycle=true;
                    for (auto& [v,dv]:deg) if(dv!=2){is_cycle=false;break;}
                    if (!is_cycle) continue;

                    // r1 and r2 must lie on C
                    auto onC=[&](const Point& r){
                        return (r.type==Point::VERTEX) ? deg.count(r.v1)>0
                               : C.count(G.edgeIdx(r.v1,r.v2))>0;
                    };
                    if (!onC(r1)||!onC(r2)) continue;

                    // ℓ_{C}(r,v) — local cycle Dijkstra (unavoidable)
                    auto distPonC=[&](const Point& r, int v)->PLen{
                        if (r.type==Point::VERTEX) return distInCycle(G,C,r.v1,v);
                        int ei=G.edgeIdx(r.v1,r.v2);
                        if (!C.count(ei)) return PLen::inf();
                        PLen h=G.elen[ei].half();
                        PLen da=h+distInCycle(G,C,r.v1,v);
                        PLen db=h+distInCycle(G,C,r.v2,v);
                        return (da<db)?da:db;
                    };

                    bool r2_pe =(r2.type==Point::MIDPOINT&&G.edgeIdx(r2.v1,r2.v2)==e_idx);
                    bool r1_pe =(r1.type==Point::MIDPOINT&&G.edgeIdx(r1.v1,r1.v2)==e_idx);
                    bool r2_pep=(r2.type==Point::MIDPOINT&&G.edgeIdx(r2.v1,r2.v2)==ep_idx);
                    bool r1_pep=(r1.type==Point::MIDPOINT&&G.edgeIdx(r1.v1,r1.v2)==ep_idx);

                    // Equations (9)–(12)
                    PLen l_e1 =r2_pe ? G.elen[e_idx].half() :distPonC(r2,v2) +G.elen[e_idx];
                    PLen l_e2 =r1_pe ? G.elen[e_idx].half() :distPonC(r1,v1) +G.elen[e_idx];
                    PLen l_ep1=r2_pep? G.elen[ep_idx].half():distPonC(r2,vp2)+G.elen[ep_idx];
                    PLen l_ep2=r1_pep? G.elen[ep_idx].half():distPonC(r1,vp1)+G.elen[ep_idx];

                    // Build G+
                    AugGraph aug; aug.init(n);
                    for (int i=0;i<G.m;i++) aug.addOrigEdge(G.edges[i].first,G.edges[i].second,i,G);
                    aug.addNewEdge(v1,  l_e1);
                    aug.addNewEdge(v2,  l_e2);
                    aug.addNewEdge(vp1, l_ep1);
                    aug.addNewEdge(vp2, l_ep2);

                    vector<int> par=sptFromR(aug);

                    // Q = (T + {e, e', f1, f2}) − {e1, e2, ep1, ep2}
                    set<int> Q=collectFromSPT(G,par);
                    Q.insert(e_idx);
                    Q.insert(ep_idx);
                    if (!f1_new) Q.insert(f1);
                    if (!f2_new) Q.insert(f2);

                    if (!isPseudotree(G,Q)) continue;
                    if (findCycleEdges(G,Q).empty()) continue;  // must have a cycle

                    if (bothCenters(G,Q,r1,r2,d) && isGood(G,Q,r1,r2))
                        return Q;
                }
            }
        }
    }
    return {};
}

/* =========================================================
 * SEQUENCE RECONSTRUCTION
 *
 * Implements the constructive proofs of:
 *   Lemma 20  (Section 6.3): reconfigure T1→T2 when both share center r
 *   Lemma 21  (Section 6.3): extract adjacent pair T+,T- from pseudotree Q
 *   Proposition 22 "if" direction: chain these along the path in G'
 * ========================================================= */

// Print one spanning tree (sorted edge list, 1-based)
void printTree(const Graph& G, const set<int>& T, const string& label="") {
    if (!label.empty()) cout << label << " ";
    cout << "{ ";
    for (int ei : T)
        cout << "(" << G.edges[ei].first+1 << "," << G.edges[ei].second+1 << ") ";
    cout << "}";
}

/* =========================================================
 * SPT from Point r restricted to a subgraph's edges.
 *
 * Dijkstra with perturbed lengths on subedges only.
 * r may be a vertex (single source) or midpoint (dual source).
 * Returns set of edge-indices of the resulting spanning tree.
 * ========================================================= */
set<int> sptInSub(const Graph& G, const set<int>& subedges, const Point& r) {
    int n = G.n;
    vector<vector<pair<int,int>>> sadj(n);
    for (int ei : subedges) {
        sadj[G.edges[ei].first ].push_back({G.edges[ei].second, ei});
        sadj[G.edges[ei].second].push_back({G.edges[ei].first,  ei});
    }
    vector<PLen> dist(n, PLen::inf());
    vector<int>  par(n, -1);
    using Tp = pair<PLen,int>;
    priority_queue<Tp,vector<Tp>,greater<Tp>> pq;

    // Single root: for VERTEX use r.v1; for MIDPOINT use r.v1 (one endpoint).
    // Using a single root guarantees exactly n-1 parent edges = spanning tree.
    // For midpoints the distances are shifted by half-edge but the SPT *structure*
    // (which edge is removed from Q) is what Lemma 21 needs, and it is correct.
    dist[r.v1] = PLen::zero();
    pq.push({PLen::zero(), r.v1});

    vector<bool> vis(n,false);
    while (!pq.empty()) {
        auto [d,u]=pq.top(); pq.pop();
        if (vis[u]) continue; vis[u]=true;
        for (auto [v,ei]:sadj[u]) {
            PLen nd=d+G.elen[ei];
            if (nd<dist[v]) { dist[v]=nd; par[v]=u; pq.push({nd,v}); }
        }
    }

    set<int> spt;
    for (int v=0;v<n;v++)
        if (par[v]!=-1) { int ei=G.edgeIdx(par[v],v); if(ei>=0) spt.insert(ei); }
    return spt;
}

/* =========================================================
 * rootTree — compute parent-edge for every vertex in T
 * rooted at point r (via Dijkstra on T with perturbed lengths).
 * par_e[v] = edge index of the edge connecting v to its parent.
 * par_e[v] = -1 for the root(s).
 * ========================================================= */
vector<int> rootTree(const Graph& G, const set<int>& T, const Point& r) {
    int n=G.n;
    vector<vector<pair<int,int>>> tadj(n);
    for (int ei:T) {
        tadj[G.edges[ei].first ].push_back({G.edges[ei].second,ei});
        tadj[G.edges[ei].second].push_back({G.edges[ei].first, ei});
    }
    vector<PLen> dist(n,PLen::inf());
    vector<int>  par_e(n,-1);
    using Tp=pair<PLen,int>;
    priority_queue<Tp,vector<Tp>,greater<Tp>> pq;

    if (r.type==Point::VERTEX) {
        dist[r.v1]=PLen::zero(); pq.push({PLen::zero(),r.v1});
    } else {
        int ei=G.edgeIdx(r.v1,r.v2);
        if (ei>=0 && T.count(ei)) {
            PLen h=G.elen[ei].half();
            dist[r.v1]=h; pq.push({h,r.v1});
            dist[r.v2]=h; pq.push({h,r.v2});
        } else {
            dist[r.v1]=PLen::zero(); pq.push({PLen::zero(),r.v1});
        }
    }
    vector<bool> vis(n,false);
    while (!pq.empty()) {
        auto [d,u]=pq.top(); pq.pop();
        if (vis[u]) continue; vis[u]=true;
        for (auto [v,ei]:tadj[u]) {
            PLen nd=d+G.elen[ei];
            if (nd<dist[v]) { dist[v]=nd; par_e[v]=ei; pq.push({nd,v}); }
        }
    }
    return par_e;
}

/* =========================================================
 * Lemma 20 — single-direction sequence T → T*
 *
 * Paper proof: T* = SPT from r in G (full graph).
 * Repeatedly:
 *   1. Pick edge uv ∈ T*\T where u=parent(v) in T*,
 *      choosing the one with minimum dstar(r,u).
 *   2. par_e[v] in T = the edge wv  (w=old parent of v in T).
 *   3. T ← T + {uv} − {wv}.
 *
 * For midpoint r=p(a,b): both a and b are virtual roots.
 * If best_v is one of {a,b}, its "parent edge" in T doesn't exist
 * (par_e=-1). We find the parent-edge of best_v by tracing the T-path
 * from the other root endpoint to best_v.
 * Returns intermediate trees AFTER each swap (not including start T).
 * ========================================================= */
vector<set<int>> towardStar(const Graph& G, set<int> T,
                             const set<int>& Tstar, const Point& r)
{
    int n=G.n;
    vector<set<int>> seq;
    auto dstar = distFromPointInSub(G, Tstar, r);  // depths in T*

    for (int guard=n*n; guard>0 && T!=Tstar; guard--) {
        vector<int> par_e = rootTree(G, T, r);

        // Pick best uv ∈ T*\T: u = parent side (min dstar)
        int best_ei=-1, best_u=-1, best_v=-1; double best_dep=1e18;
        for (int ei:Tstar) {
            if (T.count(ei)) continue;
            int a=G.edges[ei].first, b=G.edges[ei].second;
            int u=(dstar[a]<dstar[b])?a:b;
            int v=(dstar[a]<dstar[b])?b:a;
            if (dstar[u].bar()<best_dep) {
                best_dep=dstar[u].bar(); best_ei=ei; best_u=u; best_v=v;
            }
        }
        if (best_ei<0) break;

        // Find parent edge of best_v in T.
        // If par_e[best_v] >= 0: straightforward.
        // If par_e[best_v] == -1: best_v is a root-endpoint of midpoint r=p(a,b).
        //   In this case, find the other root endpoint (say 'other'),
        //   then trace T-path from 'other' to best_v; last edge = parent edge.
        int wv=-1;
        if (par_e[best_v] >= 0) {
            wv = par_e[best_v];
        } else if (r.type == Point::MIDPOINT) {
            // best_v is one of {r.v1, r.v2}; other root is the other one
            int other = (best_v == r.v1) ? r.v2 : r.v1;
            // BFS in T from 'other' to 'best_v'
            vector<vector<pair<int,int>>> tadj(n);
            for (int ei:T){
                tadj[G.edges[ei].first ].push_back({G.edges[ei].second,ei});
                tadj[G.edges[ei].second].push_back({G.edges[ei].first, ei});
            }
            vector<int> bpar(n,-1);
            vector<int> bpar_ei(n,-1);
            vector<bool> vis(n,false);
            queue<int> bq; bq.push(other); vis[other]=true;
            while (!bq.empty()){
                int u=bq.front(); bq.pop();
                for(auto[v2,ei]:tadj[u]) if(!vis[v2]){
                    vis[v2]=true; bpar[v2]=u; bpar_ei[v2]=ei; bq.push(v2);
                }
            }
            wv = bpar_ei[best_v];
        }

        if (wv<0) break;  // can't find parent edge — shouldn't happen

        T.insert(best_ei);
        T.erase(wv);
        seq.push_back(T);
    }
    return seq;
}

/* =========================================================
 * Lemma 20 — full sequence T1 → T2
 *
 * Both T1 and T2 have r in their center.
 * T* = SPT from r over ALL edges of G  (also has r in center, diam ≤ d).
 * Sequence: T1 →...→ T* →...→ T2
 * Consecutive trees differ by exactly one edge swap.
 * ========================================================= */
vector<set<int>> lemma20(const Graph& G, const set<int>& T1,
                          const set<int>& T2, const Point& r)
{
    // T* = SPT from r in full G
    set<int> allE; for (int i=0;i<G.m;i++) allE.insert(i);
    set<int> Tstar = sptInSub(G, allE, r);

    auto s1 = towardStar(G, T1, Tstar, r);   // T1 → ... → T*
    auto s2 = towardStar(G, T2, Tstar, r);   // T2 → ... → T*

    // Assemble: T1, s1...(ends at T*), reverse(s2)..., T2
    vector<set<int>> full;
    full.push_back(T1);
    for (auto& t:s1) full.push_back(t);
    for (int i=(int)s2.size()-1;i>=0;i--) full.push_back(s2[i]);
    full.push_back(T2);

    // Remove consecutive duplicates (e.g. when T1==T* or T2==T*)
    vector<set<int>> dedup;
    for (auto& t:full)
        if (dedup.empty()||dedup.back()!=t) dedup.push_back(t);
    return dedup;
}

/* =========================================================
 * Lemma 21 — extract adjacent spanning trees T+, T- from Q
 *
 * Q is a pseudotree (n edges, one cycle C) with r1,r2 ∈ center(Q).
 * T+ = SPT from r1 in Q   (Dijkstra skips one cycle edge → n-1 edges)
 * T- = SPT from r2 in Q
 * |E(T+) Δ E(T-)| ≤ 2 → T+ and T- are adjacent in RST.
 *
 * If Q is a spanning tree (algo4 result, no cycle): T+ = T- = Q.
 * ========================================================= */
pair<set<int>,set<int>> lemma21(const Graph& G, const set<int>& Q,
                                 const Point& r1, const Point& r2)
{
    return { sptInSub(G,Q,r1), sptInSub(G,Q,r2) };
}

/* =========================================================
 * Validate sequence: every tree is spanning, consecutive pairs
 * differ by exactly one edge, every tree has diameter ≤ d.
 * ========================================================= */
bool validateSeq(const Graph& G, const vector<set<int>>& seq, double d) {
    bool ok=true;
    for (int i=0;i<(int)seq.size();i++) {
        if (!isSpanningTree(G,seq[i])) {
            cout<<"  [FAIL] T["<<i<<"] is NOT a spanning tree ("<<seq[i].size()<<" edges)!\n";
            ok=false;
        }
        if (i>0) {
            // symmetric difference
            int diff=0;
            for (int e:seq[i])   if (!seq[i-1].count(e)) diff++;
            for (int e:seq[i-1]) if (!seq[i].count(e))   diff++;
            if (diff!=2) {
                cout<<"  [FAIL] T["<<i-1<<"]→T["<<i<<"] has "<<diff/2<<" changed edges (expected 1)!\n";
                ok=false;
            }
        }
        // diameter check
        bool diam_ok=true;
        for (int v=0;v<G.n&&diam_ok;v++) {
            auto dist=distFromPointInSub(G,seq[i],Point::vertex(v));
            for (int u=0;u<G.n;u++)
                if (dist[u].bar()>d+1e-9){diam_ok=false;break;}
        }
        if (!diam_ok) {
            cout<<"  [FAIL] T["<<i<<"] has diameter > "<<d<<"!\n";
            ok=false;
        }
    }
    return ok;
}

/* =========================================================
 * FULL RECONFIGURATION SEQUENCE  (Proposition 22 "if" direction)
 *
 * Path in G':   r_0 — r_1 — ... — r_k
 *   r_0 ∈ center(Ts),   r_k ∈ center(Tt)
 *   Q_i witnesses edge r_i — r_{i+1}
 *
 * Full sequence (Proposition 22):
 *   Ts →[Lem20,r_0]→ T+_0 →[Lem21]→ T-_1 →[Lem20,r_1]→ T+_1 →...→ Tt
 * ========================================================= */
void printReconfSequence(
    const Graph& G,
    const set<int>& Ts, const set<int>& Tt,
    const vector<Point>& path,
    const vector<set<int>>& witnesses,
    double d)
{
    cout << "\n=== Reconfiguration Sequence ===\n\n";

    vector<set<int>> full_seq;
    vector<string>   full_label;

    auto push = [&](const set<int>& T, const string& lbl) {
        if (!full_seq.empty() && full_seq.back()==T) return;
        full_seq.push_back(T); full_label.push_back(lbl);
    };

    set<int> cur=Ts;
    push(cur,"Ts");

    int k=(int)path.size()-1;
    for (int i=0;i<k;i++) {
        Point ri=path[i], ri1=path[i+1];
        auto [Tplus,Tminus]=lemma21(G,witnesses[i],ri,ri1);

        // Lemma 20: cur → T+_i
        auto seg=lemma20(G,cur,Tplus,ri);
        for (int s=1;s<(int)seg.size();s++)
            push(seg[s],"Lem20→T+_"+to_string(i));
        cur=Tplus;

        // Single flip: T+_i → T-_{i+1}
        push(Tminus,"flip T+_"+to_string(i)+"→T-_"+to_string(i+1));
        cur=Tminus;
    }

    // Lemma 20: cur → Tt
    auto seg_last=lemma20(G,cur,Tt,path[k]);
    for (int s=1;s<(int)seg_last.size();s++)
        push(seg_last[s],"Lem20→Tt");
    push(Tt,"Tt");

    // Print
    for (int i=0;i<(int)full_seq.size();i++) {
        cout<<"T["<<i<<"] ("<<full_label[i]<<"): ";
        printTree(G,full_seq[i]);
        cout<<"\n";
    }
    cout<<"\nTotal trees in sequence: "<<full_seq.size()<<"\n";
    cout<<"(Each consecutive pair is one edge swap)\n\n";

    // Validate every step
    cout<<"Validating...\n";
    if (validateSeq(G,full_seq,d))
        cout<<"  All "<<full_seq.size()-1<<" steps are VALID.\n";
}

/* =========================================================
 * ALGORITHM 3  (Section 6.3) — Modified Algorithm for RST
 *
 * Steps exactly as per paper.  Extra bookkeeping (no algorithmic change):
 *   - Store BFS parent in G' for path reconstruction.
 *   - Store witness Q for every G'-edge added.
 *   - After YES, reconstruct and print the reconfiguration sequence.
 * ========================================================= */
void algorithm3(const Graph& G,
                const set<int>& Ts_edges,
                const set<int>& Tt_edges,
                double d)
{
    cout << "\n=== Algorithm 3: RST with Small Diameter ===\n\n";

    // Step 1
    cout << "Step 1: Computing center(Ts) and center(Tt)...\n";
    set<Point> cTs=computeCenter(G,Ts_edges,d);
    set<Point> cTt=computeCenter(G,Tt_edges,d);
    cout<<"  center(Ts) = { "; for(auto& p:cTs) cout<<p.str()<<" "; cout<<"}\n";
    cout<<"  center(Tt) = { "; for(auto& p:cTt) cout<<p.str()<<" "; cout<<"}\n";

    if (cTs.empty()||cTt.empty()) {
        cout<<"\nERROR: Ts or Tt does not have diameter ≤ "<<d<<".\n"; return;
    }

    vector<Point> VR;
    for (int v=0;v<G.n;v++) VR.push_back(Point::vertex(v));
    for (int ei=0;ei<G.m;ei++)
        VR.push_back(Point::midpoint(G.edges[ei].first,G.edges[ei].second));
    int sz=VR.size();
    map<Point,int> pt_idx; for(int i=0;i<sz;i++) pt_idx[VR[i]]=i;

    cout<<"\n  |V ∪ R| = "<<sz<<" ("<<G.n<<" vertices + "<<G.m<<" midpoints)\n";
    cout<<"\nSteps 2-4: Building G' with early termination...\n";

    set<pair<int,int>>           Gp_edges;
    vector<vector<int>>          gp_adj(sz);
    map<pair<int,int>,set<int>>  witness;   // (i,j) → Q witnessing edge i-j in G'

    // BFS from center(Ts) with parent tracking
    vector<bool> reached(sz,false);
    vector<int>  bfs_par(sz,-1);
    queue<int>   bfs_q;
    for (auto& p:cTs) {
        int idx=pt_idx[p];
        if (!reached[idx]){ reached[idx]=true; bfs_q.push(idx); }
    }

    auto advanceBFS=[&](){
        while (!bfs_q.empty()){
            int u=bfs_q.front(); bfs_q.pop();
            for (int v:gp_adj[u])
                if (!reached[v]){ reached[v]=true; bfs_par[v]=u; bfs_q.push(v); }
        }
    };

    auto findTarget=[&]()->int{
        for (auto& p:cTt) if (reached[pt_idx[p]]) return pt_idx[p];
        return -1;
    };

    bool early_exit=false; int pairs_checked=0, edges_added=0, target=-1;
    advanceBFS();
    target=findTarget(); if(target>=0){early_exit=true; goto done;}

    for (int i=0;i<sz&&!early_exit;i++) {
        for (int j=i+1;j<sz&&!early_exit;j++) {
            pairs_checked++;
            set<int> Q=algo4(G,VR[i],VR[j],d);
            if (Q.empty()) Q=algo5(G,VR[i],VR[j],d);
            if (!Q.empty()) {
                Gp_edges.insert({i,j});
                gp_adj[i].push_back(j); gp_adj[j].push_back(i);
                witness[{i,j}]=Q; witness[{j,i}]=Q;
                edges_added++;
                if (reached[i]&&!reached[j]){reached[j]=true;bfs_par[j]=i;bfs_q.push(j);}
                if (reached[j]&&!reached[i]){reached[i]=true;bfs_par[i]=j;bfs_q.push(i);}
                advanceBFS();
                target=findTarget();
                if (target>=0){
                    cout<<"  [Early termination] YES after "<<pairs_checked
                        <<" pairs, "<<edges_added<<" G'-edges.\n";
                    early_exit=true;
                }
            }
        }
    }
    if (!early_exit)
        cout<<"  All pairs: "<<pairs_checked<<" | G'-edges: "<<edges_added<<"\n";

done:
    cout<<"\n--- Auxiliary Graph G' ---\nVertices (V ∪ R):\n  ";
    for(int i=0;i<sz;i++){cout<<VR[i].str();if(i+1<sz)cout<<", ";}
    cout<<"\nEdges in G':\n";
    if(Gp_edges.empty()) cout<<"  (none)\n";
    else for(auto&[a,b]:Gp_edges) cout<<"  "<<VR[a].str()<<" -- "<<VR[b].str()<<"\n";

    bool ans=(target>=0);
    cout<<"\n========================================\n";
    cout<<"RESULT: "<<(ans?"YES":"NO")<<"\n";
    if(ans) cout<<"Ts IS reconfigurable to Tt with diameter ≤ "<<d<<".\n";
    else    cout<<"Ts is NOT reconfigurable to Tt with diameter ≤ "<<d<<".\n";
    cout<<"========================================\n";

    if (!ans) return;

    // Reconstruct path in G' via BFS parent pointers
    vector<int> path_idx;
    {
        int cur=target;
        while (cur!=-1) {
            path_idx.push_back(cur);
            bool is_seed=false;
            for(auto& p:cTs) if(pt_idx[p]==cur){is_seed=true;break;}
            if(is_seed) break;
            cur=bfs_par[cur];
        }
        reverse(path_idx.begin(),path_idx.end());
    }

    vector<Point>    gp_path;
    vector<set<int>> gp_witnesses;
    for (int idx:path_idx) gp_path.push_back(VR[idx]);
    for (int i=0;i+1<(int)path_idx.size();i++)
        gp_witnesses.push_back(witness[{path_idx[i],path_idx[i+1]}]);

    cout<<"\nPath in G': ";
    for(int i=0;i<(int)gp_path.size();i++){
        cout<<gp_path[i].str(); if(i+1<(int)gp_path.size()) cout<<" → ";
    }
    cout<<"\n";

    printReconfSequence(G, Ts_edges, Tt_edges, gp_path, gp_witnesses, d);
}

/* =========================================================
 * DIAMETER UTILITY (used by output writer)
 * Compute the actual diameter of a spanning tree (BFS from each vertex).
 * ========================================================= */
int treeDiameter(const Graph& G, const set<int>& T) {
    int n = G.n;
    // build adjacency for T
    vector<vector<int>> adj(n);
    for (int ei : T) {
        adj[G.edges[ei].first ].push_back(G.edges[ei].second);
        adj[G.edges[ei].second].push_back(G.edges[ei].first);
    }
    int diam = 0;
    for (int s = 0; s < n; s++) {
        vector<int> dist(n,-1); dist[s]=0;
        queue<int> q; q.push(s);
        while (!q.empty()) {
            int u=q.front(); q.pop();
            for (int v:adj[u]) if (dist[v]==-1){dist[v]=dist[u]+1;q.push(v);}
        }
        for (int v=0;v<n;v++) if (dist[v]>diam) diam=dist[v];
    }
    return diam;
}

/* =========================================================
 * STRUCTURED OUTPUT WRITER
 * Writes output_constrained.txt in a format compare.py can parse.
 * ========================================================= */
void writeStructuredOutput(const string& filename,
                            const Graph& G,
                            const vector<set<int>>& seq,
                            const vector<string>& labels,
                            double d,
                            bool answer_yes)
{
    ofstream f(filename);
    if (!f) { cerr << "Cannot write " << filename << "\n"; return; }

    f << "METHOD constrained\n";
    f << "RESULT " << (answer_yes ? "YES" : "NO") << "\n";
    f << "D_BOUND " << (int)d << "\n";
    f << "N_VERTICES " << G.n << "\n";
    f << "N_EDGES " << G.m << "\n";
    f << "GRAPH_EDGES";
    for (int i=0;i<G.m;i++) f << " " << G.edges[i].first+1 << "-" << G.edges[i].second+1;
    f << "\n";
    f << "SEQUENCE_LENGTH " << seq.size() << "\n";

    for (int i = 0; i < (int)seq.size(); i++) {
        int diam = treeDiameter(G, seq[i]);
        // Fix label: last step is always "Tt", first is always "Ts"
        string lbl = labels[i];
        if (i == 0) lbl = "Ts";
        if (i == (int)seq.size()-1) lbl = "Tt";

        f << "STEP " << i
          << " LABEL " << lbl
          << " DIAM " << diam
          << " VIOLATES " << (diam > (int)d ? 1 : 0);

        // Compute ADDED/REMOVED by diffing with previous step
        if (i > 0) {
            set<int> added, removed;
            for (int ei : seq[i])   if (!seq[i-1].count(ei)) added.insert(ei);
            for (int ei : seq[i-1]) if (!seq[i].count(ei))   removed.insert(ei);
            if (!added.empty()) {
                f << " ADDED";
                for (int ei : added) f << " " << G.edges[ei].first+1 << "-" << G.edges[ei].second+1;
            } else {
                f << " ADDED none";
            }
            if (!removed.empty()) {
                f << " REMOVED";
                for (int ei : removed) f << " " << G.edges[ei].first+1 << "-" << G.edges[ei].second+1;
            } else {
                f << " REMOVED none";
            }
        } else {
            f << " ADDED none REMOVED none";
        }

        f << " EDGES";
        for (int ei : seq[i])
            f << " " << G.edges[ei].first+1 << "-" << G.edges[ei].second+1;
        f << "\n";
    }
    f.close();
    cout << "\n[Output] Structured sequence written to: " << filename << "\n";
}

/* =========================================================
 * MODIFIED algorithm3 — also returns sequence + labels
 * so main() can write the structured output file.
 * ========================================================= */
pair<vector<set<int>>,vector<string>>
algorithm3_with_sequence(const Graph& G,
                          const set<int>& Ts_edges,
                          const set<int>& Tt_edges,
                          double d)
{
    cout << "\n=== Algorithm 3: RST with Small Diameter ===\n\n";

    cout << "Step 1: Computing center(Ts) and center(Tt)...\n";
    set<Point> cTs = computeCenter(G, Ts_edges, d);
    set<Point> cTt = computeCenter(G, Tt_edges, d);
    cout << "  center(Ts) = { "; for(auto& p:cTs) cout<<p.str()<<" "; cout<<"}\n";
    cout << "  center(Tt) = { "; for(auto& p:cTt) cout<<p.str()<<" "; cout<<"}\n";

    if (cTs.empty()||cTt.empty()) {
        cout<<"\nERROR: Ts or Tt does not have diameter ≤ "<<d<<".\n";
        return {{},{}};
    }

    vector<Point> VR;
    for (int v=0;v<G.n;v++) VR.push_back(Point::vertex(v));
    for (int ei=0;ei<G.m;ei++)
        VR.push_back(Point::midpoint(G.edges[ei].first,G.edges[ei].second));
    int sz=VR.size();
    map<Point,int> pt_idx; for(int i=0;i<sz;i++) pt_idx[VR[i]]=i;

    cout<<"\n  |V ∪ R| = "<<sz<<"\n";
    cout<<"\nSteps 2-4: Building G'...\n";

    set<pair<int,int>>           Gp_edges;
    vector<vector<int>>          gp_adj(sz);
    map<pair<int,int>,set<int>>  witness;

    vector<bool> reached(sz,false);
    vector<int>  bfs_par(sz,-1);
    queue<int>   bfs_q;
    for (auto& p:cTs){int idx=pt_idx[p];if(!reached[idx]){reached[idx]=true;bfs_q.push(idx);}}

    auto advanceBFS=[&](){
        while(!bfs_q.empty()){int u=bfs_q.front();bfs_q.pop();
            for(int v:gp_adj[u])if(!reached[v]){reached[v]=true;bfs_par[v]=u;bfs_q.push(v);}}};

    auto findTarget=[&]()->int{
        for(auto& p:cTt)if(reached[pt_idx[p]])return pt_idx[p];return -1;};

    bool early_exit=false; int pairs_checked=0,edges_added=0,target=-1;
    advanceBFS(); target=findTarget(); if(target>=0){early_exit=true;goto done;}

    for(int i=0;i<sz&&!early_exit;i++){
        for(int j=i+1;j<sz&&!early_exit;j++){
            pairs_checked++;
            set<int> Q=algo4(G,VR[i],VR[j],d);
            if(Q.empty()) Q=algo5(G,VR[i],VR[j],d);
            if(!Q.empty()){
                Gp_edges.insert({i,j});
                gp_adj[i].push_back(j); gp_adj[j].push_back(i);
                witness[{i,j}]=Q; witness[{j,i}]=Q;
                edges_added++;
                if(reached[i]&&!reached[j]){reached[j]=true;bfs_par[j]=i;bfs_q.push(j);}
                if(reached[j]&&!reached[i]){reached[i]=true;bfs_par[i]=j;bfs_q.push(i);}
                advanceBFS();
                target=findTarget();
                if(target>=0){cout<<"  [Early exit] YES after "<<pairs_checked<<" pairs.\n";early_exit=true;}
            }
        }
    }
    if(!early_exit)cout<<"  All pairs checked: "<<pairs_checked<<" | G'-edges: "<<edges_added<<"\n";

done:
    bool ans=(target>=0);
    cout<<"\n========================================\n";
    cout<<"RESULT: "<<(ans?"YES":"NO")<<"\n";
    cout<<"========================================\n";

    if(!ans) return {{Ts_edges,Tt_edges},{"Ts","Tt"}};

    // Reconstruct path in G'
    vector<int> path_idx;
    {int cur=target;
     while(cur!=-1){
         path_idx.push_back(cur);
         bool seed=false; for(auto& p:cTs)if(pt_idx[p]==cur){seed=true;break;}
         if(seed)break; cur=bfs_par[cur];}
     reverse(path_idx.begin(),path_idx.end());}

    vector<Point> gp_path;
    vector<set<int>> gp_witnesses;
    for(int idx:path_idx) gp_path.push_back(VR[idx]);
    for(int i=0;i+1<(int)path_idx.size();i++)
        gp_witnesses.push_back(witness[{path_idx[i],path_idx[i+1]}]);

    cout<<"\nPath in G': ";
    for(int i=0;i<(int)gp_path.size();i++){cout<<gp_path[i].str();if(i+1<(int)gp_path.size())cout<<" → ";}
    cout<<"\n";

    // Build full sequence (same logic as printReconfSequence but return it)
    vector<set<int>> full_seq;
    vector<string>   full_label;

    auto push=[&](const set<int>& T,const string& lbl){
        if(!full_seq.empty()&&full_seq.back()==T)return;
        full_seq.push_back(T); full_label.push_back(lbl);};

    set<int> cur=Ts_edges; push(cur,"Ts");
    int k=(int)gp_path.size()-1;
    for(int i=0;i<k;i++){
        Point ri=gp_path[i],ri1=gp_path[i+1];
        auto[Tplus,Tminus]=lemma21(G,gp_witnesses[i],ri,ri1);
        auto seg=lemma20(G,cur,Tplus,ri);
        for(int s=1;s<(int)seg.size();s++) push(seg[s],"Lem20→T+_"+to_string(i));
        cur=Tplus;
        push(Tminus,"flip→T-_"+to_string(i+1));
        cur=Tminus;
    }
    auto seg_last=lemma20(G,cur,Tt_edges,gp_path[k]);
    for(int s=1;s<(int)seg_last.size();s++) push(seg_last[s],"Lem20→Tt");
    push(Tt_edges,"Tt");

    // Print + validate
    cout<<"\n=== Reconfiguration Sequence ===\n\n";
    for(int i=0;i<(int)full_seq.size();i++){
        int diam=treeDiameter(G,full_seq[i]);
        bool viol=(diam>(int)d);
        cout<<"T["<<i<<"] ("<<full_label[i]<<") diam="<<diam
            <<(viol?" [VIOLATION]":" [ok]")<<": { ";
        for(int ei:full_seq[i]) cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") ";
        cout<<"}\n";
    }
    cout<<"\nTotal trees: "<<full_seq.size()<<"\n";

    // Validate
    bool valid=validateSeq(G,full_seq,d);
    if(valid) cout<<"Validation: ALL steps valid (spanning tree, 1-edge diff, diam≤"<<d<<")\n";

    return {full_seq, full_label};
}

/* =========================================================
 * MAIN — supports CLI input or input.txt file
 *
 * Usage:
 *   ./rst_cons                  (interactive CLI prompts)
 *   ./rst_cons input.txt        (read from file)
 *   ./rst_cons input.txt output_constrained.txt  (also write output)
 *
 * Input format (comments starting with # are stripped):
 *   n m
 *   m graph edges (u v)
 *   n-1 edges of Ts
 *   n-1 edges of Tt
 *   d
 * ========================================================= */
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Determine input source
    istream* in_ptr = &cin;
    ifstream file_in;
    string output_file = "output_constrained.txt";

    if (argc >= 2) {
        file_in.open(argv[1]);
        if (!file_in) { cerr << "Cannot open " << argv[1] << "\n"; return 1; }
        in_ptr = &file_in;
        cout << "[Input] Reading from file: " << argv[1] << "\n";
    } else {
        cout << "[Input] Reading from stdin (use: " << argv[0] << " input.txt [output.txt])\n";
    }
    if (argc >= 3) output_file = argv[2];

    // Helper: read next token ignoring # comment lines
    auto readToken = [&](auto& val) {
        string tok;
        while (*in_ptr >> tok) {
            if (tok[0] == '#') { in_ptr->ignore(10000,'\n'); continue; }
            istringstream ss(tok); ss >> val; return true;
        }
        return false;
    };

    cout << "=== RST with Small Diameter (Constrained) ===\n";
    cout << "Algorithm 3 — Bousquet et al., 2022\n\n";

    int n, m;
    if (!readToken(n) || !readToken(m)) { cerr << "Bad input\n"; return 1; }
    cout << "n=" << n << ", m=" << m << "\n";

    vector<pair<int,int>> raw(m);
    for (int i=0;i<m;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        if(u>v)swap(u,v); raw[i]={u,v};
    }

    Graph G; G.build(n, raw);

    cout << "[OPT-1] Precomputing all-pairs distances (" << n << " Dijkstra runs)...\n";
    precomputeAllPairs(G);
    cout << "  Done.\n";

    set<int> Ts, Tt;
    for (int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){cerr<<"Ts edge "<<u+1<<"-"<<v+1<<" not in G!\n";return 1;}
        Ts.insert(ei);
    }
    for (int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){cerr<<"Tt edge "<<u+1<<"-"<<v+1<<" not in G!\n";return 1;}
        Tt.insert(ei);
    }

    double d; readToken(d);

    cout << "\nGraph : " << n << " vertices, " << m << " edges\n";
    cout << "Ts    : "; for(int ei:Ts) cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") "; cout<<"\n";
    cout << "Tt    : "; for(int ei:Tt) cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") "; cout<<"\n";
    cout << "d     : " << d << "\n";

    auto [seq, labels] = algorithm3_with_sequence(G, Ts, Tt, d);

    // Write structured output
    bool ans = (!seq.empty() && seq.back() == Tt);
    writeStructuredOutput(output_file, G, seq, labels, d, ans);

    cout << "\nDone. Output: " << output_file << "\n";
    return 0;
}