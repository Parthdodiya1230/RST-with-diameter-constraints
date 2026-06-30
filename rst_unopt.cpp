/*
 * Reconfiguration of Spanning Trees with Small Diameter
 * UNOPTIMIZED version — OPT-3 (early termination) REMOVED.
 *
 * OPT-1 (all-pairs Dijkstra precompute) and OPT-2 (center reuse) are kept.
 * OPT-3 removed: the G' construction loop always checks ALL (i,j) pairs
 * before deciding YES/NO, instead of stopping as soon as BFS from center(Ts)
 * reaches center(Tt).
 *
 * Based on: Bousquet et al., "Reconfiguration of Spanning Trees with Degree
 *           Constraint or Diameter Constraint", arXiv:2201.04354v1, 2022.
 *
 * INPUT (1-based vertex indexing):
 *   n m
 *   u_1 v_1  ...  (m edges of G)
 *   u_1 v_1  ...  (n-1 edges of Ts)
 *   u_1 v_1  ...  (n-1 edges of Tt)
 *   d
 *
 * OUTPUT:
 *   YES/NO  <time_ms>
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
#include <chrono>
using namespace std;

/* =========================================================
 * PERTURBED LENGTH
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

    double bar() const { return base; }
};

/* =========================================================
 * POINT in V ∪ R
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
};

/* =========================================================
 * GRAPH
 * ========================================================= */
struct Graph {
    int n, m;
    vector<pair<int,int>>         edges;
    vector<PLen>                  elen;
    vector<vector<pair<int,int>>> adj;
    map<pair<int,int>,int>        edge_map;

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
 * OPT-1: ALL-PAIRS SHORTEST DISTANCES (kept)
 * ========================================================= */
vector<vector<PLen>> gdist;
vector<vector<int>>  gpar;

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

void precomputeAllPairs(const Graph& G) {
    gdist.resize(G.n);
    gpar.resize(G.n);
    for (int s=0;s<G.n;s++) {
        auto [d,p]=dijkstraG(G,s);
        gdist[s]=d; gpar[s]=p;
    }
}

PLen distVV(int s, int v) { return gdist[s][v]; }

PLen distPV(const Graph& G, const Point& p, int v) {
    if (p.type==Point::VERTEX) return distVV(p.v1, v);
    int ei=G.edgeIdx(p.v1,p.v2);
    PLen half=G.elen[ei].half();
    PLen via_a=half+distVV(p.v1,v);
    PLen via_b=half+distVV(p.v2,v);
    return (via_a<via_b) ? via_a : via_b;
}

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
 * Sub-graph utilities
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

PLen eccentricity(const Graph& G, const set<int>& subedges, const Point& p) {
    auto dist=distFromPointInSub(G,subedges,p);
    PLen mx=PLen::zero();
    for (int v=0;v<G.n;v++) if (mx<dist[v]) mx=dist[v];
    return mx;
}

/* OPT-2: center computed once per tree (kept) */
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

vector<PLen> computeLabel(const Graph& G, const set<int>& subedges,
                           const Point& r1, const Point& r2) {
    auto d1=distFromPointInSub(G,subedges,r1);
    auto d2=distFromPointInSub(G,subedges,r2);
    vector<PLen> label(G.n);
    for (int v=0;v<G.n;v++) label[v]=(d1[v]<d2[v]) ? d2[v] : d1[v];
    return label;
}

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

bool isGood(const Graph& G, const set<int>& subedges,
            const Point& r1, const Point& r2) {
    auto label=computeLabel(G,subedges,r1,r2);
    const double EPS=1e-9;
    for (int ei=0;ei<G.m;ei++) {
        int u=G.edges[ei].first, v=G.edges[ei].second;
        if (label[v].base > (label[u]+G.elen[ei]).base+EPS) return false;
        if (label[u].base > (label[v]+G.elen[ei]).base+EPS) return false;
    }
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

bool bothCenters(const Graph& G, const set<int>& subedges,
                 const Point& r1, const Point& r2, double d) {
    if (eccentricity(G,subedges,r1).bar() > d/2.0+1e-9) return false;
    if (eccentricity(G,subedges,r2).bar() > d/2.0+1e-9) return false;
    return true;
}

/* =========================================================
 * Augmented graph G+
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
 * Algorithm 4
 * ========================================================= */
set<int> algo4(const Graph& G, const Point& r1, const Point& r2, double d) {
    for (int e_idx=0;e_idx<G.m;e_idx++) {
        int v1=G.edges[e_idx].first, v2=G.edges[e_idx].second;
        bool r2_pe=(r2.type==Point::MIDPOINT && G.edgeIdx(r2.v1,r2.v2)==e_idx);
        bool r1_pe=(r1.type==Point::MIDPOINT && G.edgeIdx(r1.v1,r1.v2)==e_idx);
        PLen l_e1 = r2_pe ? G.elen[e_idx].half() : distPV(G,r2,v2) + G.elen[e_idx];
        PLen l_e2 = r1_pe ? G.elen[e_idx].half() : distPV(G,r1,v1) + G.elen[e_idx];
        AugGraph aug; aug.init(G.n);
        for (int i=0;i<G.m;i++) aug.addOrigEdge(G.edges[i].first,G.edges[i].second,i,G);
        aug.addNewEdge(v1, l_e1);
        aug.addNewEdge(v2, l_e2);
        vector<int> par=sptFromR(aug);
        set<int> Q=collectFromSPT(G,par);
        Q.insert(e_idx);
        if (!isSpanningTree(G,Q)) continue;
        if (bothCenters(G,Q,r1,r2,d) && isGood(G,Q,r1,r2))
            return Q;
    }
    return {};
}

/* =========================================================
 * Algorithm 5
 * ========================================================= */
set<int> algo5(const Graph& G, const Point& r1, const Point& r2, double d) {
    int n=G.n;
    for (int e_idx=0;e_idx<G.m;e_idx++) {
        int v1=G.edges[e_idx].first, v2=G.edges[e_idx].second;
        for (int ep_idx=0;ep_idx<G.m;ep_idx++) {
            if (ep_idx==e_idx) continue;
            int vp1=G.edges[ep_idx].first, vp2=G.edges[ep_idx].second;
            int f_total=G.m+2;
            for (int f1=0;f1<f_total;f1++) {
                bool f1_new=(f1>=G.m);
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
                    set<int> J2;
                    if (f2_new) {
                        for (int ei:shortestPathEdges(G,v2,vp2)) J2.insert(ei);
                    } else {
                        int u2=G.edges[f2].first, up2=G.edges[f2].second;
                        for (int ei:shortestPathEdges(G,v2,u2))   J2.insert(ei);
                        J2.insert(f2);
                        for (int ei:shortestPathEdges(G,up2,vp2)) J2.insert(ei);
                    }
                    set<int> C;
                    for (int ei:J1) C.insert(ei);
                    C.insert(ep_idx);
                    for (int ei:J2) C.insert(ei);
                    C.insert(e_idx);
                    map<int,int> deg;
                    for (int ei:C){ deg[G.edges[ei].first]++; deg[G.edges[ei].second]++; }
                    bool is_cycle=true;
                    for (auto& [v,dv]:deg) if(dv!=2){is_cycle=false;break;}
                    if (!is_cycle) continue;
                    auto onC=[&](const Point& r){
                        return (r.type==Point::VERTEX) ? deg.count(r.v1)>0
                               : C.count(G.edgeIdx(r.v1,r.v2))>0;
                    };
                    if (!onC(r1)||!onC(r2)) continue;
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
                    PLen l_e1 =r2_pe ? G.elen[e_idx].half() :distPonC(r2,v2) +G.elen[e_idx];
                    PLen l_e2 =r1_pe ? G.elen[e_idx].half() :distPonC(r1,v1) +G.elen[e_idx];
                    PLen l_ep1=r2_pep? G.elen[ep_idx].half():distPonC(r2,vp2)+G.elen[ep_idx];
                    PLen l_ep2=r1_pep? G.elen[ep_idx].half():distPonC(r1,vp1)+G.elen[ep_idx];
                    AugGraph aug; aug.init(n);
                    for (int i=0;i<G.m;i++) aug.addOrigEdge(G.edges[i].first,G.edges[i].second,i,G);
                    aug.addNewEdge(v1,  l_e1);
                    aug.addNewEdge(v2,  l_e2);
                    aug.addNewEdge(vp1, l_ep1);
                    aug.addNewEdge(vp2, l_ep2);
                    vector<int> par=sptFromR(aug);
                    set<int> Q=collectFromSPT(G,par);
                    Q.insert(e_idx);
                    Q.insert(ep_idx);
                    if (!f1_new) Q.insert(f1);
                    if (!f2_new) Q.insert(f2);
                    if (!isPseudotree(G,Q)) continue;
                    if (findCycleEdges(G,Q).empty()) continue;
                    if (bothCenters(G,Q,r1,r2,d) && isGood(G,Q,r1,r2))
                        return Q;
                }
            }
        }
    }
    return {};
}

/* =========================================================
 * Algorithm 3 — NO OPT-3
 *
 * Difference from the optimized version:
 *   - The inner pair loop does NOT call return "YES" early when BFS
 *     from center(Ts) first reaches center(Tt).
 *   - After ALL (i,j) pairs have been processed, a single final
 *     BFS reachability check determines the answer.
 * ========================================================= */
string runAlgo3(const Graph& G,
                const set<int>& Ts_edges,
                const set<int>& Tt_edges,
                double d)
{
    // Step 1: centers (OPT-2 kept)
    set<Point> cTs = computeCenter(G, Ts_edges, d);
    set<Point> cTt = computeCenter(G, Tt_edges, d);
    if (cTs.empty() || cTt.empty()) return "NO";

    vector<Point> VR;
    for (int v=0;v<G.n;v++) VR.push_back(Point::vertex(v));
    for (int ei=0;ei<G.m;ei++)
        VR.push_back(Point::midpoint(G.edges[ei].first,G.edges[ei].second));
    int sz = VR.size();
    map<Point,int> pt_idx;
    for (int i=0;i<sz;i++) pt_idx[VR[i]]=i;

    // G' adjacency
    vector<vector<int>> gp_adj(sz);

    // BFS state from center(Ts) — updated incrementally but no early exit
    vector<bool> reached(sz, false);
    queue<int>   bfs_q;
    for (auto& p : cTs) {
        int idx = pt_idx[p];
        if (!reached[idx]) { reached[idx]=true; bfs_q.push(idx); }
    }

    auto advanceBFS = [&]() {
        while (!bfs_q.empty()) {
            int u = bfs_q.front(); bfs_q.pop();
            for (int v : gp_adj[u])
                if (!reached[v]) { reached[v]=true; bfs_q.push(v); }
        }
    };

    // Seed BFS with initial center(Ts) nodes
    advanceBFS();

    // Steps 2–4: enumerate ALL pairs — NO early termination (OPT-3 removed)
    for (int i = 0; i < sz; i++) {
        for (int j = i+1; j < sz; j++) {
            set<int> Q = algo4(G, VR[i], VR[j], d);
            if (Q.empty()) Q = algo5(G, VR[i], VR[j], d);
            if (!Q.empty()) {
                gp_adj[i].push_back(j);
                gp_adj[j].push_back(i);
                if (reached[i] && !reached[j]) { reached[j]=true; bfs_q.push(j); }
                if (reached[j] && !reached[i]) { reached[i]=true; bfs_q.push(i); }
                advanceBFS();
                // OPT-3 REMOVED: do NOT return "YES" here even if target reached
            }
        }
    }

    // Final reachability check — only after all pairs processed
    for (auto& p : cTt)
        if (reached[pt_idx[p]]) return "YES";
    return "NO";
}

/* =========================================================
 * MAIN
 * ========================================================= */
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    istream* in_ptr = &cin;
    ifstream file_in;
    if (argc >= 2) {
        file_in.open(argv[1]);
        if (!file_in) { cerr << "Cannot open " << argv[1] << "\n"; return 1; }
        in_ptr = &file_in;
    }

    auto readToken = [&](auto& val) {
        string tok;
        while (*in_ptr >> tok) {
            if (tok[0] == '#') { in_ptr->ignore(10000,'\n'); continue; }
            istringstream ss(tok); ss >> val; return true;
        }
        return false;
    };

    int n, m;
    if (!readToken(n) || !readToken(m)) { cerr << "Bad input\n"; return 1; }

    vector<pair<int,int>> raw(m);
    for (int i=0;i<m;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        if(u>v) swap(u,v); raw[i]={u,v};
    }

    Graph G; G.build(n, raw);
    precomputeAllPairs(G);

    set<int> Ts, Tt;
    for (int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){ cerr<<"Ts edge not in G\n"; return 1; }
        Ts.insert(ei);
    }
    for (int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){ cerr<<"Tt edge not in G\n"; return 1; }
        Tt.insert(ei);
    }

    double d; readToken(d);

    auto t0 = chrono::high_resolution_clock::now();
    string ans = runAlgo3(G, Ts, Tt, d);
    auto t1 = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(t1 - t0).count();

    cout << ans << " " << fixed << ms << "\n";
    return 0;
}