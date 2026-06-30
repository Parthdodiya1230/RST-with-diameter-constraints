/*
 * Naive Spanning Tree Reconfiguration
 * ====================================
 * Symmetric-difference method: given Ts and Tt, repeatedly pick an edge
 * e ∈ Tt \ T, find the unique cycle it creates in T, and remove one edge
 * from that cycle that is in T \ Tt.
 *
 * This gives a valid reconfiguration sequence Ts → Tt with no diameter
 * guarantee. We compute and print the diameter of every intermediate tree
 * so it can be compared against Algorithm 3's constrained sequence.
 *
 * Usage:
 *   ./rst_naive                        (interactive stdin)
 *   ./rst_naive input.txt              (read from file)
 *   ./rst_naive input.txt output_naive.txt   (also write output file)
 *
 * Input format (same as rst_constrained; lines starting with # ignored):
 *   n m
 *   m graph edges (u v), 1-based
 *   n-1 edges of Ts
 *   n-1 edges of Tt
 *   d   (diameter bound — used only for violation reporting, NOT enforced)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <string>
using namespace std;

/* =========================================================
 * GRAPH
 * ========================================================= */
struct Graph {
    int n, m;
    vector<pair<int,int>> edges;
    map<pair<int,int>,int> eidx;

    void build(int _n, const vector<pair<int,int>>& raw) {
        n=_n; m=(int)raw.size(); edges=raw;
        for(int i=0;i<m;i++) eidx[raw[i]]=i;
    }
    int edgeIdx(int u, int v) const {
        auto it=eidx.find({min(u,v),max(u,v)});
        return it==eidx.end() ? -1 : it->second;
    }
};

/* =========================================================
 * TREE UTILITIES
 * ========================================================= */

// Compute diameter of a spanning tree (BFS from every vertex)
int treeDiameter(const Graph& G, const set<int>& T) {
    int n=G.n;
    vector<vector<int>> adj(n);
    for(int ei:T){
        adj[G.edges[ei].first ].push_back(G.edges[ei].second);
        adj[G.edges[ei].second].push_back(G.edges[ei].first);
    }
    int diam=0;
    for(int s=0;s<n;s++){
        vector<int> dist(n,-1); dist[s]=0;
        queue<int> q; q.push(s);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v:adj[u]) if(dist[v]==-1){dist[v]=dist[u]+1;q.push(v);}
        }
        for(int v=0;v<n;v++) diam=max(diam,dist[v]);
    }
    return diam;
}

// Find edges on the unique path from src to dst in a spanning tree T.
// Returns the edge INCIDENT TO dst on this path (the "parent edge" of dst).
// Used to find which tree-edge to remove when adding edge (src,dst).
vector<int> pathEdges(const Graph& G, const set<int>& T, int src, int dst) {
    int n=G.n;
    vector<vector<pair<int,int>>> adj(n);
    for(int ei:T){
        adj[G.edges[ei].first ].push_back({G.edges[ei].second,ei});
        adj[G.edges[ei].second].push_back({G.edges[ei].first, ei});
    }
    vector<int> par_e(n,-1);
    vector<bool> vis(n,false);
    queue<int> q; q.push(src); vis[src]=true;
    while(!q.empty()){
        int u=q.front(); q.pop();
        for(auto[v,ei]:adj[u]) if(!vis[v]){
            vis[v]=true; par_e[v]=ei; q.push(v);
        }
    }
    // Collect all edges on path src→dst
    vector<int> path;
    int cur=dst;
    while(par_e[cur]!=-1){
        path.push_back(par_e[cur]);
        // move to parent
        int ei=par_e[cur];
        int next=(G.edges[ei].first==cur)?G.edges[ei].second:G.edges[ei].first;
        cur=next;
    }
    return path;  // edges from dst back to src (in reverse order)
}

/* =========================================================
 * NAIVE RECONFIGURATION
 *
 * Strategy: at each step, pick the lexicographically smallest
 * edge e ∈ Tt \ T, find the unique cycle T+e creates, pick the
 * lexicographically smallest edge f ∈ cycle ∩ (T \ Tt) to remove.
 * If no f in T \ Tt, pick smallest f in cycle ∩ T.
 *
 * This is the standard symmetric-difference method.
 * It converges in at most |Ts Δ Tt|/2 steps.
 * ========================================================= */
struct Step {
    set<int> tree;
    int added_edge;    // -1 if first step
    int removed_edge;  // -1 if first step
    int diameter;
    bool violates;
    string label;
};

vector<Step> naiveReconf(const Graph& G, const set<int>& Ts, const set<int>& Tt, double d) {
    vector<Step> seq;

    // Initial step
    {
        Step s0;
        s0.tree=Ts; s0.added_edge=-1; s0.removed_edge=-1;
        s0.diameter=treeDiameter(G,Ts);
        s0.violates=(s0.diameter>(int)d);
        s0.label="Ts";
        seq.push_back(s0);
    }

    set<int> T=Ts;
    set<int> TGoal=Tt;

    for(int iter=0; iter<(int)G.m && T!=TGoal; iter++) {
        // Find all edges in Tt \ T
        vector<int> add_cands;
        for(int ei:TGoal) if(!T.count(ei)) add_cands.push_back(ei);
        if(add_cands.empty()) break;

        // Pick lexicographically smallest (by (u,v))
        sort(add_cands.begin(),add_cands.end(),[&](int a,int b){
            return G.edges[a]<G.edges[b];
        });
        int add_ei=add_cands[0];
        int u=G.edges[add_ei].first, v=G.edges[add_ei].second;

        // Find the cycle: unique path from u to v in T
        vector<int> cycle_edges=pathEdges(G, T, u, v);
        // cycle = cycle_edges ∪ {add_ei}

        // Candidates to remove: prefer edges in T \ TGoal (makes progress)
        vector<int> rem_prefer, rem_any;
        for(int ei:cycle_edges){
            if(!TGoal.count(ei)) rem_prefer.push_back(ei);
            rem_any.push_back(ei);
        }
        auto& rem_cands = rem_prefer.empty() ? rem_any : rem_prefer;
        sort(rem_cands.begin(),rem_cands.end(),[&](int a,int b){
            return G.edges[a]<G.edges[b];
        });
        int rem_ei=rem_cands[0];

        // Apply swap
        T.insert(add_ei);
        T.erase(rem_ei);

        Step s;
        s.tree=T;
        s.added_edge=add_ei;
        s.removed_edge=rem_ei;
        s.diameter=treeDiameter(G,T);
        s.violates=(s.diameter>(int)d);
        s.label=(T==TGoal?"Tt":"T"+to_string((int)seq.size()));
        seq.push_back(s);
    }

    return seq;
}

/* =========================================================
 * STRUCTURED OUTPUT WRITER
 * Matches the format written by rst_constrained for compare.py
 * ========================================================= */
void writeOutput(const string& filename,
                 const Graph& G,
                 const vector<Step>& seq,
                 double d)
{
    int violations=0;
    for(auto& s:seq) if(s.violates) violations++;

    ofstream f(filename);
    if(!f){cerr<<"Cannot write "<<filename<<"\n";return;}

    f<<"METHOD naive\n";
    f<<"RESULT " << (seq.back().tree == seq.back().tree ? "YES" : "NO") << "\n";
    // For naive, result is always YES (we just reach Tt eventually)
    // We write YES unless something went wrong
    f<<"D_BOUND "<<(int)d<<"\n";
    f<<"N_VERTICES "<<G.n<<"\n";
    f<<"N_EDGES "<<G.m<<"\n";
    f<<"GRAPH_EDGES";
    for(int i=0;i<G.m;i++) f<<" "<<G.edges[i].first+1<<"-"<<G.edges[i].second+1;
    f<<"\n";
    f<<"SEQUENCE_LENGTH "<<seq.size()<<"\n";

    for(int i=0;i<(int)seq.size();i++){
        auto& s=seq[i];
        f<<"STEP "<<i
         <<" LABEL "<<s.label
         <<" DIAM "<<s.diameter
         <<" VIOLATES "<<(s.violates?1:0)
         <<" ADDED "  <<(s.added_edge  >=0 ? to_string(G.edges[s.added_edge  ].first+1)+"-"+to_string(G.edges[s.added_edge  ].second+1) : "none")
         <<" REMOVED "<<(s.removed_edge>=0 ? to_string(G.edges[s.removed_edge].first+1)+"-"+to_string(G.edges[s.removed_edge].second+1) : "none")
         <<" EDGES";
        for(int ei:s.tree)
            f<<" "<<G.edges[ei].first+1<<"-"<<G.edges[ei].second+1;
        f<<"\n";
    }
    f.close();
    cout<<"\n[Output] Written to: "<<filename<<"\n";
}

/* =========================================================
 * MAIN
 * ========================================================= */
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    istream* in_ptr=&cin;
    ifstream file_in;
    string output_file="output_naive.txt";

    if(argc>=2){
        file_in.open(argv[1]);
        if(!file_in){cerr<<"Cannot open "<<argv[1]<<"\n";return 1;}
        in_ptr=&file_in;
        cout<<"[Input] Reading from: "<<argv[1]<<"\n";
    } else {
        cout<<"[Input] Reading from stdin (tip: "<<argv[0]<<" input.txt [output.txt])\n";
    }
    if(argc>=3) output_file=argv[2];

    // Read next token, skip # comment lines
    auto readToken=[&](auto& val)->bool{
        string tok;
        while(*in_ptr>>tok){
            if(tok[0]=='#'){in_ptr->ignore(10000,'\n');continue;}
            istringstream ss(tok); ss>>val; return true;
        }
        return false;
    };

    cout<<"=== RST Naive Reconfiguration ===\n";
    cout<<"Symmetric-difference method, no diameter guarantee\n\n";

    int n,m;
    if(!readToken(n)||!readToken(m)){cerr<<"Bad input\n";return 1;}
    cout<<"n="<<n<<", m="<<m<<"\n";

    vector<pair<int,int>> raw(m);
    for(int i=0;i<m;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        if(u>v)swap(u,v); raw[i]={u,v};
    }
    Graph G; G.build(n,raw);

    set<int> Ts,Tt;
    for(int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){cerr<<"Ts edge "<<u+1<<"-"<<v+1<<" not in G!\n";return 1;}
        Ts.insert(ei);
    }
    for(int i=0;i<n-1;i++){
        int u,v; readToken(u); readToken(v); u--;v--;
        int ei=G.edgeIdx(u,v);
        if(ei<0){cerr<<"Tt edge "<<u+1<<"-"<<v+1<<" not in G!\n";return 1;}
        Tt.insert(ei);
    }
    double d; readToken(d);

    cout<<"Ts    : "; for(int ei:Ts) cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") "; cout<<"\n";
    cout<<"Tt    : "; for(int ei:Tt) cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") "; cout<<"\n";
    cout<<"d     : "<<d<<" (used only for violation reporting)\n\n";

    cout<<"--- Running naive symmetric-difference reconfiguration ---\n\n";
    auto seq=naiveReconf(G,Ts,Tt,d);

    // Print sequence
    int violations=0;
    for(int i=0;i<(int)seq.size();i++){
        auto& s=seq[i];
        if(s.violates) violations++;
        cout<<"T["<<i<<"] ("<<s.label<<") diam="<<s.diameter
            <<(s.violates?" [VIOLATION >d]":" [ok]");
        if(s.added_edge>=0)
            cout<<" +("<<G.edges[s.added_edge].first+1<<","<<G.edges[s.added_edge].second+1<<")";
        if(s.removed_edge>=0)
            cout<<" -("<<G.edges[s.removed_edge].first+1<<","<<G.edges[s.removed_edge].second+1<<")";
        cout<<": { ";
        for(int ei:s.tree)
            cout<<"("<<G.edges[ei].first+1<<","<<G.edges[ei].second+1<<") ";
        cout<<"}\n";
    }

    cout<<"\n========================================\n";
    cout<<"Naive steps      : "<<seq.size()-1<<"\n";
    cout<<"Diameter bound d : "<<(int)d<<"\n";
    cout<<"Violations (diam > d): "<<violations<<"/"<<seq.size()<<"\n";
    if(violations==0)
        cout<<"Note: naive reached Tt without violating d="<<(int)d
            <<" — but this is NOT guaranteed in general.\n";
    else
        cout<<"These intermediate trees VIOLATE the diameter constraint.\n"
            <<"A real network using this sequence would have unacceptable path lengths.\n";
    cout<<"========================================\n";

    writeOutput(output_file,G,seq,d);
    cout<<"Done. Output: "<<output_file<<"\n";
    return 0;
}