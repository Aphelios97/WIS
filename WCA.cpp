#include <iostream>
#include <vector>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "Psapi.lib")


#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;


struct TemporalWedge {
    int a;          // position in chronological order
    int b;          // position in chronological order
    int inter_ab;   // |e_a ∩ e_b|
    bool active;
};

static int triple_intersection_size(
        const unordered_set<int>& A,
        const unordered_set<int>& B,
        const unordered_set<int>& C) {

    const unordered_set<int>* smallest = &A;
    if (B.size() < smallest->size()) smallest = &B;
    if (C.size() < smallest->size()) smallest = &C;

    int cnt = 0;
    for (int x : *smallest) {
        if (A.count(x) && B.count(x) && C.count(x)) {
            ++cnt;
        }
    }
    return cnt;
}


static double peak_memory_mb() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(
            GetCurrentProcess(),
            &pmc,
            sizeof(pmc))) {
        return static_cast<double>(pmc.PeakWorkingSetSize)
               / 1024.0;
    }
    return -1.0;
}

int main(int argc, char* argv[]) {
    clock_t run_start = clock();

    if (argc != 3) {
        cerr << "Usage: " << argv[0]
             << " <dataset> <delta>\\n";
        return 1;
    }

    string dataset = argv[1];
    double delta = stod(argv[2]);

    string graphFile = dataset + ".txt";

    cout << "====================================================\n";
    cout << "WCA Exact Temporal Hypergraph Motif Counting\n";
    cout << "====================================================\n";
    cout << "Dataset: " << dataset << "\n";
    cout << "delta: " << fixed << setprecision(0) << delta << "\n";


    vector<vector<int>> node2hyperedge;
    vector<vector<int>> hyperedge2node;
    vector<unordered_set<int>> hyperedge2node_set;
    vector<double> hyperedge2time;

    read_data(
        graphFile,
        node2hyperedge,
        hyperedge2node,
        hyperedge2node_set,
        hyperedge2time
    );

    const int V = (int)node2hyperedge.size();
    const int E = (int)hyperedge2node.size();

    if (E == 0) {
        cerr << "ERROR: empty temporal hypergraph.\n";
        return 1;
    }

    cout << "# Nodes: " << V << "\n";
    cout << "# Temporal hyperedges: " << E << "\n";

    vector<int> order(E);
    iota(order.begin(), order.end(), 0);

    sort(order.begin(), order.end(), [&](int x, int y) {
        if (hyperedge2time[x] == hyperedge2time[y]) {
            return x < y;
        }
        return hyperedge2time[x] < hyperedge2time[y];
    });

    vector<unsigned long long> motif_count(96, 0ULL);

    vector<unordered_set<int>> node2active(V);

    vector<unordered_set<int>> H_C(E);

    vector<TemporalWedge> Lambda_C;
    Lambda_C.reserve((size_t)E * 4);

    deque<int> active_queue;
    vector<char> active(E, 0);

    unsigned long long total_wedges_created = 0ULL;
    unsigned long long total_candidates = 0ULL;
    unsigned long long total_motifs = 0ULL;

    auto remove_expired_edge = [&](int p) {

        if (!active[p]) {
            return;
        }

        int global_id = order[p];


        vector<int> incident_wedges(
            H_C[p].begin(),
            H_C[p].end()
        );

        for (int wedge_id : incident_wedges) {

            if (wedge_id < 0 ||
                wedge_id >= (int)Lambda_C.size()) {
                continue;
            }

            TemporalWedge& w = Lambda_C[wedge_id];

            if (!w.active) {
                continue;
            }

            w.active = false;

            H_C[w.a].erase(wedge_id);
            H_C[w.b].erase(wedge_id);
        }

        H_C[p].clear();


        for (int node : hyperedge2node[global_id]) {
            if (node >= 0 && node < V) {
                node2active[node].erase(p);
            }
        }

        active[p] = 0;
    };


    for (int i = 0; i < E; ++i) {

        int global_i = order[i];
        double t_i = hyperedge2time[global_i];



        while (!active_queue.empty()) {

            int r = active_queue.front();
            int global_r = order[r];

            double t_r = hyperedge2time[global_r];

            if (t_i - t_r <= delta) {
                break;
            }

            active_queue.pop_front();
            remove_expired_edge(r);
        }


        unordered_map<int, int> neighbor_inter;

        for (int node : hyperedge2node[global_i]) {

            if (node < 0 || node >= V) {
                continue;
            }

            for (int m : node2active[node]) {

                if (active[m]) {
                    ++neighbor_inter[m];
                }
            }
        }

        vector<int> N_i;
        N_i.reserve(neighbor_inter.size());

        for (const auto& kv : neighbor_inter) {

            if (kv.second > 0 && active[kv.first]) {
                N_i.push_back(kv.first);
            }
        }


        unordered_set<int> Omega_i;

        for (int m : N_i) {

            for (int wedge_id : H_C[m]) {

                if (wedge_id >= 0 &&
                    wedge_id < (int)Lambda_C.size() &&
                    Lambda_C[wedge_id].active) {

                    Omega_i.insert(wedge_id);
                }
            }
        }

        total_candidates +=
            (unsigned long long)Omega_i.size();

        unordered_set<unsigned long long> omega_pairs;
        omega_pairs.reserve(Omega_i.size() * 2 + 1);

        for (int wedge_id : Omega_i) {
            const TemporalWedge& w = Lambda_C[wedge_id];
            int x = min(w.a, w.b);
            int y = max(w.a, w.b);
            unsigned long long key =
                (static_cast<unsigned long long>(
                    static_cast<unsigned int>(x)) << 32)
                | static_cast<unsigned int>(y);
            omega_pairs.insert(key);
        }


        for (int wedge_id : Omega_i) {

            const TemporalWedge& wedge =
                Lambda_C[wedge_id];

            if (!wedge.active) {
                continue;
            }

            int j = wedge.a;
            int k = wedge.b;

            if (!active[j] || !active[k]) {
                continue;
            }

            int global_j = order[j];
            int global_k = order[k];



            auto it_ij = neighbor_inter.find(j);
            auto it_ik = neighbor_inter.find(k);

            int C_ij =
                (it_ij == neighbor_inter.end())
                ? 0
                : it_ij->second;

            int C_ik =
                (it_ik == neighbor_inter.end())
                ? 0
                : it_ik->second;

            if (C_ij == 0 && C_ik == 0) {
                continue;
            }


            double t_j = hyperedge2time[global_j];
            double t_k = hyperedge2time[global_k];

            double t_min = min(t_i, min(t_j, t_k));
            double t_max = max(t_i, max(t_j, t_k));

            if (t_max - t_min > delta) {
                continue;
            }


            int C_jk = wedge.inter_ab;


            int g_ijk =
                triple_intersection_size(
                    hyperedge2node_set[global_i],
                    hyperedge2node_set[global_j],
                    hyperedge2node_set[global_k]
                );

            int motif_idx =
                get_motif_index(
                    (int)hyperedge2node[global_i].size(),
                    (int)hyperedge2node[global_j].size(),
                    (int)hyperedge2node[global_k].size(),

                    C_ij,
                    C_jk,
                    C_ik,

                    g_ijk,

                    t_i,
                    t_j,
                    t_k
                );

            if (motif_idx < 0 || motif_idx >= 96) {

                cerr
                    << "WARNING: invalid motif index "
                    << motif_idx
                    << " for temporal edges "
                    << global_j << ", "
                    << global_k << ", "
                    << global_i << "\n";

                continue;
            }

            ++motif_count[motif_idx];
            ++total_motifs;
        }

        for (size_t x = 0; x < N_i.size(); ++x) {
            int j = N_i[x];

            for (size_t y = x + 1; y < N_i.size(); ++y) {
                int k = N_i[y];

                int p = min(j, k);
                int q = max(j, k);
                unsigned long long key =
                    (static_cast<unsigned long long>(
                        static_cast<unsigned int>(p)) << 32)
                    | static_cast<unsigned int>(q);

                // If (j,k) is already a wedge, it was counted above.
                if (omega_pairs.find(key) != omega_pairs.end()) {
                    continue;
                }

                int C_ij = neighbor_inter[j];
                int C_ik = neighbor_inter[k];

                // Both j and k are in N_i, so these should be positive.
                if (C_ij <= 0 || C_ik <= 0) {
                    continue;
                }

                int global_j = order[j];
                int global_k = order[k];

                double t_j = hyperedge2time[global_j];
                double t_k = hyperedge2time[global_k];

                double t_min = min(t_i, min(t_j, t_k));
                double t_max = max(t_i, max(t_j, t_k));

                if (t_max - t_min > delta) {
                    continue;
                }

                /*
                 * Since j and k do not intersect:
                 *     C_jk = 0
                 * and therefore the triple intersection is also 0.
                 */
                int motif_idx =
                    get_motif_index(
                        (int)hyperedge2node[global_i].size(),
                        (int)hyperedge2node[global_j].size(),
                        (int)hyperedge2node[global_k].size(),

                        C_ij,
                        0,
                        C_ik,

                        0,

                        t_i,
                        t_j,
                        t_k
                    );

                if (motif_idx < 0 || motif_idx >= 96) {
                    cerr
                        << "WARNING: invalid motif index "
                        << motif_idx
                        << " for temporal edges "
                        << global_j << ", "
                        << global_k << ", "
                        << global_i << "\n";
                    continue;
                }

                ++total_candidates;
                ++motif_count[motif_idx];
                ++total_motifs;
            }
        }

        for (int m : N_i) {

            TemporalWedge new_wedge;

            new_wedge.a = m;
            new_wedge.b = i;

            // |e_m ∩ e_i|
            new_wedge.inter_ab =
                neighbor_inter[m];

            new_wedge.active = true;

            int wedge_id =
                (int)Lambda_C.size();

            Lambda_C.push_back(new_wedge);

            H_C[m].insert(wedge_id);
            H_C[i].insert(wedge_id);

            ++total_wedges_created;
        }

        for (int node : hyperedge2node[global_i]) {

            if (node >= 0 && node < V) {
                node2active[node].insert(i);
            }
        }

        active[i] = 1;
        active_queue.push_back(i);

        // Optional progress information.
        if ((i + 1) % 100000 == 0) {

            cout
                << "Processed "
                << (i + 1)
                << " / "
                << E
                << " temporal hyperedges...\n";
        }
    }
    double runtime =
        (double)(clock() - run_start)
        / CLOCKS_PER_SEC;

    cout << "\n====================================================\n";
    cout << "Exact WCA Result\n";
    cout << "====================================================\n";

    cout << "ID | Count\n";
    cout << "------------------------------\n";

    for (int i = 0; i < 96; ++i) {

        cout
            << setw(2)
            << (i + 1)
            << " | "
            << motif_count[i]
            << "\n";
    }

    cout << "# Motif instances counted: "
         << total_motifs << "\n";

    cout << "WCA Running Time: "
         << fixed
         << setprecision(6)
         << runtime
         << " sec\n";

    return 0;
}

