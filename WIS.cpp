#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <numeric>
#include <ctime>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include <set>
#include <cmath>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;


string hyperedge2hash(vector<int>& nodes) {
    stringstream ss;

    sort(nodes.begin(), nodes.end());

    for (size_t i = 0; i < nodes.size(); i++) {
        ss << nodes[i];

        if (i < nodes.size() - 1)
            ss << ",";
    }

    return ss.str();
}



struct ActiveTemporalWedge {
    int u;
    int v;
    int inter_uv;
    bool active;
};

static int temporal_intersection_size(
        const unordered_set<int>& A,
        const unordered_set<int>& B) {

    const unordered_set<int>* small = &A;
    const unordered_set<int>* large = &B;

    if (A.size() > B.size()) {
        small = &B;
        large = &A;
    }

    int cnt = 0;

    for (int x : *small) {
        if (large->count(x)) {
            ++cnt;
        }
    }

    return cnt;
}

static int temporal_triple_intersection_size(
        const unordered_set<int>& A,
        const unordered_set<int>& B,
        const unordered_set<int>& C) {

    const unordered_set<int>* smallest = &A;

    if (B.size() < smallest->size()) {
        smallest = &B;
    }

    if (C.size() < smallest->size()) {
        smallest = &C;
    }

    int cnt = 0;

    for (int x : *smallest) {
        if (A.count(x) && B.count(x) && C.count(x)) {
            ++cnt;
        }
    }

    return cnt;
}

static void count_window_algorithm1(
        const vector<int>& edges_in_window,
        const vector<vector<int>>& hyperedge2node,
        const vector<unordered_set<int>>& hyperedge2node_set,
        const vector<double>& hyperedge2time,
        int V,
        double delta,
        double window_width,
        double q_j,
        vector<double>& sample_counts) {

    vector<unordered_set<int>> node2active(V);
    unordered_map<int, unordered_set<int>> H_C;
    unordered_set<int> active_edges;
    vector<ActiveTemporalWedge> Lambda_C;

    size_t left = 0;

    auto remove_expired_edge = [&](int expired) {

        auto hit = H_C.find(expired);

        if (hit != H_C.end()) {

            vector<int> incident(
                hit->second.begin(),
                hit->second.end()
            );

            for (int wedge_id : incident) {

                if (wedge_id < 0 ||
                    wedge_id >= (int)Lambda_C.size()) {
                    continue;
                }

                ActiveTemporalWedge& w =
                    Lambda_C[wedge_id];

                if (!w.active) {
                    continue;
                }

                w.active = false;

                int other =
                    (w.u == expired)
                    ? w.v
                    : w.u;

                auto other_it =
                    H_C.find(other);

                if (other_it != H_C.end()) {
                    other_it->second.erase(wedge_id);
                }
            }

            H_C.erase(hit);
        }

        for (int node : hyperedge2node[expired]) {

            if (node >= 0 && node < V) {
                node2active[node].erase(expired);
            }
        }

        active_edges.erase(expired);
    };

    for (size_t pos = 0;
         pos < edges_in_window.size();
         ++pos) {

        int i = edges_in_window[pos];
        double t_i = hyperedge2time[i];

        while (left < pos) {

            int r = edges_in_window[left];
            double t_r = hyperedge2time[r];

            if (t_i - t_r <= delta) {
                break;
            }

            remove_expired_edge(r);
            ++left;
        }

        unordered_map<int, int> neighbor_inter;

        for (int node : hyperedge2node[i]) {

            if (node < 0 || node >= V) {
                continue;
            }

            for (int m : node2active[node]) {

                if (active_edges.count(m)) {
                    ++neighbor_inter[m];
                }
            }
        }

        vector<int> N_i;
        N_i.reserve(neighbor_inter.size());

        for (const auto& kv : neighbor_inter) {

            if (kv.second > 0 &&
                active_edges.count(kv.first)) {

                N_i.push_back(kv.first);
            }
        }

        unordered_set<int> C_i;

        for (int m : N_i) {

            auto hit = H_C.find(m);

            if (hit == H_C.end()) {
                continue;
            }

            for (int wedge_id : hit->second) {

                if (wedge_id >= 0 &&
                    wedge_id < (int)Lambda_C.size() &&
                    Lambda_C[wedge_id].active) {

                    C_i.insert(wedge_id);
                }
            }
        }

        for (int wedge_id : C_i) {

            const ActiveTemporalWedge& wedge =
                Lambda_C[wedge_id];

            if (!wedge.active) {
                continue;
            }

            int j = wedge.u;
            int k = wedge.v;

            if (!active_edges.count(j) ||
                !active_edges.count(k)) {

                continue;
            }

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

            double t_j = hyperedge2time[j];
            double t_k = hyperedge2time[k];

            double t_min =
                min(t_i, min(t_j, t_k));

            double t_max =
                max(t_i, max(t_j, t_k));

            double duration =
                t_max - t_min;

            if (duration > delta) {
                continue;
            }

            int C_jk =
                wedge.inter_uv;

            int g_ijk =
                temporal_triple_intersection_size(
                    hyperedge2node_set[i],
                    hyperedge2node_set[j],
                    hyperedge2node_set[k]
                );

            int motif_idx =
                get_motif_index(
                    (int)hyperedge2node[i].size(),
                    (int)hyperedge2node[j].size(),
                    (int)hyperedge2node[k].size(),

                    C_ij,
                    C_jk,
                    C_ik,

                    g_ijk,

                    t_i,
                    t_j,
                    t_k
                );

            double p_m =
                1.0
                - (
                    duration
                    / window_width
                );

            if (p_m > 1e-9) {

                sample_counts[motif_idx]
                    += 1.0
                       / (
                           p_m
                           * q_j
                       );
            }
        }

        for (int m : N_i) {

            ActiveTemporalWedge new_wedge;

            new_wedge.u = m;
            new_wedge.v = i;
            new_wedge.inter_uv =
                neighbor_inter[m];
            new_wedge.active = true;

            int wedge_id =
                (int)Lambda_C.size();

            Lambda_C.push_back(
                new_wedge
            );

            H_C[m].insert(
                wedge_id
            );

            H_C[i].insert(
                wedge_id
            );
        }

        for (int node : hyperedge2node[i]) {

            if (node >= 0 && node < V) {

                node2active[node]
                    .insert(i);
            }
        }

        active_edges.insert(i);
    }
}

int main(int argc, char* argv[]) {

    clock_t run_start;

    if (argc != 6) {
        cerr << "Usage: " << argv[0]
             << " <dataset> <sample_times(S)> <b> <c> <delta>\n";
        return 1;
    }

    string dataset = argv[1];
    int S = stoi(argv[2]);
    int b_shifts = stoi(argv[3]);
    double c_param = stod(argv[4]);
    double delta = stod(argv[5]);

    int SEED = static_cast<unsigned int>(time(0));
    mt19937 init_gen(SEED);

    cout << "Dataset: " << dataset << endl;
    cout << "Sample Times (S): " << S << endl;
    cout << "b: " << b_shifts << endl;
    cout << "c: " << c_param << endl;
    cout << "delta: " << delta << endl;


    string graphFile =
            dataset + ".txt";


    run_start = clock();


    vector<vector<int>> node2hyperedge,
                        hyperedge2node;

    vector<unordered_set<int>>
            hyperedge2node_set;

    vector<double>
            hyperedge2time;


    read_data(
            graphFile,
            node2hyperedge,
            hyperedge2node,
            hyperedge2node_set,
            hyperedge2time
    );


    int V =
            (int)node2hyperedge.size();

    int E =
            (int)hyperedge2node.size();


    cout << "Dataset: "
         << dataset
         << endl;

    cout << "# Nodes: "
         << V
         << ", # Edges: "
         << E
         << endl;

    cout << "Data Load Time: "
         << (double)(clock() - run_start)
            / CLOCKS_PER_SEC
         << " sec\n";


    run_start = clock();


    unordered_map<string, int>
            hash2index;

    vector<int>
            hyperedge2index;

    vector<vector<int>>
            E2V;

    vector<unordered_set<int>>
            E2V_set;


    for (int i = 0; i < E; i++) {

        string hash =
                hyperedge2hash(
                        hyperedge2node[i]
                );


        if (hash2index.find(hash)
            == hash2index.end()) {

            hash2index[hash] =
                    (int)hash2index.size();

            E2V.push_back(
                    hyperedge2node[i]
            );

            E2V_set.push_back(
                    hyperedge2node_set[i]
            );
        }


        hyperedge2index.push_back(
                hash2index[hash]
        );
    }


    int E_static =
            (int)E2V.size();


    cout << "# Induced Static Hyperedges: "
         << E_static
         << endl;


    int e_max = 0;

    long long total_nodes_sum = 0;


    for (const auto& nodes :
            hyperedge2node) {

        int current_size =
                (int)nodes.size();


        if (current_size > e_max) {
            e_max = current_size;
        }


        total_nodes_sum +=
                current_size;
    }


    double avg_nodes =
            (E > 0)
            ? (double)total_nodes_sum / E
            : 0;


    cout << "Max nodes in a temporal hyperedge (Emax): "
         << e_max
         << endl;

    cout << "Average nodes per temporal hyperedge: "
         << fixed
         << setprecision(2)
         << avg_nodes
         << endl;



    vector<int> order(E);

    iota(
            order.begin(),
            order.end(),
            0
    );


    sort(
            order.begin(),
            order.end(),
            [&](int i, int j) {

                if (hyperedge2time[i]
                    == hyperedge2time[j])
                    return i < j;

                return hyperedge2time[i]
                       < hyperedge2time[j];
            }
    );


    double window_width =
            c_param * delta;


    cout << "Window Width: "
         << window_width
         << " (c="
         << c_param
         << " * delta)"
         << endl;


    cout << "Preprocessing Done: "
         << (double)(clock() - run_start)
            / CLOCKS_PER_SEC
         << " sec\n";


    cout << "-----------------------------------------------------------\n";


    run_start = clock();


    vector<double>
            global_weighted_counts(
                    96,
                    0.0
            );


    mt19937 gen(SEED);


    vector<vector<int>>
            local_adj(E_static);

    vector<unordered_map<int, int>>
            local_inter(E_static);

    vector<vector<int>>
            static2temporal(E_static);

    vector<vector<int>>
            node2static_in_window(V);


    for (int s = 0;
         s < S;
         s++) {


        if (s % 10 == 0
            && s > 0) {

            cout << "Sample "
                 << s
                 << " / "
                 << S
                 << "..."
                 << endl;
        }


        vector<double>
                sample_counts(
                        96,
                        0.0
                );


        for (int shift_k = 0;
             shift_k < b_shifts;
             shift_k++) {


            uniform_real_distribution<>
                    shift_dist(
                            -window_width,
                            0
                    );


            double s_shift =
                    shift_dist(gen);


            uniform_int_distribution<>
                    edge_dist(
                            0,
                            E - 1
                    );


            double t_pivot =
                    hyperedge2time[
                            order[
                                    edge_dist(gen)
                            ]
                    ];


            double T_start =
                    floor(
                            (t_pivot - s_shift)
                            / window_width
                    )
                    * window_width
                    + s_shift;


            double T_end =
                    T_start
                    + window_width;


            auto it_start =
                    lower_bound(
                            order.begin(),
                            order.end(),
                            T_start,

                            [&](int idx,
                                double val) {

                                return hyperedge2time[idx]
                                       < val;
                            }
                    );


            int start_pos =
                    distance(
                            order.begin(),
                            it_start
                    );


            vector<int>
                    edges_in_window;


            for (int i = start_pos;
                 i < E;
                 i++) {

                int e_idx =
                        order[i];

                double t =
                        hyperedge2time[e_idx];


                if (t >= T_end)
                    break;


                edges_in_window.push_back(
                        e_idx
                );
            }


            int n_j =
                    edges_in_window.size();


            if (n_j < 3)
                continue;


            double q_j =
                    (double)n_j
                    / (double)E;


            count_window_algorithm1(
                    edges_in_window,
                    hyperedge2node,
                    hyperedge2node_set,
                    hyperedge2time,
                    V,
                    delta,
                    window_width,
                    q_j,
                    sample_counts
            );


        } // End Shifts


        for (int i = 0;
             i < 96;
             i++) {

            global_weighted_counts[i]
                    += sample_counts[i]
                       / (double)b_shifts;
        }

    } // End Samples



    vector<double>
            final_counts(96);


    for (int i = 0;
         i < 96;
         i++) {

        final_counts[i] =
                global_weighted_counts[i]
                / (double)S;
    }


    cout << "Counting Time: "
         << fixed
         << setprecision(2)
         << 1000.0
            * (double)(
                    clock()
                    - run_start
            )
            / CLOCKS_PER_SEC
         << " ms"
         << endl;


    return 0;
}