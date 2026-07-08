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
#include <string>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;

string hyperedge2hash(vector<int>& nodes) {
    stringstream ss;
    sort(nodes.begin(), nodes.end());
    for (size_t i = 0; i < nodes.size(); i++) {
        ss << nodes[i];
        if (i < nodes.size() - 1) ss << ",";
    }
    return ss.str();
}

int main(int argc, char* argv[]) {
    clock_t run_start;

    // ================== 参数配置 ==================
    string dataset = argv[1];
    int S = atoi(argv[2]);
    double delta = atof(argv[3]);

    int SEED = static_cast<unsigned int>(time(0));
    mt19937 init_gen(SEED);
    uniform_int_distribution<> param_dist(1, 5);
    // 窗口参数
    double c_param = static_cast<double>(param_dist(init_gen));
    // 移位次数
    int b_shifts = param_dist(init_gen);

    cout << "Dataset: " << dataset << endl;
    cout << "S: " << S << endl;
    cout << "Delta: " << delta << endl;
    cout << "Randomly selected c_param (c): " << c_param << endl;
    cout << "Randomly selected b_shifts (b): " << b_shifts << endl;

    string graphFile = dataset + ".txt";

    // ================== 1. 读取数据 ==================
    run_start = clock();
    vector<vector<int>> node2hyperedge, hyperedge2node;
    vector<unordered_set<int>> hyperedge2node_set;
    vector<double> hyperedge2time;
    read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set, hyperedge2time);

    int V = (int)node2hyperedge.size();
    int E = (int)hyperedge2node.size();

    cout << "Dataset: " << dataset << endl;
    cout << "# Nodes: " << V << ", # Edges: " << E << endl;
    cout << "Data Load Time: " << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec\n";

    // ================== 2. 静态图预处理 ==================
    run_start = clock();
    unordered_map<string, int> hash2index;
    vector<int> hyperedge2index;
    vector<vector<int>> E2V;
    vector<unordered_set<int>> E2V_set;

    for (int i = 0; i < E; i++) {
        string hash = hyperedge2hash(hyperedge2node[i]);
        if (hash2index.find(hash) == hash2index.end()) {
            hash2index[hash] = (int)hash2index.size();
            E2V.push_back(hyperedge2node[i]);
            E2V_set.push_back(hyperedge2node_set[i]);
        }
        hyperedge2index.push_back(hash2index[hash]);
    }
    int E_static = (int)E2V.size();
    cout << "# Induced Static Hyperedges: " << E_static << endl;

    int e_max = 0;
    long long total_nodes_sum = 0;
    for (const auto& nodes : hyperedge2node) {
        int current_size = (int)nodes.size();
        if (current_size > e_max) {
            e_max = current_size;
        }
        total_nodes_sum += current_size;
    }
    double avg_nodes = (E > 0) ? (double)total_nodes_sum / E : 0;

    cout << "Max nodes in a temporal hyperedge (Emax): " << e_max << endl;
    cout << "Average nodes per temporal hyperedge: " << fixed << setprecision(2) << avg_nodes << endl;

    // ================== 3. 时间排序 ==================
    vector<int> order(E);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int i, int j) {
        if (hyperedge2time[i] == hyperedge2time[j]) return i < j;
        return hyperedge2time[i] < hyperedge2time[j];
    });

    double window_width = c_param * delta;
    cout << "Window Width: " << window_width << " (c=" << c_param << " * delta)" << endl;
    cout << "Preprocessing Done: " << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec\n";
    cout << "-----------------------------------------------------------\n";

    // ================== 4. 核心计数逻辑 ==================
    run_start = clock();
    vector<double> global_weighted_counts(96, 0.0);
    mt19937 gen(SEED);

    vector<vector<int>> local_adj(E_static);
    vector<unordered_map<int, int>> local_inter(E_static);
    vector<vector<int>> static2temporal(E_static);
    vector<vector<int>> node2static_in_window(V);

    for (int s = 0; s < S; s++) {
        if (s % 10 == 0 && s > 0) cout << "Sample " << s << " / " << S << "..." << endl;

        vector<double> sample_counts(96, 0.0);

        for (int shift_k = 0; shift_k < b_shifts; shift_k++) {
            uniform_real_distribution<> shift_dist(-window_width, 0);
            double s_shift = shift_dist(gen);

            uniform_int_distribution<> edge_dist(0, E - 1);
            double t_pivot = hyperedge2time[order[edge_dist(gen)]];

            double T_start = floor((t_pivot - s_shift) / window_width) * window_width + s_shift;
            double T_end = T_start + window_width;

            auto it_start = lower_bound(order.begin(), order.end(), T_start, [&](int idx, double val){
                return hyperedge2time[idx] < val;
            });
            int start_pos = distance(order.begin(), it_start);

            vector<int> edges_in_window;
            for (int i = start_pos; i < E; i++) {
                int e_idx = order[i];
                double t = hyperedge2time[e_idx];
                if (t >= T_end) break;
                edges_in_window.push_back(e_idx);
            }

            int n_j = edges_in_window.size();
            if (n_j < 3) continue;


            double q_j = (double)n_j / (double)E;

            set<int> active_static_indices;
            for(int e : edges_in_window) {
                int s_idx = hyperedge2index[e];
                active_static_indices.insert(s_idx);
                static2temporal[s_idx].clear();
                local_adj[s_idx].clear();
                local_inter[s_idx].clear();
            }

            for(int e : edges_in_window) {
                int s_idx = hyperedge2index[e];
                static2temporal[s_idx].push_back(e);
            }
            auto C2 = [](long double n) -> long double {
                return n >= 2 ? n * (n - 1) / 2.0L : 0.0L;
            };

            auto C3 = [](long double n) -> long double {
                return n >= 3 ? n * (n - 1) * (n - 2) / 6.0L : 0.0L;
            };

            long double est96 = 0.0L;
            long double estDup = 0.0L;
            long double estABC = 0.0L;

            for (int u : active_static_indices) {
                long double ku = static2temporal[u].size();
                est96 += C3(ku);
            }

            for (int u : active_static_indices) {
                for (int v : local_adj[u]) {
                    if (u >= v) continue;

                    long double ku = static2temporal[u].size();
                    long double kv = static2temporal[v].size();

                    estDup += C2(ku) * kv + C2(kv) * ku;

                    for (int c : local_adj[u]) {
                        if (v >= c) continue;
                        long double kc = static2temporal[c].size();
                        estABC += ku * kv * kc;
                    }

                    for (int c : local_adj[v]) {
                        if (u >= c) continue;
                        if (local_inter[u].count(c)) continue;
                        long double kc = static2temporal[c].size();
                        estABC += ku * kv * kc;
                    }
                }
            }

            for(int s_idx : active_static_indices) {
                for(int node : E2V[s_idx]) node2static_in_window[node].clear();
            }
            for(int s_idx : active_static_indices) {
                for(int node : E2V[s_idx]) node2static_in_window[node].push_back(s_idx);
            }

            for(int u : active_static_indices) {
                unordered_map<int, int> neighbors_count;
                for(int node : E2V[u]) {
                    for(int v : node2static_in_window[node]) {
                        if(u != v) neighbors_count[v]++;
                    }
                }
                for(auto& p : neighbors_count) {
                    local_adj[u].push_back(p.first);
                    local_inter[u][p.first] = p.second;
                }
            }

            for (int u : active_static_indices) {
                const auto& times_u = static2temporal[u];
                if (times_u.size() < 3) continue;

                for (int i = 0; i < times_u.size(); i++) {
                    for (int j = i + 1; j < times_u.size(); j++) {
                        for (int k = j + 1; k < times_u.size(); k++) {
                            double t1 = hyperedge2time[times_u[i]];
                            double t2 = hyperedge2time[times_u[j]];
                            double t3 = hyperedge2time[times_u[k]];

                            double t_min = min({t1, t2, t3});
                            double t_max = max({t1, t2, t3});
                            double duration = t_max - t_min;

                            if (duration > delta) continue;

                            int sz = E2V[u].size();
                            int m_idx = get_motif_index(sz, sz, sz, sz, sz, sz, sz, t1, t2, t3);

                            double p_m = 1.0 - (duration / window_width);
                            if (p_m > 1e-9) {
                                sample_counts[m_idx] += 1.0 / (p_m * q_j);
                            }
                        }
                    }
                }
            }
            // ==============================================================================

            for(int u : active_static_indices) {
                for(int v : local_adj[u]) {
                    if(u >= v) continue;

                    int C_uv = local_inter[u][v];
                    const auto& times_u = static2temporal[u];
                    const auto& times_v = static2temporal[v];

                    // ================= 新增逻辑：处理 87-95 号模体 (双重复边) =================
                    // 情况 A: (u, u, v)
                    if (times_u.size() >= 2 && times_v.size() >= 1) {
                        for (int i = 0; i < times_u.size(); i++) {
                            for (int j = i + 1; j < times_u.size(); j++) {
                                for (int t_v_id : times_v) {
                                    double tu1 = hyperedge2time[times_u[i]];
                                    double tu2 = hyperedge2time[times_u[j]];
                                    double tv = hyperedge2time[t_v_id];

                                    double t_min = min({tu1, tu2, tv}), t_max = max({tu1, tu2, tv});
                                    double duration = t_max - t_min;
                                    if (duration > delta) continue;

                                    int m_idx = get_motif_index(E2V[u].size(), E2V[u].size(), E2V[v].size(),
                                                                E2V[u].size(), C_uv, C_uv, C_uv, tu1, tu2, tv);
                                    double p_m = 1.0 - (duration / window_width);
                                    if (p_m > 1e-9) sample_counts[m_idx] += 1.0 / (p_m * q_j);
                                }
                            }
                        }
                    }

                    // 情况 B: (v, v, u)
                    if (times_v.size() >= 2 && times_u.size() >= 1) {
                        for (int i = 0; i < times_v.size(); i++) {
                            for (int j = i + 1; j < times_v.size(); j++) {
                                for (int t_u_id : times_u) {
                                    double tv1 = hyperedge2time[times_v[i]];
                                    double tv2 = hyperedge2time[times_v[j]];
                                    double tu = hyperedge2time[t_u_id];

                                    double t_min = min({tv1, tv2, tu}), t_max = max({tv1, tv2, tu});
                                    double duration = t_max - t_min;
                                    if (duration > delta) continue;

                                    int m_idx = get_motif_index(E2V[v].size(), E2V[v].size(), E2V[u].size(),
                                                                E2V[v].size(), C_uv, C_uv, C_uv, tv1, tv2, tu);
                                    double p_m = 1.0 - (duration / window_width);
                                    if (p_m > 1e-9) sample_counts[m_idx] += 1.0 / (p_m * q_j);
                                }
                            }
                        }
                    }
                    // ==============================================================================

                    for(int c : local_adj[u]) {
                        if(v >= c) continue;

                        int C_uc = local_inter[u][c];
                        int C_vc = 0;
                        if(local_inter[v].count(c)) C_vc = local_inter[v][c];

                        int g_abc = 0;
                        if(C_vc > 0) {
                            for(int node : E2V[u]) {
                                if(E2V_set[v].count(node) && E2V_set[c].count(node)) g_abc++;
                            }
                        }

                        for(int t_u_id : static2temporal[u]) {
                            double tu = hyperedge2time[t_u_id];
                            for(int t_v_id : static2temporal[v]) {
                                double tv = hyperedge2time[t_v_id];
                                for(int t_c_id : static2temporal[c]) {
                                    double tc = hyperedge2time[t_c_id];

                                    double t_min = min({tu, tv, tc});
                                    double t_max = max({tu, tv, tc});
                                    double duration = t_max - t_min;

                                    if(duration > delta) continue;

                                    int m_idx = get_motif_index(
                                            E2V[u].size(), E2V[v].size(), E2V[c].size(),
                                            C_uv, C_vc, C_uc, g_abc,
                                            tu, tv, tc
                                    );

                                    double p_m = 1.0 - (duration / window_width);
                                    if(p_m > 1e-9) {
                                        double weight = 1.0 / (p_m * q_j);
                                        sample_counts[m_idx] += weight;
                                    }
                                }
                            }
                        }
                    } // End Case A

                    for(int c : local_adj[v]) {
                        if(u >= c) continue;
                        if(local_inter[u].count(c)) continue;

                        int C_vc = local_inter[v][c];
                        int C_uc = 0;
                        int g_abc = 0;

                        for(int t_u_id : static2temporal[u]) {
                            double tu = hyperedge2time[t_u_id];
                            for(int t_v_id : static2temporal[v]) {
                                double tv = hyperedge2time[t_v_id];
                                for(int t_c_id : static2temporal[c]) {
                                    double tc = hyperedge2time[t_c_id];

                                    double t_min = min({tu, tv, tc});
                                    double t_max = max({tu, tv, tc});
                                    double duration = t_max - t_min;
                                    if(duration > delta) continue;

                                    int m_idx = get_motif_index(
                                            E2V[u].size(), E2V[v].size(), E2V[c].size(),
                                            C_uv, C_vc, C_uc, g_abc,
                                            tu, tv, tc
                                    );

                                    double p_m = 1.0 - (duration / window_width);
                                    if(p_m > 1e-9) {
                                        sample_counts[m_idx] += 1.0 / (p_m * q_j);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for(int i=0; i<96; i++) {
            global_weighted_counts[i] += sample_counts[i] / (double)b_shifts;
        }

    }

    // ================== 5. 输出结果 ==================

    vector<double> final_counts(96);
    for(int i=0; i<96; i++) final_counts[i] = global_weighted_counts[i] / (double)S;
    vector<double> truth;
    for (int i = 0; i < 96; i++) {
        double val = final_counts[i];
        cout << i + 1 << " "
             << fixed << setprecision(1) << val
             << endl;
    }
    cout << "Counting Done: " << 1000 * (double)(clock() - run_start) / CLOCKS_PER_SEC << " \n";
    return 0;
}