#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

int main() {
    // --- 解决控制台乱码：设置输出编码为 UTF-8 ---
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    // 1. 核心配置
    string dataset_name;
    cout << "Dataset_name: ";
    cin >> dataset_name;
    string folder = dataset_name + "/";
    string nverts_file    = folder + dataset_name + "-nverts.txt";
    string simplices_file = folder + dataset_name + "-simplices.txt";
    string times_file     = folder + dataset_name + "-times.txt";
    string output_file    = folder + dataset_name + ".txt";

    // 2. 打开文件流
    ifstream fnverts(nverts_file);
    ifstream fsimplices(simplices_file);
    ifstream ftimes(times_file);
    ofstream foutput(output_file);

    // 检查文件是否成功打开
    if (!fnverts.is_open() || !fsimplices.is_open() || !ftimes.is_open()) {
        cerr << "Error!" << endl;
        return 1;
    }

    int num_vertices;     // 对应每个单纯形的顶点数
    long long timestamp;  // 对应时间戳
    int node_id;          // 单个节点ID

    cout << "start deal with " << folder << " dataset" << endl;

    // 3. 逐个单纯形进行处理
    // 逻辑：nverts.txt 和 times.txt 的行数是一致的
    while (fnverts >> num_vertices && ftimes >> timestamp) {
        vector<int> current_simplex;

        // 根据 nverts 指定的数量，从 simplices.txt 中读取连续的节点
        for (int i = 0; i < num_vertices; ++i) {
            if (fsimplices >> node_id) {
                current_simplex.push_back(node_id);
            }
        }

        // 4. 按照目标格式输出：节点1,节点2,节点3\t时间戳
        for (size_t i = 0; i < current_simplex.size(); ++i) {
            foutput << current_simplex[i];
            if (i < current_simplex.size() - 1) {
                foutput << ","; // 节点间用逗号分隔
            }
        }
        foutput << "\t" << timestamp << "\n"; // 接制表符和时间戳
    }

    // 5. 关闭文件
    fnverts.close();
    fsimplices.close();
    ftimes.close();
    foutput.close();

    cout << "Over: " << output_file << endl;

    return 0;
}