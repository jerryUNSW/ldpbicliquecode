#ifndef BICLIQUE_DDS_ELDP_V17_H
#define BICLIQUE_DDS_ELDP_V17_H

#include "bigraph.h"
#include <vector>
#include <string>

struct DDSResultV17 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
    
    // 新增：详细的分析数据
    std::string trace_file;           // greedy trace 文件路径
    std::string engagement_file;      // engagement 数据文件路径
    std::string variance_file;        // variance 分析文件路径
};

// v17 = v14 的核心 + 详细的数据收集
// 目标：验证 OPTIMIZATION_ROADMAP 中的假设
DDSResultV17 eldp_v17(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed,
                      const std::string& dataset_name, int round_id);

#endif
