#ifndef BICLIQUE_DDS_ELDP_V16_H
#define BICLIQUE_DDS_ELDP_V16_H

#include "bigraph.h"
#include <vector>

struct DDSResultV16 {
    std::vector<int> vertex_set;
    long double estimated_density;
    long double real_density;
    
    // 新增：优化统计信息
    int stopped_by_variance;      // 是否因方差停止
    int stopped_by_confidence;    // 是否因置信区间停止
    int stopped_by_precision;     // 是否因精度停止
    int stopped_by_size;          // 是否因大小停止
    int final_size;               // 最终集合大小
};

// v16 优化算法
DDSResultV16 eldp_v16(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed,
                      // 优化参数
                      int size_threshold = 100,
                      long double variance_threshold = 0.5,
                      long double precision_threshold = 0.5,
                      int max_size = 500,
                      int top_k = 3);

#endif
