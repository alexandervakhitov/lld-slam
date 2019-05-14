//
// Created by alexander on 04.07.18.
//

#ifndef ORB_SLAM2_LINEOPTIMIZER_H
#define ORB_SLAM2_LINEOPTIMIZER_H

#include <Thirdparty/g2o/g2o/core/sparse_optimizer.h>
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"

class LineOptimizer {
public:
    LineOptimizer(g2o::SparseOptimizer& optimizer, int maxPtId, double thHuberLinesStereo, double thHuberLinesMono, double gamma);

    void AddLineMinimal(int line_id, const Eigen::Matrix3d& K_eig, double b, const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const std::map<int, std::pair<KeyLine, KeyLine>>& proj_map);

    void DisableOutliers();


    bool GetLineData(int line_id, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir, std::vector<int>* outlier_projs);

private:
    int ln_filter = 4;

    g2o::SparseOptimizer& optimizer;
    int maxPtId;
    double thHuberLinesStereo, thHuberLinesMono;
    double gamma;
    double infoLines;
    std::map<int, std::vector<g2o::EdgeSE3ProjectLine*>> line2edge;
    std::map<g2o::EdgeSE3ProjectLine*, g2o::VertexSBALine *> line_edge2vert;
    std::map<g2o::EdgeSE3ProjectLine*, bool> line_edge2type;
};


#endif //ORB_SLAM2_LINEOPTIMIZER_H
