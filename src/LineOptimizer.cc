/**
* This file is part of LLD-SLAM.
*
* Copyright (C) 2018 Alexander Vakhitov <alexander.vakhitov at gmail dot com> (Skoltech)
* For more information see <https://github.com/alexandervakhitov/lld-slam>
*
* LLD-SLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* LLD-SLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with LLD-SLAM. If not, see <http://www.gnu.org/licenses/>.
*/

#include <utility>
#include <Eigen/Dense>
#include <Thirdparty/g2o/g2o/types/types_sba.h>
#include <Thirdparty/g2o/g2o/types/types_six_dof_expmap.h>
#include <Thirdparty/g2o/g2o/core/robust_kernel_impl.h>
#include <include/LineMatching.h>
#include "LineOptimizer.h"

LineOptimizer::LineOptimizer(g2o::SparseOptimizer &optimizer, int maxPtId, double thHuberLinesStereo,
                             double thHuberLinesMono, double gamma):
        optimizer(optimizer), maxPtId(maxPtId), thHuberLinesStereo(thHuberLinesStereo), thHuberLinesMono(thHuberLinesMono), gamma(gamma)
{
    infoLines = 1.0;
    this->thHuberLinesStereo *= gamma;
    this->thHuberLinesMono *= gamma;
    infoLines *= gamma*gamma;
}

void LineOptimizer::AddLineMinimal(int line_id, const Eigen::Matrix3d& K_eig, double stereo_b, const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const std::map<int, std::pair<KeyLine, KeyLine>>& proj_map)
{
    g2o::VertexSBALine *v_line = new g2o::VertexSBALine();
    int id = line_id + maxPtId + 1;
    v_line->setId(id);
    Eigen::Matrix3d R_line;
    R_line.col(0) = line_dir;
    R_line.col(1) = X0 / X0.norm();
    R_line.col(2) = line_dir.cross(X0) / X0.norm();
    Eigen::Quaterniond q_line(R_line);
    double alpha = X0.norm();
    v_line->setEstimate(g2o::LineParams(q_line, alpha));
    v_line->setFixed(false);
    v_line->setMarginalized(true);
    optimizer.addVertex(v_line);

    for (std::map<int, std::pair<KeyLine, KeyLine>>::const_iterator mit = proj_map.begin(), mend = proj_map.end(); mit != mend; mit++)
    {
        int kf_id = mit->first;

        for (int si = 0; si < 2; si++)
        {
            if (si == 1 &&  mit->second.second.startPointX < 0)
            {
                continue;
            }
            g2o::EdgeSE3ProjectLine* e = new g2o::EdgeSE3ProjectLine();
            e->f = K_eig(0,0);
            e->cx = K_eig(0, 2);
            e->cy = K_eig(1, 2);
//            std::cout << e->f << " " << e->cx << " " << e->cy << std::endl;
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(v_line));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(kf_id)));
            double b = 0;
            if (si == 1)
            {
                b = -stereo_b;
            }
            e->b = Eigen::Vector3d(b,0,0);
            Eigen::Vector2d obs;
            obs.setZero();
            e->setMeasurement(obs);

            g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
            double thHuberLines = thHuberLinesStereo;
            if (mit->second.second.startPointX < 0)
            {
                thHuberLines = thHuberLinesMono;
            }
            rk->setDelta(thHuberLines);
            e->setRobustKernel(rk);
            KeyLine kl;
            if (si == 0) {
                kl = mit->second.first;
            } else {
                kl = mit->second.second;
            }

            double info = infoLines;
            double thrReproj = GetReprojThrPyramid(1.0, kl.octave);

            info /= thrReproj*thrReproj;
            e->setInformation(Eigen::Matrix2d::Identity()*info);

            Eigen::Vector3d xs, xe;
            xs << kl.startPointX, kl.startPointY, 1.0;
            xe << kl.endPointX, kl.endPointY, 1.0;
            Eigen::Matrix3d Ki = K_eig.inverse();

            Ki.setIdentity();

//            for (int g = 0; g < 3; g++)
//                std::cout << Ki.row(g) << std::endl;
            e->x1 =  Ki * xs;
            e->x2 = Ki * xe;
            e->computeError();

            optimizer.addEdge(e);

            if (line2edge.count(line_id) == 0)
            {
                line2edge.insert(std::make_pair(line_id, std::vector<g2o::EdgeSE3ProjectLine*>()));
            }
            line2edge[line_id].push_back(e);
            line_edge2vert[e] = v_line;
            line_edge2type[e] = (mit->second.second.startPointX >= 0);
        }
    }
}

void LineOptimizer::DisableOutliers()
{

    std::map<g2o::VertexSBALine*, int> line_proj_cnt;

    int line_edge_cnt = 0;
    for (auto& ev: line_edge2vert)
    {
        auto vp = ev.second;
        auto ep = ev.first;
        if (line_proj_cnt.find(vp) == line_proj_cnt.end())
        {
            line_proj_cnt[vp] = 0;
        }
        double thr = thHuberLinesStereo*thHuberLinesStereo;
        if (!line_edge2type[ep])
        {
            thr = thHuberLinesMono*thHuberLinesMono;
        }
        if (ep->chi2() > thr || !ep->IsDepthPositive())
        {
            ep->setLevel(1);
        } else {
            line_proj_cnt[vp] += 2;
            line_edge_cnt++;
        }
        ep->setRobustKernel(0);
    }

    for (auto& vc: line_proj_cnt)
    {
        if (vc.second <= ln_filter)
        {
            g2o::HyperGraph::EdgeSet e_set(vc.first->edges());
            for (auto& e: e_set)
            {
                optimizer.removeEdge(e);
            }
            optimizer.removeVertex(vc.first);
        }
    }
}

bool LineOptimizer::GetLineData(int line_id, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir, std::vector<int>* outlier_projs)
{
    int id = line_id + maxPtId + 1;
    g2o::VertexSBALine *v_line = static_cast<g2o::VertexSBALine *>(optimizer.vertex(id));
    if (!v_line)
    {
        return false;
    }
    Eigen::Matrix3d R_l(v_line->estimate().GetR());
    *line_dir = R_l.col(0);
    *X0 = v_line->estimate().GetAlpha() * R_l.col(1);
    int cnt = 0;
    for (auto& e: line2edge[line_id])
    {
        bool depth_pos = e->IsDepthPositive();
        e->computeError();
        double thr = thHuberLinesStereo*thHuberLinesStereo;
        if (!line_edge2type[e])
        {
            thr = thHuberLinesMono*thHuberLinesMono;
        }
        if (e->chi2() > thr || !depth_pos)
        {
            int kf_id = e->vertex(1)->id();
            outlier_projs->push_back(kf_id);
        }
        cnt++;
    }
    return true;
}

