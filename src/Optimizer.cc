/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*
* Modified by Alexander Vakhitov (2018): added MapLine-related containers and methods
*/

#include "Optimizer.h"

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"

#include "include/vgl.h"
#include "include/LineMatching.h"

#include<Eigen/StdVector>

#include "Converter.h"

#include<mutex>
#include <opencv/cxeigen.hpp>
#include <include/LineOptimizer.h>

using namespace LineMatching;

namespace ORB_SLAM2
{

void AddLineMinimal(const int maxPtId, const double infoLines, const double thHuberLines, g2o::SparseOptimizer& optimizer,
                    MapLine* pML, KeyFrame* pKF, std::vector<g2o::EdgeSE3ProjectLine*>& vpEdgesLines, std::vector<MapLine*>& vpMapLine,
                    std::vector<KeyFrame*>& vpEdgeKFLines, std::map<MapLine*, std::vector<g2o::EdgeSE3ProjectLine*>>& mapline_edges)
{
    g2o::VertexSBALine *v_line = new g2o::VertexSBALine();
    int id = pML->mnId + maxPtId + 1;
    v_line->setId(id);
    Eigen::Vector3d X0, line_dir;
    pML->GetMinimalPos(&X0, &line_dir);
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

    const map<KeyFrame *, size_t> observations = pML->GetObservations();

    Eigen::Matrix3d K_eig;
    cv::cv2eigen(pKF->mK, K_eig);

    for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
         mit != mend; mit++)
    {
        KeyFrame* pKFi = mit->first;
        if (!pKFi)
        {
            std::cout << " error pointer to KF" << std::endl;
        }

        if(pKFi->isBad()) {
            continue;
        }

        cv::Mat T_cv = pKFi->GetPose();
        Eigen::Matrix<double, 4, 4> T;
        cv::cv2eigen(T_cv, T);

        int line_id = mit->second;

        for (int si = 0; si < 2; si++) {
            g2o::EdgeSE3ProjectLine* e = new g2o::EdgeSE3ProjectLine();
            e->f = pKFi->fx;
            e->cx = pKFi->cx;
            e->cy = pKFi->cy;
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(v_line));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
            double b = 0;
            if (si == 1)
            {
                b = -pKFi->mbf / pKFi->mK.at<float>(0,0);
            }
            e->b = Eigen::Vector3d(b,0,0);
            Eigen::Vector2d obs;
            obs.setZero();
            e->setMeasurement(obs);


            g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
            rk->setDelta(thHuberLines );
            e->setRobustKernel(rk);
            KeyLine kl;
            if (si == 0) {
                kl = pKFi->mvLinesLeft[line_id];
            } else {
                if (pKFi->line_matches[line_id] < 0)
                {
                    std::cout <<"BIG ERROR no stereo line" << std::endl;
                }
                kl = pKFi->mvLinesRight[pKFi->line_matches[line_id]];
            }

            double info = infoLines;
            double thrReproj = GetReprojThrPyramid(1.0, kl.octave);
            info /= thrReproj*thrReproj;
            e->setInformation(Eigen::Matrix2d::Identity()*info);

            Eigen::Vector3d xs, xe;
            xs << kl.startPointX, kl.startPointY, 1.0;
            xe << kl.endPointX, kl.endPointY, 1.0;
            e->x1 = K_eig.inverse() * xs;
            e->x2 = K_eig.inverse() * xe;

            optimizer.addEdge(e);
            vpEdgesLines.push_back(e);
            vpMapLine.push_back(pML);
            vpEdgeKFLines.push_back(pKFi);
            if (mapline_edges.count(pML) == 0)
            {
                mapline_edges.insert(std::make_pair(pML, std::vector<g2o::EdgeSE3ProjectLine*>()));
            }
            mapline_edges[pML].push_back(e);
        }
    }
}

int AddLineMinimalGlobal(const int maxPtId, const double infoLines, const double thHuberLines, g2o::SparseOptimizer& optimizer,
                          MapLine* pML)
{
    g2o::VertexSBALine *v_line = new g2o::VertexSBALine();
    int id = pML->mnId + maxPtId + 1;
    v_line->setId(id);
    Eigen::Vector3d X0, line_dir;
    pML->GetMinimalPos(&X0, &line_dir);
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

    const map<KeyFrame *, size_t> observations = pML->GetObservations();

    int edge_cnt = 0;
//        pML->total_obs = 0;
//        pML->active_obs = 0;
//        pML->X1init = X1;
//        pML->X2init = X2;

    for (map<KeyFrame *, size_t>::const_iterator mit = observations.begin(), mend = observations.end();
         mit != mend; mit++)
    {
        KeyFrame* pKFi = mit->first;

        if (!pKFi)
        {
            std::cout << " error pointer to KF" << std::endl;
        }

        if(pKFi->isBad()) {
            continue;
        }

        Eigen::Matrix3d K_eig;
        cv::cv2eigen(pKFi->mK, K_eig);


        cv::Mat T_cv = pKFi->GetPose();
        Eigen::Matrix<double, 4, 4> T;
        cv::cv2eigen(T_cv, T);

        int line_id = mit->second;

        for (int si = 0; si < 2; si++) {
            g2o::EdgeSE3ProjectLine* e = new g2o::EdgeSE3ProjectLine();
            e->f = pKFi->fx;
            e->cx = pKFi->cx;
            e->cy = pKFi->cy;
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(v_line));
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(pKFi->mnId)));
            double b = 0;
            if (si == 1)
            {
                b = -pKFi->mbf / pKFi->mK.at<float>(0,0);
            }
            e->b = Eigen::Vector3d(b,0,0);
            Eigen::Vector2d obs;
            obs.setZero();
            e->setMeasurement(obs);
            e->setInformation(Eigen::Matrix2d::Identity());
//            e->setInformation(Eigen::Matrix3d::Identity());//temporary
            g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
            rk->setDelta(thHuberLines );
            e->setRobustKernel(rk);
            KeyLine kl;
            if (si == 0) {
                kl = pKFi->mvLinesLeft[line_id];
            } else {
                if (pKFi->line_matches[line_id] < 0)
                {
                    std::cout <<"BIG ERROR no stereo line" << std::endl;
                }
                kl = pKFi->mvLinesRight[pKFi->line_matches[line_id]];
            }
            Eigen::Vector3d xs, xe;
            xs << kl.startPointX, kl.startPointY, 1.0;
            xe << kl.endPointX, kl.endPointY, 1.0;
            e->x1 = K_eig.inverse() * xs;
            e->x2 = K_eig.inverse() * xe;

            optimizer.addEdge(e);
            edge_cnt++;
        }
    }

    if (edge_cnt == 0)
    {
        optimizer.removeVertex(v_line);
    }
    return edge_cnt;
}

void AddMainEdgeVertex(KeyFrame* refKF, MapLine* pML, int id, bool is_start, Eigen::Vector3d& X1, g2o::SparseOptimizer& optimizer, int maxPtId, double thHuberStereo,
                       std::vector<g2o::EdgeStereoSE3ProjectXYZ*>& lineMainEdges,
                       std::vector<MapLine*>& mapLines, g2o::VertexSBAPointXYZ **v_line_1_p)
{

    g2o::VertexSBAPointXYZ *v_line_1 = new g2o::VertexSBAPointXYZ();

    v_line_1->setId(id);
    v_line_1->setEstimate(X1);
    v_line_1->setFixed(false);
    v_line_1->setMarginalized(true);
    optimizer.addVertex(v_line_1);

    Eigen::Matrix<double,3,1> obs;
    KeyLine kl_left = refKF->mvLinesLeft[pML->mObservations[refKF]];
    if (is_start) {
        const float kp_ur = refKF->mvLinesRight[refKF->line_matches[pML->mObservations[refKF]]].startPointX;
        obs << kl_left.startPointX, kl_left.startPointY, kp_ur;
    } else {
        const float kp_ur = refKF->mvLinesRight[refKF->line_matches[pML->mObservations[refKF]]].endPointX;
        obs << kl_left.endPointX, kl_left.endPointY, kp_ur;
    }

//
    g2o::EdgeStereoSE3ProjectXYZ* e1 = new g2o::EdgeStereoSE3ProjectXYZ();
//
    e1->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
    e1->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(refKF->mnId)));
    e1->setMeasurement(obs);
    float invSigma2 = 1.0;
    int cur_oct = 0;
    while (cur_oct < kl_left.octave)//refKF->mvInvLevelSigma2[kl_left.octave];
    {
        invSigma2 /= LinePyrFactor*LinePyrFactor;
        cur_oct+=1;
    }
    Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
    e1->setInformation(Info);

    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
    e1->setRobustKernel(rk);
    rk->setDelta(thHuberStereo);

    e1->fx = refKF->fx;
    e1->fy = refKF->fy;
    e1->cx = refKF->cx;
    e1->cy = refKF->cy;
    e1->bf = refKF->mbf;

    optimizer.addEdge(e1);
    lineMainEdges.push_back(e1);
    mapLines.push_back(pML);
    *v_line_1_p = v_line_1;

    e1->computeError();

    g2o::VertexSE3Expmap* v_frame = dynamic_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(refKF->mnId));
    Eigen::Vector3d X1c = v_frame->estimate().map(X1);
    X1c = e1->fx * X1c / X1c(2);
    X1c(0) += e1->cx;
    X1c(1) += e1->cy;
}

void Optimizer::GlobalBundleAdjustment(Map* pMap, int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    vector<MapPoint*> vpMP = pMap->GetAllMapPoints();
    vector<MapLine*> vpML = pMap->GetAllMapLines();
    BundleAdjustment(vpKFs,vpMP, vpML, nIterations,pbStopFlag, nLoopKF, bRobust);
}


void Optimizer::BundleAdjustment(const vector<KeyFrame *> &vpKFs, const vector<MapPoint *> &vpMP, const vector<MapLine *> &vpML,
                                 int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<bool> vbNotIncludedMP;
    vbNotIncludedMP.resize(vpMP.size());

    vector<bool> vbNotIncludedML;
    vbNotIncludedML.resize(vpML.size());

    g2o::SparseOptimizer optimizer;

    g2o::BlockSolverX::LinearSolverType *linearSolver;
    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>();
    g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX(linearSolver);

    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    long unsigned int maxKFid = 0;

    // Set KeyFrame vertices
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }

    const float thHuber2D = sqrt(5.99);
    const float thHuber3D = sqrt(7.815);
    double thHuberLines = thHuber3D/2.0; //for ICRA / 4.0!!!

    int maxPtId = -1;

    // Set MapPoint vertices
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];
        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        const int id = pMP->mnId+maxKFid+1;
        if (id >= maxPtId)
        {
            maxPtId = id;
        }
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        int nEdges = 0;
        //SET EDGES
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKF = mit->first;
            if(pKF->isBad() || pKF->mnId>maxKFid)
                continue;

            nEdges++;

            const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

            if(pKF->mvuRight[mit->second]<0)
            {
                Eigen::Matrix<double,2,1> obs;
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave*2];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber2D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;

                optimizer.addEdge(e);
            }
            else
            {
                Eigen::Matrix<double,3,1> obs;
                const float kp_ur = pKF->mvuRight[mit->second];
                obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave*2];
                Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                e->setInformation(Info);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber3D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;
                e->bf = pKF->mbf;

                optimizer.addEdge(e);
            }
        }

        if(nEdges==0)
        {
            optimizer.removeVertex(vPoint);
            vbNotIncludedMP[i]=true;
        }
        else
        {
            vbNotIncludedMP[i]=false;
        }
    }

    double infoLines = 1.0;

    for (size_t i = 0; i < vpML.size(); i++)
    {
        MapLine* pML = vpML[i];

        if (!pML || pML->isBad() || pML->Observations() < 4)
        {
            vbNotIncludedML[i] = true;
            continue;
        }

//        pML->UseMinimalAsMainParam();
        int nEdges = 0;
        nEdges = AddLineMinimalGlobal(maxPtId, infoLines, thHuberLines, optimizer, pML);
        if (nEdges == 0)
        {
            vbNotIncludedML[i] = true;
        } else {
            vbNotIncludedML[i] = false;
        }
    }


    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations);

    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        if(nLoopKF==0)
        {
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }
        else
        {
            pKF->mTcwGBA.create(4,4,CV_32F);
            Converter::toCvMat(SE3quat).copyTo(pKF->mTcwGBA);
            pKF->mnBAGlobalForKF = nLoopKF;
        }
    }

    //Points
    for(size_t i=0; i<vpMP.size(); i++)
    {
        if(vbNotIncludedMP[i])
            continue;

        MapPoint* pMP = vpMP[i];

        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));

        if(nLoopKF==0)
        {
            pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
            pMP->UpdateNormalAndDepth();
        }
        else
        {
            pMP->mPosGBA.create(3,1,CV_32F);
            Converter::toCvMat(vPoint->estimate()).copyTo(pMP->mPosGBA);
            pMP->mnBAGlobalForKF = nLoopKF;
        }
    }

    for(size_t i=0; i<vpML.size(); i++) {
        if (vbNotIncludedMP[i])
            continue;

        MapLine *pML = vpML[i];

        if (pML->isBad())
            continue;

        g2o::VertexSBALine *v_line = static_cast<g2o::VertexSBALine *>(optimizer.vertex(pML->mnId + maxPtId + 1));
        Eigen::Matrix3d R_l(v_line->estimate().GetR());
        Eigen::Vector3d line_dir = R_l.col(0);
        Eigen::Vector3d X0 = v_line->estimate().GetAlpha() * R_l.col(1);
        pML->SetMinimalPos(X0, line_dir);

    }

}


void AddLineMinOnlyPose(double info_lines, double deltaLinesStereo, double deltaLinesMono, const Eigen::Matrix3d& Kti, int i, MapLine* pML, g2o::SparseOptimizer& optimizer,
                    Frame* pFrame, std::vector<g2o::EdgeSE3ProjectLineOnlyPose*>& vpEdgesLines, std::vector<size_t>& vnIndexLines, std::vector<bool>& vnStereoLines)
{
    Eigen::Vector3d X0, line_dir;
    pML->GetMinimalPos(&X0, &line_dir);

    Eigen::Matrix3d K;
    K.setIdentity();
    K(0,0) = pFrame->fx;
    K(1,1) = pFrame->fx;
    K(0,2) = pFrame->cx;
    K(1,2) = pFrame->cy;

    double deltaLines = deltaLinesStereo;
    if (pFrame->line_matches[i] < 0)
    {
        deltaLines = deltaLinesMono;
    }

    for (int si = 0; si < 2; si++) {

        if (si == 1 && pFrame->line_matches[i] < 0)
        {
            continue;
        }

        g2o::EdgeSE3ProjectLineOnlyPose *e = new g2o::EdgeSE3ProjectLineOnlyPose ();

        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));

        Eigen::Vector2d obs;
        obs.setZero();
        e->setMeasurement(obs);


        g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
        rk->setDelta(deltaLines);

        e->f = pFrame->fx;
        e->cx = pFrame->cx;
        e->cy = pFrame->cy;

        KeyLine kl;
        if (si == 0) {
            kl = pFrame->mvLinesLeft[i];
        } else {
            if (pFrame->line_matches[i] < 0) {
//                std::cout << "BIG ERROR no stereo line" << std::endl;
                continue;
            }
            kl = pFrame->mvLinesRight[pFrame->line_matches[i]];
        }

        double info = info_lines;

        double thrReproj = GetReprojThrPyramid(1.0, kl.octave);

        info /= thrReproj*thrReproj;

        e->setInformation(Eigen::Matrix2d::Identity() * info);
        Eigen::Vector3d x1, x2;
        x1 << kl.startPointX, kl.startPointY, 1.0;
        x2 << kl.endPointX, kl.endPointY, 1.0;
        e->x1 = x1;
        e->x2 = x2;

        e->X1 = X0;
        e->X2 = X0 + line_dir;

        double b = 0;
        if (si == 1)
        {
            b = -pFrame->mbf / pFrame->mK.at<float>(0,0);
        }
        e->b = Eigen::Vector3d(b,0,0);

        optimizer.addEdge(e);

        vpEdgesLines.push_back(e);
        vnIndexLines.push_back((size_t)i);
        if (pFrame->line_matches[i] < 0)
        {
            vnStereoLines.push_back(false);
        } else {
            vnStereoLines.push_back(true);
        }
    }
}


int Optimizer::PoseOptimization(Frame *pFrame, double gamma)
{

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    const int N = pFrame->N;

    vector<g2o::EdgeSE3ProjectXYZOnlyPose*> vpEdgesMono;
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(N);
    vnIndexEdgeMono.reserve(N);

    vector<g2o::EdgeStereoSE3ProjectXYZOnlyPose*> vpEdgesStereo;
    vector<size_t> vnIndexEdgeStereo;
    vpEdgesStereo.reserve(N);
    vnIndexEdgeStereo.reserve(N);

    vector<g2o::EdgeSE3ProjectLineICRABinOnlyPose*> vpEdgesLines_icra;
    vector<g2o::EdgeSE3ProjectLineOnlyPose*> vpEdgesLines_min;

    vpEdgesLines_min.reserve(pFrame->mvLinesLeft.size());

    vector<size_t> vnIndexLines;
    vector<bool> vnStereoLines;
    vnIndexLines.reserve(pFrame->mvLinesLeft.size());

    const float deltaMono = sqrt(5.991);
    const float deltaStereo = sqrt(7.815);
//    float deltaLines = deltaStereo/4.0;

//    float gamma = 0.25;
    float deltaLinesStereo = deltaStereo;///2.0;
    float deltaLinesMono = deltaMono;///2.0;
    double info_lines = 1.0;

    deltaLinesStereo *= gamma;
    deltaLinesMono *= gamma;
    info_lines *= gamma*gamma;

    Eigen::Matrix3d K_eig;
    cv::cv2eigen(pFrame->mK, K_eig);
    Eigen::Matrix3d Kti = K_eig.transpose().inverse();

    {
        unique_lock<mutex> lock(MapPoint::mGlobalMutex);
//        unique_lock<mutex> lock_line(MapLine::mGlobalMutex);

        for (int i = 0; i < N; i++) {
            MapPoint *pMP = pFrame->mvpMapPoints[i];
            if (pMP) {
                // Monocular observation
                if (pFrame->mvuRight[i] < 0) {
                    nInitialCorrespondences++;
                    pFrame->mvbOutlier[i] = false;

                    Eigen::Matrix<double, 2, 1> obs;
                    const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZOnlyPose *e = new g2o::EdgeSE3ProjectXYZOnlyPose();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));
                    e->setMeasurement(obs);
                    const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity() * invSigma2);

                    g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(deltaMono);

                    e->fx = pFrame->fx;
                    e->fy = pFrame->fy;
                    e->cx = pFrame->cx;
                    e->cy = pFrame->cy;
                    cv::Mat Xw = pMP->GetWorldPos();
                    e->Xw[0] = Xw.at<float>(0);
                    e->Xw[1] = Xw.at<float>(1);
                    e->Xw[2] = Xw.at<float>(2);

                    optimizer.addEdge(e);

                    vpEdgesMono.push_back(e);
                    vnIndexEdgeMono.push_back(i);
                } else  // Stereo observation
                {
                    nInitialCorrespondences++;
                    pFrame->mvbOutlier[i] = false;

                    //SET EDGE
                    Eigen::Matrix<double, 3, 1> obs;
                    const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                    const float &kp_ur = pFrame->mvuRight[i];
                    obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                    g2o::EdgeStereoSE3ProjectXYZOnlyPose *e = new g2o::EdgeStereoSE3ProjectXYZOnlyPose();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex *>(optimizer.vertex(0)));
                    e->setMeasurement(obs);
                    const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                    Eigen::Matrix3d Info = Eigen::Matrix3d::Identity() * invSigma2;
                    e->setInformation(Info);

                    g2o::RobustKernelHuber *rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(deltaStereo);

                    e->fx = pFrame->fx;
                    e->fy = pFrame->fy;
                    e->cx = pFrame->cx;
                    e->cy = pFrame->cy;
                    e->bf = pFrame->mbf;
                    cv::Mat Xw = pMP->GetWorldPos();
                    e->Xw[0] = Xw.at<float>(0);
                    e->Xw[1] = Xw.at<float>(1);
                    e->Xw[2] = Xw.at<float>(2);

                    optimizer.addEdge(e);

                    vpEdgesStereo.push_back(e);
                    vnIndexEdgeStereo.push_back(i);
                }
            }

        }

        for (size_t i = 0; i < pFrame->mvpMapLines.size(); i++) {
            MapLine *pML = pFrame->mvpMapLines[i];
            if (!pML) {
                continue;
            }

            AddLineMinOnlyPose(info_lines, deltaLinesStereo, deltaLinesMono, Kti, i, pML, optimizer, pFrame,
                               vpEdgesLines_min, vnIndexLines, vnStereoLines);
        }
    }

//    std::cout << " added " << vnIndexLines.size()/4 << " lines" << std::endl;

    if(nInitialCorrespondences<3)
        return 0;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2Mono[4]={5.991,5.991,5.991,5.991};
    const float chi2Stereo[4]={7.815,7.815,7.815, 7.815};
    //const int its[4]={10,10,10,10};
    const int its[4]={10,10,10,10};

    int nBad=0;
    for(size_t it=0; it<4; it++)
    {

        vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
        optimizer.initializeOptimization(0);
        optimizer.optimize(its[it]);

        nBad=0;
        for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
        {
            g2o::EdgeSE3ProjectXYZOnlyPose* e = vpEdgesMono[i];

            const size_t idx = vnIndexEdgeMono[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Mono[it])
            {                
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {
                pFrame->mvbOutlier[idx]=false;
                e->setLevel(0);
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        for(size_t i=0, iend=vpEdgesStereo.size(); i<iend; i++)
        {
            g2o::EdgeStereoSE3ProjectXYZOnlyPose* e = vpEdgesStereo[i];

            const size_t idx = vnIndexEdgeStereo[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Stereo[it])
            {
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {                
                e->setLevel(0);
                pFrame->mvbOutlier[idx]=false;
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        if(optimizer.edges().size()<10)
            break;


        for (size_t i = 0; i < vnIndexLines.size(); i++)
        {
            g2o::OptimizableGraph::Edge* e;
            e = static_cast<g2o::OptimizableGraph::Edge*>(vpEdgesLines_min[i]);
            int idx = vnIndexLines[i];
            e->computeError();
            const float chi2 = e->chi2();
            double thr = deltaLinesStereo*deltaLinesStereo;
            if (!vnStereoLines[idx])
            {
                thr = deltaLinesMono*deltaLinesMono;
            }
            if (chi2 > thr)//chi2Stereo[it]*info_lines*0.5*0.5
            {
                pFrame->mvbOutlierLines[idx] = true;
                e->setLevel(1);
            } else {
                pFrame->mvbOutlierLines[idx] = false;
                e->setLevel(0);
            }
            if(it==2)
                e->setRobustKernel(0);
        }
    }    

    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPose(pose);

    int out_cnt = 0;
    for (size_t j = 0; j < vnIndexLines.size(); j++)
    {
        int idx = vnIndexLines[j];
        if (pFrame->mvbOutlierLines[idx])
        {
            out_cnt++;
        }
    }

    return nInitialCorrespondences-nBad;
}



void Optimizer::LocalBundleAdjustment(KeyFrame *pKF, bool* pbStopFlag, Map* pMap, double gamma)
{
    std::list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKF = pKF->mnId;

    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKF = pKF->mnId;
        if(!pKFi->isBad())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapPoints seen in Local KeyFrames
    std::list<MapPoint*> lLocalMapPoints;
    std::list<MapLine*> lLocalMapLines;
    for(std::list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        vector<MapLine *> vpMLs = (*lit)->GetMapLineMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKF=pKF->mnId;
                    }
        }
        for (vector<MapLine *>::iterator vit = vpMLs.begin(), vend = vpMLs.end(); vit != vend; vit++) {
            MapLine *pML = *vit;
            if (pML)
                if (!pML->isBad()) {
                    if (pML->Observations() < 4)
                    {
                        continue;
                    }
                    if (pML->mnBALocalForKF != pKF->mnId) {
                        lLocalMapLines.push_back(pML);
                        pML->mnBALocalForKF = pKF->mnId;
                    }
                }
        }
    }



    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    std::list<KeyFrame*> lFixedCameras;
    for(std::list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKF!=pKF->mnId && pKFi->mnBAFixedForKF!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKF=pKF->mnId;
                if(!pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }
    for (std::list<MapLine *>::iterator lit = lLocalMapLines.begin(), lend = lLocalMapLines.end(); lit != lend; lit++) {
        map<KeyFrame *, size_t> observations = (*lit)->GetObservations();
        for (map<KeyFrame *, size_t>::iterator mit = observations.begin(), mend = observations.end();
             mit != mend; mit++) {
            KeyFrame *pKFi = mit->first;


            if (pKFi->mnBALocalForKF != pKF->mnId && pKFi->mnBAFixedForKF != pKF->mnId) {
                pKFi->mnBAFixedForKF = pKF->mnId;
                if (!pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }


    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType *linearSolver;
    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolverX::PoseMatrixType>();
    g2o::BlockSolverX *solver_ptr = new g2o::BlockSolverX(linearSolver);
    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);


    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    // Set Local KeyFrame vertices
    std::map<int, KeyFrame*> id2kf;
    for(std::list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
        id2kf.insert(std::make_pair(pKFi->mnId, pKFi));
    }

    // Set Fixed KeyFrame vertices
    for(std::list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;

        id2kf.insert(std::make_pair(pKFi->mnId, pKFi));
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMono;
    vpEdgesMono.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSize);

    vector<g2o::EdgeStereoSE3ProjectXYZ*> vpEdgesStereo;
    vpEdgesStereo.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFStereo;
    vpEdgeKFStereo.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeStereo;
    vpMapPointEdgeStereo.reserve(nExpectedSize);

    std::vector<g2o::EdgeStereoSE3ProjectXYZ*> lineMainEdges;

    const float thHuberMono = sqrt(5.991);
    const float thHuberStereo = sqrt(7.815);

    unsigned long maxPtId = 0;

    for(std::list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        unsigned int id = pMP->mnId+maxKFid+1;
        if (id > maxPtId)
        {
            maxPtId = id;
        }
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMono);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMono.push_back(e);
                    vpEdgeKFMono.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
                }
                else // Stereo observation
                {
                    Eigen::Matrix<double,3,1> obs;
                    const float kp_ur = pKFi->mvuRight[mit->second];
                    obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                    g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                    e->setInformation(Info);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberStereo);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;
                    e->bf = pKFi->mbf;

                    optimizer.addEdge(e);
                    vpEdgesStereo.push_back(e);
                    vpEdgeKFStereo.push_back(pKFi);
                    vpMapPointEdgeStereo.push_back(pMP);
                }
            }
        }
    }

    double thHuberLinesStereo = thHuberStereo;///2.0 //for ICRA / 4.0!!!
    double thHuberLinesMono = thHuberMono;///2.0 //for ICRA / 4.0!!!
    LineOptimizer line_optimizer(optimizer, maxPtId, thHuberLinesStereo, thHuberLinesMono, gamma);
    std::vector<std::map<int, std::pair<KeyLine, KeyLine>>> proj_maps;

    for(std::list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {

        MapLine* pML = *lit;
        const map<KeyFrame *, size_t> observations = pML->GetObservations();
        std::map<int, std::pair<KeyLine, KeyLine>> proj_map;
        for (const auto& obs_pair: observations)
        {
            KeyFrame* pKFcurr = obs_pair.first;
            if(pKFcurr->isBad()) {
                continue;
            }
            const KeyLine& kl = pKFcurr->mvLinesLeft[obs_pair.second];
            if (pKFcurr->line_matches[obs_pair.second] >= 0) {
                const KeyLine &kl_right = pKFcurr->mvLinesRight[pKFcurr->line_matches[obs_pair.second]];
                proj_map.insert(std::make_pair(pKFcurr->mnId,
                                               std::make_pair(kl, kl_right)));
            } else {
                KeyLine kl_empty;
                kl_empty.startPointX = -1;
                kl_empty.endPointX = -1;
                kl_empty.startPointY = -1;
                kl_empty.endPointY = -1;
                proj_map.insert(std::make_pair(pKFcurr->mnId,
                                               std::make_pair(kl, kl_empty)));
            }
        }
        Eigen::Matrix3d K_eig;
        cv::cv2eigen(pKF->mK, K_eig);
        Eigen::Vector3d X0, line_dir;
        pML->GetMinimalPos(&X0, &line_dir);
        line_optimizer.AddLineMinimal(pML->mnId, K_eig, pKF->mbf / pKF->mK.at<float>(0,0), X0, line_dir, proj_map);
        proj_maps.push_back(proj_map);
    }

    if(pbStopFlag)
        if(*pbStopFlag)
            return;
    optimizer.initializeOptimization();
    optimizer.optimize(5);

//    std::cout << " init opt done " << std::endl;

    bool bDoMore= true;

    if(pbStopFlag)
        if(*pbStopFlag)
            bDoMore = false;

    if(bDoMore)
    {
        for (int it = 0; it < 1; it++)
        {
            // Check inlier observations
            for (size_t i = 0, iend = vpEdgesMono.size(); i < iend; i++) {
                g2o::EdgeSE3ProjectXYZ *e = vpEdgesMono[i];
                MapPoint *pMP = vpMapPointEdgeMono[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > 5.991 || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }

            for (size_t i = 0, iend = vpEdgesStereo.size(); i < iend; i++) {
                g2o::EdgeStereoSE3ProjectXYZ *e = vpEdgesStereo[i];
                MapPoint *pMP = vpMapPointEdgeStereo[i];

                if (pMP->isBad())
                    continue;

                if (e->chi2() > 7.815 || !e->isDepthPositive()) {
                    e->setLevel(1);
                }

                e->setRobustKernel(0);
            }

            line_optimizer.DisableOutliers();

            // Optimize again without the outliers

            optimizer.initializeOptimization(0);

            optimizer.optimize(15);
        }
    }


    vector<pair<KeyFrame*,MapPoint*> > vToErase;
    vToErase.reserve(vpEdgesMono.size()+vpEdgesStereo.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMono[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    for(size_t i=0, iend=vpEdgesStereo.size(); i<iend;i++)
    {
        g2o::EdgeStereoSE3ProjectXYZ* e = vpEdgesStereo[i];
        MapPoint* pMP = vpMapPointEdgeStereo[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>7.815 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFStereo[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    vector<pair<KeyFrame*,MapLine*> > vToEraseLines;

    std::map<MapLine*, std::pair<Eigen::Vector3d, Eigen::Vector3d>> vToUpdateLines;
    for(std::list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine *pML = *lit;
        Eigen::Vector3d X0, line_dir;
        std::vector<int> bad_obs;
        if (!line_optimizer.GetLineData(pML->mnId, &X0, &line_dir, &bad_obs))
        {
            continue;
        }
        for (size_t boi = 0; boi < bad_obs.size(); boi++)
        {
            vToEraseLines.push_back(std::make_pair(id2kf[bad_obs[boi]], pML));
        }
        vToUpdateLines.insert(std::make_pair(pML, std::make_pair(X0, line_dir)));
    }


//    std::cout << " before map mutex " << std::endl;
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapPoint* pMPi = vToErase[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }


    if(!vToEraseLines.empty())
    {
        for(size_t i=0;i<vToEraseLines.size();i++)
        {

            KeyFrame* pKFi = vToEraseLines[i].first;
            MapLine* pMLi = vToEraseLines[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }


    // Recover optimized data

    //Keyframes
    for(std::list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }

    //Points
    for(std::list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }

    //Lines
    for(auto lit=vToUpdateLines.begin(), lend=vToUpdateLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = lit->first;
        pML->SetMinimalPos(lit->second.first, lit->second.second);
    }

}


void Optimizer::OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections, const bool &bFixScale)
{
    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    g2o::BlockSolver_7_3::LinearSolverType * linearSolver =
           new g2o::LinearSolverEigen<g2o::BlockSolver_7_3::PoseMatrixType>();
    g2o::BlockSolver_7_3 * solver_ptr= new g2o::BlockSolver_7_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

    solver->setUserLambdaInit(1e-16);
    optimizer.setAlgorithm(solver);

    const vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    const vector<MapPoint*> vpMPs = pMap->GetAllMapPoints();

    const unsigned int nMaxKFid = pMap->GetMaxKFid();

    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vScw(nMaxKFid+1);
    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vCorrectedSwc(nMaxKFid+1);
    vector<g2o::VertexSim3Expmap*> vpVertices(nMaxKFid+1);

    const int minFeat = 100;

    // Set KeyFrame vertices
    for(size_t i=0, iend=vpKFs.size(); i<iend;i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSim3Expmap* VSim3 = new g2o::VertexSim3Expmap();

        const int nIDi = pKF->mnId;

        LoopClosing::KeyFrameAndPose::const_iterator it = CorrectedSim3.find(pKF);

        if(it!=CorrectedSim3.end())
        {
            vScw[nIDi] = it->second;
            VSim3->setEstimate(it->second);
        }
        else
        {
            Eigen::Matrix<double,3,3> Rcw = Converter::toMatrix3d(pKF->GetRotation());
            Eigen::Matrix<double,3,1> tcw = Converter::toVector3d(pKF->GetTranslation());
            g2o::Sim3 Siw(Rcw,tcw,1.0);
            vScw[nIDi] = Siw;
            VSim3->setEstimate(Siw);
        }

        if(pKF==pLoopKF)
            VSim3->setFixed(true);

        VSim3->setId(nIDi);
        VSim3->setMarginalized(false);
        VSim3->_fix_scale = bFixScale;

        optimizer.addVertex(VSim3);

        vpVertices[nIDi]=VSim3;
    }


    set<pair<long unsigned int,long unsigned int> > sInsertedEdges;

    const Eigen::Matrix<double,7,7> matLambda = Eigen::Matrix<double,7,7>::Identity();

    // Set Loop edges
    for(map<KeyFrame *, set<KeyFrame *> >::const_iterator mit = LoopConnections.begin(), mend=LoopConnections.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        const long unsigned int nIDi = pKF->mnId;
        const set<KeyFrame*> &spConnections = mit->second;
        const g2o::Sim3 Siw = vScw[nIDi];
        const g2o::Sim3 Swi = Siw.inverse();

        for(set<KeyFrame*>::const_iterator sit=spConnections.begin(), send=spConnections.end(); sit!=send; sit++)
        {
            const long unsigned int nIDj = (*sit)->mnId;
            if((nIDi!=pCurKF->mnId || nIDj!=pLoopKF->mnId) && pKF->GetWeight(*sit)<minFeat)
                continue;

            const g2o::Sim3 Sjw = vScw[nIDj];
            const g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;

            optimizer.addEdge(e);

            sInsertedEdges.insert(make_pair(min(nIDi,nIDj),max(nIDi,nIDj)));
        }
    }

    // Set normal edges
    for(size_t i=0, iend=vpKFs.size(); i<iend; i++)
    {
        KeyFrame* pKF = vpKFs[i];

        const int nIDi = pKF->mnId;

        g2o::Sim3 Swi;

        LoopClosing::KeyFrameAndPose::const_iterator iti = NonCorrectedSim3.find(pKF);

        if(iti!=NonCorrectedSim3.end())
            Swi = (iti->second).inverse();
        else
            Swi = vScw[nIDi].inverse();

        KeyFrame* pParentKF = pKF->GetParent();

        // Spanning tree edge
        if(pParentKF)
        {
            int nIDj = pParentKF->mnId;

            g2o::Sim3 Sjw;

            LoopClosing::KeyFrameAndPose::const_iterator itj = NonCorrectedSim3.find(pParentKF);

            if(itj!=NonCorrectedSim3.end())
                Sjw = itj->second;
            else
                Sjw = vScw[nIDj];

            g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;
            optimizer.addEdge(e);
        }

        // Loop edges
        const set<KeyFrame*> sLoopEdges = pKF->GetLoopEdges();
        for(set<KeyFrame*>::const_iterator sit=sLoopEdges.begin(), send=sLoopEdges.end(); sit!=send; sit++)
        {
            KeyFrame* pLKF = *sit;
            if(pLKF->mnId<pKF->mnId)
            {
                g2o::Sim3 Slw;

                LoopClosing::KeyFrameAndPose::const_iterator itl = NonCorrectedSim3.find(pLKF);

                if(itl!=NonCorrectedSim3.end())
                    Slw = itl->second;
                else
                    Slw = vScw[pLKF->mnId];

                g2o::Sim3 Sli = Slw * Swi;
                g2o::EdgeSim3* el = new g2o::EdgeSim3();
                el->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pLKF->mnId)));
                el->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                el->setMeasurement(Sli);
                el->information() = matLambda;
                optimizer.addEdge(el);
            }
        }

        // Covisibility graph edges
        const vector<KeyFrame*> vpConnectedKFs = pKF->GetCovisiblesByWeight(minFeat);
        for(vector<KeyFrame*>::const_iterator vit=vpConnectedKFs.begin(); vit!=vpConnectedKFs.end(); vit++)
        {
            KeyFrame* pKFn = *vit;
            if(pKFn && pKFn!=pParentKF && !pKF->hasChild(pKFn) && !sLoopEdges.count(pKFn))
            {
                if(!pKFn->isBad() && pKFn->mnId<pKF->mnId)
                {
                    if(sInsertedEdges.count(make_pair(min(pKF->mnId,pKFn->mnId),max(pKF->mnId,pKFn->mnId))))
                        continue;

                    g2o::Sim3 Snw;

                    LoopClosing::KeyFrameAndPose::const_iterator itn = NonCorrectedSim3.find(pKFn);

                    if(itn!=NonCorrectedSim3.end())
                        Snw = itn->second;
                    else
                        Snw = vScw[pKFn->mnId];

                    g2o::Sim3 Sni = Snw * Swi;

                    g2o::EdgeSim3* en = new g2o::EdgeSim3();
                    en->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFn->mnId)));
                    en->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                    en->setMeasurement(Sni);
                    en->information() = matLambda;
                    optimizer.addEdge(en);
                }
            }
        }
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(15);

    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    // SE3 Pose Recovering. Sim3:[sR t;0 1] -> SE3:[R t/s;0 1]
    for(size_t i=0;i<vpKFs.size();i++)
    {
        KeyFrame* pKFi = vpKFs[i];

        const int nIDi = pKFi->mnId;

        g2o::VertexSim3Expmap* VSim3 = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(nIDi));
        g2o::Sim3 CorrectedSiw =  VSim3->estimate();
        vCorrectedSwc[nIDi]=CorrectedSiw.inverse();
        Eigen::Matrix3d eigR = CorrectedSiw.rotation().toRotationMatrix();
        Eigen::Vector3d eigt = CorrectedSiw.translation();
        double s = CorrectedSiw.scale();

        eigt *=(1./s); //[R t/s;0 1]

        cv::Mat Tiw = Converter::toCvSE3(eigR,eigt);

        pKFi->SetPose(Tiw);
    }

    // Correct points. Transform to "non-optimized" reference keyframe pose and transform back with optimized pose
    for(size_t i=0, iend=vpMPs.size(); i<iend; i++)
    {
        MapPoint* pMP = vpMPs[i];

        if(pMP->isBad())
            continue;

        int nIDr;
        if(pMP->mnCorrectedByKF==pCurKF->mnId)
        {
            nIDr = pMP->mnCorrectedReference;
        }
        else
        {
            KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();
            nIDr = pRefKF->mnId;
        }


        g2o::Sim3 Srw = vScw[nIDr];
        g2o::Sim3 correctedSwr = vCorrectedSwc[nIDr];

        cv::Mat P3Dw = pMP->GetWorldPos();
        Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
        Eigen::Matrix<double,3,1> eigCorrectedP3Dw = correctedSwr.map(Srw.map(eigP3Dw));

        cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);
        pMP->SetWorldPos(cvCorrectedP3Dw);

        pMP->UpdateNormalAndDepth();
    }
}

int Optimizer::OptimizeSim3(KeyFrame *pKF1, KeyFrame *pKF2, vector<MapPoint *> &vpMatches1, g2o::Sim3 &g2oS12, const float th2, const bool bFixScale)
{
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();

    g2o::BlockSolverX * solver_ptr = new g2o::BlockSolverX(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // Calibration
    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;

    // Camera poses
    const cv::Mat R1w = pKF1->GetRotation();
    const cv::Mat t1w = pKF1->GetTranslation();
    const cv::Mat R2w = pKF2->GetRotation();
    const cv::Mat t2w = pKF2->GetTranslation();

    // Set Sim3 vertex
    g2o::VertexSim3Expmap * vSim3 = new g2o::VertexSim3Expmap();    
    vSim3->_fix_scale=bFixScale;
    vSim3->setEstimate(g2oS12);
    vSim3->setId(0);
    vSim3->setFixed(false);
    vSim3->_principle_point1[0] = K1.at<float>(0,2);
    vSim3->_principle_point1[1] = K1.at<float>(1,2);
    vSim3->_focal_length1[0] = K1.at<float>(0,0);
    vSim3->_focal_length1[1] = K1.at<float>(1,1);
    vSim3->_principle_point2[0] = K2.at<float>(0,2);
    vSim3->_principle_point2[1] = K2.at<float>(1,2);
    vSim3->_focal_length2[0] = K2.at<float>(0,0);
    vSim3->_focal_length2[1] = K2.at<float>(1,1);
    optimizer.addVertex(vSim3);

    // Set MapPoint vertices
    const int N = vpMatches1.size();
    const vector<MapPoint*> vpMapPoints1 = pKF1->GetMapPointMatches();
    vector<g2o::EdgeSim3ProjectXYZ*> vpEdges12;
    vector<g2o::EdgeInverseSim3ProjectXYZ*> vpEdges21;
    vector<size_t> vnIndexEdge;

    vnIndexEdge.reserve(2*N);
    vpEdges12.reserve(2*N);
    vpEdges21.reserve(2*N);

    const float deltaHuber = sqrt(th2);

    int nCorrespondences = 0;

    for(int i=0; i<N; i++)
    {
        if(!vpMatches1[i])
            continue;

        MapPoint* pMP1 = vpMapPoints1[i];
        MapPoint* pMP2 = vpMatches1[i];

        const int id1 = 2*i+1;
        const int id2 = 2*(i+1);

        const int i2 = pMP2->GetIndexInKeyFrame(pKF2);

        if(pMP1 && pMP2)
        {
            if(!pMP1->isBad() && !pMP2->isBad() && i2>=0)
            {
                g2o::VertexSBAPointXYZ* vPoint1 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D1w = pMP1->GetWorldPos();
                cv::Mat P3D1c = R1w*P3D1w + t1w;
                vPoint1->setEstimate(Converter::toVector3d(P3D1c));
                vPoint1->setId(id1);
                vPoint1->setFixed(true);
                optimizer.addVertex(vPoint1);

                g2o::VertexSBAPointXYZ* vPoint2 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D2w = pMP2->GetWorldPos();
                cv::Mat P3D2c = R2w*P3D2w + t2w;
                vPoint2->setEstimate(Converter::toVector3d(P3D2c));
                vPoint2->setId(id2);
                vPoint2->setFixed(true);
                optimizer.addVertex(vPoint2);
            }
            else
                continue;
        }
        else
            continue;

        nCorrespondences++;

        // Set edge x1 = S12*X2
        Eigen::Matrix<double,2,1> obs1;
        const cv::KeyPoint &kpUn1 = pKF1->mvKeysUn[i];
        obs1 << kpUn1.pt.x, kpUn1.pt.y;

        g2o::EdgeSim3ProjectXYZ* e12 = new g2o::EdgeSim3ProjectXYZ();
        e12->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id2)));
        e12->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e12->setMeasurement(obs1);
        const float &invSigmaSquare1 = pKF1->mvInvLevelSigma2[kpUn1.octave];
        e12->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare1);

        g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
        e12->setRobustKernel(rk1);
        rk1->setDelta(deltaHuber);
        optimizer.addEdge(e12);

        // Set edge x2 = S21*X1
        Eigen::Matrix<double,2,1> obs2;
        const cv::KeyPoint &kpUn2 = pKF2->mvKeysUn[i2];
        obs2 << kpUn2.pt.x, kpUn2.pt.y;

        g2o::EdgeInverseSim3ProjectXYZ* e21 = new g2o::EdgeInverseSim3ProjectXYZ();

        e21->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id1)));
        e21->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e21->setMeasurement(obs2);
        float invSigmaSquare2 = pKF2->mvInvLevelSigma2[kpUn2.octave];
        e21->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare2);

        g2o::RobustKernelHuber* rk2 = new g2o::RobustKernelHuber;
        e21->setRobustKernel(rk2);
        rk2->setDelta(deltaHuber);
        optimizer.addEdge(e21);

        vpEdges12.push_back(e12);
        vpEdges21.push_back(e21);
        vnIndexEdge.push_back(i);
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(5);

    // Check inliers
    int nBad=0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
            optimizer.removeEdge(e12);
            optimizer.removeEdge(e21);
            vpEdges12[i]=static_cast<g2o::EdgeSim3ProjectXYZ*>(NULL);
            vpEdges21[i]=static_cast<g2o::EdgeInverseSim3ProjectXYZ*>(NULL);
            nBad++;
        }
    }

    int nMoreIterations;
    if(nBad>0)
        nMoreIterations=10;
    else
        nMoreIterations=5;

    if(nCorrespondences-nBad<10)
        return 0;

    // Optimize again only with inliers

    optimizer.initializeOptimization();
    optimizer.optimize(nMoreIterations);

    int nIn = 0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
        }
        else
            nIn++;
    }

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2oS12= vSim3_recov->estimate();

    return nIn;
}


} //namespace ORB_SLAM




