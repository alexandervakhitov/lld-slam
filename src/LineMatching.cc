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

#include <Eigen/Dense>
#include <opencv/cxeigen.hpp>
#include "LineMatching.h"
#include "vgl.h"
#include "KeyFrame.h"

double LineMatching::LinePyrFactor = 1.44;

bool EndpointsCloseness(const Eigen::Matrix3d& K, const Eigen::Matrix<double, 3, 4>& T_right, const Eigen::Vector3d& p10, const Eigen::Vector3d& p20, const KeyLine& kl2)
{
    Eigen::Vector3d p1 = K*T_right.block<3,3>(0,0).transpose() * (p10 - T_right.block<3,1>(0,3));
    Eigen::Vector3d p2 = K*T_right.block<3,3>(0,0).transpose() * (p20 - T_right.block<3,1>(0,3));
    p1 = p1/p1(2);
    p2 = p2/p2(2);
    Eigen::Vector3d xs2, xe2;
    xs2 << kl2.startPointX, kl2.startPointY, 1.0;
    xe2 << kl2.endPointX, kl2.endPointY, 1.0;
    double d11 = (xs2-p1).norm();
    double d12 = (xs2-p2).norm();
    double d21 = (xe2-p1).norm();
    double d22 = (xe2-p2).norm();
    double thr = LineLength(kl2);
    if ((d11<thr && d22 < thr) || (d12 < thr && d21 < thr))
    {
        return true;
    }
    return false;
}

double LineLength(const KeyLine& kl1)
{

    Eigen::Vector2d pt_s;
    pt_s << kl1.startPointX , kl1.startPointY;
    Eigen::Vector2d pt_e;
    pt_e << kl1.endPointX , kl1.endPointY;

    return (pt_s-pt_e).norm();
}

#define PI 3.14159265

void GetHoughCoordinates(const Eigen::Vector3d& leq0,
                         double sx, double sy,
                         std::vector<int>* dist_inds, std::vector<int>* ang_inds,
                         int step_dist, int step_ang)
{

    dist_inds->clear();
    ang_inds->clear();
    Eigen::Vector3d leq = leq0;
    leq(0) /= sx;
    leq(1) /= sy;
    leq = leq/leq.segment<2>(0).norm();
//    std::cout << leq.transpose() << std::endl;
    //to make y's coefficient leq(1) >= 0
    if (leq(1) < 0)
    {
        leq = -leq;
    }

    int dist_cell_num = FRAME_DIST_CELLS; //lines_grid.size();
    int ang_cell_num = FRAME_ANG_CELLS; //lines_grid[0].size();
    double dist_level = fabs(leq(2)/(sqrt(2.0)))*dist_cell_num;
    int dist_ind = floor(dist_level+0.5);
    dist_ind = std::min(dist_ind, dist_cell_num-1);
    dist_ind = std::max(dist_ind, 0);

    int shift_dist = -1;
    if (dist_level - dist_ind < 0)
    {
        shift_dist = 1;
    }

    double ang = atan2(leq(1), leq(0));//from 0 to pi
    double ang_level = ang/PI*ang_cell_num;
    int ang_ind = floor(ang_level+0.5);
    ang_ind = std::min(ang_ind, ang_cell_num-1);
    ang_ind = std::max(ang_ind, 0);

    int shift_ang = -1;
    if (ang_level - ang_ind < 0)
    {
        shift_ang = 1;
    }

    int ang_max = std::max(ang_ind, ang_ind+shift_ang);
    for (int i = ang_max; i < ang_max+step_ang; i++)
    {
        int ang_corr = i;
        if (i < 0)
        {
            ang_corr += ang_cell_num;
        }
        ang_corr = ang_corr % ang_cell_num;
        ang_inds->push_back(ang_corr);
    }

    int ang_min = std::min(ang_ind, ang_ind+shift_ang);
    for (int i = ang_min; i > ang_min-step_ang; i--)
    {
        int ang_corr = i;
        if (i < 0)
        {
            ang_corr += ang_cell_num;
        }
        ang_corr = ang_corr % ang_cell_num;

        ang_inds->push_back(ang_corr);
    }

    int dist_max = std::max(dist_ind, dist_ind + shift_dist);
    for (int i = dist_max; i < dist_max+step_dist; i++)
    {
        if (i>=0 && i < dist_cell_num-1)
        {
            dist_inds->push_back(i);
        }
    }

    int dist_min = std::min(dist_ind, dist_ind + shift_dist);
    for (int i = dist_min; i > dist_min-step_dist; i--)
    {
        if (i>=0 && i < dist_cell_num-1)
        {
            dist_inds->push_back(i);
        }
    }

    return;

}

void SubselectWithGrid(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir,
                       const Eigen::Matrix<double,4,4>& T_last, const Eigen::Matrix3d& K, double sx, double sy,
                       const std::vector<std::vector<std::vector<int>>>& lines_grid,
                       std::vector<int>& other_inds)
{
    Eigen::Vector3d Xl1 = K*T_last.block<3,3>(0,0).transpose() * (X0 - T_last.block<3,1>(0,3));
    Eigen::Vector3d Xl2 = K*T_last.block<3,3>(0,0).transpose() * (X0 + line_dir - T_last.block<3,1>(0,3));
    Eigen::Vector3d leq = Xl1.cross(Xl2);
    leq = leq / leq.segment<2>(0).norm();
//    std::cout << " obtained leq " << leq.transpose() <<std::endl;
    std::vector<int> dis, ais;
//    GetHoughCoordinates(leq, sx, sy, &dis, &ais, 10, 10, max_dist);
    GetHoughCoordinates(leq, sx, sy, &dis, &ais, 3, 3);
    other_inds.clear();
    std::set<int> inds;
    for (int ai: ais) {
        for (int di: dis) {
            for (int oi: lines_grid[di][ai])
            {
                inds.insert(oi);
            }
        }
    }
    std::copy(inds.begin(), inds.end(), std::back_inserter(other_inds));

    return;

    cv::Mat debug_grid = cv::Mat::zeros(10*FRAME_DIST_CELLS , 10*FRAME_ANG_CELLS, CV_8UC3 );
    size_t max_val = 0;
    for (int i = 0; i < FRAME_DIST_CELLS; i++)
    {
        for (int j = 0; j < FRAME_ANG_CELLS; j++)
        {
            if (lines_grid[i][j].size() > max_val)
            {
                max_val = lines_grid[i][j].size();
            }
        }
    }

    for (int i = 0; i < FRAME_DIST_CELLS; i++)
    {
        for (int j = 0; j < FRAME_ANG_CELLS; j++)
        {
            for (int k = 0; k < 10; k++)
            {
                for (int l = 0; l < 10; l++)
                {
                    debug_grid.at<cv::Vec3b>(10*i+k, 10*j+l)[0] = (255*lines_grid[i][j].size())/max_val;
                }
            }
        }
    }


    for (int ai: ais)
    {
        for (int di: dis)
        {
            for (int k = 0; k < 10; k++)
            {
                for (int l = 0; l < 10; l++)
                {
                    debug_grid.at<cv::Vec3b>(10*di+k, 10*ai+l)[1] = 255;
                }
            }
        }
    }
//
    cv::imshow("grid", debug_grid);
    cv::waitKey(0);
}


Eigen::Matrix<double, 4, 4> GetTForRight(const Eigen::Matrix<double, 4, 4>& T, double b_val)
{
    Eigen::Matrix<double, 4, 4> T_right;
    T_right = T;
    Eigen::Vector3d b;
    b.setZero();
    b(0) = b_val;
    T_right.block<3, 1>(0, 3) = T_right.block<3, 1>(0, 3) + T_right.block<3, 3>(0, 0) * b;
    return T_right;
};

double GetReprojThrPyramid(double thrReprojLineBase, int pyr_lev)
{
    double thrReprojLine = thrReprojLineBase;
    for (int oi = 0; oi < pyr_lev; oi++)
    {
        thrReprojLine *= LineMatching::LinePyrFactor;
    }
    return thrReprojLine;
}

Eigen::Vector3d GetNormalizedLineEq(const KeyLine& kl, const Eigen::Matrix3d& K )
{
    Eigen::Vector3d leq;
    vgl::NormalizedLineEquation(kl.startPointX, kl.startPointY, kl.endPointX, kl.endPointY, K, &leq);
    return leq;
}

Eigen::Vector3d GetLineEq(const KeyLine& kl)
{
    Eigen::Vector3d Xs(kl.startPointX, kl.startPointY, 1.0);
    Eigen::Vector3d Xe(kl.endPointX, kl.endPointY, 1.0);
//    std::cout << Xs.transpose() << std::endl;
//    std::cout << Xe.transpose() << std::endl;

    Eigen::Vector3d leq = Xs.cross(Xe);
//    std::cout << leq.transpose() << std::endl;

    return leq;
}


double GetReprojErrPixelsL1(const KeyLine& kl, const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Matrix4d& T, const Eigen::Matrix3d& K)
{
    Eigen::Vector2d xs(kl.startPointX, kl.startPointY);
    Eigen::Vector2d xe(kl.endPointX, kl.endPointY);
    return vgl::LineReprojErrorL1(xs, xe, T, X0, line_dir, K);
}

void ReprojectKeyLineTo3D(const KeyLine& kl, const Eigen::Matrix4d& T, const Eigen::Matrix3d& K, const Eigen::Vector3d& X0, const Eigen::Vector3d& lineDir, Eigen::Vector3d* X1, Eigen::Vector3d* X2)
{
//    std::cout << "x0 " << X0 << " d " << lineDir << std::endl;
    Eigen::Vector3d X0rot = vgl::MapPoint(T, X0);
    Eigen::Vector3d dirRot = T.block<3,3>(0,0).transpose() * lineDir;
    double d,p;
    Eigen::Vector2d x_s(kl.startPointX, kl.startPointY);
//    std::cout << "x0r" << X0rot << " dr " << dirRot << std::endl;
    vgl::ReprojectLinePointTo3D(X0rot, dirRot, x_s, K, &d, &p);
//    std::cout << " reprojected " << x_s << " obtained " << p << " " << d << std::endl;
    *X1 = X0+ p*lineDir;
    Eigen::Vector2d x_e(kl.endPointX, kl.endPointY);
    vgl::ReprojectLinePointTo3D(X0rot, dirRot, x_e, K, &d, &p);
//    std::cout << " reprojected " << x_e << " obtained " << p << std::endl;
    *X2 = X0 + p*lineDir;
}