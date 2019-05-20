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

#ifndef ORB_SLAM2_LINEMATCHING_H
#define ORB_SLAM2_LINEMATCHING_H

#include <vector>
#include "KeyFrame.h"
#include "vgl.h"
typedef cv::line_descriptor::KeyLine KeyLine;

void GetHoughCoordinates(const Eigen::Vector3d& leq0, double sx, double sy,
                         std::vector<int>* dist_inds, std::vector<int>* ang_inds,
                         int step_dist, int step_ang);

void SubselectWithGrid(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir,
                       const Eigen::Matrix<double,4,4>& T_last,
                       const Eigen::Matrix3d& K,
                       double sx, double sy,
                       const std::vector<std::vector<std::vector<int>>>& lines_grid,
                       std::vector<int>& other_inds);

double LineLength(const KeyLine& kl1);

bool EndpointsCloseness(const Eigen::Matrix3d& K, const Eigen::Matrix<double, 3, 4>& T_right, const Eigen::Vector3d& p10,
                        const Eigen::Vector3d& p20, const KeyLine& kl2);

Eigen::Matrix<double, 4, 4> GetTForRight(const Eigen::Matrix<double, 4, 4>& T, double b_val);

Eigen::Vector3d GetNormalizedLineEq(const KeyLine& kl, const Eigen::Matrix3d& K);

Eigen::Vector3d GetLineEq(const KeyLine& kl);

namespace LineMatching {
    extern double LinePyrFactor;
}

double GetReprojThrPyramid(double thrReprojLineBase, int pyr_lev);

double GetReprojErrPixelsL1(const KeyLine& kl, const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Matrix4d& T, const Eigen::Matrix3d& K);

void ReprojectKeyLineTo3D(const KeyLine& kl, const Eigen::Matrix4d& T, const Eigen::Matrix3d& K, const Eigen::Vector3d& X0, const Eigen::Vector3d& lineDir, Eigen::Vector3d* X1, Eigen::Vector3d* X2);



#endif //ORB_SLAM2_LINEMATCHING_H
