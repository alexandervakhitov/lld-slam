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

#ifndef LINEDESCDATACOLLECTOR_VGL_H
#define LINEDESCDATACOLLECTOR_VGL_H

#include <Eigen/Dense>
#include<Eigen/StdVector>
#include <opencv2/core/types.hpp>
#include <opencv2/line_descriptor.hpp>

typedef std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d> > posevector;
typedef cv::line_descriptor::KeyLine KeyLine;

#include <opencv2/core/types.hpp>

namespace vgl {
    bool MultiTriangulateLine(const posevector& Ts,
                                   const std::vector<Eigen::Vector3d>& lines,
                                   Eigen::Vector3d* X0_p, Eigen::Vector3d* line_dir_p);

    bool TriangulateLine(const Eigen::Matrix<double, 3, 4> &T1,
                         const Eigen::Matrix<double, 3, 4> &T2, const Eigen::Vector3d &line2d_n_1,
                         const Eigen::Vector3d &line2d_n_2, Eigen::Vector3d *X0, Eigen::Vector3d *line_dir);

    void NormalizedLineEquation(double sx, double sy, double ex, double ey, const Eigen::Matrix3d& K, Eigen::Vector3d* lineEq);

    //line is represented as [a,b,x,y]: n = [cos(a)cos(b), cos(a)sin(b), sin(a)] is a normal to the plane Pi orthogonal
    // to a line and passing through the origin, (x,y) is coordinate of intersection point between the line and the plane,
    // coords in the plane are defined as orthogonal with axis X being a projection of the z axis (if zero, then y, then x),
    // Y axis is the cross product between normal and X axis
    void Line3DFromPluecker(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, int* ax_code, double* coords);
    void PlueckerFromLine3D(const double coords[4], int ax_code, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir);

    void RLine3DFromPluecker(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, double* coords);
    void PlueckerFromRLine3D(const double coords[4], Eigen::Vector3d* X0, Eigen::Vector3d* line_dir, Eigen::Matrix<double, 3, 4>* J_X0 = NULL, Eigen::Matrix<double, 3, 4>* J_line_dir = NULL);

    bool IsLineCrossingRect(const cv::Point2f& lu, const cv::Point2f& rd, const Eigen::Vector3d& line_eq, std::vector<Eigen::Vector3d>* pps, bool segment_priority=false);

    void ReprojectLinePointTo3D(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Vector2d& pp, const Eigen::Matrix3d& K,
                                double* depth, double* line_param);

    void ProjectLine(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Matrix<double, 3, 4>& T, Eigen::Vector3d* line_eq);

    double LineReprojErrorL1(const Eigen::Vector2d& xs, const Eigen::Vector2d& xe, const Eigen::Matrix<double, 4, 4>& T_other, const Eigen::Vector3d& X0,
                                  const Eigen::Vector3d& line_dir, const Eigen::Matrix3d& K);

    double LineEptReprojError(const Eigen::Vector3d& l, const Eigen::Matrix<double, 4, 4>& T, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, double f);

    void ReprojectEndpointTo3D(const Eigen::Vector3d& endpoint, const Eigen::Vector3d& X0, const Eigen::Vector3d& d_rot, double* p);


    void LineFrom3DEndpoints(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir);

    Eigen::Vector3d MapPoint(const Eigen::Matrix4d& T_c2w, const Eigen::Vector3d& X);

    Eigen::Vector2d ProjectPoint(const Eigen::Matrix3d& K, Eigen::Vector3d& X);
}

#endif //LINEDESCDATACOLLECTOR_VGL_H
