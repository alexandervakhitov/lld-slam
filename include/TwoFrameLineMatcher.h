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

#ifndef ORB_SLAM2_LRLINEMATCHER_H
#define ORB_SLAM2_LRLINEMATCHER_H

#include "LineMatching.h"
#include "MapLine.h"

using namespace ORB_SLAM2;

class TwoFrameLineMatcher {
public:
    TwoFrameLineMatcher(const Eigen::Matrix3d& K,
                   double b,
                   double tau, int minLineLength, LineMatcher* lineMatcher):
            K(K), tau(tau), is_stereo(true), lineMatcher(lineMatcher), minLineLength(minLineLength)
    {
        T.setIdentity();
        T_right = GetTForRight(T, b);
    };

    void MatchLines(const std::vector<KeyLine>& lines, const std::vector<KeyLine>& other_lines, const cv::Mat& descsLeft,
                        const cv::Mat& descsRight, std::vector<int>* desc_matches);
    bool CheckLinePair(const KeyLine& lf, const KeyLine& lf2, const cv::Mat& desc1, const cv::Mat& desc2, double* min_d_p, double* sec_min_d_p);
private:
    Eigen::Matrix3d K;
    Eigen::Matrix4d T;
    Eigen::Matrix4d T_right;
    double tau;
    bool is_descdist_passed, is_triangang2_passed, is_ovlap_passed, is_farclose_passed, is_X0_passed, is_triangang_passed, is_3dlen_passed, is_len_passed;
    bool is_stereo;
    LineMatcher* lineMatcher;
    int minLineLength;
    std::vector<bool> is_matched, is_other_matched;
};


#endif //ORB_SLAM2_LRLINEMATCHER_H
