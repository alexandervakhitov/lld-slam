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

#include <opencv2/core/mat.hpp>
#include <include/vgl.h>
#include "LineMatching.h"
#include "TwoFrameLineMatcher.h"

void TwoFrameLineMatcher::MatchLines(const std::vector<KeyLine>& lines, const std::vector<KeyLine>& other_lines, const cv::Mat& descsLeft,
                                   const cv::Mat& descsRight,
                                   std::vector<int>* descMatches)
{
    if (is_stereo)
    {
        is_matched = std::vector<bool>(lines.size(), false);
        is_other_matched = std::vector<bool>(other_lines.size(), false);
    }

    descMatches->resize(lines.size(), -1);


    for (size_t j = 0; j < lines.size(); j++)
    {
        if (is_matched[j])
        {
            continue;
        }
        double min_d = std::numeric_limits<double>::max();
        double sec_min_d = std::numeric_limits<double>::max();
        int min_j = -1;

        for (size_t oi = 0; oi < other_lines.size(); oi++)
        {
            if (is_other_matched[oi]) {
                continue;
            }

            if (CheckLinePair(lines[j], other_lines[oi], descsLeft.row(j), descsRight.row(oi), &min_d, &sec_min_d))
            {
                min_j = oi;
            }

        }

        if (min_j >= 0)
        {
            is_other_matched[min_j] = true;
        }
        (*descMatches)[j] = min_j;
    }

    int match_cnt = 0;
    for (size_t i = 0; i < is_other_matched.size(); i++)
    {
        if (is_other_matched[i])
        {
            match_cnt++;
        }
    }
}

bool TwoFrameLineMatcher::CheckLinePair(const KeyLine &kl1, const KeyLine &kl2, const cv::Mat& desc, const cv::Mat& otherDesc, double *min_d_p, double *sec_min_d_p)
{
    if (is_stereo && kl1.octave != kl2.octave)
    {
        return false;
    }

    double len_thr = minLineLength;

    if (LineLength(kl1) < len_thr || LineLength(kl2) < len_thr)
    {
        return false;
    }

    double min_d = *min_d_p;

    Eigen::Vector3d X0;
    Eigen::Vector3d line_dir;

    Eigen::Vector3d leftEq = GetNormalizedLineEq(kl1, K);
    Eigen::Vector3d rightEq = GetNormalizedLineEq(kl2, K);
    if (!vgl::TriangulateLine(T.block<3,4>(0,0), T_right.block<3,4>(0,0), leftEq, rightEq, &X0, &line_dir) || X0.norm() < 0.5)
    {
        return false;
    }
    Eigen::Vector3d p1, p2;
    ReprojectKeyLineTo3D(kl1, T, K, X0, line_dir, &p1, &p2);
    if (p1(2)<0 || p2(2)<0)
    {
        return false;
    }


    double d = lineMatcher->MatchLineDescriptors(desc, otherDesc);

    if (d < min_d)
    {
        min_d = d;
    }
    if (min_d < *min_d_p && min_d < tau)
    {
        *sec_min_d_p = *min_d_p;
        *min_d_p = min_d;
        return true;
    }
    return false;
}