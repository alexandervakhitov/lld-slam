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



#ifndef ORB_SLAM2_MAPLINE_H
#define ORB_SLAM2_MAPLINE_H

#include <Eigen/Dense>
#include <opencv2/core/mat.hpp>
#include <mutex>
#include <map>
#include <include/LineExtractor.h>

namespace ORB_SLAM2 {

    class KeyFrame;
    class Map;
    class Frame;


    class MapLine
    {
    public:
        MapLine(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, KeyFrame* pRefKF, Map* pMap, int i);
        void AddObservation(KeyFrame* pKF, size_t idx);
        void EraseObservation(KeyFrame* pKF);
        int Observations();
        std::map<KeyFrame*,size_t> GetObservations();
        void ComputeDistinctiveDescriptors();

        void SetPos(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2);

        bool CheckOctave(int oct, const Eigen::Vector3d& c);

        int GetIndexInKeyFrame(KeyFrame *pKF);
        bool IsInKeyFrame(KeyFrame *pKF);

        void IncreaseVisible(int n=1);

        void SetBadFlag();
        bool isBad();

        float GetFoundRatio();

        void IncreaseFound(int n=1);

        void Replace(MapLine* pMP);

        void OutputProjectionErrors();

        void GetPos(Eigen::Vector3d* X1p, Eigen::Vector3d* X2p);
        void GetMinimalPos(Eigen::Vector3d* X0p, Eigen::Vector3d* line_dir_p);

        void GetDescriptor(cv::Mat* desc);

        void SetMinimalPos(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir);

        bool GetMainPoints3D(Eigen::Vector3d* p1, Eigen::Vector3d* p2);

        bool CheckEndpointsConsistency(const Eigen::Matrix<double, 3, 4>& T_other, const Eigen::Matrix3d& K, const KeyLine& kl2);


        KeyFrame* mpRefKF;
        Map* mpMap;


        std::map<KeyFrame*,size_t> mObservations;
        int nObs;

        long unsigned int mnBALocalForKF;

        std::mutex mMutexPos;
        std::mutex mMutexFeatures;
        bool mbBad;

        static std::mutex mGlobalMutex;

        long unsigned int mnId;
        long unsigned int mnSaveId;
        static long unsigned int nNextId;

        int mnVisible;
        int mnFound;

        bool mbTrackInView;
        Eigen::Vector3d track_leq, track_leq_right;

        int mnLastFrameSeen;
        int mnFirstKFid;

        MapLine* mpReplaced;

        int mnFuseCandidateForKF;

        int mnTrackReferenceForFrame;

        cv::Scalar color;

        int tracked_last_id;

    protected:
        Eigen::Vector3d X0, line_dir;
        Eigen::Vector3d X1, X2;
        cv::Mat mDescriptor;
        bool is_ept;

    };
}

#endif //ORB_SLAM2_MAPLINE_H
