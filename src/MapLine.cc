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

#include <include/KeyFrame.h>
#include "MapLine.h"
#include <opencv/cxeigen.hpp>
#include <include/vgl.h>
#include <include/LineMatching.h>

namespace ORB_SLAM2 {
    long unsigned int MapLine::nNextId=0;

    void GetMinimalFromNonMinimal(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir)
    {
        *line_dir = X2-X1;
        *line_dir = (*line_dir) / line_dir->norm();
        *X0 = X1;
        *X0 = *X0 - ((*X0).dot(*line_dir))*(*line_dir);
    }

    mutex MapLine::mGlobalMutex;

    MapLine::MapLine(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, KeyFrame *pRefKF, Map *pMap, int i) :
            mpRefKF(pRefKF), mpMap(pMap), nObs(0), mnBALocalForKF(0), mbBad(false), mnVisible(1),
            mnFound(1), mnLastFrameSeen(0), mnFirstKFid(pRefKF->mnId), mnFuseCandidateForKF(-1), mnTrackReferenceForFrame(-1), tracked_last_id(-1)
    {
        this->X0 = X0;
        this->line_dir = line_dir;
        AddObservation(pRefKF, i);
        pRefKF->AddMapLine(this, i);

        unique_lock<mutex> lock(mpMap->mMutexPointCreation);
        ComputeDistinctiveDescriptors();

        const KeyLine& kl = pRefKF->mvLinesLeft[i];
        mnSaveId = pRefKF->mnFrameId * ((((int)kl.startPointX) * ((int)kl.startPointY))^(((int)kl.endPointX) * ((int)kl.endPointY)));
        mnId=nNextId++;

        int R = (rand() % (int) (255 + 1));
        int G = (rand() % (int) (255 + 1));
        int B = (rand() % (int) (255 + 1));
        color = cv::Scalar(R,G,B);

    }

    void MapLine::AddObservation(KeyFrame* pKF, size_t idx)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        if(mObservations.count(pKF))
            return;
        mObservations[pKF]=idx;

        if(pKF->line_matches[idx]>=0)
            nObs+=2;
        else {
            nObs++;
//            std::cout << " adding line without right" << std::endl;
        }
    }
    void MapLine::EraseObservation(KeyFrame* pKF)
    {
        bool bBad=false;
        {
            unique_lock<mutex> lock(mMutexFeatures);
            if(mObservations.count(pKF))
            {
                int idx = mObservations[pKF];
                if(pKF->line_matches[idx]>=0)
                    nObs-=2;
                else {
                    nObs--;
                }

                mObservations.erase(pKF);

                if(mpRefKF==pKF)
                    mpRefKF=mObservations.begin()->first;

                // If only 3 observations or less, discard point
                if(nObs<=4)
                    bBad=true;
            }
        }

        if(bBad)
            SetBadFlag();
    }

    void MapLine::SetBadFlag()
    {
        map<KeyFrame*,size_t> obs;
        {
            unique_lock<mutex> lock1(mMutexFeatures);
            unique_lock<mutex> lock2(mMutexPos);
            mbBad=true;
            obs = mObservations;
            mObservations.clear();
        }
        for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
        {
            KeyFrame* pKF = mit->first;
            pKF->EraseMapLineMatch(this);
        }

        mpMap->EraseMapLine(this);
    }


    map<KeyFrame*, size_t> MapLine::GetObservations()
    {
        unique_lock<mutex> lock(mMutexFeatures);
        return mObservations;
    }


    void MapLine::ComputeDistinctiveDescriptors()
    {
        // Retrieve all observed descriptors
        vector<cv::Mat> vDescriptors;

        map<KeyFrame*,size_t> observations;

        {
            unique_lock<mutex> lock1(mMutexFeatures);
            if(mbBad)
                return;
            observations=mObservations;
        }

        if(observations.empty())
            return;


        vDescriptors.reserve(observations.size());

        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKF = mit->first;

            if(!pKF->isBad())
            {
                vDescriptors.push_back(pKF->mDescriptorsLines.row(mit->second));
            }
        }

        if(vDescriptors.empty())
            return;

        // Compute distances between them
        const size_t N = vDescriptors.size();

        float Distances[N][N];
        for(size_t i=0;i<N;i++)
        {
            Distances[i][i]=0;
            for(size_t j=i+1;j<N;j++)
            {
                double distij = cv::norm(vDescriptors[i] - vDescriptors[j]);
                Distances[i][j]=distij;
                Distances[j][i]=distij;
            }
        }

        // Take the descriptor with least median distance to the rest
        int BestMedian = INT_MAX;
        int BestIdx = 0;
        for(size_t i=0;i<N;i++)
        {
            vector<float> vDists(Distances[i],Distances[i]+N);
            sort(vDists.begin(),vDists.end());
            int median = vDists[0.5*(N-1)];

            if(median<BestMedian)
            {
                BestMedian = median;
                BestIdx = i;
            }
        }

        {
            unique_lock<mutex> lock(mMutexFeatures);
            mDescriptor = vDescriptors[BestIdx].clone();
        }
    }

    int MapLine::GetIndexInKeyFrame(KeyFrame *pKF)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        if(mObservations.count(pKF))
            return mObservations[pKF];
        else
            return -1;
    }

    bool MapLine::isBad()
    {
        unique_lock<mutex> lock(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        return mbBad;
    }

    void MapLine::IncreaseVisible(int n)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        mnVisible+=n;
    }

    int MapLine::Observations()
    {
        unique_lock<mutex> lock(mMutexFeatures);
        int cnt  = 0;
        for (auto obs: mObservations)
        {
            if (obs.first->line_matches[obs.second]>=0)
            {
                cnt += 2;
            } else {
                cnt ++;
            }
        }
        return cnt;
    }

    bool MapLine::IsInKeyFrame(KeyFrame *pKF)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        return (mObservations.count(pKF));
    }


    float MapLine::GetFoundRatio()
    {
        unique_lock<mutex> lock(mMutexFeatures);
        return static_cast<float>(mnFound)/mnVisible;
    }

    void MapLine::IncreaseFound(int n)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        mnFound+=n;
    }

//    void MapLine::UseMinimalAsMainParam()
//    {
//        unique_lock<mutex> lock2(mGlobalMutex);
//        unique_lock<mutex> lock(mMutexPos);
//
////        Eigen::Vector3d dX1 = X1-X0;
////        Eigen::Vector3d dX2 = X2-X0;
////        double deps = 1e-6;
////        if (fabs(dX1.transpose() * line_dir) > deps  || fabs(dX2.transpose() * line_dir ) > deps) {
//            X1 = X0;
//            X2 = X0 + line_dir;
////        }
//    }

//    void MapLine::UseMinimalAsMainParam(const Eigen::Matrix<double, 3, 4>& T, KeyFrame* KF)
//    {
//        unique_lock<mutex> lock2(mGlobalMutex);
//        unique_lock<mutex> lock(mMutexPos);
//        Eigen::Vector3d dX1 = X1-X0;
//        Eigen::Vector3d dX2 = X2-X0;
//        double deps = 1e-6;
//        if (fabs(dX1.transpose() * line_dir) > deps  || fabs(dX2.transpose() * line_dir ) > deps) {
//            int idx = mObservations[KF];
//            std::vector<double> params;
//            KF->lines_left[idx].GetPointParams(X0, line_dir, T, &params);
//            X1 = X0 + params[0]*line_dir;
//            X2 = X0 + params[1]*line_dir;
//        }
//    }

//    void MapLine::UseNONMinimalAsMainParam()
//    {
//        unique_lock<mutex> lock2(mGlobalMutex);
//        unique_lock<mutex> lock(mMutexPos);
//        if (isnanf(X1(0)) || isnanf(X2(0)))
//        {
//            X1 = X0;
//            X2 = X0 + line_dir;
//            return;
//        }
//        Eigen::Vector3d old_line_dir = line_dir;
//        line_dir = X2-X1;
//        if (line_dir.norm() == 0)
//        {
//            line_dir = old_line_dir;
//            X1 = X0;
//            X2 = X0 + line_dir;
//            return;
//        }
//        line_dir = line_dir / line_dir.norm();
//        X0 = X1 + ((X0-X1).dot( line_dir))*line_dir;
////        double pc = X0.transpose() * line_dir;
////        X0 = X0 - pc*line_dir;
//    }

    void MapLine::OutputProjectionErrors()
    {
        std::cout << "errors for line " << mnId<< " ";
        for (std::map<KeyFrame*, size_t>::iterator it = mObservations.begin(); it != mObservations.end(); it++)
        {
            KeyFrame* pKF = it->first;
            int idx = it->second;
            Eigen::Vector3d leq = GetNormalizedLineEq(pKF->mvLinesLeft[idx], pKF->Ke);
            cv::Mat Tcw = pKF->GetPose();
            Eigen::Matrix3d R;
            Eigen::Vector3d t;
            cv::cv2eigen(Tcw.colRange(0,3).rowRange(0,3), R);
            cv::cv2eigen(Tcw.rowRange(0,3).col(3), t);
            Eigen::Vector3d X1c = R*X1 + t;
            Eigen::Vector3d X2c = R*X2 + t;
            Eigen::Vector3d X1cr = X1c;
            Eigen::Vector3d X2cr = X2c;
            X1cr(0) = X1cr(0) - pKF->mbf / pKF->mK.at<float>(0,0);
            X2cr(0) = X2cr(0) - pKF->mbf / pKF->mK.at<float>(0,0);
            X1c = X1c / X1c(2);
            X2c = X2c / X2c(2);
            X1cr = X1cr / X1cr(2);
            X2cr = X2cr / X2cr(2);
            leq = leq/leq.segment<2>(0).norm();
            double projerr = fabs(X1c.transpose() * leq) + fabs(X2c.transpose() * leq);
            double projerr_right = -1;
            if (pKF->line_matches[idx] >= 0) {
                Eigen::Vector3d leq_right = GetNormalizedLineEq(pKF->mvLinesRight[pKF->line_matches[idx]], pKF->Ke);
                leq_right = leq_right / leq_right.segment<2>(0).norm();
                projerr_right = fabs(X1cr.transpose() * leq_right) + fabs(X2cr.transpose() * leq_right);
            }
            std::cout << " " << pKF->fx * projerr << " ( " << pKF->fx * projerr_right << ")";
            if (isnanf(pKF->fx * projerr))
            {

//                std::cout << leq.transpose() << std::endl;
//                std::cout << X1.transpose() << std::endl;
//                std::cout << X2.transpose() << std::endl;
            }
        }
        std::cout << std::endl;
    }


    void MapLine::Replace(MapLine* pMP)
    {
        if(pMP->mnId==this->mnId)
            return;

        int nvisible, nfound;
        map<KeyFrame*,size_t> obs;
        {
            unique_lock<mutex> lock1(mMutexFeatures);
            unique_lock<mutex> lock2(mMutexPos);
            obs=mObservations;
            mObservations.clear();
            mbBad=true;
            nvisible = mnVisible;
            nfound = mnFound;
            mpReplaced = pMP;
        }

        for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
        {
            // Replace measurement in keyframe
            KeyFrame* pKF = mit->first;

            if(!pMP->IsInKeyFrame(pKF))
            {
                pKF->ReplaceMapLineMatch(mit->second, pMP);
                pMP->AddObservation(pKF,mit->second);
            }
            else
            {
                pKF->EraseMapLineMatch(mit->second);
            }
        }
        pMP->IncreaseFound(nfound);
        pMP->IncreaseVisible(nvisible);
        pMP->ComputeDistinctiveDescriptors();

        mpMap->EraseMapLine(this);
    }

    void MapLine::SetPos(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2)
    {
        unique_lock<mutex> lock2(mGlobalMutex);
        unique_lock<mutex> lock(mMutexPos);
        this->X1 = X1;
        this->X2 = X2;
    }

    void MapLine::SetMinimalPos(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir)
    {
        unique_lock<mutex> lock2(mGlobalMutex);
        unique_lock<mutex> lock(mMutexPos);
        this->X0 = X0;
        this->line_dir = line_dir;
    }

    void MapLine::GetPos(Eigen::Vector3d* X1p, Eigen::Vector3d* X2p)
    {
        unique_lock<mutex> lock2(mGlobalMutex);
        unique_lock<mutex> lock(mMutexPos);
//        LineData det_left = mpRefKF->lines_left[mObservations[mpRefKF]].GetMainDetection();
//        cv::Mat T = mpRefKF->GetPose();
//        Eigen::Matrix<double, 4, 4> T_eig;
//        cv::cv2eigen(T, T_eig);
//        Eigen::Vector3d X0r = T_eig.block<3,3>(0,0)*X0 + T_eig.block<3,1>(0,3);
//        Eigen::Vector3d ldr = T_eig.block<3,3>(0,0)*line_dir;
//        double p_s, p_e;
//        vgl::ReprojectEndpointTo3D(det_left.x_s, X0r, ldr, &p_s);
//        vgl::ReprojectEndpointTo3D(det_left.x_e, X0r, ldr, &p_e);
//        *X1p = X0 + line_dir * p_s;
//        *X2p = X0 + line_dir * p_e;
        *X1p = X1;
        *X2p = X2;
    }

    void MapLine::GetMinimalPos(Eigen::Vector3d* X0p, Eigen::Vector3d* line_dir_p)
    {
        unique_lock<mutex> lock2(mGlobalMutex);
        unique_lock<mutex> lock(mMutexPos);
        if (is_ept)
        {
            GetMinimalFromNonMinimal(X1, X2, X0p, line_dir_p);
        } else {
            *X0p = X0;
            *line_dir_p = line_dir;
        }
    }

    void MapLine::GetDescriptor(cv::Mat* desc)
    {
        unique_lock<mutex> lock(mMutexFeatures);
        mDescriptor.copyTo(*desc);
    }

    bool MapLine::GetMainPoints3D(Eigen::Vector3d* p1, Eigen::Vector3d* p2)
    {
        Eigen::Vector3d X0, ld;
        GetMinimalPos(&X0, &ld);
        int i0 = GetIndexInKeyFrame(mpRefKF);
        if (i0 < 0)
        {
            return false;
        }
        cv::Mat T0 = mpRefKF->GetPoseInverse();
        std::vector<double> params;
        Eigen::Matrix<double, 4, 4> T_eig;
        cv::cv2eigen(T0, T_eig);
        ReprojectKeyLineTo3D(mpRefKF->mvLinesLeft[i0], T_eig, mpRefKF->Ke, X0, ld, p1, p2);
        return true;
    }

    bool MapLine::CheckEndpointsConsistency(const Eigen::Matrix<double, 3, 4>& T_other, const Eigen::Matrix3d& K, const KeyLine& kl2)
    {
        Eigen::Vector3d p1, p2;
        if (GetMainPoints3D(&p1, &p2))
        {
//            p1 = K*T_other.block<3,3>(0,0).transpose() * (p1 - T_other.block<3,1>(0,3));
//            p2 = K*T_other.block<3,3>(0,0).transpose() * (p2 - T_other.block<3,1>(0,3));
//            p1 = p1/p1(2);
//            p2 = p2/p2(2);
            return EndpointsCloseness(K, T_other, p1, p2, kl2);
        } else {
            return false;
        }
    }

//    const LineFragment& MapLine::GetMainDetection(Eigen::Matrix<double, 4, 4>* T_main)
//    {
//        cv::Mat T_cv = mpRefKF->GetPoseInverse();
//        cv::cv2eigen(T_cv, *T_main);
//        return mpRefKF->lines_left[mObservations[mpRefKF]];
//    }

    bool MapLine::CheckOctave(int oct, const Eigen::Vector3d& c)
    {

        int best_oct = -1;
        double best_dist = 1e3;
        double signed_dist = 0;
        Eigen::Vector3d p1, p2;
        GetMainPoints3D(&p1, &p2);
        double dist = (c-p1).norm();
        for (auto& obs: mObservations)
        {
            int c_oct = obs.first->mvLinesLeft[obs.second].octave;
            cv::Mat c_T = obs.first->GetPoseInverse();
            Eigen::Matrix<double, 3, 4> c_T_eig;
            cv::cv2eigen(c_T.rowRange(0, 3).colRange(0, 4), c_T_eig);
            double c_dist = (p1-c_T_eig.col(3)).norm();
            double dist_diff = fabs(dist-c_dist);
            if (dist_diff < best_dist)
            {
                best_dist = dist_diff;
                signed_dist = dist-c_dist;
                best_oct = c_oct;
            }
        }
        if (best_oct<0)
        {
            return false;
        }
        if ((signed_dist>0 && oct <= best_oct) || (signed_dist<=0 && oct >= best_oct))
        {
            return true;
        }
        return false;
    }

}