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

#include "vgl.h"
#include <iostream>
#include <opencv/cv.hpp>
#include <opencv/cxeigen.hpp>


//
bool vgl::MultiTriangulateLine(const posevector& Ts,
                          const std::vector<Eigen::Vector3d>& lines,
                          Eigen::Vector3d* X0_p, Eigen::Vector3d* line_dir_p)
{
    if (lines.size() < 3 || lines.size() != Ts.size())
    {
        return false;
    }
    std::vector<Eigen::Vector3d> normals;
    for (size_t i = 0; i < lines.size(); i++)
    {
        Eigen::Vector3d leq = lines[i];
        leq = leq/leq.norm();
        normals.push_back(Ts[i].block<3,3>(0,0)*leq);
    }

    for (size_t i = 1; i < normals.size(); i++)
    {
        if (fabs(normals[0].dot(normals[i]))/normals[0].norm() / normals[i].norm()>0.975)
        {
            return false;
        }
    }

    Eigen::MatrixXd M2(lines.size(), 3);
    for (size_t i = 0; i < lines.size(); i++)
    {
        M2.row(i) = normals[i].transpose();
    }
    Eigen::Matrix3d Vm = M2.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).matrixV();
    Eigen::Vector3d line_dir = Vm.col(2);

    Eigen::MatrixXd M1(lines.size(), 3);
    Eigen::VectorXd b1(lines.size());
    for (size_t i = 0; i < lines.size(); i++)
    {
        M1.row(i) = normals[i].transpose();
        b1(i) = normals[i].dot(Ts[i].block<3,1>(0,3));
    }
//    M1.row(lines.size()) = line_dir.transpose();
//    b1(lines.size()) = 0;
    Eigen::Vector3d X0 = M1.colPivHouseholderQr().solve(b1);

    X0 = X0 - (X0.dot(line_dir))*line_dir;
    *line_dir_p = line_dir;
    *X0_p = X0;
    return true;
}


bool vgl::TriangulateLine(const Eigen::Matrix<double, 3, 4>& T1,
                     const Eigen::Matrix<double, 3, 4>& T2, const Eigen::Vector3d& line2d_n_1,
                     const Eigen::Vector3d& line2d_n_2, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir)
{
    Eigen::Vector3d normal_1 = T1.block<3,3>(0,0) * line2d_n_1;
    Eigen::Vector3d normal_2 = T2.block<3,3>(0,0) * line2d_n_2;

    if (fabs(normal_1.dot(normal_2)) / normal_1.norm() / normal_2.norm() > 0.975)
    {
        return false;
    }
    *line_dir = normal_1.cross(normal_2);
    *line_dir = (*line_dir) / line_dir->norm();

    Eigen::Matrix3d M;
    Eigen::Vector3d b;
    M.row(0) = normal_1.transpose();
    M.row(1) = normal_2.transpose();
    b(0) = normal_1.transpose() * T1.block<3,1>(0, 3);
    b(1) = normal_2.transpose() * T2.block<3,1>(0, 3);
    b(2) = 0;
    M.row(2) = line_dir->transpose();
    auto qr = M.colPivHouseholderQr();
    if (qr.rank() < 3)
    {
        return false;
    }
    *X0 = qr.solve(b);

//    double res = line2d_n_2.transpose() * T2.block<3,3>(0,0).transpose() * (*X0 + *line_dir - T2.block<3,1>(0,3) );
//    std::cout << res << std::endl;
//    std::cout << line2d_n_1.transpose() * T1.block<3,3>(0,0).transpose() * (*X0 + *line_dir - T1.block<3,1>(0,3 ) << std::endl;
//
//    std::cout << line2d_n_2.transpose() * T2.block<3,3>(0,0).transpose() * (*X0 - T2.block<3,1>(0,3) ) << std::endl;
//    std::cout << line2d_n_2.transpose() * T2.block<3,3>(0,0).transpose() * (*X0 + *line_dir - T2.block<3,1>(0,3) ) << std::endl;

    return true;
//    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullV);
//    *X0 = svd.matrixV().col(2);


//    Eigen::Vector3d test0 = M * (*X0);

//    std::cout << test0.norm() << std::endl;

}

int GetAxisCode(const Eigen::Vector3d& line_dir)
{
    Eigen::Matrix3d rem = Eigen::Matrix3d::Identity() - line_dir * line_dir.transpose();
    int ax_code = -1;
    double abs_proj = 0;
    Eigen::Vector3d ax;
    do
    {
        ax_code++;
        ax = Eigen::Vector3d::Zero();
        ax(2-ax_code) = 1.0;
        abs_proj = (rem * ax).norm();
    } while (abs_proj < 1e-8);
    return ax_code;
}

void GetAxes(int ax_code, const Eigen::Vector3d& line_dir, Eigen::Vector3d* x_axis, Eigen::Vector3d* y_axis)
{
    Eigen::Vector3d ax;
    ax.setZero();
    ax(ax_code) = 1.0;
    *x_axis = ax - line_dir * (line_dir.transpose() * ax);
    *x_axis = *x_axis / (x_axis->norm());
    *y_axis = line_dir.cross(*x_axis);
}

void vgl::Line3DFromPluecker(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, int* ax_code_p, double* coords_p)
{
    double b = atan2(line_dir(1), line_dir(0));
    double a = asin(line_dir(2)); // cos(a) always positive

    int ax_code = GetAxisCode(line_dir);
    Eigen::Vector3d x_axis, y_axis;
    GetAxes(ax_code, line_dir, &x_axis, &y_axis);
    double x = x_axis.dot(X0);
    double y = y_axis.dot(X0);

    coords_p[0] = a;
    coords_p[1] = b;
    coords_p[2] = x;
    coords_p[3] = y;
    *ax_code_p = ax_code;
}

void vgl::PlueckerFromLine3D(const double coords[4], int ax_code, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir)
{
    double a = coords[0];
    double b = coords[1];
    *line_dir << cos(a)*cos(b), cos(a)*sin(b), sin(a);
    double x = coords[2];
    double y = coords[3];
    Eigen::Vector3d x_axis, y_axis;
    GetAxes(ax_code, *line_dir, &x_axis, &y_axis);
    *X0 = x*x_axis + y*y_axis;
}

void vgl::RLine3DFromPluecker(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, double* coords)
{
    cv::Mat R(3,3,CV_64FC1);
    cv::Mat r1(3,1,CV_64FC1);
    cv::eigen2cv(line_dir, r1);
    cv::Mat r2(3,1,CV_64FC1);
    double alpha = X0.norm();
    Eigen::Vector3d X0n = X0/alpha;
    cv::eigen2cv(X0n, r2);
    cv::Mat r3 = r1.cross(r2);

    r1.copyTo(R.col(0));
    r2.copyTo(R.col(1));
    r3.copyTo(R.col(2));
    cv::Mat rvec;
//    std::cout << R << std::endl;
    cv::Rodrigues(R, rvec);
    cv::Mat Rtest;
    cv::Rodrigues(rvec, Rtest);
//    std::cout << cv::norm(R-Rtest) << std::endl;
    for (int i = 0; i < 3; i++)
    {
        coords[i] = rvec.at<double>(i,0);
    }
    coords[3] = alpha;
}

void vgl::PlueckerFromRLine3D(const double coords[4], Eigen::Vector3d* X0, Eigen::Vector3d* line_dir, Eigen::Matrix<double, 3, 4>* J_X0, Eigen::Matrix<double, 3, 4>* J_line_dir)
{
    cv::Mat rvec(3, 1, CV_64FC1);
    for (int i = 0; i < 3; i++)
    {
        rvec.at<double>(i, 0) = coords[i];
    }
    cv::Mat R, JR;
    cv::Rodrigues(rvec, R, JR);
    cv::Mat r1 = R.col(0);
    cv::cv2eigen(r1, *line_dir);
    cv::Mat r2 = R.col(1);
    cv::cv2eigen(r2, *X0);
    *X0 = *X0 * coords[3];
    if (J_X0)
    {
        //assume column-major order of R
//        cv::Mat J_X0_theta_cv = coords[3]*JR(cv::Rect(3,0,3,3)).t();
//        Eigen::Matrix3d J_X0_theta;
//        cv::cv2eigen(J_X0_theta_cv, J_X0_theta);
//        J_X0->block<3,3>(0,0) = J_X0_theta;
//        Eigen::Vector3d r2_eig;
//        cv::cv2eigen(r2, r2_eig);
//        J_X0->block<3,1>(0,3) = r2_eig;
//
//        cv::Mat J_d_theta_cv = JR(cv::Rect(0,0,3,3)).t();
//        Eigen::Matrix3d J_d_theta;
//        cv::cv2eigen(J_d_theta_cv, J_d_theta);
//        J_line_dir->block<3,3>(0,0) = J_d_theta;
//        J_line_dir->block<3,1>(0,3) = Eigen::Vector3d::Zero();
    }
}

bool vgl::IsLineCrossingRect(const cv::Point2f& lu, const cv::Point2f& rd, const Eigen::Vector3d& line_eq,
                             std::vector<Eigen::Vector3d>* pps, bool segment_priority)
{
    std::vector<Eigen::Vector3d> corners;
    int ci = 0;
    int cj = 0;
    for (int ti = 0; ti < 4; ti++)
    {
        if (ti % 2 == 0)
        {
            ci += 1;
        } else {
            cj += 1;
        }
        int ri = ci % 2;
        int rj = cj % 2;
        Eigen::Vector3d v(lu.x + (rd.x-lu.x)*(ri+0.0), lu.y + (rd.y-lu.y)*(rj+0.0), 1.0);
        corners.push_back(v );
    }
    if (pps != NULL && !segment_priority)
    {
        pps->clear();
    }

    bool rv = false;
    Eigen::Vector3d center_pt(0.5*(lu.x+rd.x), 0.5*(lu.y+rd.y), 1.0);
    for (int i = 0; i < 4; i++)
    {
        Eigen::Vector3d v1 = corners[i];
        int j = i+1;
        if (j == 4)
        {
            j = 0;
        }
        Eigen::Vector3d v2 = corners[j];
        double p1 = (line_eq.transpose() * v1);
        double p2 = (line_eq.transpose() * v2);
        if (p1 * p2 < 0)
        {
            Eigen::Vector3d side_vec = v1.cross(v2);
            Eigen::Vector3d pp = side_vec.cross(line_eq);
            if (pps != NULL)
            {
                if (segment_priority)
                {
                    int pt_id = -1;
                    for (int k = 0; k < 2; k++)
                    {
                        double p1 = (side_vec.transpose() * (*pps)[k]);
                        double p2 = (side_vec.transpose() * center_pt);
                        if (p1 * p2 < 0)
                        {
                            if (pt_id < 0)
                            {
                                pt_id = k;
                            } else {
                                rv = false;
                                return rv;
                            }
                        }
                    }
                    if (pt_id >= 0)
                    {
                        (*pps)[pt_id] = pp/pp(2);
                    }
                } else {
                    pps->push_back(pp / pp(2));
                }
            }
            rv = true;
        }
    }

    return rv;
}



//bool vgl::IsLineProjectedToCamera(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Matrix3d& K_inv,
//                                  const Eigen::Matrix<double, 3, 4>& T, const cv::Size& im_size, Eigen::Vector3d* line_eq, std::vector<Eigen::Vector3d>* pps)
//{
//
//    ProjectLine(X0, line_dir, T, line_eq);
//
//    Eigen::Vector3d p1 = T.block<3,3>(0,0).transpose() * (X0 - T.block<3,1>(0,3));
//    Eigen::Vector3d dir = T.block<3,3>(0,0).transpose() * line_dir;
//
//    Eigen::Vector3d lu = K_inv * Eigen::Vector3d(0,0,1);
//    Eigen::Vector3d rd = K_inv * Eigen::Vector3d(im_size.width,im_size.height,1);
//    bool rv = IsLineCrossingRect(cv::Point2f(lu(0), lu(1)), cv::Point2f(rd(0), rd(1)), *line_eq, pps);
//    if (pps == NULL)
//    {
//        return rv;
//    }
//    rv = false;
//    for (int i = 0; i < pps->size(); i++)
//    {
//
//        double depth;
//        double line_param;
//        Eigen::Vector3d pp = (*pps)[i];
//        ReprojectLinePointTo3D(p1, dir, pp, &depth, &line_param);
//        if (pp(2)*depth > 0)
//        {
//            rv = true;
//        }
//    }
//    return rv;
//}

void vgl::ReprojectLinePointTo3D(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Vector2d& pp, const Eigen::Matrix3d& K,
                            double* depth, double* line_param)
{
    Eigen::Matrix<double, 3, 2> M;
    M.col(1) = - K*line_dir;
    Eigen::Vector3d pp_h (pp(0), pp(1), 1.0);
    M.col(0) = pp_h;
    Eigen::Matrix<double, 2, 1> sol = M.colPivHouseholderQr().solve(K*X0);
    *depth = sol(0);
    *line_param = sol(1);
}

void vgl::ProjectLine(const Eigen::Vector3d& X0, const Eigen::Vector3d& line_dir, const Eigen::Matrix<double, 3, 4>& T, Eigen::Vector3d* line_eq)
{
    Eigen::Vector3d p1 = T.block<3,3>(0,0).transpose() * (X0 - T.block<3,1>(0,3));
    Eigen::Vector3d p2 = T.block<3,3>(0,0).transpose() * (X0+line_dir - T.block<3,1>(0,3));
    *line_eq = p1.cross(p2);
    *line_eq = (*line_eq)/line_eq->segment<2>(0).norm();
}

void GetMinMax(const int& a, const int& b, int* min_p, int* max_p)
{
    int min_val = std::min(a,b);
    int max_val = std::max(a,b);
    *max_p = max_val;
    *min_p = min_val;
}


//bool vgl::ReprojBehind(const Eigen::Vector3d& xs3, const Eigen::Vector3d& xe3, const Eigen::Matrix<double, 4, 4>& T_other, const Eigen::Vector3d& X0,
//                            const Eigen::Vector3d& line_dir)
//{
//    Eigen::Vector3d X0c = T_other.block<3, 3>(0, 0).transpose() * (X0 - T_other.block<3, 1>(0, 3));
//    Eigen::Vector3d line_dir_c = T_other.block<3, 3>(0, 0).transpose() * line_dir;
//    double d3, p3;
//    ReprojectLinePointTo3D(X0c, line_dir_c, xs3, &d3, &p3);
//    double d4, p4;
//    ReprojectLinePointTo3D(X0c, line_dir_c, xe3, &d4, &p4);
//}
void GetYSortedEndPoints(const KeyLine& kl1, int* x1_s, int* x1_e, int* y1_s, int* y1_e)
{
    if (kl1.startPointY < kl1.endPointY)
    {
        *x1_s = kl1.startPointX;
        *y1_s = kl1.startPointY;
        *x1_e = kl1.endPointX;
        *y1_e = kl1.endPointY;
    } else {
        *x1_s = kl1.endPointX;
        *y1_s = kl1.endPointY;
        *x1_e = kl1.startPointX;
        *y1_e = kl1.startPointY;
    }
}

double GetXByY(const int& x_s, const int& x_e, const int& y_s, const int& y_e, int y_ref)
{
    return (double(y_ref-y_s) * x_e + double (y_e - y_ref) * x_s)/ double(y_e-y_s);
}

bool RefinePoint(const double& x_left, const double& x_right_init, const int& y, const cv::Mat& im_left,
                 const cv::Mat& im_right, const int& w, const int& L, double* x_right_fin)
{
//    const int w = 5;
//    const int L = 10;
    if (y < w || y+w >= im_left.rows)
    {
        return false;
    }

    int x_left_i = int(floor(x_left+0.5));
    if (x_left_i < w || x_left_i+w >= im_left.cols)
    {
        return false;
    }
    int x_right_i = int(floor(x_right_init + 0.5));

    if (x_right_i-L < w || x_right_i+L+w >= im_left.cols)
    {
        return false;
    }
    cv::Mat IL = im_left.rowRange(y-w,y+w+1).colRange(x_left_i - w, x_left_i + w + 1);
    IL.convertTo(IL,CV_32F);
    IL = IL - IL.at<float>(w,w) *cv::Mat::ones(IL.rows,IL.cols,CV_32F);

    int bestDist = INT_MAX;
    int bestincR = 0;

    std::vector<float> vDists;
    vDists.resize(2*L+1);


    const float iniu = x_right_i + L - w;
    const float endu = x_right_i + L + w + 1;

    bool rv = true;
    if(iniu<0 || endu >= im_right.cols) {
//        std::cout << " case 2 " << std::endl;
        rv = false;
    }

    for(int incR=-L; incR<=+L; incR++)
    {
        cv::Mat IR = im_right.rowRange(y-w,y+w+1).colRange(x_right_i + incR - w,x_right_i + incR + w + 1);
        IR.convertTo(IR,CV_32F);
        IR = IR - IR.at<float>(w,w) *cv::Mat::ones(IR.rows,IR.cols,CV_32F);

        float dist = cv::norm(IL,IR,cv::NORM_L1);
        if(dist<bestDist)
        {
            bestDist =  dist;
            bestincR = incR;
        }

        vDists[L+incR] = dist;
    }

    if(bestincR==-L || bestincR==L) {
//        std::cout << " case 1 " << std::endl;
        rv = false;
    }

    const float dist1 = vDists[L+bestincR-1];
    const float dist2 = vDists[L+bestincR];
    const float dist3 = vDists[L+bestincR+1];

    const float deltaR = (dist1-dist3)/(2.0f*(dist1+dist3-2.0f*dist2));

    if(deltaR<-1 || deltaR>1) {
//        std::cout << " case 3 " << std::endl;
        rv = false;
    }

    if (!rv) {
//        std::cout << " costs:" << std::endl;
//        for (int i = 0; i < vDists.size(); i++)
//        {
//            std::cout << vDists[i] << " ";
//        }
//        std::cout << std::endl;
//        cv::Mat iml_d, imr_d;
//        cv::cvtColor(im_left, iml_d, CV_GRAY2BGR);
//        cv::cvtColor(im_right, imr_d, CV_GRAY2BGR);
//        cv::rectangle(iml_d, cv::Point(x_left - w, y - w), cv::Point(x_left + w + 1, y + w + 1), cv::Scalar(0, 0, 255));
//        cv::rectangle(imr_d, cv::Point(x_right - w - L, y - w), cv::Point(x_right + L + w + 1, y + w + 1),
//                      cv::Scalar(0, 255, 0));
//        cv::imshow("iml", iml_d);
//        cv::imshow("imr", imr_d);
//        cv::waitKey(0);
    }

    double dx = x_left - x_left_i;
    double x_for_i = x_right_i + deltaR + bestincR;
    *x_right_fin = x_for_i + dx;
    return rv;
}


//works only for integer endpoint coordinates
//bool vgl::RefineLineStereo(const KeyLine& kl1, const KeyLine& kl2, const cv::Mat& frame1, const cv::Mat& frame2, KeyLine* kl1_ref_p, KeyLine* kl2_ref_p)
//{
//    int x1_s, y1_s, x1_e, y1_e;
//    GetYSortedEndPoints(kl1, &x1_s, &x1_e, &y1_s, &y1_e);
//    int x2_s, y2_s, x2_e, y2_e;
//    GetYSortedEndPoints(kl2, &x2_s, &x2_e, &y2_s, &y2_e);
//
//    int y_l = std::max(y1_s, y2_s);
//    int y_h = std::min(y1_e, y2_e);
//
//    if (y_l > y_h)
//    {
//        return false;
//    }
//
//    if (y_h-y_l < 10)
//    {
//        return false;
//    }
//
//    const int w = 5;
//    const int w_half = (w-1)/2;
//    const int L = 3;
//
//    int y_ref_start = y_l + w_half;
//    int y_ref_end = y_h-w_half;
//
//    double x1_ref_start = GetXByY(x1_s, x1_e, y1_s, y1_e, y_ref_start);
//    double x1_ref_end = GetXByY(x1_s, x1_e, y1_s, y1_e, y_ref_end);
//    double x2_ref_start = GetXByY(x2_s, x2_e, y2_s, y2_e, y_ref_start);
//    double x2_ref_end = GetXByY(x2_s, x2_e, y2_s, y2_e, y_ref_end);
//
//
//
//    double x2_start, x2_end;
//    bool rv_start = RefinePoint(x1_ref_start, x2_ref_start, y_ref_start, frame1, frame2, w, L, &x2_start);
//    bool rv_end = RefinePoint(x1_ref_end, x2_ref_end, y_ref_end, frame1, frame2, w, L, &x2_end);
//    KeyLine kl_f = kl2;
//    kl_f.startPointY = y_ref_start;
//    kl_f.endPointY = y_ref_end;
////    kl_f.startPointX = x2_ref_start;//x2_start;
////    kl_f.endPointX = x2_ref_end;//x2_end;
//    kl_f.startPointX = x2_start;
//    kl_f.endPointX = x2_end;
//
//    KeyLine kl_f_1 = kl1;
//    kl_f_1.startPointY = y_ref_start;
//    kl_f_1.endPointY = y_ref_end;
//    kl_f_1.startPointX = x1_ref_start;
//    kl_f_1.endPointX = x1_ref_end;
//    *kl2_ref_p = kl_f;
//    *kl1_ref_p = kl_f_1;
//
////    return true;
//
//    return rv_start && rv_end;
//}


double vgl::LineEptReprojError(const Eigen::Vector3d& l, const Eigen::Matrix<double, 4, 4>& T, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, double f)
{
    Eigen::Vector3d p1c = T.block<3,3>(0,0).transpose() * (X1 - T.block<3,1>(0,3));
    Eigen::Vector3d p2c = T.block<3,3>(0,0).transpose() * (X2 - T.block<3,1>(0,3));
    return f * (fabs(1.0/p1c(2) * l.dot(p1c)) + fabs(1.0/p2c(2) * l.dot(p2c)));
}



double vgl::LineReprojErrorL1(const Eigen::Vector2d& xs, const Eigen::Vector2d& xe, const Eigen::Matrix<double, 4, 4>& T_other, const Eigen::Vector3d& X0,
                       const Eigen::Vector3d& line_dir, const Eigen::Matrix3d& K)
{
    Eigen::Vector3d Xc1 = K*MapPoint(T_other, X0);
    Eigen::Vector3d Xc2 = K*MapPoint(T_other, X0 + line_dir);
    Eigen::Vector3d leq3 = Xc1.cross(Xc2);
    leq3 = leq3 / leq3.segment<2>(0).norm();
    double err1 = fabs(xs.dot(leq3.segment<2>(0))+leq3(2));
    double err2 = fabs(xe.dot(leq3.segment<2>(0))+leq3(2));
    double se = err1 + err2;
    return se;
}

void vgl::ReprojectEndpointTo3D(const Eigen::Vector3d& endpoint, const Eigen::Vector3d& X0, const Eigen::Vector3d& d_rot, double* p)
{
    Eigen::Matrix<double, 3, 2> M;
    M.col(1) = - d_rot;
    M.col(0) = endpoint;
    Eigen::Matrix<double, 2, 1> sol = M.colPivHouseholderQr().solve(X0);
    *p = sol(1);
}

void vgl::LineFrom3DEndpoints(const Eigen::Vector3d& X1, const Eigen::Vector3d& X2, Eigen::Vector3d* X0, Eigen::Vector3d* line_dir)
{
    *line_dir = X2-X1;
    *line_dir = *line_dir / (line_dir->norm());
    double a = (X1.dot(*line_dir));
    *X0 = X1 - a*(*line_dir);
}

void vgl::NormalizedLineEquation(double sx, double sy, double ex, double ey, const Eigen::Matrix3d& K, Eigen::Vector3d* lineEq)
{
    Eigen::Vector3d Xs(sx, sy, 1.0);
    Eigen::Vector3d Xe(ex, ey, 1.0);
    Eigen::Vector3d lineImg = Xs.cross(Xe);
    *lineEq = K.transpose() * lineImg;
    *lineEq = *lineEq / lineEq->segment<2>(0).norm();
}

Eigen::Vector3d vgl::MapPoint(const Eigen::Matrix4d& T_c2w, const Eigen::Vector3d& X)
{
    return T_c2w.block<3,3>(0,0).transpose() * (X - T_c2w.block<3,1>(0,3));
}

Eigen::Vector2d vgl::ProjectPoint(const Eigen::Matrix3d& K, Eigen::Vector3d& X)
{
    Eigen::Vector3d x_h = K*X;
    return x_h.segment<2>(0)/x_h(2);
}