// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "types_six_dof_expmap.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {

    using namespace std;

    Eigen::Matrix3d cpmat(const Eigen::Vector3d& t)
    {
        Eigen::Matrix3d t_hat;
        t_hat << 0, -t(2), t(1),
                t(2), 0, -t(0),
                -t(1), t(0), 0;
        return t_hat;
    }

    Vector2d project2d(const Vector3d& v)  {
        Vector2d res;
        res(0) = v(0)/v(2);
        res(1) = v(1)/v(2);
        return res;
    }

    Vector3d unproject2d(const Vector2d& v)  {
        Vector3d res;
        res(0) = v(0);
        res(1) = v(1);
        res(2) = 1;
        return res;
    }

    VertexSE3Expmap::VertexSE3Expmap() : BaseVertex<6, SE3Quat>() {
    }

    bool VertexSE3Expmap::read(std::istream& is) {
        Vector7d est;
        for (int i=0; i<7; i++)
            is  >> est[i];
        SE3Quat cam2world;
        cam2world.fromVector(est);
        setEstimate(cam2world.inverse());
        return true;
    }

    bool VertexSE3Expmap::write(std::ostream& os) const {
        SE3Quat cam2world(estimate().inverse());
        for (int i=0; i<7; i++)
            os << cam2world[i] << " ";
        return os.good();
    }


    EdgeSE3ProjectXYZ::EdgeSE3ProjectXYZ() : BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>() {
    }

    bool EdgeSE3ProjectXYZ::read(std::istream& is){
        for (int i=0; i<2; i++){
            is >> _measurement[i];
        }
        for (int i=0; i<2; i++)
            for (int j=i; j<2; j++) {
                is >> information()(i,j);
                if (i!=j)
                    information()(j,i)=information()(i,j);
            }
        return true;
    }

    bool EdgeSE3ProjectXYZ::write(std::ostream& os) const {

        for (int i=0; i<2; i++){
            os << measurement()[i] << " ";
        }

        for (int i=0; i<2; i++)
            for (int j=i; j<2; j++){
                os << " " <<  information()(i,j);
            }
        return os.good();
    }


    void EdgeSE3ProjectXYZ::linearizeOplus() {
        VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
        SE3Quat T(vj->estimate());
        VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
        Vector3d xyz = vi->estimate();
        Vector3d xyz_trans = T.map(xyz);

        double x = xyz_trans[0];
        double y = xyz_trans[1];
        double z = xyz_trans[2];
        double z_2 = z*z;

        Matrix<double,2,3> tmp;
        tmp(0,0) = fx;
        tmp(0,1) = 0;
        tmp(0,2) = -x/z*fx;

        tmp(1,0) = 0;
        tmp(1,1) = fy;
        tmp(1,2) = -y/z*fy;

        _jacobianOplusXi =  -1./z * tmp * T.rotation().toRotationMatrix();

        _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
        _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
        _jacobianOplusXj(0,2) = y/z *fx;
        _jacobianOplusXj(0,3) = -1./z *fx;
        _jacobianOplusXj(0,4) = 0;
        _jacobianOplusXj(0,5) = x/z_2 *fx;

        _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
        _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
        _jacobianOplusXj(1,2) = -x/z *fy;
        _jacobianOplusXj(1,3) = 0;
        _jacobianOplusXj(1,4) = -1./z *fy;
        _jacobianOplusXj(1,5) = y/z_2 *fy;
    }

    Vector2d EdgeSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz) const{
        Vector2d proj = project2d(trans_xyz);
        Vector2d res;
        res[0] = proj[0]*fx + cx;
        res[1] = proj[1]*fy + cy;
        return res;
    }


    Vector3d EdgeStereoSE3ProjectXYZ::cam_project(const Vector3d & trans_xyz, const float &bf) const{
        const float invz = 1.0f/trans_xyz[2];
        Vector3d res;
        res[0] = trans_xyz[0]*invz*fx + cx;
        res[1] = trans_xyz[1]*invz*fy + cy;
        res[2] = res[0] - bf*invz;
        return res;
    }

    EdgeStereoSE3ProjectXYZ::EdgeStereoSE3ProjectXYZ() : BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>() {
    }

    bool EdgeStereoSE3ProjectXYZ::read(std::istream& is){
        for (int i=0; i<=3; i++){
            is >> _measurement[i];
        }
        for (int i=0; i<=2; i++)
            for (int j=i; j<=2; j++) {
                is >> information()(i,j);
                if (i!=j)
                    information()(j,i)=information()(i,j);
            }
        return true;
    }

    bool EdgeStereoSE3ProjectXYZ::write(std::ostream& os) const {

        for (int i=0; i<=3; i++){
            os << measurement()[i] << " ";
        }

        for (int i=0; i<=2; i++)
            for (int j=i; j<=2; j++){
                os << " " <<  information()(i,j);
            }
        return os.good();
    }

    void EdgeStereoSE3ProjectXYZ::linearizeOplus() {
        VertexSE3Expmap * vj = static_cast<VertexSE3Expmap *>(_vertices[1]);
        SE3Quat T(vj->estimate());
        VertexSBAPointXYZ* vi = static_cast<VertexSBAPointXYZ*>(_vertices[0]);
        Vector3d xyz = vi->estimate();
        Vector3d xyz_trans = T.map(xyz);

        const Matrix3d R =  T.rotation().toRotationMatrix();

        double x = xyz_trans[0];
        double y = xyz_trans[1];
        double z = xyz_trans[2];
        double z_2 = z*z;

        _jacobianOplusXi(0,0) = -fx*R(0,0)/z+fx*x*R(2,0)/z_2;
        _jacobianOplusXi(0,1) = -fx*R(0,1)/z+fx*x*R(2,1)/z_2;
        _jacobianOplusXi(0,2) = -fx*R(0,2)/z+fx*x*R(2,2)/z_2;

        _jacobianOplusXi(1,0) = -fy*R(1,0)/z+fy*y*R(2,0)/z_2;
        _jacobianOplusXi(1,1) = -fy*R(1,1)/z+fy*y*R(2,1)/z_2;
        _jacobianOplusXi(1,2) = -fy*R(1,2)/z+fy*y*R(2,2)/z_2;

        _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*R(2,0)/z_2;
        _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)-bf*R(2,1)/z_2;
        _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2)-bf*R(2,2)/z_2;

        _jacobianOplusXj(0,0) =  x*y/z_2 *fx;
        _jacobianOplusXj(0,1) = -(1+(x*x/z_2)) *fx;
        _jacobianOplusXj(0,2) = y/z *fx;
        _jacobianOplusXj(0,3) = -1./z *fx;
        _jacobianOplusXj(0,4) = 0;
        _jacobianOplusXj(0,5) = x/z_2 *fx;

        _jacobianOplusXj(1,0) = (1+y*y/z_2) *fy;
        _jacobianOplusXj(1,1) = -x*y/z_2 *fy;
        _jacobianOplusXj(1,2) = -x/z *fy;
        _jacobianOplusXj(1,3) = 0;
        _jacobianOplusXj(1,4) = -1./z *fy;
        _jacobianOplusXj(1,5) = y/z_2 *fy;

        _jacobianOplusXj(2,0) = _jacobianOplusXj(0,0)-bf*y/z_2;
        _jacobianOplusXj(2,1) = _jacobianOplusXj(0,1)+bf*x/z_2;
        _jacobianOplusXj(2,2) = _jacobianOplusXj(0,2);
        _jacobianOplusXj(2,3) = _jacobianOplusXj(0,3);
        _jacobianOplusXj(2,4) = 0;
        _jacobianOplusXj(2,5) = _jacobianOplusXj(0,5)-bf/z_2;
    }


//Only Pose

    bool EdgeSE3ProjectXYZOnlyPose::read(std::istream& is){
        for (int i=0; i<2; i++){
            is >> _measurement[i];
        }
        for (int i=0; i<2; i++)
            for (int j=i; j<2; j++) {
                is >> information()(i,j);
                if (i!=j)
                    information()(j,i)=information()(i,j);
            }
        return true;
    }

    bool EdgeSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

        for (int i=0; i<2; i++){
            os << measurement()[i] << " ";
        }

        for (int i=0; i<2; i++)
            for (int j=i; j<2; j++){
                os << " " <<  information()(i,j);
            }
        return os.good();
    }


    void EdgeSE3ProjectXYZOnlyPose::linearizeOplus() {
        VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
        Vector3d xyz_trans = vi->estimate().map(Xw);

        double x = xyz_trans[0];
        double y = xyz_trans[1];
        double invz = 1.0/xyz_trans[2];
        double invz_2 = invz*invz;

        _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
        _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
        _jacobianOplusXi(0,2) = y*invz *fx;
        _jacobianOplusXi(0,3) = -invz *fx;
        _jacobianOplusXi(0,4) = 0;
        _jacobianOplusXi(0,5) = x*invz_2 *fx;

        _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
        _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
        _jacobianOplusXi(1,2) = -x*invz *fy;
        _jacobianOplusXi(1,3) = 0;
        _jacobianOplusXi(1,4) = -invz *fy;
        _jacobianOplusXi(1,5) = y*invz_2 *fy;
    }

    Vector2d EdgeSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
        Vector2d proj = project2d(trans_xyz);
        Vector2d res;
        res[0] = proj[0]*fx + cx;
        res[1] = proj[1]*fy + cy;
        return res;
    }


    Vector3d EdgeStereoSE3ProjectXYZOnlyPose::cam_project(const Vector3d & trans_xyz) const{
        const float invz = 1.0f/trans_xyz[2];
        Vector3d res;
        res[0] = trans_xyz[0]*invz*fx + cx;
        res[1] = trans_xyz[1]*invz*fy + cy;
        res[2] = res[0] - bf*invz;
        return res;
    }


    bool EdgeStereoSE3ProjectXYZOnlyPose::read(std::istream& is){
        for (int i=0; i<=3; i++){
            is >> _measurement[i];
        }
        for (int i=0; i<=2; i++)
            for (int j=i; j<=2; j++) {
                is >> information()(i,j);
                if (i!=j)
                    information()(j,i)=information()(i,j);
            }
        return true;
    }

    bool EdgeStereoSE3ProjectXYZOnlyPose::write(std::ostream& os) const {

        for (int i=0; i<=3; i++){
            os << measurement()[i] << " ";
        }

        for (int i=0; i<=2; i++)
            for (int j=i; j<=2; j++){
                os << " " <<  information()(i,j);
            }
        return os.good();
    }

    void EdgeStereoSE3ProjectXYZOnlyPose::linearizeOplus() {
        VertexSE3Expmap * vi = static_cast<VertexSE3Expmap *>(_vertices[0]);
        Vector3d xyz_trans = vi->estimate().map(Xw);

        double x = xyz_trans[0];
        double y = xyz_trans[1];
        double invz = 1.0/xyz_trans[2];
        double invz_2 = invz*invz;

        _jacobianOplusXi(0,0) =  x*y*invz_2 *fx;
        _jacobianOplusXi(0,1) = -(1+(x*x*invz_2)) *fx;
        _jacobianOplusXi(0,2) = y*invz *fx;
        _jacobianOplusXi(0,3) = -invz *fx;
        _jacobianOplusXi(0,4) = 0;
        _jacobianOplusXi(0,5) = x*invz_2 *fx;

        _jacobianOplusXi(1,0) = (1+y*y*invz_2) *fy;
        _jacobianOplusXi(1,1) = -x*y*invz_2 *fy;
        _jacobianOplusXi(1,2) = -x*invz *fy;
        _jacobianOplusXi(1,3) = 0;
        _jacobianOplusXi(1,4) = -invz *fy;
        _jacobianOplusXi(1,5) = y*invz_2 *fy;

        _jacobianOplusXi(2,0) = _jacobianOplusXi(0,0)-bf*y*invz_2;
        _jacobianOplusXi(2,1) = _jacobianOplusXi(0,1)+bf*x*invz_2;
        _jacobianOplusXi(2,2) = _jacobianOplusXi(0,2);
        _jacobianOplusXi(2,3) = _jacobianOplusXi(0,3);
        _jacobianOplusXi(2,4) = 0;
        _jacobianOplusXi(2,5) = _jacobianOplusXi(0,5)-bf*invz_2;
    }

    void EdgeSE3ProjectLineICRABin::linearize(JacobianXiOplusType& _jacobianOplusXi, JacobianXjOplusType& _jacobianOplusXj)
    {
        const VertexSBAPointXYZ* v0 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
        const VertexSE3Expmap *v2 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;

        SE3Quat T(v2->estimate());
        const Matrix3d R = T.rotation().toRotationMatrix();
        Vector3d xyz_trans = T.map(v0->estimate());
        Vector3d xyz_trans_b = xyz_trans + b;
        Eigen::Vector3d xp = K*xyz_trans_b;

        Eigen::Matrix<double, 3, 6> J_l_2;
        J_l_2.block<3,3>(0,0) = - K*cpmat(xyz_trans);
        J_l_2.block<3,3>(0,3) = K*Eigen::Matrix3d::Identity();

        Eigen::Matrix<double, 2, 3> J_proj;
        J_proj.setZero();
        J_proj(0,0) = 1.0/xp(2);
        J_proj(0,2) = -1.0/(xp(2)*xp(2))*xp(0);
        J_proj(1,1) = 1.0/xp(2);
        J_proj(1,2) = -1.0/(xp(2)*xp(2))*xp(1);

        Eigen::Matrix<double, 2, 6> J2 = J_proj * J_l_2;

        Eigen::Vector3d l_proj = x1.cross(x2);
        l_proj = l_proj / l_proj.segment<2>(0).norm();

        _jacobianOplusXj = l_proj.segment<2>(0).transpose() * J2;

        double x = xyz_trans_b[0];
        double y = xyz_trans_b[1];
        double z = xyz_trans_b[2];
        double z_2 = z * z;

        MatrixXd J;
        J.resize(2, 3);
        J(0, 0) = -f * R(0, 0) / z + f * x * R(2, 0) / z_2;
        J(0, 1) = -f * R(0, 1) / z + f * x * R(2, 1) / z_2;
        J(0, 2) = -f * R(0, 2) / z + f * x * R(2, 2) / z_2;

        J(1, 0) = -f * R(1, 0) / z + f * y * R(2, 0) / z_2;
        J(1, 1) = -f * R(1, 1) / z + f * y * R(2, 1) / z_2;
        J(1, 2) = -f * R(1, 2) / z + f * y * R(2, 2) / z_2;

        _jacobianOplusXi =  -l_proj.segment<2>(0).transpose() * J;
    }

    void EdgeSE3ProjectLineICRABin::linearizeOplus()
    {
        linearize(_jacobianOplusXi, _jacobianOplusXj);
    }

    void EdgeSE3ProjectLineICRABinOnlyPose::linearize(JacobianXiOplusType& _jacobianOplusXi)
    {
        const VertexSE3Expmap *v2 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;

        SE3Quat T(v2->estimate());
        Vector3d xyz_trans = T.map(X);
        Vector3d xyz_trans_b = xyz_trans + b;
        Eigen::Vector3d xp = K*xyz_trans_b;

        Eigen::Matrix<double, 3, 6> J_l_2;
        J_l_2.block<3,3>(0,0) = - K*cpmat(xyz_trans);
        J_l_2.block<3,3>(0,3) = K*Eigen::Matrix3d::Identity();

        Eigen::Matrix<double, 2, 3> J_proj;
        J_proj.setZero();
        J_proj(0,0) = 1.0/xp(2);
        J_proj(0,2) = -1.0/(xp(2)*xp(2))*xp(0);
        J_proj(1,1) = 1.0/xp(2);
        J_proj(1,2) = -1.0/(xp(2)*xp(2))*xp(1);

        Eigen::Matrix<double, 2, 6> J2 = J_proj * J_l_2;

        Eigen::Vector3d l_proj = x1.cross(x2);
        l_proj = l_proj / l_proj.segment<2>(0).norm();

        _jacobianOplusXi = l_proj.segment<2>(0).transpose() * J2;
    }

    void EdgeSE3ProjectLineICRABinOnlyPose::linearizeOplus()
    {
        linearize(_jacobianOplusXi);
    }


    void FormJacobianLineWRTCam(const Eigen::Vector3d& X1m, const Eigen::Vector3d& X2m, const Eigen::Vector3d& b, const Eigen::Matrix3d& K, Eigen::Matrix<double, 3, 6>* J_lp, Eigen::Matrix3d* D_l_ltilde_p)
    {
        Eigen::Matrix3d D_l_ltilde;

        Eigen::Vector3d l_tilde = (K*(X1m+b)).cross(K*(X2m+b));
        double n = l_tilde.segment<2>(0).norm();

        Eigen::Matrix<double, 1, 3> dn;
        dn(0) = -l_tilde(0) / (n*n*n);
        dn(1) = -l_tilde(1) / (n*n*n);
        dn(2) = 0;
        D_l_ltilde = l_tilde * dn + 1.0/n * Eigen::Matrix3d::Identity();
        *D_l_ltilde_p = D_l_ltilde;

        Eigen::Matrix<double, 3, 6> J_l;

        Eigen::Matrix<double, 3, 6> J_l_2;
        J_l_2.block<3,3>(0,0) = - K*cpmat(X2m);
        J_l_2.block<3,3>(0,3) = K*Eigen::Matrix3d::Identity();
        Eigen::Matrix<double, 3, 6> J_l_1;
        J_l_1.block<3,3>(0,0) = - K*cpmat(X1m);
        J_l_1.block<3,3>(0,3) = K*Eigen::Matrix3d::Identity();

        J_l = cpmat(K*(X1m+b)) * J_l_2 - cpmat(K*(X2m+b)) * J_l_1;

        J_l = D_l_ltilde * J_l;
        *J_lp = J_l;
    }


    void EdgeSE3ProjectLine::linearizeOplus()
    {
        linearize(_jacobianOplusXi, _jacobianOplusXj);
    }
//
    void EdgeSE3ProjectLine::linearize(JacobianXiOplusType& _jacobianOplusXi, JacobianXjOplusType& _jacobianOplusXj)
    {
        const VertexSBALine *v0 = static_cast<const VertexSBALine *>(_vertices[0]);
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        Eigen::Matrix3d R = v0->estimate().GetR();
        Eigen::Vector3d  X1 = R.col(1) * v0->estimate().GetAlpha();
        Eigen::Vector3d  X2 = X1 + R.col(0);
        Eigen::Vector3d  X1m = v1->estimate().map(X1);
        Eigen::Vector3d  X2m = v1->estimate().map(X2);

        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;

        Eigen::Matrix<double, 3, 6> J_l;
        Eigen::Matrix3d D_l_lt;
        FormJacobianLineWRTCam(X1m, X2m, b, K, &J_l, &D_l_lt);

        Eigen::Matrix<double, 1, 6> r1 = x1.transpose() * J_l;
        Eigen::Matrix<double, 1, 6> r2 = x2.transpose() * J_l;

        _jacobianOplusXj.row(0) = r1;
        _jacobianOplusXj.row(1) = r2;

        Eigen::Matrix3d dX1dr = -cpmat(v0->estimate().GetAlpha() * R.col(1));
        Eigen::Vector3d dX1da = R.col(1);
        Eigen::Matrix<double, 3, 4> dX1;
        dX1.setZero();
        dX1.block<3,3>(0, 0) = 2*dX1dr;
        dX1.block<3,1>(0, 3) = dX1da;

        Eigen::Matrix<double, 3, 4> dX2 = dX1;
        dX2.block<3,3>(0,0) = dX2.block<3,3>(0,0) - 2*cpmat(R.col(0));

        Eigen::Matrix4d T_se3 = v1->estimate().to_homogeneous_matrix();
        Eigen::Matrix3d R_cam = T_se3.block<3,3>(0,0);
        Eigen::Matrix<double, 3, 4> d_l_tilde = cpmat(K*(X1m+b)) * K*R_cam * dX2 - cpmat(K*(X2m+b))*K*R_cam*dX1;
        Eigen::Matrix<double, 3, 4> D_l = D_l_lt * d_l_tilde;

        Eigen::Matrix<double, 1, 4> r1l = x1.transpose() * D_l;
        Eigen::Matrix<double, 1, 4> r2l = x2.transpose() * D_l;

        _jacobianOplusXi.row(0) = r1l;
        _jacobianOplusXi.row(1) = r2l;

        for (int i = 0; i < _jacobianOplusXj.rows(); i++)
        {
            for (int j = 0; j < _jacobianOplusXj.cols(); j++)
            {
                if (isnanf(_jacobianOplusXj(i, j)))
                {
                    std::cout << " EdgeSBALine: nan in Xj " << std::endl;
                }
            }
        }

        for (int i = 0; i < _jacobianOplusXi.rows(); i++)
        {
            for (int j = 0; j < _jacobianOplusXi.cols(); j++)
            {
                if (isnanf(_jacobianOplusXi(i, j)))
                {
                    std::cout << " EdgeSBALine: nan in Xi " << std::endl;
                }
            }
        }
//  for (int i = 0; i < 4; i++)
//  {
//    _jacobianOplusXi(0, i) = f*r1l(0, i);
//    _jacobianOplusXi(1, i) = f*r2l(0, i);
//  }
    }

    void EdgeSE3ProjectLineOnlyPose::linearize(JacobianXiOplusType& _jacobianOplusXi)
    {

        Eigen::Matrix<double, 3, 6> J_l;
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
        Eigen::Vector3d X1m = v1->estimate().map(X1);
        Eigen::Vector3d X2m = v1->estimate().map(X2);

        Eigen::Matrix3d K;
        K.setIdentity();
        K(0,0) = f;
        K(1,1) = f;
        K(0,2) = cx;
        K(1,2) = cy;
        Eigen::Matrix3d D_l_lt;
        FormJacobianLineWRTCam(X1m, X2m, b, K, &J_l, &D_l_lt);

        Eigen::Matrix<double, 1, 6> r1 = x1.transpose() * J_l;
        Eigen::Matrix<double, 1, 6> r2 = x2.transpose() * J_l;

        _jacobianOplusXi.row(0) = r1;
//    for (int i = 0; i < 6; i++)
//    {
//        _jacobianOplusXi(0, i) = f*r1(0, i);
//    }
        _jacobianOplusXi.row(1) = r2;
//    for (int i = 0; i < 6; i++)
//    {
//        _jacobianOplusXi(1, i) = f*r2(0, i);
//    }
    }

    void EdgeSE3ProjectLineOnlyPose::linearizeOplus()
    {
        linearize(_jacobianOplusXi);
    }


} // end namespace
