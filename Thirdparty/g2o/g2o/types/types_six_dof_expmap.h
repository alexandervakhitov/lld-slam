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

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Modified by Alexander Vakhitov (2018)
// Added EdgeSE3ProjectLine, EdgeSE3ProjectLineOnlyPose

#ifndef G2O_SIX_DOF_TYPES_EXPMAP
#define G2O_SIX_DOF_TYPES_EXPMAP

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>
#include "../../../../include/vgl.h"

namespace g2o {
namespace types_six_dof_expmap {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;


/**
 * \brief SE3 Vertex parameterized internally with a transformation matrix
 and externally with its exponential map
 */
class  VertexSE3Expmap : public BaseVertex<6, SE3Quat>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexSE3Expmap();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate = SE3Quat();
  }

  virtual void oplusImpl(const double* update_)  {
    Eigen::Map<const Vector6d> update(update_);
    setEstimate(SE3Quat::exp(update)*estimate());
  }
};


class  EdgeSE3ProjectXYZ: public  BaseBinaryEdge<2, Vector2d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(v2->estimate()));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }
    

  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZ: public  BaseBinaryEdge<3, Vector3d, VertexSBAPointXYZ, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZ();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(v2->estimate()),bf);
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[1]);
    const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz, const float &bf) const;

  double fx, fy, cx, cy, bf;
};

class  EdgeSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<2, Vector2d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector2d obs(_measurement);
    _error = obs-cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy;
};


class  EdgeStereoSE3ProjectXYZOnlyPose: public  BaseUnaryEdge<3, Vector3d, VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeStereoSE3ProjectXYZOnlyPose(){}

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  void computeError()  {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    Vector3d obs(_measurement);
    _error = obs - cam_project(v1->estimate().map(Xw));
  }

  bool isDepthPositive() {
    const VertexSE3Expmap* v1 = static_cast<const VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }


  virtual void linearizeOplus();

  Vector3d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy, bf;
};

class EdgeSE3ProjectLineICRABin: public BaseBinaryEdge<1, double, VertexSBAPointXYZ, VertexSE3Expmap> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectLineICRABin() : BaseBinaryEdge<1, double, VertexSBAPointXYZ, VertexSE3Expmap>() {
    }

    bool read(std::istream &) { return true;};

    bool write(std::ostream &) const { return true; } ;

    bool IsDepthPositive()
    {
        const VertexSBAPointXYZ* v0 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
        const VertexSE3Expmap *v2 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        Eigen::Vector3d Xm = v2->estimate().map(v0->estimate())+b;
        return (Xm(2) > 0);
    }

//    Eigen::Vector2d  cam_project(const Eigen::Vector3d &trans_xyz) const;

    void linearize(JacobianXiOplusType& _jacobianOplusXi, JacobianXjOplusType& _jacobianOplusXj);

    void computeError()
    {
        const VertexSBAPointXYZ* v0 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
        const VertexSE3Expmap *v2 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;
        Eigen::Vector3d Xm = K*(v2->estimate().map(v0->estimate()) + b);
        Eigen::Vector2d xmp = Xm.segment<2>(0)/Xm(2);
        Eigen::Vector3d l_proj = x1.cross(x2);
        l_proj = l_proj / l_proj.segment<2>(0).norm();
        _error(0) = l_proj.segment<2>(0).transpose() * xmp + l_proj(2);
    }

    virtual void linearizeOplus();

    Eigen::Vector3d x1;
    Eigen::Vector3d x2;
//
    Eigen::Vector3d b;

    double f, cx, cy;
};


class EdgeSE3ProjectLineICRABinOnlyPose: public BaseUnaryEdge<1, double, VertexSE3Expmap> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectLineICRABinOnlyPose() : BaseUnaryEdge<1, double, VertexSE3Expmap>() {
    }

    bool read(std::istream &) { return true;};

    bool write(std::ostream &) const { return true; } ;

    void linearize(JacobianXiOplusType& _jacobianOplusXi);

    void computeError()
    {
        const VertexSE3Expmap *v2 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;
        Eigen::Vector3d Xm = K*(v2->estimate().map(X) + b);
        Eigen::Vector2d xmp = Xm.segment<2>(0)/Xm(2);
        Eigen::Vector3d l_proj = x1.cross(x2);
        l_proj = l_proj / l_proj.segment<2>(0).norm();
        _error(0) = l_proj.segment<2>(0).transpose() * xmp + l_proj(2);
    }

    virtual void linearizeOplus();

    Eigen::Vector3d X;
    Eigen::Vector3d x1;
    Eigen::Vector3d x2;
//
    Eigen::Vector3d b;

    double f, cx, cy;
};


class EdgeSE3ProjectLine: public BaseBinaryEdge<2, Eigen::Vector2d, VertexSBALine, VertexSE3Expmap> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectLine() : BaseBinaryEdge<2, Vector2d, VertexSBALine, VertexSE3Expmap>()
    {}

    bool read(std::istream &) { return true;};

    bool write(std::ostream &) const { return true; } ;

    void linearize(JacobianXiOplusType& _jacobianOplusXi, JacobianXjOplusType& _jacobianOplusXj);

    bool IsDepthPositive()
    {
        const VertexSBALine *v0 = static_cast<const VertexSBALine *>(_vertices[0]);
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        const Eigen::Matrix3d R = v0->estimate().GetR();
        Eigen::Vector3d  X0 = R.col(1) * v0->estimate().GetAlpha();
        Eigen::Vector3d  ld = R.col(0);
        Eigen::Vector3d X0l = v1->estimate().map(X0) + b;
        Eigen::Vector3d ldl = v1->estimate().map(X0+ld) + b - X0l;
        double depth1, depth2, line_param;
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;


        vgl::ReprojectLinePointTo3D(X0l, ldl, x1.segment<2>(0), K,
                                    &depth1, &line_param);

        vgl::ReprojectLinePointTo3D(X0l, ldl, x2.segment<2>(0), K,
                                    &depth2, &line_param);

        if (depth1<0 || depth2<0)
        {
            return false;
        } else {
            return true;
        }
    }

    void computeError()
    {
//        std::cout << " line error compute start" << std::endl;
//        std::cout << "class fields: " << f << " " << cx << " " << cy << std::endl;
//        _error(0) = f;
//        _error(1) = cx;
        const VertexSBALine *v0 = static_cast<const VertexSBALine *>(_vertices[0]);
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[1]);
        Eigen::Matrix3d R = v0->estimate().GetR();
        Eigen::Vector3d  X1 = R.col(1) * v0->estimate().GetAlpha();
        Eigen::Vector3d  X2 = X1 + R.col(0);

//        std::cout << "decode X0: " << X1.transpose() << std::endl;
//        std::cout << "decode ld: " << R.col(0).transpose() << std::endl;
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0, 0) = f;
        K(1, 1) = f;
        K(0, 2) = cx;
        K(1, 2) = cy;

        Eigen::Vector3d  X1m = K*(v1->estimate().map(X1) + b);
        Eigen::Vector3d  X2m = K*(v1->estimate().map(X2) + b);

        Eigen::Vector3d  l_tilde = X1m.cross(X2m);
        Eigen::Vector3d l = l_tilde / l_tilde.segment<2>(0).norm();
//        std::cout << " f " << f << " x1 " << x1 << " l " << l << std::endl;
        _error(0) = x1.transpose() * l;
        _error(1) = x2.transpose() * l;
////        _error(2) = 0.5*(_error(0)+_error(1));
//        std::cout << " line error compute done" << std::endl;
    }

    virtual void linearizeOplus();
//
    Eigen::Vector3d x1;
    Eigen::Vector3d x2;
//
    Eigen::Vector3d b;

    double f, cx, cy;
};


class EdgeSE3ProjectLineOnlyPose : public BaseUnaryEdge<2, Eigen::Vector2d, VertexSE3Expmap>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeSE3ProjectLineOnlyPose() {
        b.setZero();
    }

    bool read(std::istream &) { return true;};

    bool write(std::ostream &) const { return true; } ;

    void linearize(JacobianXiOplusType& _jacobianOplusXi);

    void computeError()
    {
        const VertexSE3Expmap *v1 = static_cast<const VertexSE3Expmap *>(_vertices[0]);
        Eigen::Matrix3d K;
        K.setIdentity();
        K(0,0) = f;
        K(1,1) = f;
        K(0,2) = cx;
        K(1,2) = cy;
        Eigen::Vector3d  X1m = K*(v1->estimate().map(X1) + b);
        Eigen::Vector3d  X2m = K*(v1->estimate().map(X2) + b);
        Eigen::Vector3d  l_tilde = X1m.cross(X2m);
        Eigen::Vector3d l = l_tilde / l_tilde.segment<2>(0).norm();
        _error(0) = x1.transpose() * l;
        _error(1) = x2.transpose() * l;
    }

    virtual void linearizeOplus();

    Eigen::Vector3d X1, X2;
    Eigen::Vector3d x1;
    Eigen::Vector3d x2;

    Eigen::Vector3d b;

    double f, cx, cy;
};

} // end namespace

#endif
