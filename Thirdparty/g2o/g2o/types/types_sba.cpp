// g2o - General Graph Optimization
// Copyright (C) 2011 Kurt Konolige
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

//Modified by Alexander Vakhitov (2018)

#include "types_sba.h"
#include <iostream>

namespace g2o {

    using namespace std;


    VertexSBAPointXYZ::VertexSBAPointXYZ() : BaseVertex<3, Vector3d>()
    {
    }

    bool VertexSBAPointXYZ::read(std::istream& is)
    {
        Vector3d lv;
        for (int i=0; i<3; i++)
            is >> _estimate[i];
        return true;
    }

    bool VertexSBAPointXYZ::write(std::ostream& os) const
    {
        Vector3d lv=estimate();
        for (int i=0; i<3; i++){
            os << lv[i] << " ";
        }
        return os.good();
    }

    LineParams::LineParams() : q(Eigen::Quaterniond::Identity()), alpha (0.0)
    {

    }

    LineParams::LineParams(const LineParams &line_params)
    {
        this->SetQ(line_params.GetQ());
        this->SetAlpha(line_params.GetAlpha());
    }

    LineParams::LineParams(const Eigen::Quaterniond &q, double alpha)  : q(q), alpha(alpha){

    }

    double LineParams::GetAlpha() const {
        return alpha;
    }

    Eigen::Quaterniond LineParams::GetQ() const {
        return q.normalized();
    }

    void LineParams::SetQ(const Eigen::Quaterniond &q) {
        this->q = q;
    }

    void LineParams::SetAlpha(double alpha) {
        this->alpha = alpha;
    }

    Eigen::Matrix3d LineParams::GetR() const {
        Eigen::Quaterniond qn = GetQ();
        return qn.toRotationMatrix();
    }

} // end namespace
