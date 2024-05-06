#pragma once
#include <Math/Vector4D.h>
#include <Math/GenVector/VectorUtil.h>

namespace WWTools
{
    using namespace ROOT::Math;
    template <class CoordSystem>
    LorentzVector<CoordSystem>
    starVector(const LorentzVector<CoordSystem>& W_lvec, const LorentzVector<CoordSystem>& lep_lvec, const XYZVector& e_minus_vec = XYZVector(0, 0, 1))
    {
        // z* points along the W flight direction in the lab frame
        const auto z = W_lvec.Vect().Unit();
        // y* points in the direction of e- x W
        const auto y = e_minus_vec.Cross(z).Unit();
        // the obvious third dimension
        const auto x = y.Cross(z);

        // .Inverse() to obtain passive rotation (rotate coord system not vector)
        const auto rot = Rotation3D(x, y, z).Inverse();

        // get the boost to boost something into the W CMS, boost and then rotate
        const auto W_boost = W_lvec.BoostToCM();
        const auto boosted_lvec = VectorUtil::boost(lep_lvec, W_boost);
        const auto transformed_lvec = rot(boosted_lvec);

        return transformed_lvec;
    }
}