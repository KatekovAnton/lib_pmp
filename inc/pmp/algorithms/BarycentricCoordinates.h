// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/mat_vec.h"

namespace pmp {

template <typename Scalar>
const Vector<Scalar, 3> barycentric_coordinates(const Vector<Scalar, 3>& p,
                                                const Vector<Scalar, 3>& u,
                                                const Vector<Scalar, 3>& v,
                                                const Vector<Scalar, 3>& w)
{
    Vector<Scalar, 3> result(Scalar(1.0 / 3.0)); // default: barycenter

    Vector<Scalar, 3> vu = v - u, wu = w - u, pu = p - u;

    // find largest absolute coodinate of normal
    Scalar nx = vu[1] * wu[2] - vu[2] * wu[1],
           ny = vu[2] * wu[0] - vu[0] * wu[2],
           nz = vu[0] * wu[1] - vu[1] * wu[0], ax = std::fabs(nx), ay = std::fabs(ny),
           az = std::fabs(nz);

    unsigned char maxCoord;

    if (ax > ay)
    {
        if (ax > az)
        {
            maxCoord = 0;
        }
        else
        {
            maxCoord = 2;
        }
    }
    else
    {
        if (ay > az)
        {
            maxCoord = 1;
        }
        else
        {
            maxCoord = 2;
        }
    }

    // solve 2D problem
    switch (maxCoord)
    {
        case 0:
        {
            if (Scalar(1.0) + ax != Scalar(1.0))
            {
                result[1] = Scalar(1.0) + (pu[1] * wu[2] - pu[2] * wu[1]) / nx - Scalar(1.0);
                result[2] = Scalar(1.0) + (vu[1] * pu[2] - vu[2] * pu[1]) / nx - Scalar(1.0);
                result[0] = Scalar(1.0) - result[1] - result[2];
            }
            break;
        }

        case 1:
        {
            if (Scalar(1.0) + ay != Scalar(1.0))
            {
                result[1] = Scalar(1.0) + (pu[2] * wu[0] - pu[0] * wu[2]) / ny - Scalar(1.0);
                result[2] = Scalar(1.0) + (vu[2] * pu[0] - vu[0] * pu[2]) / ny - Scalar(1.0);
                result[0] = Scalar(1.0) - result[1] - result[2];
            }
            break;
        }

        case 2:
        {
            if (Scalar(1.0) + az != Scalar(1.0))
            {
                result[1] = Scalar(1.0) + (pu[0] * wu[1] - pu[1] * wu[0]) / nz - Scalar(1.0);
                result[2] = Scalar(1.0) + (vu[0] * pu[1] - vu[1] * pu[0]) / nz - Scalar(1.0);
                result[0] = Scalar(1.0) - result[1] - result[2];
            }
            break;
        }
    }

    return result;
}

} // namespace pmp
