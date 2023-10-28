// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include "pmp/types.h"
#include "pmp/surface_mesh.h"

namespace pmp {

//! \addtogroup algorithms
//! @{

//! clamp cotangent values as if angles are in [1, 179]
inline double clamp_cot(const double v)
{
    const double bound = 19.1; // 3 degrees
    return (v < -bound ? -bound : (v > bound ? bound : v));
}

//! clamp cosine values as if angles are in [1, 179]
inline double clamp_cos(const double v)
{
    const double bound = 0.9986; // 3 degrees
    return (v < -bound ? -bound : (v > bound ? bound : v));
}

//! compute area of triangle f
Scalar triangle_area(const SurfaceMesh& mesh, Face f);

//! compute barycentric Voronoi area of vertex v
double voronoi_area_barycentric(const SurfaceMesh& mesh, Vertex v);

//! compute the sum of angles around vertex v (used for Gaussian curvature)
Scalar angle_sum(const SurfaceMesh& mesh, Vertex v);

//! discrete curvature information for a vertex. used for vertex_curvature()
struct VertexCurvature
{
    VertexCurvature() : mean(0.0), gauss(0.0), max(0.0), min(0.0) {}

    Scalar mean;
    Scalar gauss;
    Scalar max;
    Scalar min;
};

//! compute min, max, mean, and Gaussian curvature for vertex v. this will not
//! give reliable values for boundary vertices.
VertexCurvature vertex_curvature(const SurfaceMesh& mesh, Vertex v);

//! @}

} // namespace pmp
