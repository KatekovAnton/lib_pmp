// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/DifferentialGeometry.h"
#include "pmp/algorithms/differential_geometry.h"

#include <cmath>

#include <limits>

#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

namespace pmp {

Scalar triangle_area(const SurfaceMesh& mesh, Face f)
{
    assert(mesh.valence(f) == 3);

    auto fv = mesh.vertices(f);
    const auto& p0 = mesh.position(*fv);
    const auto& p1 = mesh.position(*(++fv));
    const auto& p2 = mesh.position(*(++fv));

    return triangle_area(p0, p1, p2);
}

double voronoi_area_barycentric(const SurfaceMesh& mesh, Vertex v)
{
    double area(0.0);

    if (!mesh.is_isolated(v))
    {
        const Point p = mesh.position(v);
        Halfedge h0, h1;
        Point q, r, pq, pr;

        for (auto h : mesh.halfedges(v))
        {
            if (mesh.is_boundary(h))
                continue;

            h0 = h;
            h1 = mesh.next_halfedge(h0);

            pq = mesh.position(mesh.to_vertex(h0));
            pq -= p;

            pr = mesh.position(mesh.to_vertex(h1));
            pr -= p;

            area += double(norm(cross(pq, pr)) / Scalar(3.0));
        }
    }

    return area;
}

Scalar angle_sum(const SurfaceMesh& mesh, Vertex v)
{
    Scalar angles(0.0);

    if (!mesh.is_boundary(v))
    {
        const Point& p0 = mesh.position(v);

        for (auto h : mesh.halfedges(v))
        {
            const Point& p1 = mesh.position(mesh.to_vertex(h));
            const Point& p2 =
                mesh.position(mesh.to_vertex(mesh.ccw_rotated_halfedge(h)));

            const Point p01 = normalize(p1 - p0);
            const Point p02 = normalize(p2 - p0);

            Scalar cos_angle = Scalar(clamp_cos(double(dot(p01, p02))));

            angles += std::acos(cos_angle);
        }
    }

    return angles;
}

VertexCurvature vertex_curvature(const SurfaceMesh& mesh, Vertex v)
{
    VertexCurvature c;

    const Scalar area = Scalar(voronoi_area(mesh, v));
    if (area > std::numeric_limits<Scalar>::min())
    {
        c.mean = Scalar(0.5) * norm(laplace(mesh, v));
        c.gauss = (Scalar(2.0) * Scalar(M_PI) - angle_sum(mesh, v)) / area;

        const Scalar s = sqrt(std::max(Scalar(0.0), c.mean * c.mean - c.gauss));
        c.min = c.mean - s;
        c.max = c.mean + s;

        assert(!std::isnan(c.mean));
        assert(!std::isnan(c.gauss));
        assert(!std::isinf(c.mean));
        assert(!std::isinf(c.gauss));

        assert(c.min <= c.mean);
        assert(c.mean <= c.max);
    }

    return c;
}

} // namespace pmp
