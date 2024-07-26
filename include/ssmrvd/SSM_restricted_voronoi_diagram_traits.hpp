#pragma once

#include <ssmrvd/Parametric_line_3.hpp>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Default.h>

namespace CGAL::SSM_restricted_voronoi_diagram {
template <class K, class SurfaceMesh, class MetricPolyhedron, class VoronoiDiagram, class MetricAABBTree = Default>
class SSM_restricted_voronoi_diagram_traits : public K {
   public:
    using Kernel = K;
    using FT = Kernel::FT;
    using T = double;

    using Point_3 = Kernel::Point_3;
    using Vector_3 = Kernel::Vector_3;
    using Ray_3 = Kernel::Ray_3;
    using Plane_3 = Kernel::Plane_3;
    using Line_3 = Kernel::Line_3;
    using Pline_3 = Parametric_line_3<Kernel, T>;

    using Surface_mesh = SurfaceMesh;
    using Metric_polyhedron = MetricPolyhedron;
    using Voronoi_diagram = VoronoiDiagram;

    using Metric_AABB_tree =
        Default::Get<MetricAABBTree,
                     AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<Metric_polyhedron>>>>::type;
};
}  // namespace CGAL::SSM_restricted_voronoi_diagram
