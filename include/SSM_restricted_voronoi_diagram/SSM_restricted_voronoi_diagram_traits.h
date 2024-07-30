#ifndef SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H
#define SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H

#include <Parametric_line/Parametric_line_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Default.h>

namespace CGAL::SSM_restricted_voronoi_diagram {
template <class K, class SurfaceMesh, class MetricPolyhedron, class VoronoiDiagram, class MetricAABBTree = Default>
class SSM_restricted_voronoi_diagram_traits : public K, public Parametric_line_traits_3<K> {
   public:
    using Kernel = K;
    using FT = Kernel::FT;
    using T = double;

    using Point_3 = Kernel::Point_3;
    using Vector_3 = Kernel::Vector_3;
    using Ray_3 = Kernel::Ray_3;
    using Plane_3 = Kernel::Plane_3;
    using Line_3 = Kernel::Line_3;
    using Pline_traits = Parametric_line_traits_3<Kernel, T>;
    using Pline_3 = Pline_traits::Parametric_line_3;

    using Surface_mesh = SurfaceMesh;
    using Metric_polyhedron = MetricPolyhedron;
    using Voronoi_diagram = VoronoiDiagram;

    using Metric_AABB_tree =
        Default::Get<MetricAABBTree,
                     AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<Metric_polyhedron>>>>::type;

    using Pline_traits::construct_point_on_3_object;
    using Pline_traits::construct_vector_3_object;
    using Pline_traits::intersect_3_object;
};
}  // namespace CGAL::SSM_restricted_voronoi_diagram

#endif  // SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H
