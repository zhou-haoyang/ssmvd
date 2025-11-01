#ifndef SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H
#define SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H

#include <Parametric_line/Parametric_line_3.h>
#include <CGAL/Default.h>

namespace CGAL::SSM_restricted_voronoi_diagram {
template <class K, class SurfaceMesh, class MetricPolyhedron, class VoronoiDiagram, class MetricTraits>
class SSM_restricted_voronoi_diagram_traits : public K, public Parametric_line_traits_3<K>, public MetricTraits {
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

    using Metric_traits = MetricTraits;
    using Metric_traits_data = typename Metric_traits::Data;

    using Pline_traits::construct_point_3_object;
    using Pline_traits::construct_point_on_3_object;
    using Pline_traits::construct_vector_3_object;
    using Pline_traits::intersect_3_object;

    using MetricTraits::construct_metric_data_object;
    using MetricTraits::metric_any_intersected_face_object;
    using MetricTraits::metric_any_intersection_object;
};
}  // namespace CGAL::SSM_restricted_voronoi_diagram

#endif  // SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_TRAITS_H
