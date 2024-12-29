#ifndef SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_TRAITS_H
#define SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_TRAITS_H

#include <Parametric_line/Parametric_line_2.h>

#include <CGAL/Polygon_2.h>

namespace CGAL::SSM_voronoi_diagram {
template <class K, class VoronoiDiagram, class MetricTraits>
class SSM_voronoi_diagram_traits : public K, public Parametric_line_traits_2<K>, public MetricTraits {
   public:
    using Kernel = K;
    using FT = typename Kernel::FT;

    using Point_2 = typename Kernel::Point_2;
    using Vector_2 = typename Kernel::Vector_2;
    using Ray_2 = typename Kernel::Ray_2;
    using Line_2 = typename Kernel::Line_2;
    using Parametric_line_traits = Parametric_line_traits_2<Kernel, FT>;
    using Parametric_line_2 = typename Parametric_line_traits::Parametric_line_2;

    using Polygon_2 = Polygon_2<K>;
    using Metric_traits = MetricTraits;

    using Voronoi_diagram = VoronoiDiagram;

    using Parametric_line_traits::construct_point_2_object;
    using Parametric_line_traits::construct_point_on_2_object;
    using Parametric_line_traits::construct_vector_2_object;
    using Parametric_line_traits::intersect_2_object;
};
}  // namespace CGAL::SSM_voronoi_diagram

#endif  // SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_TRAITS_H
