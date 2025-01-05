#ifndef SSM_VORONOI_DIAGRAM_POLYGON_METRIC_TRAITS_H
#define SSM_VORONOI_DIAGRAM_POLYGON_METRIC_TRAITS_H

#include <SSM_voronoi_diagram/Polygon_2_vertex_pair_circulator.h>

#include <CGAL/Origin.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <optional>
#include <utility>

namespace CGAL::SSM_voronoi_diagram {
template <class K>
class Polygon_metric_traits {
   public:
    using Kernel = K;
    using Polygon_2 = Polygon_2<K>;
    using Point_2 = typename Kernel::Point_2;
    using Vector_2 = typename Kernel::Vector_2;
    using Ray_2 = typename Kernel::Ray_2;
    using Segment_2 = typename Kernel::Segment_2;

    using Metric = Polygon_2;
    using Metric_vertex_circulator = typename Polygon_2::Vertex_const_circulator;
    using Metric_edge_circulator = Polygon_2_vertex_pair_circulator<Metric_vertex_circulator>;

    struct Metric_any_intersection {
        std::optional<std::pair<Point_2, Metric_edge_circulator>> operator()(const Polygon_2 &p,
                                                                             const Vector_2 &d) const {
            auto it0 = Metric_edge_circulator(p.vertices_circulator());
            auto it = it0;
            do {
                auto [m0, m1] = *it;
                auto segment = K().construct_segment_2_object()(*m0, *m1);
                auto isect = CGAL::intersection(Ray_2(ORIGIN, d), segment);
                if (isect) {
                    // The intersection can not be segment here for a valid metric
                    return std::make_pair(std::get<Point_2>(*isect), it);
                }
                ++it;
            } while (it != it0);
            return std::nullopt;
        }
    };

    struct Metric_any_intersected_edge {
        std::optional<Metric_edge_circulator> operator()(const Polygon_2 &p, const Vector_2 &d) const {
            auto it0 = Metric_edge_circulator(p.vertices_circulator());
            auto it = it0;
            do {
                auto [m0, m1] = *it;
                auto segment = K().construct_segment_2_object()(*m0, *m1);
                if (CGAL::do_intersect(Ray_2(ORIGIN, d), segment)) {
                    return it;
                }
                ++it;
            } while (it != it0);
            return std::nullopt;
        }
    };

    auto metric_any_intersection_object() const { return Metric_any_intersection(); }

    auto metric_any_intersected_edge_object() const { return Metric_any_intersected_edge(); }
};
}  // namespace CGAL::SSM_voronoi_diagram

#endif  // SSM_VORONOI_DIAGRAM_POLYGON_METRIC_TRAITS_H
