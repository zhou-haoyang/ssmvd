#ifndef SSM_RESTRICTED_VORONOI_DIAGRAM_AABB_METRIC_TRAITS_H
#define SSM_RESTRICTED_VORONOI_DIAGRAM_AABB_METRIC_TRAITS_H

#include <CGAL/Origin.h>
#include <optional>

namespace CGAL::SSM_restricted_voronoi_diagram {
template <class K, class Tree>
class AABB_metric_traits {
   public:
    using Kernel = K;
    using Data = Tree;

    using Id = typename Tree::Primitive_id;

    using Vector_3 = typename Kernel::Vector_3;
    using Point_3 = typename Kernel::Point_3;
    using Ray_3 = typename Kernel::Ray_3;

    struct Construct_metric_data {
        template <class Graph>
        Data operator()(const Graph& g) const {
            auto [first, last] = faces(g);
            return Tree(first, last, g);
        }
    };

    struct Metric_any_intersection {
        std::optional<std::pair<Point_3, Id>> operator()(const Data& data, const Vector_3& d) const {
            Ray_3 ray(ORIGIN, d);
            auto res = data.any_intersection(ray);
            if (!res) return std::nullopt;

            auto [obj, fd] = *res;
            auto pm = std::get<Point_3>(&obj);
            if (!pm) return std::nullopt;

            return std::make_pair(*pm, fd);
        }
    };

    struct Metric_any_intersected_face {
        std::optional<Id> operator()(const Data& data, const Vector_3& d) const {
            Ray_3 ray(ORIGIN, d);
            auto res = data.any_intersected_primitive(ray);
            if (!res) return std::nullopt;
            return *res;
        }
    };

    auto construct_metric_data_object() const { return Construct_metric_data(); }
    auto metric_any_intersection_object() const { return Metric_any_intersection(); }
    auto metric_any_intersected_face_object() const { return Metric_any_intersected_face(); }
};
}  // namespace CGAL::SSM_restricted_voronoi_diagram
#endif  // SSM_RESTRICTED_VORONOI_DIAGRAM_AABB_METRIC_TRAITS_H
