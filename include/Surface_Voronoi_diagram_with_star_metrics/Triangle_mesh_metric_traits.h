#ifndef SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRIANGLE_MESH_METRIC_TRAITS_H
#define SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRIANGLE_MESH_METRIC_TRAITS_H

#include <CGAL/Origin.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/graph_traits.hpp>
#include <optional>

namespace CGAL::Surface_Voronoi_diagram_with_star_metrics {
/*!
 * \brief Traits class using a triangle mesh as the star-metric. Provides trivial ray-mesh intersection.
 *
 * \tparam K Kernel
 * \tparam G Graph type for the metric polyhedron. Should model the FaceGraph concept.
 */
template <class K, class G> class Triangle_mesh_metric_traits
{
public:
  using Kernel = K;
  using Graph = G;

  using Vector_3 = typename Kernel::Vector_3;
  using Point_3 = typename Kernel::Point_3;
  using Ray_3 = typename Kernel::Ray_3;
  using Triangle_3 = typename Kernel::Triangle_3;

  using graph_traits = typename boost::graph_traits<Graph>;
  using face_descriptor = typename graph_traits::face_descriptor;

  using Data = std::vector<std::pair<face_descriptor, Triangle_3>>;

  struct Construct_metric_data
  {
    template <class Graph> Data operator()(const Graph& g) const {
      auto construct_triangle = Kernel().construct_triangle_3_object();

      auto vpm = get(vertex_point, g);
      Data res;
      res.reserve(num_faces(g));
      for(auto fd : faces(g)) {
        auto hd = halfedge(fd, g);
        auto v0 = source(hd, g), v1 = target(hd, g), v2 = target(next(hd, g), g);
        auto p0 = get(vpm, v0), p1 = get(vpm, v1), p2 = get(vpm, v2);
        auto t = construct_triangle(p0, p1, p2);
        res.emplace_back(fd, t);
      }
      return res;
    }
  };

  struct Metric_any_intersection
  {
    std::optional<std::pair<Point_3, face_descriptor>> operator()(const Data& data, const Vector_3& d) const {
      auto intersect = Kernel().intersect_3_object();

      Ray_3 ray(ORIGIN, d);
      for(const auto& [fd, t] : data) {
        auto res = intersection(ray, t);
        if(!res)
          continue;
        auto pm = std::get_if<Point_3>(&*res);
        if(!pm)
          continue;
        return std::make_pair(*pm, fd);
      }
      return std::nullopt;
    }
  };

  struct Metric_any_intersected_face
  {
    std::optional<face_descriptor> operator()(const Data& data, const Vector_3& d) const {
      auto do_intersect = Kernel().do_intersect_3_object();

      Ray_3 ray(ORIGIN, d);
      for(const auto& [fd, t] : data) {
        if(do_intersect(ray, t))
          return fd;
      }
      return std::nullopt;
    }
  };

  auto construct_metric_data_object() const { return Construct_metric_data(); }
  auto metric_any_intersection_object() const { return Metric_any_intersection(); }
  auto metric_any_intersected_face_object() const { return Metric_any_intersected_face(); }
};
} // namespace CGAL::Surface_Voronoi_diagram_with_star_metrics
#endif // SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRIANGLE_MESH_METRIC_TRAITS_H
