#ifndef SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRAITS_H
#define SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRAITS_H

#include <CGAL/Default.h>
#include <Parametric_line/Parametric_line_3.h>

namespace CGAL::Surface_Voronoi_diagram_with_star_metrics {
/*!
 * \brief Traits class for Surface Voronoi diagrams with star metrics.
 *
 * \tparam K Kernel
 * \tparam SurfaceMesh, MetricPolyhedron Graph types for the surface mesh and metric polyhedron. Both should
 *   model the FaceGraph concept.
 * \tparam VoronoiDiagram Voronoi diagram graph type. Should model the MutableFaceGraph concept.
 * \tparam MetricTraits Metric traits class providing metric-related types and ray intersection functors. Built-in
 * implementations include AABB_metric_traits and Triangle_mesh_metric_traits.
 */
template <class K, class SurfaceMesh, class MetricPolyhedron, class VoronoiDiagram, class MetricTraits>
class Surface_Voronoi_diagram_with_star_metrics_traits : public K,
                                                         public Parametric_line_traits_3<K, typename K::FT>,
                                                         public MetricTraits
{
public:
  using Kernel = K;
  using FT = Kernel::FT;
  using T = FT;

  using Point_3 = Kernel::Point_3;
  using Vector_3 = Kernel::Vector_3;
  using Ray_3 = Kernel::Ray_3;
  using Plane_3 = Kernel::Plane_3;
  using Line_3 = Kernel::Line_3;
  using Parametric_line_traits = Parametric_line_traits_3<Kernel, T>;
  using Parametric_line_3 = Parametric_line_traits::Parametric_line_3;

  using Surface_mesh = SurfaceMesh;
  using Metric_polyhedron = MetricPolyhedron;
  using Voronoi_diagram = VoronoiDiagram;

  using Metric_traits = MetricTraits;
  using Metric_traits_data = typename Metric_traits::Data;

  using Parametric_line_traits::construct_point_3_object;
  using Parametric_line_traits::construct_point_on_3_object;
  using Parametric_line_traits::construct_vector_3_object;
  using Parametric_line_traits::intersect_3_object;

  using MetricTraits::construct_metric_data_object;
  using MetricTraits::metric_any_intersected_face_object;
  using MetricTraits::metric_any_intersection_object;
};
} // namespace CGAL::Surface_Voronoi_diagram_with_star_metrics

#endif // SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_TRAITS_H
