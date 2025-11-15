#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <Voronoi_diagram_with_star_metrics_2.h>
#include <Voronoi_diagram_with_star_metrics_2/Polygon_metric_traits.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/point_generators_2.h>

namespace VD = CGAL::Voronoi_diagram_with_star_metrics_2;

using Kernel = CGAL::Simple_cartesian<double>;
using Voronoi_diagram = CGAL::Surface_mesh<Kernel::Point_2>;
using Metric_traits = VD::Polygon_metric_traits<Kernel>;
using Traits = VD::Voronoi_diagram_with_star_metrics_2_traits<Kernel, Voronoi_diagram, Metric_traits>;
using Voronoi_diagram_with_star_metrics_2 = VD::Voronoi_diagram_with_star_metrics_2<Traits>;
using Point_2 = typename Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

int main(int argc, char* argv[]) {
  const size_t n_sites = argc > 1 ? std::stoi(argv[1]) : 30;
  Polygon_2 boundary;
  boundary.push_back(Point_2(-1, -1));
  boundary.push_back(Point_2(1, -1));
  boundary.push_back(Point_2(1, 1));
  boundary.push_back(Point_2(-1, 1));

  Voronoi_diagram_with_star_metrics_2 ssm_vd(boundary);

  Polygon_2 metric;
  metric.push_back(Point_2(0, -1));
  metric.push_back(Point_2(1, 0));
  metric.push_back(Point_2(0, 1));
  metric.push_back(Point_2(-1, 0));

  std::vector<Point_2> sites;
  CGAL::Random_points_in_square_2<Point_2> gen;
  std::copy_n(gen, n_sites, std::back_inserter(sites));

  auto m = ssm_vd.add_metric(metric);
  for(auto& site : sites) {
    ssm_vd.add_site(site, m);
  }

  ssm_vd.build();

  CGAL_assertion(CGAL::is_valid_face_graph(ssm_vd.voronoi_diagram().graph, true));

  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(boundary, scene);

  auto& voronoi = ssm_vd.voronoi_diagram();

  for(const auto& ed : CGAL::edges(voronoi.graph)) {
    auto v0 = CGAL::source(ed, voronoi.graph);
    auto v1 = CGAL::target(ed, voronoi.graph);
    scene.add_segment(voronoi.vertex_point_map[v0], voronoi.vertex_point_map[v1], CGAL::IO::green());
  }

  int i = 0;
  for(auto& p : sites) {
    scene.add_point(p, CGAL::IO::yellow());
  }

  CGAL::draw_graphics_scene(scene);
}
