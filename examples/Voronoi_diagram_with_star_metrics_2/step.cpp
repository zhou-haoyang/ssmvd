#include <CGAL/Dynamic_property_map.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <Voronoi_diagram_with_star_metrics_2.h>
#include <Voronoi_diagram_with_star_metrics_2/Polygon_metric_traits.h>

#include <CGAL/Point_set_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>

#include <CGAL/Graphics_scene_options.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_surface_mesh.h>
#include <qcolor.h>
#include <qnamespace.h>

#include <algorithm>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

namespace VD = CGAL::Voronoi_diagram_with_star_metrics_2;

typedef CGAL::Simple_cartesian<double> Kernel;
using GT = typename CGAL::Projection_traits_xy_3<Kernel>;

typedef CGAL::Surface_mesh<Kernel::Point_2> Voronoi_diagram;
typedef CGAL::Surface_mesh<Kernel::Point_3> Voronoi_diagram_3;

typedef VD::Polygon_metric_traits<Kernel> Metric_traits;
typedef VD::Voronoi_diagram_with_star_metrics_2_traits<Kernel, Voronoi_diagram, Metric_traits> Traits;
typedef VD::Voronoi_diagram_with_star_metrics_2<Traits> Voronoi_diagram_with_star_metrics_2;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
typedef CGAL::Creator_uniform_2<double, Point_2> Creator;
using PS2 = typename CGAL::Point_set_2<Kernel>;

using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram_3>;
using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
using vd_face_descriptor = typename vd_graph_traits::face_descriptor;

struct Graphics_scene_options_vd
    : public CGAL::
          Graphics_scene_options<Voronoi_diagram_3, vd_vertex_descriptor, vd_edge_descriptor, vd_face_descriptor>
{
  Graphics_scene_options_vd() { disable_faces(); }

  bool colored_vertex(const Voronoi_diagram_3&, vd_vertex_descriptor) const { return true; }
  CGAL::IO::Color vertex_color(const Voronoi_diagram_3&, vd_vertex_descriptor) const {
    return CGAL::IO::Color(220, 0, 0);
  }

  bool colored_edge(const Voronoi_diagram_3&, vd_edge_descriptor) const { return true; }
  CGAL::IO::Color edge_color(const Voronoi_diagram_3&, vd_edge_descriptor) const { return CGAL::IO::Color(220, 0, 0); }
};

int main(int argc, char** argv) {
  CGAL::get_default_random() = CGAL::Random(0);

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
  auto m = ssm_vd.add_metric(metric);
  // CGAL::perturb_points_2(metric.vertices_begin(), metric.vertices_end(), 0.1);

  auto random_sites = [&](size_t n) {
    std::vector<Point_2> sites;
    CGAL::Random_points_in_square_2<Point_2> gen(2);
    std::copy_n(gen, n, std::back_inserter(sites));

    ssm_vd.clear_sites();
    for(auto& site : sites) {
      ssm_vd.add_site(site, m);
    }
  };

  CGAL::Graphics_scene scene;

  auto update = [&]() {
    scene.clear();

    CGAL::add_to_graphics_scene(boundary, scene);

    auto& voronoi = ssm_vd.voronoi_diagram();
    for(auto ed : CGAL::edges(voronoi.graph)) {
      auto v0 = CGAL::source(ed, voronoi.graph);
      auto v1 = CGAL::target(ed, voronoi.graph);
      scene.add_segment(voronoi.vertex_point_map[v0], voronoi.vertex_point_map[v1], CGAL::IO::red());
    }

    size_t i = 0;
    for(auto& c : ssm_vd.sites()) {
      scene.add_point(c.point(), CGAL::IO::red());
      scene.add_text(c.point(), std::format("c{}", i++));

      for(auto m : c.metric()->polygon().vertices()) {
        scene.add_ray(c.point(), m - CGAL::ORIGIN, CGAL::IO::gray());
      }
    }

    for(auto vd : CGAL::vertices(voronoi.graph)) {
      scene.add_text(voronoi.graph.point(vd), std::format("v{}", vd.id()));
    }

    for(auto& tr : ssm_vd.i_traces()) {
      scene.add_ray(tr.bisector.p(), tr.bisector.d(), CGAL::IO::green());
    }
  };

  size_t n = argc > 1 ? std::stoi(argv[1]) : 30;
  random_sites(n);
  update();

  CGAL::Qt::QApplication_and_basic_viewer app(scene, "Step Building");
  if(app) {
    app.basic_viewer().on_key_pressed = [&](QKeyEvent* e, CGAL::Qt::Basic_viewer* basic_viewer) -> bool {
      const ::Qt::KeyboardModifiers modifiers = e->modifiers();
      if(modifiers == Qt::NoButton) {
        switch(e->key()) {
        case Qt::Key_C:
          random_sites(n);
          update();
          basic_viewer->redraw();
          break;
        case Qt::Key_B:
          ssm_vd.trace_all_boundaries(true, true);
          update();
          basic_viewer->redraw();
          break;
        case Qt::Key_Space:
          ssm_vd.step();
          update();
          basic_viewer->redraw();
          break;
        case Qt::Key_R:
          ssm_vd.reset();
          update();
          basic_viewer->redraw();
          break;
        case Qt::Key_X:
          ssm_vd.build();
          update();
          basic_viewer->redraw();
          break;
        default:
          return false;
        }
      }

      return true;
    };

    // Here we add shortcut descriptions
    app.basic_viewer().setKeyDescription(::Qt::Key_C, "Random sites");
    app.basic_viewer().setKeyDescription(::Qt::Key_B, "Trace all boundaries");
    app.basic_viewer().setKeyDescription(::Qt::Key_Space, "Step");
    app.basic_viewer().setKeyDescription(::Qt::Key_R, "Reset");
    app.basic_viewer().setKeyDescription(::Qt::Key_X, "Build");

    // Then we run the app
    app.run();
  }
}
