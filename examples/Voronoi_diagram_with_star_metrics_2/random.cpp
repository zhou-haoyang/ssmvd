#include <CGAL/Dynamic_property_map.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <Voronoi_diagram_with_star_metrics_2.h>
#include <Voronoi_diagram_with_star_metrics_2/Polygon_metric_traits.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Point_set_2.h>

#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/IO/Color.h>

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

struct Graphics_scene_options_vd : public CGAL::Graphics_scene_options<Voronoi_diagram_3, vd_vertex_descriptor,
                                                                       vd_edge_descriptor, vd_face_descriptor> {
    Graphics_scene_options_vd() { disable_faces(); }

    bool colored_vertex(const Voronoi_diagram_3 &, vd_vertex_descriptor) const { return true; }
    CGAL::IO::Color vertex_color(const Voronoi_diagram_3 &, vd_vertex_descriptor) const {
        return CGAL::IO::Color(220, 0, 0);
    }

    bool colored_edge(const Voronoi_diagram_3 &, vd_edge_descriptor) const { return true; }
    CGAL::IO::Color edge_color(const Voronoi_diagram_3 &, vd_edge_descriptor) const {
        return CGAL::IO::Color(220, 0, 0);
    }
};

int main() {
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
    CGAL::perturb_points_2(metric.vertices_begin(), metric.vertices_end(), 0.1);

    std::vector<Point_2> sites;
    CGAL::Random_points_in_square_2<Point_2> gen;
    std::copy_n(gen, 30, std::back_inserter(sites));
    // CGAL::points_on_square_grid_2(1, 100, std::back_inserter(sites), Creator());
    CGAL::perturb_points_2(sites.begin(), sites.end(), 0.1);
    // for (auto c : sites) {
    //     std::cout << c << std::endl;
    // }

    auto m = ssm_vd.add_metric(metric);
    for (auto &site : sites) {
        ssm_vd.add_site(site, m);
    }

    // vd.trace_all_boundaries(boundary.edges_begin(), ++boundary.edges_begin(), true, true);
    ssm_vd.trace_all_boundaries(true, true);
    ssm_vd.trace_faces();
    // ssm_vd.build();

    CGAL_assertion(CGAL::is_valid_face_graph(ssm_vd.voronoi_diagram().graph, true));

    CGAL::Graphics_scene scene;
    CGAL::add_to_graphics_scene(boundary, scene);

    // /// convert Voronoi_diagram to Voronoi_diagram_3
    // Voronoi_diagram_3 vd3;
    // using Vertex_point_3_property = CGAL::dynamic_vertex_property_t<Point_3>;
    // using Vertex_point_3 = typename boost::property_map<Voronoi_diagram, Vertex_point_3_property>::const_type;
    auto &voronoi = ssm_vd.voronoi_diagram();
    // Vertex_point_3 vpm3 = CGAL::get(Vertex_point_3_property{}, voronoi.graph);
    // for (auto vd : CGAL::vertices(voronoi.graph)) {
    //     Point_2 p = voronoi.vpm[vd];
    //     // vpm3[vd] = Point_3(p.x(), p.y(), 0);
    //     put(vpm3, vd, Point_3(p.x(), p.y(), 0));
    // }
    // CGAL::copy_face_graph(voronoi.graph, vd3, CGAL::parameters::vertex_point_map(vpm3));
    // ///

    // CGAL::add_to_graphics_scene(vd3, scene, Graphics_scene_options_vd());

    for (auto ed : CGAL::edges(voronoi.graph)) {
        auto v0 = CGAL::source(ed, voronoi.graph);
        auto v1 = CGAL::target(ed, voronoi.graph);
        scene.add_segment(voronoi.vpm[v0], voronoi.vpm[v1], CGAL::IO::red());
    }

    int i = 0;
    for (auto &p : sites) {
        scene.add_point(p, CGAL::IO::red());
        scene.add_text(p, std::format("c{}", i++));
    }

    for (auto vd : CGAL::vertices(voronoi.graph)) {
        scene.add_text(voronoi.graph.point(vd), std::format("v{}", vd.id()));
    }

    CGAL::draw_graphics_scene(scene);
}
