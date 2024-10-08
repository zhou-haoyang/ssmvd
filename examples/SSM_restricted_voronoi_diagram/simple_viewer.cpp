#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <SSM_restricted_voronoi_diagram/Triangle_mesh_metric_traits.h>
#include <SSM_restricted_voronoi_diagram/AABB_metric_traits.h>
#include <SSM_restricted_voronoi_diagram.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/draw_point_set_3.h>

#include <iostream>

namespace RVD = CGAL::SSM_restricted_voronoi_diagram;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef CGAL::Surface_mesh<Point_3> Metric_polyhedron;
typedef CGAL::Surface_mesh<Point_3> Voronoi_diagram;

typedef CGAL::AABB_face_graph_triangle_primitive<Metric_polyhedron> Primitive;
typedef CGAL::AABB_traits_3<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Tree;

typedef RVD::AABB_metric_traits<Kernel, Tree> Metric_traits;
// typedef RVD::Triangle_mesh_metric_traits<Kernel, Metric_polyhedron> Metric_traits;
typedef RVD::SSM_restricted_voronoi_diagram_traits<Kernel, Surface_mesh, Metric_polyhedron, Voronoi_diagram,
                                                   Metric_traits>
    Traits;
typedef RVD::SSM_restricted_voronoi_diagram<Traits> SSM_restricted_voronoi_diagram;
typedef SSM_restricted_voronoi_diagram::Voronoi_diagram_data Voronoi_diagram_data;

using PS3 = typename CGAL::Point_set_3<Point_3>;

using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram>;
using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
using vd_face_descriptor = typename vd_graph_traits::face_descriptor;

struct Graphics_scene_options_sites
    : public CGAL::Graphics_scene_options<PS3, typename PS3::const_iterator, typename PS3::const_iterator,
                                          typename PS3::const_iterator> {
    bool colored_vertex(const PS3&, typename PS3::const_iterator) const { return true; }
    CGAL::IO::Color vertex_color(const PS3&, typename PS3::const_iterator) const { return CGAL::IO::Color(0, 220, 0); }
};

struct Graphics_scene_options_vd : public CGAL::Graphics_scene_options<Voronoi_diagram, vd_vertex_descriptor,
                                                                       vd_edge_descriptor, vd_face_descriptor> {
    Graphics_scene_options_vd() { disable_faces(); }

    bool colored_vertex(const Voronoi_diagram&, vd_vertex_descriptor) const { return true; }
    CGAL::IO::Color vertex_color(const Voronoi_diagram&, vd_vertex_descriptor) const {
        return CGAL::IO::Color(220, 0, 0);
    }

    bool colored_edge(const Voronoi_diagram&, vd_edge_descriptor) const { return true; }
    CGAL::IO::Color edge_color(const Voronoi_diagram&, vd_edge_descriptor) const { return CGAL::IO::Color(220, 0, 0); }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input mesh> <site file>" << std::endl;
        return EXIT_FAILURE;
    }
    const std::string mesh_file = argv[1], site_file = argv[2];

    Surface_mesh sm;
    if (!CGAL::IO::read_polygon_mesh(mesh_file, sm)) {
        std::cerr << "Invalid input file: " << mesh_file << std::endl;
        return EXIT_FAILURE;
    }

    SSM_restricted_voronoi_diagram ssm_vd(sm);

    std::ifstream is(site_file);
    if (!ssm_vd.read_sites(is)) {
        std::cerr << "Invalid site file: " << site_file << std::endl;
        return EXIT_FAILURE;
    }

    ssm_vd.build();

    PS3 ps_sites;
    for (auto it = ssm_vd.site_cbegin(); it != ssm_vd.site_cend(); ++it) {
        ps_sites.insert(it->point);
    }

    CGAL::Graphics_scene scene;
    CGAL::add_to_graphics_scene(sm, scene);
    CGAL::add_to_graphics_scene(ps_sites, scene, Graphics_scene_options_sites());
    CGAL::add_to_graphics_scene(ssm_vd.voronoi_diagram().graph, scene, Graphics_scene_options_vd());

    CGAL::draw_graphics_scene(scene);

    return EXIT_SUCCESS;
}
