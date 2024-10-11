#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/point_generators_3.h>
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
using Voronoi_diagram_data = SSM_restricted_voronoi_diagram::Voronoi_diagram_data;

struct Graphics_scene_options_sites
    : public CGAL::Graphics_scene_options<PS3, typename PS3::const_iterator, typename PS3::const_iterator,
                                          typename PS3::const_iterator> {
    bool colored_vertex(const PS3&, typename PS3::const_iterator) const { return true; }
    CGAL::IO::Color vertex_color(const PS3&, typename PS3::const_iterator) const { return CGAL::IO::Color(0, 220, 0); }
};

struct Graphics_scene_options_vd : public CGAL::Graphics_scene_options<Voronoi_diagram, vd_vertex_descriptor,
                                                                       vd_edge_descriptor, vd_face_descriptor> {
    const Voronoi_diagram_data& voronoi;
    Graphics_scene_options_vd(const Voronoi_diagram_data& voronoi) : voronoi(voronoi) {
        //  disable_faces();
    }

    bool colored_vertex(const Voronoi_diagram&, vd_vertex_descriptor) const { return true; }
    CGAL::IO::Color vertex_color(const Voronoi_diagram&, vd_vertex_descriptor) const {
        return CGAL::IO::Color(220, 0, 0);
    }

    bool colored_edge(const Voronoi_diagram&, vd_edge_descriptor ed) const { return true; }
    CGAL::IO::Color edge_color(const Voronoi_diagram&, vd_edge_descriptor ed) const {
        return get(voronoi.edge_bisector_map, ed) ? CGAL::IO::red() : CGAL::IO::black();
    }
};

int main(int argc, char* argv[]) {
    const std::string mesh_file = argc > 1 ? argv[1] : CGAL::data_file_path("data/meshes/elephant.off");
    const size_t n_sites = argc > 2 ? std::stoi(argv[2]) : 10;

    Surface_mesh sm;
    if (!CGAL::IO::read_polygon_mesh(mesh_file, sm)) {
        std::cerr << "Invalid input file: " << mesh_file << std::endl;
        return EXIT_FAILURE;
    }

    SSM_restricted_voronoi_diagram ssm_vd(sm);

    Metric_polyhedron mp;
    CGAL::make_tetrahedron(Point_3(-1, -1, -1), Point_3(1, 1, -1), Point_3(-1, 1, 1), Point_3(1, -1, 1), mp);
    auto m = ssm_vd.add_metric(mp);

    CGAL::Random_points_in_triangle_mesh_3<Surface_mesh> gen(sm);
    std::vector<Point_3> sites;
    std::copy_n(gen, n_sites, std::back_inserter(sites));
    for (const auto& site : sites) {
        ssm_vd.add_site(site, m);
    }

    ssm_vd.build();

    PS3 ps_sites;
    for (auto it = ssm_vd.site_cbegin(); it != ssm_vd.site_cend(); ++it) {
        ps_sites.insert(it->point);
    }

    CGAL::Graphics_scene scene;
    // CGAL::add_to_graphics_scene(sm, scene);
    CGAL::add_to_graphics_scene(ps_sites, scene, Graphics_scene_options_sites());
    CGAL::add_to_graphics_scene(ssm_vd.voronoi_diagram().graph, scene,
                                Graphics_scene_options_vd(ssm_vd.voronoi_diagram()));

    CGAL::draw_graphics_scene(scene);

    return EXIT_SUCCESS;
}
