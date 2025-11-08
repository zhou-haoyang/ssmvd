#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Random.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/generators.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/point_generators_3.h>
#include <Surface_Voronoi_diagram_with_star_metrics/Triangle_mesh_metric_traits.h>
#include <Surface_Voronoi_diagram_with_star_metrics/AABB_metric_traits.h>
#include <Surface_Voronoi_diagram_with_star_metrics.h>

#include <iostream>

namespace RVD = CGAL::Surface_Voronoi_diagram_with_star_metrics;

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
typedef RVD::Surface_Voronoi_diagram_with_star_metrics_traits<Kernel, Surface_mesh, Metric_polyhedron, Voronoi_diagram,
                                                   Metric_traits>
    Traits;
typedef RVD::Surface_Voronoi_diagram_with_star_metrics<Traits> Surface_Voronoi_diagram_with_star_metrics;
typedef Surface_Voronoi_diagram_with_star_metrics::Voronoi_diagram_data Voronoi_diagram_data;

using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram>;
using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
using vd_face_descriptor = typename vd_graph_traits::face_descriptor;
using Voronoi_diagram_data = Surface_Voronoi_diagram_with_star_metrics::Voronoi_diagram_data;

int main(int argc, char* argv[]) {
    CGAL::get_default_random() = CGAL::Random(0);

    const std::string mesh_file = argc > 1 ? argv[1] : CGAL::data_file_path("data/meshes/elephant.off");
    const size_t n_sites = argc > 2 ? std::stoi(argv[2]) : 10;
    long n_threads = argc > 3 ? std::stol(argv[3]) : -1;
    if (n_threads < 0) {
        n_threads = std::thread::hardware_concurrency();
    }

    Surface_mesh sm;
    if (!CGAL::IO::read_polygon_mesh(mesh_file, sm)) {
        std::cerr << "Invalid input file: " << mesh_file << std::endl;
        return EXIT_FAILURE;
    }

    Surface_Voronoi_diagram_with_star_metrics ssm_vd(sm, n_threads);

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

    std::cout << ssm_vd.i_timer().time() << " " << ssm_vd.b_timer().time() << std::endl;

    return EXIT_SUCCESS;
}
