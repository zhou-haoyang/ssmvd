#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <SSM_voronoi_diagram.h>
#include <SSM_voronoi_diagram/Polygon_metric_traits.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>

#include <gtest/gtest.h>
#include <algorithm>
#include <optional>
#include <variant>
#include <vector>

namespace VD = CGAL::SSM_voronoi_diagram;

typedef CGAL::Simple_cartesian<double> Kernel;

typedef CGAL::Surface_mesh<Kernel::Point_2> Voronoi_diagram;
typedef VD::Polygon_metric_traits<Kernel> Metric_traits;
typedef VD::SSM_voronoi_diagram_traits<Kernel, Voronoi_diagram, Metric_traits> Traits;
typedef VD::SSM_voronoi_diagram<Traits> SSM_voronoi_diagram;
using Point_2 = typename Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
typedef CGAL::Creator_uniform_2<double, Point_2> Creator;

TEST(Simple_cartesian_random, test) {
    Polygon_2 boundary;
    boundary.push_back(Point_2(-1, -1));
    boundary.push_back(Point_2(1, -1));
    boundary.push_back(Point_2(1, 1));
    boundary.push_back(Point_2(-1, 1));

    SSM_voronoi_diagram vd(boundary);

    Polygon_2 metric;
    metric.push_back(Point_2(0, -1));
    metric.push_back(Point_2(1, 0));
    metric.push_back(Point_2(0, 1));
    metric.push_back(Point_2(-1, 0));
    CGAL::perturb_points_2(metric.vertices_begin(), metric.vertices_end(), 0.1);

    std::vector<Point_2> sites;
    CGAL::Random_points_on_square_2<Point_2> gen;
    std::copy_n(gen, 30, std::back_inserter(sites));
    // CGAL::points_on_square_grid_2(1, 100, std::back_inserter(sites), Creator());
    CGAL::perturb_points_2(sites.begin(), sites.end(), 0.1);
    // for (auto c : sites) {
    //     std::cout << c << std::endl;
    // }

    auto m = vd.add_metric(metric);
    for (auto &site : sites) {
        vd.add_site(site, m);
    }

    // vd.trace_all_boundaries(boundary.edges_begin(), ++boundary.edges_begin(), true, true);
    vd.trace_all_boundaries(true, true);

    // vd.build();

    ASSERT_TRUE(vd.check_voronoi_diagram());
}
