#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <Voronoi_diagram_with_star_metrics_2.h>
#include <Voronoi_diagram_with_star_metrics_2/Polygon_metric_traits.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_2.h>

#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

namespace VD = CGAL::Voronoi_diagram_with_star_metrics_2;

typedef CGAL::Simple_cartesian<double> Kernel;

using FT = typename Kernel::FT;
typedef CGAL::Surface_mesh<Kernel::Point_2> Voronoi_diagram;
typedef VD::Polygon_metric_traits<Kernel> Metric_traits;
typedef VD::Voronoi_diagram_with_star_metrics_2_traits<Kernel, Voronoi_diagram, Metric_traits> Traits;
typedef VD::Voronoi_diagram_with_star_metrics_2<Traits> Voronoi_diagram_with_star_metrics_2;
using Point_2 = typename Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Metric_iterator = typename Voronoi_diagram_with_star_metrics_2::Metric_iterator;
typedef CGAL::Creator_uniform_2<double, Point_2> Creator;

void make_rect(Polygon_2& p) {
  p.push_back(Point_2(-1, -1));
  p.push_back(Point_2(1, -1));
  p.push_back(Point_2(1, 1));
  p.push_back(Point_2(-1, 1));
}

void make_cross(Polygon_2& p) {
  p.push_back(Point_2(0, -1));
  p.push_back(Point_2(1, 0));
  p.push_back(Point_2(0, 1));
  p.push_back(Point_2(-1, 0));
  // CGAL::perturb_points_2(p.begin(), p.end(), 0.1);
}

void test_b_trace(FT a, size_t n, size_t n_sites) {
  CGAL::get_default_random() = CGAL::Random(0);

  Polygon_2 boundary;
  make_rect(boundary);

  Voronoi_diagram_with_star_metrics_2 vd(boundary);

  Polygon_2 metric;
  make_cross(metric);
  auto m = vd.add_metric(metric);

  for(size_t i = 0; i < n; i++) {
    std::vector<Point_2> sites;
    CGAL::Random_points_in_square_2<Point_2> gen(a);
    std::copy_n(gen, n_sites, std::back_inserter(sites));

    vd.clear_sites();
    for(auto& site : sites) {
      vd.add_site(site, m);
    }

    vd.reset();
    vd.trace_all_boundaries();

    ASSERT_TRUE(vd.check_voronoi_diagram());
  }
}

void test_build(FT a, size_t n, size_t n_sites) {
  CGAL::get_default_random() = CGAL::Random(0);

  Polygon_2 boundary;
  make_rect(boundary);

  Voronoi_diagram_with_star_metrics_2 vd(boundary);

  Polygon_2 metric;
  make_cross(metric);
  auto m = vd.add_metric(metric);

  for(size_t i = 0; i < n; i++) {
    std::vector<Point_2> sites;
    CGAL::Random_points_in_square_2<Point_2> gen(a);
    std::copy_n(gen, n_sites, std::back_inserter(sites));

    vd.clear_sites();
    for(auto& site : sites) {
      vd.add_site(site, m);
    }

    vd.build();

    ASSERT_TRUE(vd.check_voronoi_diagram());
  }
}

TEST(Simple_cartesian_random, trace_boundary) { test_b_trace(1, 100, 30); }

TEST(Simple_cartesian_random, trace_boundary_outside) { test_b_trace(2, 100, 30); }

TEST(Simple_cartesian_random, trace_boundary_outside_2) { test_b_trace(2, 100, 2); }

TEST(Simple_cartesian_random, build) { test_build(1, 100, 30); }

TEST(Simple_cartesian_random, build_outside) { test_build(2, 100, 30); }

TEST(Simple_cartesian_random, build_outside_2) { test_build(2, 100, 2); }
