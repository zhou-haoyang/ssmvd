#include <SSM_voronoi_diagram.h>
#include <SSM_voronoi_diagram/Polygon_metric_traits.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polygon_2.h>

#include <gtest/gtest.h>
#include <optional>

namespace VD = CGAL::SSM_voronoi_diagram;

typedef CGAL::Simple_cartesian<double> Kernel;

typedef CGAL::Surface_mesh<Kernel::Point_2> Voronoi_diagram;
typedef VD::Polygon_metric_traits<Kernel> Metric_traits;
typedef VD::SSM_voronoi_diagram_traits<Kernel, Voronoi_diagram, Metric_traits> Traits;
typedef VD::SSM_voronoi_diagram<Traits> SSM_voronoi_diagram;
using Point_2 = typename Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

class SSM_voronoi_diagram_test : public SSM_voronoi_diagram, public ::testing::Test {
   protected:
    SSM_voronoi_diagram_test() : SSM_voronoi_diagram(boundary) {}
    Polygon_2 boundary;

    void SetUp() override {
        Polygon_2 metric;
        metric.push_back(Point_2(0, -1));
        metric.push_back(Point_2(1, 0));
        metric.push_back(Point_2(0, 1));
        metric.push_back(Point_2(-1, 0));
        add_site(Point_2(1, 1), add_metric(std::move(metric)));
    }
};

TEST_F(SSM_voronoi_diagram_test, Empty_intervals) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), res, 1, 0);
    EXPECT_TRUE(res.empty());
}

TEST_F(SSM_voronoi_diagram_test, One_interval) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), res, 0, 0.4);
    EXPECT_TRUE(res.num_cones() == 1);
}

TEST_F(SSM_voronoi_diagram_test, One_interval_ray) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(-1, -1)), res, 0);
    EXPECT_TRUE(res.num_cones() == 1);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), res, 0, 1);
    EXPECT_TRUE(res.num_cones() == 2);
    auto it = res.begin();
    EXPECT_TRUE((*it).second.tmin == 0);
    EXPECT_TRUE((*it++).second.tmax == 0.5);
    EXPECT_TRUE((*it).second.tmin == 0.5);
    EXPECT_TRUE((*it).second.tmax == 1);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals_ray) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), res, 0);
    EXPECT_TRUE(res.num_cones() == 2);
    EXPECT_TRUE(res.is_infinite());
    auto it = res.begin();
    EXPECT_TRUE((*it).second.tmin == 0);
    EXPECT_TRUE((*it++).second.tmax == 0.5);
    EXPECT_TRUE((*it).second.tmin == 0.5);
    EXPECT_TRUE((*it).second.tmax == std::nullopt);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals_ccw) {
    Segment_cone_intervals res;
    find_segment_cone_intervals(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(0, 2)), res, 0, 1);
    EXPECT_TRUE(res.num_cones() == 2);
    auto it = res.begin();
    EXPECT_TRUE((*it).second.tmin == 0);
    EXPECT_TRUE((*it++).second.tmax == 0.5);
    EXPECT_TRUE((*it).second.tmin == 0.5);
    EXPECT_TRUE((*it).second.tmax == 1);
}
