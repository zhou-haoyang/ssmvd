#include <SSM_voronoi_diagram.h>
#include <SSM_voronoi_diagram/Polygon_metric_traits.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polygon_2.h>

#include <gtest/gtest.h>
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

    std::vector<FT> find_interval_endpoints(Site_const_iterator site, const Parametric_line_2 &segment, FT tmin,
                                            std::optional<FT> tmax = std::nullopt) {
        if (tmax && *tmax < tmin) return {};

        auto p = construct_vector(site->point(), construct_point(segment));
        auto d = construct_vector(segment);
        std::vector<FT> res{tmin};
        for (auto edge : site->metric()->polygon().edges()) {
            auto isect = intersect_ray(p, d, construct_vector(CGAL::ORIGIN, construct_source(edge)), tmin, tmax);
            if (isect && !std::holds_alternative<Colinear>(*isect)) {
                auto [t, _] = std::get<std::pair<FT, FT>>(*isect);
                res.push_back(t);
            }
        }
        if (tmax) {
            res.push_back(*tmax);
        }
        std::sort(res.begin(), res.end());
        return res;
    }

    void check_interval_endpoints(Site_const_iterator site, const Parametric_line_2 &segment, FT tmin,
                                  std::optional<FT> tmax = std::nullopt) {
        Intervals res;
        find_intervals(site, segment, res, tmin, tmax);
        auto endpoints = find_interval_endpoints(site, segment, tmin, tmax);
        EXPECT_EQ(res.num_endpoints(), endpoints.size());
        auto it = res.endpoints().begin();
        for (auto t : endpoints) {
            EXPECT_FLOAT_EQ(*it, t);
            ++it;
        }
    }
};

TEST_F(SSM_voronoi_diagram_test, Empty_intervals) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), 1, 0);
}

TEST_F(SSM_voronoi_diagram_test, One_interval) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), 0, 0.4);
}

TEST_F(SSM_voronoi_diagram_test, One_interval_ray) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(-1, -1)), 0);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), 0, 1);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals_ray) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(2, 0)), 0);
}

TEST_F(SSM_voronoi_diagram_test, Two_intervals_ccw) {
    Intervals res;
    check_interval_endpoints(m_sites.begin(), Parametric_line_2(Point_2(0, 0), Point_2(0, 2)), 0, 1);
}
