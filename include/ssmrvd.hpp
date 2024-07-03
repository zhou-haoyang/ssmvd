#pragma once

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>

#include <limits>
#include <unordered_map>
#include <vector>

namespace CGAL::SSM_restricted_voronoi_diagram {

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel> MetricPolyhedron;
using MeshVertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type;
using MetricVertexPointPMap = typename boost::property_map<MetricPolyhedron, vertex_point_t>::const_type;
using AABBTree =
    AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<MetricPolyhedron, Default, Tag_false>>>;

template <class Kernel, class NumericLimits = typename std::numeric_limits<typename Kernel::FT>>
struct Parametric_line_3 {
    static_assert(NumericLimits::has_infinity, "FT must have infinity");
    static constexpr auto INF = NumericLimits::infinity();

    using FT = Kernel::FT;
    using Point_3 = Kernel::Point_3;
    using Vector_3 = Kernel::Vector_3;

    Point_3 p;
    Vector_3 d;
    FT t_min = -INF, t_max = INF;

    Parametric_line_3(Point_3 p = {}, Vector_3 d = {}, FT t0 = -INF, FT t1 = INF)
        : p(std::move(p)), d(std::move(d)), t_min(std::min(t0, t1)), t_max(std::max(t0, t1)) {}

    static Parametric_line_3 segment(const Point_3 &p0, const Point_3 &p1) {
        return Parametric_line_3(p0, p1 - p0, 0, 1);
    }

    static Parametric_line_3 ray(const Point_3 &p0, const Vector_3 &d) { return Parametric_line_3(p0, d, 0, INF); }

    static Parametric_line_3 line(const Point_3 &p0, const Point_3 &p1) {
        return Parametric_line_3(p0, p1 - p0, -INF, INF);
    }

    Point_3 operator()(FT t) const { return p + t * d; }

    Parametric_line_3 reverse() const { return Parametric_line_3(p, -d, -t_max, -t_min); }

    void reverse_inplace() { *this = reverse(); }

    Point_3 p_min() const { return p + t_min * d; }

    Point_3 p_max() const { return p + t_max * d; }
};

// template <class Kernel, class MetricPolyhedron, class SurfaceMesh,
//           class VertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type,
//           class AABBTree =
//               AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<MetricPolyhedron, Default,
//               Tag_false>>>>
class SSM_restricted_voronoi_diagram {
   public:
    using FT = Kernel::FT;
    using Point_3 = Kernel::Point_3;
    using Ray_3 = Kernel::Ray_3;
    using Pline_3 = Parametric_line_3<Kernel>;

    using index_t = size_t;

    using mesh_graph_traits = typename boost::graph_traits<SurfaceMesh>;
    using mesh_vertex_descriptor = mesh_graph_traits::vertex_descriptor;
    using mesh_face_descriptor = mesh_graph_traits::face_descriptor;
    using mesh_halfedge_descriptor = mesh_graph_traits::halfedge_descriptor;

    using metric_graph_traits = typename boost::graph_traits<MetricPolyhedron>;
    using metric_face_descriptor = metric_graph_traits::face_descriptor;
    using metric_halfedge_descriptor = metric_graph_traits::halfedge_descriptor;

   private:
    struct M_cone_3 {
        index_t site_idx;
        metric_face_descriptor face;
    };

    struct Metric {
        MetricPolyhedron graph;
        MetricVertexPointPMap vpm;
    };

    struct Site {
        Point_3 point;
        index_t metric_idx;
    };

   public:
    SSM_restricted_voronoi_diagram(const SurfaceMesh &mesh, MeshVertexPointPMap vpm) : mesh(mesh), vpm(vpm) {}
    SSM_restricted_voronoi_diagram(const SurfaceMesh &mesh) : mesh(mesh), vpm(get(vertex_point, mesh)) {}

    void trace_boundary(mesh_halfedge_descriptor bhd) {
        M_cone_3 k0;
        find_nearest_site(vpm[source(bhd, mesh)], k0);
        auto [bhit, bhit_end] = halfedges_around_face(bhd, mesh);

        // Trace boundary
        for (;;) {
            auto bh_line = Pline_3::segment(vpm[source(*bhit, mesh)], vpm[target(*bhit, mesh)]);
        }
        // opposite()
    }

   private:
    std::vector<Site> sites;

    std::vector<Metric> metrics;
    AABBTree tree;

    const SurfaceMesh &mesh;
    MeshVertexPointPMap vpm;

    void build_tree() {
        tree.clear();
        for (auto &metric : metrics) {
            auto [fbegin, fend] = faces(metric.graph);
            tree.insert(fbegin, fend, metric.graph);
        }
        tree.build();
    }

    FT find_nearest_site(const Point_3 &p, M_cone_3 &m_cone) {
        build_tree();  // TODO: rebuild only when necessary
        Ray_3 ray(ORIGIN, p);
        std::vector<typename AABBTree::template Intersection_and_primitive_id<Ray_3>::Type> intersections;
        tree.all_intersections(ray, std::back_inserter(intersections));
        std::unordered_map<const MetricPolyhedron *, std::pair<FT, metric_face_descriptor>> m_weights;
        for (auto &[obj, m_face] : intersections) {
            auto &[fd, m_ptr] = m_face;
            auto isect_pt = boost::get<Point_3>(obj);
            FT weight = 1.0 / (isect_pt - ORIGIN).squared_length();
            m_weights[m_ptr] = std::make_pair(weight, fd);
        }

        FT d_min;

        for (index_t i = 0; i < sites.size(); ++i) {
            auto &[c, m_idx] = sites[i];
            auto it = m_weights.find(&metrics[m_idx].graph);
            if (it == m_weights.end()) {
                continue;
            }
            auto d = approximate_sqrt((p - c).squared_length() * it->second.first);
            if (d < d_min) {
                d_min = d;
                m_cone.site_idx = i;
                m_cone.face = it->second.second;
            }
        }
        return d_min;
    }

    bool isect(const Pline_3 &l, const M_cone_3 &k, FT &t_min, FT &t_max,
               std::optional<metric_halfedge_descriptor> &h_min, std::optional<metric_halfedge_descriptor> &h_max) {
        auto &[c, m_idx] = sites[k.site_idx];
        auto &metric = metrics[m_idx];
        Pline_3 l_inf(ORIGIN + (l.p - c), l.d);
        std::vector<std::tuple<metric_halfedge_descriptor, FT, Sign>> isects;
        for (auto hd : halfedges_around_face(halfedge(k.face, metric.graph), metric.graph)) {
            Vector_3 v0 = metric.vpm[source(hd, metric.graph)] - ORIGIN;
            Vector_3 v1 = metric.vpm[target(hd, metric.graph)] - ORIGIN;
            Vector_3 n = cross_product(v0, v1);
            if (is_zero(scalar_product(l_inf.d, n))) {
                // TODO: handle the case that the line lies on the face
                continue;
            }

            FT d_proj = scalar_product(l_inf.p - ORIGIN, n);
            FT ti = -scalar_product(l_inf.p - ORIGIN, n) / d_proj;
            Point_3 pi = l_inf(ti);
            Vector_3 vi = pi - ORIGIN;
            if (!(orientation(n, v0, vi) == POSITIVE && orientation(n, vi, v1) == POSITIVE)) {
                // p_i is outside the 2D cone
                continue;
            }

            isects.emplace_back(hd, ti, sign(d_proj));
            if (isects.size() == 2) {
                break;
            }
        }
        if (isects.empty()) {
            h_min = h_max = std::nullopt;
            return false;
        } else if (isects.size() == 1) {
            auto &[hd, t, sign] = isects.front();
            if (sign == POSITIVE) {
                t_min = t;
                h_min = hd;
                t_max = l.t_max;
                h_max = std::nullopt;
            } else {
                t_min = l.t_min;
                h_min = std::nullopt;
                t_max = t;
                h_max = hd;
            }

            if (t_min > t_max) {
                return false;
            }

        } else {
            auto &[hd0, t0, sign0] = isects.front();
            auto &[hd1, t1, sign1] = isects.back();
            if (t0 > t1) {
                std::swap(t0, t1);
                std::swap(hd0, hd1);
            }
            t_min = t0;
            h_min = hd0;
            t_max = t1;
            h_max = hd1;
        }

        if (t_min < l.t_min) {
            t_min = l.t_min;
            h_min = std::nullopt;
        }

        if (t_max > l.t_max) {
            t_max = l.t_max;
            h_max = std::nullopt;
        }

        return t_min <= t_max;
    }
};
}  // namespace CGAL::SSM_restricted_voronoi_diagram
