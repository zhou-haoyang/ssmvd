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
#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>

#include <limits>
#include <unordered_map>
#include <vector>

namespace CGAL {

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel> MetricPolyhedron;
using MeshVertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type;
using MetricVertexPointPMap = typename boost::property_map<MetricPolyhedron, vertex_point_t>::const_type;
using AABBTree =
    AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<MetricPolyhedron, Default, Tag_false>>>;

template <class Kernel, class T = double>
struct Parametric_line_3 {
    static constexpr auto INF = std::numeric_limits<T>::infinity();

    using FT = Kernel::FT;
    using Point_3 = Kernel::Point_3;
    using Vector_3 = Kernel::Vector_3;

    Point_3 p;
    Vector_3 d;
    T t_min = -INF, t_max = INF;

    Parametric_line_3(Point_3 p = {}, Vector_3 d = {}, T t0 = -INF, T t1 = INF)
        : p(std::move(p)), d(std::move(d)), t_min(std::min(t0, t1)), t_max(std::max(t0, t1)) {}

    Parametric_line_3(const Kernel::Line_3 &line, T t0 = -INF, T t1 = INF)
        : Parametric_line_3(line.point(0), line.to_vector(), t0, t1) {}

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

template <class K, class T>
bool isect(const Parametric_line_3<K, T> &l, const typename K::Vector &n, const typename K::FT &D, T &t) {
    auto nd = scalar_product(n, l.d), np = scalar_product(n, l.p - ORIGIN) + D;
    if (is_zero(nd)) {
        t = std::numeric_limits<T>::infinity();
        return is_zero(np);
    }
    t = -np / nd;
    return true;
}

template <class K, class T>
bool isect(const Parametric_line_3<K, T> &l, const typename K::Plane_3 &p, T &t) {
    return isect(l, p.orthogonal_vector(), p.d(), t);
}

template <class Kernel, class FaceGraph, class VPMap>
Kernel::Plane_3 supporting_plane(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph &g,
                                 const VPMap &vpm) {
    auto i0 = source(hd, g), i1 = target(hd, g), i2 = target(next(hd, g), g);
    auto v0 = vpm[i0], v1 = vpm[i1], v2 = vpm[i2];
    return {v0, v1, v2};
}
namespace SSM_restricted_voronoi_diagram {
// template <class Kernel, class MetricPolyhedron, class SurfaceMesh,
//           class VertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type,
//           class AABBTree =
//               AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<MetricPolyhedron, Default,
//               Tag_false>>>>
class SSM_restricted_voronoi_diagram {
   public:
    using FT = Kernel::FT;
    using T = FT;
    using Point_3 = Kernel::Point_3;
    using Ray_3 = Kernel::Ray_3;
    using Plane_3 = Kernel::Plane_3;
    using Line_3 = Kernel::Line_3;
    using Pline_3 = Parametric_line_3<Kernel, T>;

    using index_t = size_t;

    using mesh_graph_traits = typename boost::graph_traits<SurfaceMesh>;
    using mesh_vertex_descriptor = mesh_graph_traits::vertex_descriptor;
    using mesh_face_descriptor = mesh_graph_traits::face_descriptor;
    using mesh_halfedge_descriptor = mesh_graph_traits::halfedge_descriptor;

    using metric_graph_traits = typename boost::graph_traits<MetricPolyhedron>;
    using metric_face_descriptor = metric_graph_traits::face_descriptor;
    using metric_face_iterator = metric_graph_traits::face_iterator;
    using metric_halfedge_descriptor = metric_graph_traits::halfedge_descriptor;

    static constexpr T INF = std::numeric_limits<T>::infinity();

   private:
    struct Metric {
        MetricPolyhedron graph;
        MetricVertexPointPMap vpm;
    };

    struct Site {
        Point_3 point;
        index_t metric_idx;
    };

    using site_iterator = std::vector<Site>::iterator;

    struct Cone_descriptor {
        index_t site_idx;
        metric_face_descriptor face;

        bool operator==(const Cone_descriptor &other) const { return site_idx == other.site_idx && face == other.face; }
    };

    struct InternalTrace {
        Plane_3 bisec;
        mesh_halfedge_descriptor he;
        Cone_descriptor k0, k1;
    };

    // class Cone_iterator {
    //     index_t site_idx;
    //     site_iterator site_it;
    //     metric_face_iterator face_it, face_end;

    //    public:
    //     using iterator_category = std::forward_iterator_tag;
    //     using value_type = Cone_descriptor;
    //     using difference_type = std::ptrdiff_t;
    //     using pointer = Cone_descriptor *;
    //     using reference = Cone_descriptor &;

    //     Cone_iterator(index_t site_idx, site_iterator site_it, metric_face_iterator face_it,
    //                   metric_face_iterator face_end)
    //         : site_idx(site_idx), site_it(site_it), face_it(face_it), face_end(face_end) {}

    //     value_type operator*() const { return {site_idx, *face_it}; }
    // };

   public:
    SSM_restricted_voronoi_diagram(const SurfaceMesh &mesh, MeshVertexPointPMap vpm) : mesh(mesh), vpm(vpm) {}
    SSM_restricted_voronoi_diagram(const SurfaceMesh &mesh) : mesh(mesh), vpm(get(vertex_point, mesh)) {}

    void trace_boundary(mesh_halfedge_descriptor bhd) {
        Cone_descriptor k0;
        find_nearest_site(vpm[source(bhd, mesh)], k0);
        // auto [bhit, bhit_end] = halfedges_around_face(bhd, mesh);

        // Loop over boundary halfedges
        for (auto bh : halfedges_around_face(bhd, mesh)) {
            auto b_line = Pline_3::segment(vpm[source(bh, mesh)], vpm[target(bh, mesh)]);
            // auto b_plane = mesh_face_plane(opposite(*bhit, mesh));
            Cone_descriptor k_prev;

            // Find all boundary-cone / boundary-bisector intersections
            for (;;) {
                Pline_3 b_segment;
                std::optional<metric_halfedge_descriptor> h_min, h_max;
                auto res = isect(k0, b_line, b_segment, h_min, h_max);
                CGAL_assertion_msg(res, "Boundary must intersect the cone");

                // Find nearest intersection of 2-site bisector plane with the boundary segment
                T dist_min = INF;
                T tb_min;
                Plane_3 bi_plane_min;
                Cone_descriptor k1_min;
                for (index_t site_idx = 0; site_idx < sites.size(); ++site_idx) {
                    if (site_idx == k0.site_idx) {
                        continue;
                    }

                    auto [c, m] = site(site_idx);
                    for (auto fd : faces(m.graph)) {
                        Cone_descriptor k1{site_idx, fd};
                        if (k1 == k_prev) continue;

                        Pline_3 b_overlap;
                        if (!isect(k1, b_segment, b_overlap)) continue;
                        Plane_3 bi_plane = get_bisect_plane(k0, k1);
                        if (bi_plane.is_degenerate()) continue;

                        T tb;
                        if (!CGAL::isect(b_overlap, bi_plane, tb)) continue;
                        T dist = tb - b_segment.t_min;
                        CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
                        if (dist < dist_min) {
                            dist_min = dist;
                            tb_min = tb;
                            bi_plane_min = bi_plane;
                            k1_min = k1;
                        }
                    }
                }

                if (dist_min < INF) {
                    // A 2-site bisector intersects the boundary segment
                    i_traces.push_back({bi_plane_min, opposite(bh, mesh), k0, k1_min});
                    // auto bi_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane.orthogonal_vector());
                    // auto ori = orientation(b_plane.orthogonal_vector(), bi_dir, b_line.d);
                    // auto bisect = Pline_3::ray(b_segment(tb_min), ori == POSITIVE ? bi_dir : -bi_dir);

                    b_line.t_min = tb_min;
                    k_prev = k0;
                    k0 = k1_min;
                } else {
                    // The bondary segment does not intersect leave the cone. Switch to the next boundary segment
                    if (!h_max) break;

                    // Otherwise switch to the next cone
                    k_prev = k0;
                    auto [c0, m0] = site(k0.site_idx);
                    k0.face = face(opposite(*h_max, m0.graph), m0.graph);
                }
            }
        }
    }

   private:
    std::vector<Site> sites;

    std::vector<Metric> metrics;
    AABBTree tree;

    const SurfaceMesh &mesh;
    MeshVertexPointPMap vpm;

    std::deque<InternalTrace> i_traces;

    std::pair<Site &, Metric &> site(index_t idx) {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    std::pair<const Site &, const Metric &> site(index_t idx) const {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    void build_tree() {
        tree.clear();
        for (auto &metric : metrics) {
            auto [fbegin, fend] = faces(metric.graph);
            tree.insert(fbegin, fend, metric.graph);
        }
        tree.build();
    }

    FT find_nearest_site(const Point_3 &p, Cone_descriptor &m_cone) {
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

    Plane_3 metric_face_plane(const Metric &m, metric_face_descriptor fd) const {
        return supporting_plane<Kernel>(halfedge(fd, m.graph), m.graph, m.vpm);
    }

    Plane_3 mesh_face_plane(mesh_halfedge_descriptor hd) const { return supporting_plane<Kernel>(hd, mesh, vpm); }

    Plane_3 get_bisect_plane(const Cone_descriptor &k0, const Cone_descriptor &k1) const {
        auto [c0, m0] = site(k0.site_idx);
        auto [c1, m1] = site(k1.site_idx);
        auto p0 = metric_face_plane(m0, k0.face);
        auto p1 = metric_face_plane(m1, k1.face);
        auto nd0 = p0.orthogonal_vector() / abs(p0.d());
        auto nd1 = p1.orthogonal_vector() / abs(p1.d());
        auto n = nd0 - nd1;
        auto d = scalar_product(nd1, c1.point - ORIGIN) - scalar_product(nd0, c0.point - ORIGIN);
        return Plane_3(n.x(), n.y(), n.z(), d);
    }

    bool isect(const Cone_descriptor &k, const Pline_3 &l, FT &t_min, FT &t_max,
               std::optional<metric_halfedge_descriptor> &h_min,
               std::optional<metric_halfedge_descriptor> &h_max) const {
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

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li,
               std::optional<metric_halfedge_descriptor> &h_min,
               std::optional<metric_halfedge_descriptor> &h_max) const {
        li = l;
        return isect(k, l, li.t_min, li.t_max, h_min, h_max);
    }

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li) const {
        li = l;
        std::optional<metric_halfedge_descriptor> h_min, h_max;
        return isect(k, l, li.t_min, li.t_max, h_min, h_max);
    }
};
}  // namespace SSM_restricted_voronoi_diagram
}  // namespace CGAL
