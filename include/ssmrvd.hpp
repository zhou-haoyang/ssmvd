#pragma once

#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Origin.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>
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

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
// typedef CGAL::Polyhedron_3<Kernel> MetricPolyhedron;
// typedef CGAL::Surface_mesh<Kernel::Point_3> VoronoiDiagramGraph;

// using MeshVertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type;
// using MetricVertexPointPMap = typename boost::property_map<MetricPolyhedron, vertex_point_t>::const_type;
// using VoronoiDiagramVertexPointPMap = typename boost::property_map<VoronoiDiagramGraph, vertex_point_t>::const_type;

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
bool isect(const Parametric_line_3<K, T> &l, const typename K::Vector_3 &n, const typename K::FT &D, T &t) {
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

template <class K, class T>
bool isect(const Parametric_line_3<K, T> &l0, const Parametric_line_3<K, T> &l1, T &t0, T &t1, bool coplanar = false) {
    auto d0 = l0.d, d1 = l1.d;
    auto p0 = l0.p, p1 = l1.p;
    auto n = cross_product(d0, d1);
    if (is_zero(n.squared_length())) {
        // parallel
        return false;
    }

    auto tn0 = cross_product((p1 - p0), d1);
    if (!coplanar && !is_zero(scalar_product(tn0, n))) {
        // not coplanar
        return false;
    }
    t0 = tn0.x() / n.x();
    auto tn1 = cross_product((p1 - p0), d0);
    t1 = tn1.x() / n.x();

    if (t0 < l0.t_min || t0 > l0.t_max || t1 < l1.t_min || t1 > l1.t_max) {
        return false;
    }
    return true;
}

template <class Kernel, class FaceGraph, class VPMap>
Kernel::Plane_3 supporting_plane(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph &g,
                                 const VPMap &vpm) {
    auto i0 = source(hd, g), i1 = target(hd, g), i2 = target(next(hd, g), g);
    auto v0 = vpm[i0], v1 = vpm[i1], v2 = vpm[i2];
    return {v0, v1, v2};
}

template <class K, class T, class FaceGraph, class VPMap>
Parametric_line_3<K, T> edge_segment(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                                     const FaceGraph &g, const VPMap &vpm) {
    auto i0 = source(hd, g), i1 = target(hd, g);
    return Parametric_line_3<K, T>::segment(vpm[i0], vpm[i1]);
}

namespace SSM_restricted_voronoi_diagram {
template <class Kernel, class MetricPolyhedron, class SurfaceMesh, class VoronoiDiagramGraph,
          class MeshVertexPointPMap = typename boost::property_map<SurfaceMesh, vertex_point_t>::const_type,
          class MetricVertexPointPMap = typename boost::property_map<MetricPolyhedron, vertex_point_t>::const_type,
          class VoronoiDiagramVertexPointPMap =
              typename boost::property_map<VoronoiDiagramGraph, vertex_point_t>::const_type>
class SSM_restricted_voronoi_diagram {
   public:
    using FT = Kernel::FT;
    using T = FT;
    using Point_3 = Kernel::Point_3;
    using Vector_3 = Kernel::Vector_3;
    using Ray_3 = Kernel::Ray_3;
    using Plane_3 = Kernel::Plane_3;
    using Line_3 = Kernel::Line_3;
    using Pline_3 = Parametric_line_3<Kernel, T>;

    using Metric_AABB_tree = AABB_tree<AABB_traits<Kernel, AABB_face_graph_triangle_primitive<MetricPolyhedron>>>;

    using index_t = std::ptrdiff_t;

    using mesh_graph_traits = typename boost::graph_traits<SurfaceMesh>;
    using mesh_vertex_descriptor = mesh_graph_traits::vertex_descriptor;
    using mesh_face_descriptor = mesh_graph_traits::face_descriptor;
    using mesh_halfedge_descriptor = mesh_graph_traits::halfedge_descriptor;

    using metric_graph_traits = typename boost::graph_traits<MetricPolyhedron>;
    using metric_face_descriptor = metric_graph_traits::face_descriptor;
    using metric_face_iterator = metric_graph_traits::face_iterator;
    using metric_halfedge_descriptor = metric_graph_traits::halfedge_descriptor;
    using metric_halfedge_descriptor_opt = std::optional<metric_halfedge_descriptor>;

    static constexpr T INF = std::numeric_limits<T>::infinity();

    //    private:
    struct Metric {
        MetricPolyhedron graph;
        MetricVertexPointPMap vpm;
        Metric_AABB_tree tree;

        auto cone_face_bases(metric_halfedge_descriptor hd) const {
            auto v0 = vpm[source(hd, graph)] - ORIGIN;
            auto v1 = vpm[target(hd, graph)] - ORIGIN;
            return std::make_pair(v0, v1);
        }

        /**
         * @brief Return the normal of the 2D cone face, pointing inward
         *
         * @param hd
         * @return Vector_3
         */
        Vector_3 cone_face_orthogonal_vector(metric_halfedge_descriptor hd) const {
            // TODO: cache the result
            auto [v0, v1] = cone_face_bases(hd);
            Vector_3 n = cross_product(v0, v1);
            return n;
        }

        // bool isect(const Pline_3 &l) {
        //     auto n = m.cone_face_orthogonal_vector(hd);
        //     FT nd = scalar_product(n, l_inf.d);
        //     if (is_zero(nd)) {
        //         // TODO: handle the case that the line lies on the face
        //         return false;
        //     }

        //     FT np = scalar_product(n, l_inf.p - ORIGIN);
        //     FT ti = -np / nd;
        //     Point_3 pi = l_inf(ti);
        //     Vector_3 vi = pi - ORIGIN;
        //     if (!(orientation(n, v0, vi) == POSITIVE && orientation(n, vi, v1) == POSITIVE)) {
        //         // p_i is outside the 2D cone
        //         continue;
        //     }
        // }
    };

    struct Voronoi_diagram {
        using graph_traits = typename boost::graph_traits<VoronoiDiagramGraph>;
        using vertex_descriptor = graph_traits::vertex_descriptor;
        using edge_descriptor = graph_traits::edge_descriptor;
        using halfedge_descriptor = graph_traits::halfedge_descriptor;

        VoronoiDiagramGraph graph;
        VoronoiDiagramVertexPointPMap vpm;

        halfedge_descriptor opposite(halfedge_descriptor hd) const { return CGAL::opposite(hd, graph); }

        vertex_descriptor add_vertex(const Point_3 &p) {
            auto vd = CGAL::add_vertex(graph);
            vpm[vd] = p;
            return vd;
        }

        halfedge_descriptor add_ray(vertex_descriptor v) {
            auto ed = CGAL::add_edge(graph);
            auto hd = halfedge(ed, graph);
            set_target(hd, v);
            return opposite(hd);
        }

        halfedge_descriptor add_ray(const Point_3 &p) { return add_ray(add_vertex(p)); }

        void set_target(halfedge_descriptor hd, vertex_descriptor v) { CGAL::set_target(hd, v, graph); }

        void set_next(halfedge_descriptor hd, halfedge_descriptor next) { CGAL::set_next(hd, next, graph); }

        void set_next_double_sided(halfedge_descriptor hd, halfedge_descriptor next) {
            set_next(hd, next);
            set_next(opposite(hd), opposite(next));
        }
    };

    struct Site {
        Point_3 point;
        index_t metric_idx;
    };

    using site_iterator = std::vector<Site>::iterator;

    struct Cone_descriptor {
        index_t site_idx = -1;
        metric_face_descriptor face;

        bool operator==(const Cone_descriptor &other) const { return site_idx == other.site_idx && face == other.face; }
    };

    struct Cone_index {};

    struct InternalTrace {
        Pline_3 bisect_line;
        Plane_3 bisect_plane, face_plane;
        mesh_halfedge_descriptor face_hd;
        Cone_descriptor k0, k1, k_prev;
        Voronoi_diagram::halfedge_descriptor v_hd;
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

    void add_site(const Point_3 &p, index_t metric_idx) { sites.push_back({p, metric_idx}); }

    index_t add_metric(const MetricPolyhedron &m, MetricVertexPointPMap vpm) {
        auto idx = metrics.size();
        auto [fbegin, fend] = faces(m);
        metrics.push_back({m, vpm, Metric_AABB_tree(fbegin, fend, m)});
        return idx;
    }

    index_t add_metric(const MetricPolyhedron &m) { return add_metric(m, get(vertex_point, m)); }

    void clear_sites() { sites.clear(); }

    void clear_metrics() { metrics.clear(); }

    void trace_boundary(mesh_halfedge_descriptor bhd) {
        Cone_descriptor k0;
        find_nearest_site(vpm[source(bhd, mesh)], k0);
        // auto [bhit, bhit_end] = halfedges_around_face(bhd, mesh);

        // Loop over boundary halfedges
        for (auto bh : halfedges_around_face(bhd, mesh)) {
            auto b_line = Pline_3::segment(vpm[source(bh, mesh)], vpm[target(bh, mesh)]);
            auto b_plane = mesh_face_plane(opposite(bh, mesh));
            Cone_descriptor k_prev;

            // Find all boundary-cone / boundary-bisector intersections
            for (;;) {
                Pline_3 b_segment;
                metric_halfedge_descriptor_opt h_min, h_max;
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
                    auto bisect_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane.orthogonal_vector());
                    auto orient = orientation(b_plane.orthogonal_vector(), bisect_dir, b_line.d);
                    auto pt_start = b_segment(tb_min);
                    auto bisect_line = Pline_3::ray(pt_start, orient == POSITIVE ? bisect_dir : -bisect_dir);
                    auto v_hd = vd.add_ray(pt_start);
                    i_traces.push_back({
                        bisect_line,
                        bi_plane_min,
                        b_plane,
                        opposite(bh, mesh),
                        k0,
                        k1_min,
                        {},
                        v_hd,
                    });
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

    //    private:
    std::vector<Site> sites;

    std::vector<Metric> metrics;
    // AABBTree tree;

    const SurfaceMesh &mesh;
    MeshVertexPointPMap vpm;

    Voronoi_diagram vd;

    std::deque<InternalTrace> i_traces;

    std::pair<Site &, Metric &> site(index_t idx) {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    std::pair<const Site &, const Metric &> site(index_t idx) const {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    // void build_tree() {
    //     tree.clear();
    //     for (auto &metric : metrics) {
    //         auto [fbegin, fend] = faces(metric.graph);
    //         tree.insert(fbegin, fend, metric.graph);
    //     }
    //     tree.build();
    // }

    T find_nearest_site(const Point_3 &p, Cone_descriptor &m_cone) const {
        // build_tree();  // TODO: rebuild only when necessary
        // Ray_3 ray(ORIGIN, p);
        // std::vector<typename AABBTree::template Intersection_and_primitive_id<Ray_3>::Type> intersections;
        // tree.all_intersections(ray, std::back_inserter(intersections));
        // std::unordered_map<const MetricPolyhedron *, std::pair<FT, metric_face_descriptor>> m_weights;
        // for (auto &[obj, m_face] : intersections) {
        //     auto &[fd, m_ptr] = m_face;
        //     auto isect_pt = boost::get<Point_3>(obj);
        //     FT weight = 1.0 / (isect_pt - ORIGIN).squared_length();
        //     m_weights[m_ptr] = std::make_pair(weight, fd);
        // }

        T d_min = INF;

        for (index_t i = 0; i < sites.size(); ++i) {
            auto [c, m] = site(i);
            Ray_3 ray(ORIGIN, p - c.point);
            auto res = m.tree.any_intersection(ray);
            if (!res) continue;

            auto [obj, fd] = *res;
            auto pm = boost::get<Point_3>(&obj);
            if (!pm) continue;

            FT weight = 1.0 / (*pm - ORIGIN).squared_length();

            // auto it = m_weights.find(&metrics[m_idx].graph);
            // if (it == m_weights.end()) {
            //     continue;
            // }
            auto d = approximate_sqrt((p - c.point).squared_length() * weight);
            if (d < d_min) {
                d_min = d;
                m_cone.site_idx = i;
                m_cone.face = fd;
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

    struct Cone_face_isect_data {
        FT t;
        metric_halfedge_descriptor hd;
        Vector_3 n;
        Sign orient;
    };

    bool isect(const Cone_descriptor &k, const Pline_3 &l, FT &t_min, FT &t_max, metric_halfedge_descriptor_opt &h_min,
               metric_halfedge_descriptor_opt &h_max) const {
        auto [c, m] = site(k.site_idx);
        Pline_3 l_inf(ORIGIN + (l.p - c.point), l.d);
        std::vector<std::tuple<metric_halfedge_descriptor, FT, Sign>> isects;
        for (auto hd : halfedges_around_face(halfedge(k.face, m.graph), m.graph)) {
            auto [v0, v1] = m.cone_face_bases(hd);
            auto n = m.cone_face_orthogonal_vector(hd);
            FT nd = scalar_product(n, l_inf.d);
            if (is_zero(nd)) {
                // TODO: handle the case that the line lies on the face
                continue;
            }

            FT np = scalar_product(n, l_inf.p - ORIGIN);
            FT ti = -np / nd;
            Point_3 pi = l_inf(ti);
            Vector_3 vi = pi - ORIGIN;
            if (!(orientation(n, v0, vi) == POSITIVE && orientation(n, vi, v1) == POSITIVE)) {
                // p_i is outside the 2D cone
                continue;
            }

            isects.emplace_back(hd, ti, sign(nd));
            if (isects.size() == 2) {
                break;
            }
        }
        if (isects.empty()) {
            h_min = h_max = std::nullopt;
            return false;
        } else if (isects.size() == 1) {
            auto &[hd, t, sign] = isects.front();
            // sign: positive: line enters the cone, negative: line leaves the cone
            if (sign == POSITIVE) {
                // positive half line [t, +inf)
                t_min = t;
                h_min = hd;
                t_max = l.t_max;
                h_max = std::nullopt;
            } else {
                // negative half line (-inf, t]
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

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li, metric_halfedge_descriptor_opt &h_min,
               metric_halfedge_descriptor_opt &h_max) const {
        li = l;
        return isect(k, l, li.t_min, li.t_max, h_min, h_max);
    }

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li) const {
        li = l;
        metric_halfedge_descriptor_opt h_min, h_max;
        return isect(k, l, li.t_min, li.t_max, h_min, h_max);
    }

    bool isect(const Cone_descriptor &k0, const Cone_descriptor &k1, const Pline_3 &l, FT &t_min, FT &t_max,
               index_t &k_min, index_t &k_max, metric_halfedge_descriptor_opt &h_min,
               metric_halfedge_descriptor_opt &h_max) const {
        T t0_min, t0_max, t1_min, t1_max;
        metric_halfedge_descriptor_opt h0_min, h0_max, h1_min, h1_max;
        if (!isect(k0, l, t0_min, t0_max, h0_min, h0_max) || !isect(k1, l, t1_min, t1_max, h1_min, h1_max)) {
            return false;
        }
        t_min = std::max(t0_min, t1_min);
        t_max = std::min(t0_max, t1_max);
        k_min = t0_min < t1_min ? 1 : 0;
        k_max = t0_max < t1_max ? 0 : 1;
        h_min = k_min == 0 ? h0_min : h1_min;
        h_max = k_max == 0 ? h0_max : h1_max;
        return t_min < t_max;
    }

    void process_i_trace(const InternalTrace &tr) {
        auto bi_ray = tr.bisect_line;

        // Clip the bisector ray with the face on mesh
        bool has_edge_isect = false;
        mesh_halfedge_descriptor edge_hd;
        Pline_3 edge;
        for (auto hd : halfedges_around_face(tr.face_hd, mesh)) {
            edge = edge_segment<Kernel, T>(hd, mesh, vpm);
            T t_edge, _;
            auto res = CGAL::isect(bi_ray, edge, t_edge, _);
            if (res && !is_zero(t_edge)) {
                has_edge_isect = true;
                edge_hd = hd;
                bi_ray.t_max = t_edge;
                break;
            }
        }
        CGAL_assertion(has_edge_isect);

        // Clip the bisector ray with the cones
        index_t k_min, k_max;
        metric_halfedge_descriptor_opt h_min, h_max;
        bool has_cone_isect = isect(tr.k0, tr.k1, bi_ray, bi_ray.t_min, bi_ray.t_max, k_min, k_max, h_min, h_max);
        CGAL_assertion(has_cone_isect);

        // Find the nearest 3-site bisector intersection
        T dist_min = INF;
        T tb_min;
        Plane_3 bi_plane_min;
        Cone_descriptor k2_min;

        for (index_t site_idx = 0; site_idx < sites.size(); ++site_idx) {
            if (site_idx == tr.k0.site_idx || site_idx == tr.k1.site_idx) {
                continue;
            }

            auto [c, m] = site(site_idx);
            for (auto fd : faces(m.graph)) {
                Cone_descriptor k2{site_idx, fd};
                if (k2 == tr.k_prev) continue;

                Pline_3 s_overlap;
                if (!isect(k2, bi_ray, s_overlap)) continue;

                Plane_3 bi_plane = get_bisect_plane(tr.k0, k2);
                if (bi_plane.is_degenerate()) continue;

                T tb;
                if (!CGAL::isect(s_overlap, bi_plane, tb)) continue;

                auto dist = tb - bi_ray.t_min;
                CGAL_assertion(dist >= 0);

                if (dist < dist_min) {
                    dist_min = dist;
                    tb_min = tb;
                    bi_plane_min = bi_plane;
                    k2_min = k2;
                }
            }
        }

        // Find intersection with the boundary
        // TODO

        if (dist_min < INF) {
        } else if (h_max) {
            // The bisector ray intersects the cone
            auto k0_next = k_max == 0 ? tr.k1 : tr.k0;
            auto k1_prev = k_max == 0 ? tr.k0 : tr.k1;
            auto [c1_prev, m1_prev] = site(k1_prev.site_idx);
            Cone_descriptor k1_next{k1_prev.site_idx, face(opposite(*h_max, m1_prev.graph), m1_prev.graph)};
            Plane_3 bi_plane = get_bisect_plane(k0_next, k1_next);
            auto bisect_dir = cross_product(bi_plane.orthogonal_vector(), tr.face_plane.orthogonal_vector());
            auto n = m1_prev.cone_face_orthogonal_vector(*h_max);
            if (sign(scalar_product(n, bi_ray.d)) != sign(scalar_product(n, bisect_dir))) {
                bisect_dir = -bisect_dir;
            }
            auto bisect_line = Pline_3::ray(bi_ray.p_max(), bisect_dir);
            auto v_vd = vd.add_vertex(bi_ray.p_max());
            auto v_hd = vd.add_ray(v_vd);
            vd.set_target(tr.v_hd, v_vd);
            vd.set_next_double_sided(tr.v_hd, v_hd);
            i_traces.push_back({
                bisect_line,
                bi_plane,
                tr.face_plane,
                edge_hd,
                k0_next,
                k1_next,
                k1_prev,
                v_hd,
            });
        } else {
            // The bisector leaves the face on mesh
            auto edge_next_hd = opposite(edge_hd, mesh);
            auto b_plane = mesh_face_plane(edge_next_hd);
            auto bisect_dir = cross_product(tr.bisect_plane.orthogonal_vector(), b_plane.orthogonal_vector());
            auto orient = orientation(b_plane.orthogonal_vector(), bisect_dir, edge.d);
            auto bisect_line = Pline_3::ray(bi_ray.p_max(), orient == POSITIVE ? bisect_dir : -bisect_dir);
            auto v_vd = vd.add_vertex(bi_ray.p_max());
            auto v_hd = vd.add_ray(v_vd);
            vd.set_target(tr.v_hd, v_vd);
            vd.set_next_double_sided(tr.v_hd, v_hd);
            i_traces.push_back({
                bisect_line,
                tr.bisect_plane,
                b_plane,
                edge_next_hd,
                tr.k0,
                tr.k1,
                tr.k_prev,
                v_hd,
            });
        }
    }
};
}  // namespace SSM_restricted_voronoi_diagram
}  // namespace CGAL
