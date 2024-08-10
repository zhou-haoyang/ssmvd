#ifndef SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_H
#define SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_H

#include <CGAL/basic.h>
#include <SSM_restricted_voronoi_diagram/SSM_restricted_voronoi_diagram_traits.h>
#include <Parametric_line/Parametric_line_3.h>

#include <CGAL/Default.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Origin.h>

#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>

#include <algorithm>
#include <deque>
#include <limits>
#include <optional>
#include <unordered_map>
#include <vector>
#include <iostream>

#define TRAIT_FUNC(ret_type, name, functor)                  \
    template <class... T>                                    \
    static ret_type name(T &&...args) {                      \
        return Traits().functor()(std::forward<T>(args)...); \
    }

namespace CGAL {
namespace SSM_restricted_voronoi_diagram {
template <class Traits, class MeshVertexPointPMap = Default, class MeshFaceIndexPMap = Default,
          class MeshEdgeIndexPMap = Default, class MetricVertexPointPMap = Default, class MetricFaceIndexPMap = Default,
          class VoronoiDiagramVertexPointPMap = Default, class VoronoiDiagramVertexIndexPMap = Default>
class SSM_restricted_voronoi_diagram {
   public:
    using T = Traits::T;
    using FT = Traits::FT;
    using Point_3 = Traits::Point_3;
    using Vector_3 = Traits::Vector_3;
    using Ray_3 = Traits::Ray_3;
    using Plane_3 = Traits::Plane_3;
    using Line_3 = Traits::Line_3;
    using Pline_3 = Traits::Pline_3;

    using Surface_mesh = Traits::Surface_mesh;
    using Metric_polyhedron = Traits::Metric_polyhedron;
    using Voronoi_diagram_graph = Traits::Voronoi_diagram;
    using Metric_AABB_tree = Traits::Metric_AABB_tree;

    using index_t = std::ptrdiff_t;

    using mesh_graph_traits = typename boost::graph_traits<Surface_mesh>;
    using mesh_vertex_descriptor = mesh_graph_traits::vertex_descriptor;
    using mesh_face_descriptor = mesh_graph_traits::face_descriptor;
    using mesh_halfedge_descriptor = mesh_graph_traits::halfedge_descriptor;

    using metric_graph_traits = typename boost::graph_traits<Metric_polyhedron>;
    using metric_face_descriptor = metric_graph_traits::face_descriptor;
    using metric_face_iterator = metric_graph_traits::face_iterator;
    using metric_halfedge_descriptor = metric_graph_traits::halfedge_descriptor;
    using metric_halfedge_descriptor_opt = std::optional<metric_halfedge_descriptor>;
    using metric_edge_descriptor = metric_graph_traits::edge_descriptor;

    using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram_graph>;
    using vd_vertex_descriptor = vd_graph_traits::vertex_descriptor;
    using vd_edge_descriptor = vd_graph_traits::edge_descriptor;
    using vd_halfedge_descriptor = vd_graph_traits::halfedge_descriptor;
    using vd_face_descriptor = vd_graph_traits::face_descriptor;

    using Mesh_vertex_point_pmap =
        Default::Get<MeshVertexPointPMap, typename boost::property_map<Surface_mesh, vertex_point_t>::const_type>::type;
    using Mesh_face_index_pmap =
        Default::Get<MeshFaceIndexPMap, typename boost::property_map<Surface_mesh, face_index_t>::const_type>::type;
    using Mesh_edge_index_pmap =
        Default::Get<MeshEdgeIndexPMap, typename boost::property_map<Surface_mesh, edge_index_t>::const_type>::type;

    using Metric_vertex_point_pmap =
        Default::Get<MetricVertexPointPMap,
                     typename boost::property_map<Metric_polyhedron, vertex_point_t>::const_type>::type;
    using Metric_face_index_pmap =
        Default::Get<MetricFaceIndexPMap,
                     typename boost::property_map<Metric_polyhedron, face_index_t>::const_type>::type;

    using Voronoi_diagram_vertex_point_pmap =
        Default::Get<VoronoiDiagramVertexPointPMap,
                     typename boost::property_map<Voronoi_diagram_graph, vertex_point_t>::const_type>::type;
    using Voronoi_diagram_vertex_index_pmap =
        Default::Get<VoronoiDiagramVertexIndexPMap,
                     typename boost::property_map<Voronoi_diagram_graph, vertex_index_t>::const_type>::type;

    static constexpr T INF = std::numeric_limits<T>::infinity();

    struct Site {
        Point_3 point;
        index_t metric_idx;
    };

    struct Cone_index {
        index_t site_idx;
        size_t face_idx;

        bool operator==(const Cone_index &other) const {
            return site_idx == other.site_idx && face_idx == other.face_idx;
        }
        bool operator>(const Cone_index &other) const {
            return site_idx > other.site_idx || (site_idx == other.site_idx && face_idx > other.face_idx);
        }

        // stream operator
        friend std::ostream &operator<<(std::ostream &os, const Cone_index &ci) {
            os << "Cone(S" << ci.site_idx << ", M" << ci.face_idx << ")";
            return os;
        }
    };

    struct Cone_index_hash {
        std::size_t operator()(const Cone_index &k) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, k.site_idx);
            boost::hash_combine(seed, k.face_idx);
            return seed;
        }
    };

    struct Internal_vertex_id {
        Cone_index k0, k1, k2;
        size_t mesh_face_idx;

        Internal_vertex_id(Cone_index ci0, Cone_index ci1, Cone_index ci2, index_t mesh_face_idx)
            : k0(ci0), k1(ci1), k2(ci2), mesh_face_idx(mesh_face_idx) {
            if (k0 > k1) std::swap(k0, k1);
            if (k1 > k2) std::swap(k1, k2);
            if (k0 > k1) std::swap(k0, k1);
        }

        bool operator==(const Internal_vertex_id &other) const {
            return k0 == other.k0 && k1 == other.k1 && k2 == other.k2 && mesh_face_idx == other.mesh_face_idx;
        }
    };

    struct Internal_vertex_id_hash {
        std::size_t operator()(const Internal_vertex_id &v) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, Cone_index_hash{}(v.k0));
            boost::hash_combine(seed, Cone_index_hash{}(v.k1));
            boost::hash_combine(seed, Cone_index_hash{}(v.k2));
            boost::hash_combine(seed, v.mesh_face_idx);
            return seed;
        }
    };

    struct Boundary_vertex_id {
        Cone_index k0, k1;
        size_t mesh_halfedge_idx;

        Boundary_vertex_id(Cone_index ci0, Cone_index ci1, index_t mesh_halfedge_idx)
            : k0(ci0), k1(ci1), mesh_halfedge_idx(mesh_halfedge_idx) {
            if (k0 > k1) std::swap(k0, k1);
        }

        bool operator==(const Boundary_vertex_id &other) const {
            return k0 == other.k0 && k1 == other.k1 && mesh_halfedge_idx == other.mesh_halfedge_idx;
        }

        friend std::ostream &operator<<(std::ostream &os, const Boundary_vertex_id &v) {
            os << "Boundary(" << v.k0 << ", " << v.k1 << ", E" << v.mesh_halfedge_idx << ")";
            return os;
        }
    };

    struct Boundary_vertex_id_hash {
        std::size_t operator()(const Boundary_vertex_id &v) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, Cone_index_hash{}(v.k0));
            boost::hash_combine(seed, Cone_index_hash{}(v.k1));
            boost::hash_combine(seed, v.mesh_halfedge_idx);
            return seed;
        }
    };

    struct Cone_descriptor {
        index_t site_idx = -1;
        metric_face_descriptor face;

        bool operator==(const Cone_descriptor &other) const { return site_idx == other.site_idx && face == other.face; }

        bool is_valid() const { return site_idx >= 0 && face != metric_graph_traits::null_face(); }
    };

    enum Vertex_type : std::size_t {
        BOUNDARY,
        BOUNDARY_CONE,
        BOUNDARY_BISECTOR,
        TWO_SITE_BISECTOR,
        THREE_SITE_BISECTOR,
    };

    struct Boundary_vertex_info {
        mesh_halfedge_descriptor hd;
        Cone_descriptor k;
    };

    struct Boundary_cone_info {
        mesh_halfedge_descriptor hd;
        Cone_descriptor k;
    };

    struct Boundary_bisector_info {
        mesh_halfedge_descriptor hd;
        Cone_descriptor k0, k1;
    };

    struct Two_site_bisector_info {
        mesh_halfedge_descriptor hd;
        Cone_descriptor k0;
        index_t c1;
        metric_halfedge_descriptor hd1;
    };

    struct Three_site_bisector_info {
        mesh_halfedge_descriptor hd;
        Cone_descriptor k0, k1, k2;
    };

    using Vertex_info = std::variant<Boundary_vertex_info, Boundary_cone_info, Boundary_bisector_info,
                                     Two_site_bisector_info, Three_site_bisector_info>;

    struct Voronoi_diagram_data {
        using Vertex_info_property = CGAL::dynamic_vertex_property_t<Vertex_info>;
        using Vertex_info_map = typename boost::property_map<Voronoi_diagram_graph, Vertex_info_property>::type;
        using Vertex_normal_property = CGAL::dynamic_vertex_property_t<Vector_3>;
        using Vertex_normal_map = typename boost::property_map<Voronoi_diagram_graph, Vertex_normal_property>::type;

        Voronoi_diagram_graph graph;
        Voronoi_diagram_vertex_point_pmap vpm;
        Voronoi_diagram_vertex_index_pmap vertex_index_map;
        Vertex_info_map vertex_info_map;
        Vertex_normal_map vertex_normal_map;
        vd_face_descriptor fd0;

        Voronoi_diagram_data()
            : vpm(get(vertex_point, graph)),
              vertex_index_map(get(vertex_index, graph)),
              vertex_info_map(get(Vertex_info_property{}, graph)),
              vertex_normal_map(get(Vertex_normal_property{}, graph)) {
            fd0 = add_face(graph);
        }

        vd_vertex_descriptor add_vertex(const Point_3 &p, const Vertex_info &info, const Vector_3 &n = {}) {
            auto vd = CGAL::add_vertex(graph);
            put(vpm, vd, p);
            // put(vertex_index_map, vd, num_vertices(graph) - 1);
            put(vertex_info_map, vd, info);
            put(vertex_normal_map, vd, n);
            return vd;
        }

        void print_halfedge_loop(vd_vertex_descriptor v) {
            for (auto hd : halfedges_around_target(v, graph)) {
                std::cerr << hd << "(" << source(hd, graph) << ", " << target(hd, graph) << ")" << " ";
            }
            std::cerr << std::endl;
            std::cerr.flush();
        }

        void insert_halfedge_loop(vd_halfedge_descriptor hd) {
            auto vt = target(hd, graph), vs = source(hd, graph);
            auto hd_cur = halfedge(vt, graph);
            if (hd_cur == vd_graph_traits::null_halfedge()) {
                set_halfedge(vt, hd, graph);
                return;
            }

            auto pt = get(vpm, vt), ps = get(vpm, vs);
            auto v = ps - pt;
            auto n = get(vertex_normal_map, vt);

            if (hd_cur != opposite(next(hd_cur, graph), graph)) {
                auto ori_cur = orientation(n, v, get(vpm, source(hd_cur, graph)) - pt);
                bool found = false;
                do {
                    auto hd_next = opposite(next(hd_cur, graph), graph);
                    auto ori_next = orientation(n, v, get(vpm, source(hd_next, graph)) - pt);
                    if (ori_cur == POSITIVE && ori_next == NEGATIVE) {
                        found = true;
                        break;
                    }
                    hd_cur = hd_next;
                    ori_cur = ori_next;
                } while (hd_cur != halfedge(vt, graph));
                CGAL_assertion(found);
            }

            set_next(hd, next(hd_cur, graph), graph);
            set_next(hd_cur, opposite(hd, graph), graph);
        }

        vd_halfedge_descriptor connect(vd_vertex_descriptor v0, vd_vertex_descriptor v1, vd_face_descriptor fd01 = {},
                                       vd_face_descriptor fd10 = {}) {
            if (halfedge(v1, graph) != vd_graph_traits::null_halfedge()) {
                for (auto hd : halfedges_around_target(v1, graph)) {
                    if (source(hd, graph) == v0) return hd;
                }
            }

            auto ed = CGAL::add_edge(graph);
            auto hd01 = halfedge(ed, graph);
            auto hd10 = opposite(hd01, graph);
            set_target(hd01, v1, graph);
            set_target(hd10, v0, graph);
            set_next(hd01, hd10, graph);
            set_next(hd10, hd01, graph);
            insert_halfedge_loop(hd10);
            insert_halfedge_loop(hd01);

            set_face(hd01, fd01, graph);
            set_face(hd10, fd10, graph);
            return hd01;
        }

        void trace_faces() {
            for (auto hd : halfedges(graph)) {
                if (face(hd, graph) != fd0) continue;
                auto fd = add_face(graph);
                set_halfedge(fd, hd, graph);
                for (auto hd_inner : halfedges_around_face(hd, graph)) {
                    set_face(hd_inner, fd, graph);
                }
            }
            remove_face(fd0, graph);
        }
    };

    struct Internal_trace {
        Pline_3 bisect_line;
        Plane_3 bisect_plane, face_plane;
        mesh_halfedge_descriptor face_hd;
        std::optional<mesh_halfedge_descriptor> prev_hd;
        Cone_descriptor k0, k1, k_prev;
        vd_vertex_descriptor v_vd;
    };

   private:
    struct Metric_data {
        Metric_polyhedron graph;
        Metric_vertex_point_pmap vpm;
        Metric_face_index_pmap face_index_map;
        Metric_AABB_tree tree;

        Metric_data(Metric_polyhedron graph, Metric_vertex_point_pmap vpm, Metric_face_index_pmap fim)
            : graph(std::move(graph)), vpm(std::move(vpm)), face_index_map(std::move(fim)) {
            auto [fbegin, fend] = faces(this->graph);
            tree.insert(fbegin, fend, this->graph);
            tree.build();
        }

        auto cone_face_bases(metric_halfedge_descriptor hd) const {
            auto v0 = get(vpm, source(hd, graph)) - ORIGIN;
            auto v1 = get(vpm, target(hd, graph)) - ORIGIN;
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
    };

   public:
    SSM_restricted_voronoi_diagram(const Surface_mesh &mesh, Mesh_vertex_point_pmap vpm,
                                   Mesh_face_index_pmap face_index_map, Mesh_edge_index_pmap edge_index_map,
                                   Traits traits = Traits())
        : mesh(mesh),
          vpm(std::move(vpm)),
          face_index_map(std::move(face_index_map)),
          edge_index_map(std::move(edge_index_map)),
          traits(traits) {}

    SSM_restricted_voronoi_diagram(const Surface_mesh &mesh)
        : SSM_restricted_voronoi_diagram(mesh, get(vertex_point, mesh), get(face_index, mesh), get(edge_index, mesh)) {}

    void add_site(const Point_3 &p, index_t metric_idx) { sites.push_back({p, metric_idx}); }

    index_t add_metric(Metric_polyhedron m, Metric_vertex_point_pmap vpm, Metric_face_index_pmap fim) {
        auto idx = metrics.size();
        metrics.emplace_back(std::move(m), std::move(vpm), std::move(fim));
        return idx;
    }

    index_t add_metric(Metric_polyhedron m) { return add_metric(m, get(vertex_point, m), get(face_index, m)); }

    void reserve_sites(std::size_t n) { sites.reserve(n); }

    void reserve_metrics(std::size_t n) { metrics.reserve(n); }

    void clear_sites() { sites.clear(); }

    void clear_metrics() { metrics.clear(); }

    std::size_t num_sites() const { return sites.size(); }

    std::size_t num_metrics() const { return metrics.size(); }

    std::size_t num_i_traces() const { return i_traces.size(); }

    auto site_cbegin() const { return sites.begin(); }

    auto site_cend() const { return sites.end(); }

    auto metric_cbegin() const { return metrics.begin(); }

    auto metric_cend() const { return metrics.end(); }

    auto i_trace_cbegin() const { return i_traces.begin(); }

    auto i_trace_cend() const { return i_traces.end(); }

    std::pair<Site &, Metric_data &> site(index_t idx) {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    std::pair<const Site &, const Metric_data &> site(index_t idx) const {
        auto &c = sites[idx];
        return {c, metrics[c.metric_idx]};
    }

    T find_nearest_site(const Point_3 &p, Cone_descriptor &m_cone) const {
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
        CGAL_assertion(d_min < INF);
        return d_min;
    }

    auto trace_boundary(mesh_halfedge_descriptor bh, Cone_descriptor k0, vd_vertex_descriptor prev_vd,
                        bool same_side = true, bool opposite_side = true, bool add_border_edges = false,
                        bool add_cone_vertices = false) {
        auto intersect_f = traits.intersect_3_object();
        auto pline_f = traits.construct_parametric_line_3_object();

        CGAL_precondition(k0.is_valid());
        auto b_line = Pline_3::segment(get(vpm, source(bh, mesh)), get(vpm, target(bh, mesh)));

        auto bh_opposite = opposite(bh, mesh);
        if (is_border(bh, mesh)) {
            same_side = false;
        }
        if (is_border(bh_opposite, mesh)) {
            opposite_side = false;
        }

        Vector_3 b_normal(NULL_VECTOR);
        Plane_3 b_plane, b_plane_opposite;
        if (same_side) {
            b_plane = mesh_face_plane(bh);
            b_normal = b_plane.orthogonal_vector();
        }
        if (opposite_side) {
            b_plane_opposite = mesh_face_plane(bh_opposite);
            b_normal += b_plane_opposite.orthogonal_vector();
        }

        auto isects = find_all_segment_cone_intersections(b_line);

        // Sub-interval of the boundary segment that intersects k0
        auto isect_iter = isects[k0.site_idx].find_face(k0.face);
        CGAL_assertion_msg(isect_iter != isects[k0.site_idx].end(), "Boundary must intersect the cone");

        // Find all boundary-cone / boundary-bisector intersections
        Cone_descriptor k_prev;
        T t_min = 0, t_max = 1;
        for (;;) {
            // Pline_3 b_segment;
            // metric_halfedge_descriptor_opt h_min, h_max;
            // auto res = isect(k0, b_line, b_segment, h_min, h_max);

            // Find nearest intersection of 2-site bisector plane with the boundary segment
            T dist_min = INF;
            T tb_min;
            Plane_3 bi_plane_min;
            Cone_descriptor k1_min;
            for (index_t site_idx = 0; site_idx < sites.size(); ++site_idx) {
                if (site_idx == k0.site_idx) {
                    continue;
                }

                auto isect = *isect_iter;
                auto isects_overlap =
                    isects[site_idx].clip(std::max(isect.second.t_min, t_min), std::min(isect.second.t_max, t_max));

                for (auto [fd, isect_overlap] : isects_overlap) {
                    Cone_descriptor k1{site_idx, fd};
                    if (k1 == k_prev) continue;

                    // Pline_3 b_overlap;
                    // if (!isect(k1, b_segment, b_overlap)) continue;
                    Plane_3 bi_plane = get_bisect_plane(k0, k1);
                    if (bi_plane.is_degenerate()) continue;

                    auto b_overlap = construct_parametric_line(b_line, isect_overlap.t_min, isect_overlap.t_max);

                    auto tb = intersect_f(b_overlap, bi_plane);
                    if (!tb) continue;
                    T dist = *tb - t_min;
                    CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
                    if (dist < dist_min) {
                        dist_min = dist;
                        tb_min = *tb;
                        bi_plane_min = bi_plane;
                        k1_min = k1;
                    }
                }
            }

            if (dist_min < INF) {
                // A 2-site bisector intersects the boundary segment
                // auto edge_hd = opposite(bh, mesh);
                Boundary_vertex_id bvid(cone_index(k0), cone_index(k1_min), get(edge_index_map, edge(bh, mesh)));
                auto pt_start = b_line(tb_min);
                vd_vertex_descriptor v_vd;
                if (auto vd_it = b_vert_map.find(bvid); vd_it != b_vert_map.cend()) {
                    v_vd = vd_it->second;
                } else {
                    v_vd = voronoi.add_vertex(pt_start, Boundary_bisector_info{bh, k0, k1_min}, b_normal);
                    b_vert_map[bvid] = v_vd;
                }

                if (debug_output) {
                    std::clog << "[b-trace] boundary 2-site bisector " << bvid << " intersects at point " << pt_start
                              << std::endl;
                }

                if (add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
                    voronoi.connect(prev_vd, v_vd);
                }
                prev_vd = v_vd;

                if (same_side) {
                    auto bisect_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane.orthogonal_vector());
                    auto orient = orientation(b_plane.orthogonal_vector(), b_line.d(), bisect_dir);
                    auto bisect_line = Pline_3::ray(pt_start, orient == POSITIVE ? bisect_dir : -bisect_dir);

                    // auto v_hd = vd.add_loop(v_vd);
                    i_traces.push_back({
                        bisect_line,
                        bi_plane_min,
                        b_plane,
                        bh,
                        bh,
                        k0,
                        k1_min,
                        {},
                        v_vd,
                    });
                }
                if (opposite_side) {
                    auto bisect_dir =
                        cross_product(bi_plane_min.orthogonal_vector(), b_plane_opposite.orthogonal_vector());
                    auto orient = orientation(b_plane_opposite.orthogonal_vector(), -b_line.d(), bisect_dir);
                    auto bisect_line = Pline_3::ray(pt_start, orient == POSITIVE ? bisect_dir : -bisect_dir);
                    i_traces.push_back({
                        bisect_line,
                        bi_plane_min,
                        b_plane_opposite,
                        bh_opposite,
                        bh_opposite,
                        k0,
                        k1_min,
                        {},
                        v_vd,
                    });
                }
                // auto bi_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane.orthogonal_vector());
                // auto ori = orientation(b_plane.orthogonal_vector(), bi_dir, b_line.d);
                // auto bisect = Pline_3::ray(b_segment(tb_min), ori == POSITIVE ? bi_dir : -bi_dir);

                // Switch current cone to new cone k1
                t_min = tb_min;
                k_prev = k0;
                k0 = k1_min;

                isect_iter = isects[k0.site_idx].find_face(k0.face);
                CGAL_assertion_msg(isect_iter != isects[k0.site_idx].end(), "Boundary must intersect the cone");
            } else {
                // The bondary segment does not intersect leave the cone. Switch to the next boundary segment
                isect_iter++;
                if (isect_iter == isects[k0.site_idx].end()) break;

                // Otherwise switch to the next cone
                auto isect = *isect_iter;
                if (add_cone_vertices) {
                    auto v_vd = voronoi.add_vertex(b_line(isect.second.t_min), Boundary_cone_info{bh, k0}, b_normal);
                    if (add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
                        voronoi.connect(prev_vd, v_vd);
                    }
                    prev_vd = v_vd;
                }

                k_prev = k0;
                // auto [c0, m0] = site(k0.site_idx);
                k0.face = isect.first;
                t_min = isect.second.t_min;

                // b_line = pline_f(b_line, b_segment.t_max());
                // if (b_line.is_point()) break;
            }
        }
        return std::make_pair(k0, prev_vd);
    }

    void trace_all_boundaries(mesh_vertex_descriptor vd, bool add_border_edges = true, bool add_cone_vertices = false) {
        Cone_descriptor k0;
        find_nearest_site(get(vpm, vd), k0);

        using edge_bool_t = CGAL::dynamic_edge_property_t<bool>;
        using edge_visited_map = typename boost::property_map<Surface_mesh, edge_bool_t>::type;

        auto edge_visited = get(edge_bool_t{}, mesh);

        std::vector<std::pair<mesh_halfedge_descriptor, Cone_descriptor>> queue;
        for (auto hd : halfedges_around_source(vd, mesh)) {
            queue.emplace_back(hd, k0);
            // put(edge_visited, edge(hd, mesh), true);
        }

        while (!queue.empty()) {
            auto [hd, k0] = queue.back();
            queue.pop_back();

            if (get(edge_visited, edge(hd, mesh))) continue;

            if (is_border(hd, mesh)) {
                auto bhd = hd;
                auto kb = k0;  // Current cone
                vd_vertex_descriptor prev_vd = vd_graph_traits::null_vertex(), vd0;
                do {
                    if (add_border_edges) {
                        auto v_vd = voronoi.add_vertex(get(vpm, source(bhd, mesh)), Boundary_vertex_info{hd, k0});
                        if (prev_vd != vd_graph_traits::null_vertex()) {
                            voronoi.connect(prev_vd, v_vd);
                        } else {
                            vd0 = v_vd;
                        }
                        prev_vd = v_vd;
                    }

                    std::tie(kb, prev_vd) =
                        trace_boundary(bhd, kb, prev_vd, false, true, add_border_edges, add_cone_vertices);
                    put(edge_visited, edge(bhd, mesh), true);
                    bhd = next(bhd, mesh);
                } while (bhd != hd);
                if (add_border_edges) voronoi.connect(prev_vd, vd0);
            } else if (is_border(opposite(hd, mesh), mesh)) {
                continue;  // We handle the border halfedge in the previous case
            } else {
                auto [k, _] = trace_boundary(hd, k0, vd_graph_traits::null_vertex(), true, true, false, false);
                put(edge_visited, edge(hd, mesh), true);
                for (auto hd_inner : halfedges_around_source(target(hd, mesh), mesh)) {
                    if (get(edge_visited, edge(hd_inner, mesh))) continue;
                    queue.emplace_back(hd_inner, k);
                    // put(edge_visited, edge(hd_inner, mesh), true);
                }
            }
        }
    }

    void trace_all_boundaries(bool add_border_edges = true, bool add_cone_vertices = false) {
        if (is_empty(mesh)) return;
        trace_all_boundaries(*vertices(mesh).first, add_border_edges, add_cone_vertices);
    }

    void reload() { vpm = get(vertex_point, mesh); }

    void reset() {
        voronoi = Voronoi_diagram_data();
        i_traces.clear();
        vert_map.clear();
        b_vert_map.clear();
    }

    bool step() {
        if (i_traces.empty()) {
            return false;
        }
        auto tr = std::move(i_traces.front());
        i_traces.pop_front();
        process_i_trace(tr);
        return true;
    }

    void build(bool add_border_edges = true, bool add_cone_vertices = false) {
        reset();
        trace_all_boundaries(add_border_edges, add_cone_vertices);
        bool stat;
        do {
            stat = step();
        } while (stat);
        voronoi.trace_faces();
        CGAL_postcondition(is_valid_face_graph(voronoi.graph, true));
    }

    const Voronoi_diagram_data &voronoi_diagram() const { return voronoi; }

   private:
    Traits traits;
    bool debug_output = false;
    std::vector<Site> sites;
    std::vector<Metric_data> metrics;

    const Surface_mesh &mesh;
    Mesh_vertex_point_pmap vpm;
    Mesh_face_index_pmap face_index_map;
    Mesh_edge_index_pmap edge_index_map;

    Voronoi_diagram_data voronoi;

    std::deque<Internal_trace> i_traces;
    std::unordered_map<Internal_vertex_id, vd_vertex_descriptor, Internal_vertex_id_hash> vert_map;

    std::unordered_map<Boundary_vertex_id, vd_vertex_descriptor, Boundary_vertex_id_hash> b_vert_map;

    Cone_index cone_index(const Cone_descriptor &k) const {
        auto [c, m] = site(k.site_idx);
        return {k.site_idx, m.face_index_map[k.face]};
    }

    template <class FaceGraph, class VPMap>
    static Plane_3 supporting_plane(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph &g,
                                    const VPMap &vpm) {
        auto i0 = source(hd, g), i1 = target(hd, g), i2 = target(next(hd, g), g);
        auto v0 = get(vpm, i0), v1 = get(vpm, i1), v2 = get(vpm, i2);
        return {v0, v1, v2};
    }

    template <class FaceGraph, class VPMap>
    static Pline_3 edge_segment(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd, const FaceGraph &g,
                                const VPMap &vpm) {
        auto i0 = source(hd, g), i1 = target(hd, g);
        return Pline_3::segment(get(vpm, i0), get(vpm, i1));
    }

    Plane_3 metric_face_plane(const Metric_data &m, metric_face_descriptor fd) const {
        return supporting_plane(halfedge(fd, m.graph), m.graph, m.vpm);
    }

    Plane_3 mesh_face_plane(mesh_halfedge_descriptor hd) const { return supporting_plane(hd, mesh, vpm); }

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
        Pline_3 l_inf(ORIGIN + (l.p() - c.point), l.d());
        std::vector<std::tuple<metric_halfedge_descriptor, FT, Sign>> isects;
        for (auto hd : halfedges_around_face(halfedge(k.face, m.graph), m.graph)) {
            auto [v0, v1] = m.cone_face_bases(hd);
            auto n = m.cone_face_orthogonal_vector(hd);
            FT nd = scalar_product(n, l_inf.d());
            if (is_zero(nd)) {
                // TODO: handle the case that the line lies on the face
                continue;
            }

            FT np = scalar_product(n, l_inf.p() - ORIGIN);
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
                t_max = l.t_max();
                h_max = std::nullopt;
            } else {
                // negative half line (-inf, t]
                t_min = l.t_min();
                h_min = std::nullopt;
                t_max = t;
                h_max = hd;
            }

            if (t_min > t_max) {
                return false;
            }

        } else {
            auto [hd0, t0, sign0] = isects.front();
            auto [hd1, t1, sign1] = isects.back();
            if (t0 > t1) {
                std::swap(t0, t1);
                std::swap(hd0, hd1);
            }
            t_min = t0;
            h_min = hd0;
            t_max = t1;
            h_max = hd1;
        }

        if (t_min < l.t_min()) {
            t_min = l.t_min();
            h_min = std::nullopt;
        }

        if (t_max > l.t_max()) {
            t_max = l.t_max();
            h_max = std::nullopt;
        }

        return t_min <= t_max;
    }

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li, metric_halfedge_descriptor_opt &h_min,
               metric_halfedge_descriptor_opt &h_max) const {
        T t_min, t_max;
        bool res = isect(k, l, t_min, t_max, h_min, h_max);
        if (res) li = traits.construct_parametric_line_3_object()(l, t_min, t_max);
        return res;
    }

    bool isect(const Cone_descriptor &k, const Pline_3 &l, Pline_3 &li) const {
        metric_halfedge_descriptor_opt h_min, h_max;
        return isect(k, l, li, h_min, h_max);
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

    TRAIT_FUNC(Ray_3, construct_ray, construct_ray_3_object)
    TRAIT_FUNC(Point_3, construct_point, construct_point_3_object)
    TRAIT_FUNC(Point_3, construct_point_on, construct_point_on_3_object)
    TRAIT_FUNC(Vector_3, construct_vector, construct_vector_3_object)
    TRAIT_FUNC(Pline_3, construct_parametric_line, construct_parametric_line_3_object)
    TRAIT_FUNC(FT, scalar_product, compute_scalar_product_3_object)
    TRAIT_FUNC(Orientation, orientation, orientation_3_object)

    struct Cone_line_intersection {
        T t_min, t_max;
        metric_halfedge_descriptor ed_min, ed_max;
    };

    struct Segment_cone_intersections;

    class Segment_cone_intersections_iterator {
       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = typename Segment_cone_intersections::value_type;
        using difference_type = std::ptrdiff_t;
        using iterator = Segment_cone_intersections_iterator;

        Segment_cone_intersections_iterator(Segment_cone_intersections *data, std::size_t idx) : data(data), idx(idx) {}

        value_type operator*() const { return (*data)[idx]; }

        iterator &operator++() {
            ++idx;
            return *this;
        }

        iterator operator++(int) {
            iterator tmp = *this;
            ++idx;
            return tmp;
        }

        bool operator==(const iterator &other) const { return data == other.data && idx == other.idx; }

       private:
        Segment_cone_intersections *data = nullptr;
        std::size_t idx;
    };

    struct Segment_cone_intersections {
        using value_type = std::pair<metric_face_descriptor, Cone_line_intersection>;
        using iterator = Segment_cone_intersections_iterator;

        value_type operator[](std::size_t idx) const {
            CGAL_assertion(idx < fds.size());
            return std::make_pair(fds[idx], Cone_line_intersection{ts[idx], ts[idx + 1], eds[idx], eds[idx + 1]});
        }

        iterator begin() { return iterator(this, 0); }

        iterator end() { return iterator(this, fds.size()); }

        bool empty() const { return fds.empty(); }

        std::size_t num_faces() const { return fds.size(); }

        T t_min() const {
            CGAL_assertion(!empty());
            return ts.front();
        }

        T t_max() const {
            CGAL_assertion(!empty());
            return ts.back();
        }

        iterator find_face(metric_face_descriptor fd) {
            auto it = idx_map.find(fd);
            if (it == idx_map.end()) {
                return end();
            }

            return iterator(this, it->second);
        }

        void add_face(metric_face_descriptor fd) {
            idx_map[fd] = fds.size();
            fds.push_back(fd);
        }

        void add_intersection(T t, metric_halfedge_descriptor ed) {
            CGAL_assertion(empty() || t >= t_max());
            ts.push_back(t);
            eds.push_back(ed);
        }

        Segment_cone_intersections clip(T tmin, T tmax) const {
            Segment_cone_intersections res;

            // Empty cases
            if (empty() || tmin > tmax || tmin > t_max() || tmax < t_min()) return res;

            auto it0 = std::lower_bound(ts.begin(), ts.end(), tmin);
            auto it1 = std::upper_bound(ts.begin(), ts.end(), tmax);

            if (tmin < *it0) {
                // tmin is in range of it0 - 1 and it0
                res.ts.push_back(tmin);
                res.eds.push_back(metric_graph_traits::null_halfedge());
                auto i0 = std::distance(ts.begin(), it0);
                CGAL_assertion(i0 > 0);
                res.fds.push_back(fds[i0 - 1]);
            }

            for (auto it = it0; it != it1; ++it) {
                res.ts.push_back(*it);
                auto i = std::distance(ts.begin(), it);
                res.eds.push_back(eds[i]);
                if (*it < tmax) res.fds.push_back(fds[i]);
            }

            if (res.ts.empty() || tmax > res.ts.back()) {
                // tmax is in range of it1 and it1 + 1
                res.ts.push_back(tmax);
                res.eds.push_back(metric_graph_traits::null_halfedge());
            }

            return res;
        }

       private:
        std::vector<T> ts;                            // N+1
        std::vector<metric_halfedge_descriptor> eds;  // N+1
        std::vector<metric_face_descriptor> fds;      // N

        mutable std::unordered_map<metric_face_descriptor, size_t> idx_map;
    };

    auto find_segment_cone_intersections(index_t site_idx, const Pline_3 &segment) const {
        // CGAL_precondition(traits.is_segment_3_object()(segment));
        auto [c, m] = site(site_idx);
        auto pmin = traits.construct_pmin_3_object()(segment);
        auto pmax = traits.construct_pmax_3_object()(segment);
        auto tmin = traits.compute_tmin_3_object()(segment);
        auto tmax = traits.compute_tmax_3_object()(segment);
        auto p = construct_vector(c.point, construct_point(segment));
        auto d = construct_vector(segment);

        if (debug_output) {
            std::clog << __func__ << ": site=" << c.point << ", segment=" << segment << std::endl;
        }

        auto ray0 = construct_ray(ORIGIN, construct_vector(c.point, pmin));
        auto ray1 = construct_ray(ORIGIN, construct_vector(c.point, pmax));

        auto isect0 = m.tree.any_intersected_primitive(ray0);
        auto isect1 = m.tree.any_intersected_primitive(ray1);

        CGAL_assertion(isect0 && isect1);
        auto fd0 = *isect0;
        auto fd1 = *isect1;

        Segment_cone_intersections res;
        res.add_intersection(tmin, metric_graph_traits::null_halfedge());
        res.add_face(fd0);

        metric_halfedge_descriptor hd_prev;
        for (auto fd = fd0; fd != fd1;) {
            bool found = false;
            auto hd0 = halfedge(fd, m.graph);
            for (auto hd : halfedges_around_face(hd0, m.graph)) {
                if (hd == hd_prev) continue;
                auto [v0, v1] = m.cone_face_bases(hd);
                auto n = m.cone_face_orthogonal_vector(hd);

                if (debug_output) {
                    std::clog << __func__ << ": test segment-cone face intersection: v0=" << v0 << ", v1=" << v1
                              << ", n=" << n << std::endl;
                }

                FT nd = scalar_product(n, d);
                if (is_zero(nd)) {
                    // TODO: handle the case that the line lies on the face
                    continue;
                }

                FT np = scalar_product(n, p);
                FT ti = -np / nd;
                if (ti < tmin || ti > tmax) continue;

                Vector_3 vi = p + ti * d;

                if (!(orientation(n, v0, vi) == POSITIVE && orientation(n, vi, v1) == POSITIVE)) {
                    // p_i is outside the 2D cone
                    continue;
                }

                // The segment leaves the cone from edge hd
                hd_prev = opposite(hd, m.graph);
                fd = face(hd_prev, m.graph);

                // res.ts.push_back(ti);
                // res.eds.push_back(hd);
                // res.fds.push_back(fd);

                res.add_face(fd);
                res.add_intersection(ti, hd);
                found = true;
                break;
            }
            CGAL_assertion_msg(found, "No intersection found");
        }
        res.add_intersection(tmax, metric_graph_traits::null_halfedge());
        return res;
    }

    auto find_all_segment_cone_intersections(const Pline_3 &segment) const {
        std::vector<Segment_cone_intersections> res;
        for (index_t site_idx = 0; site_idx < sites.size(); ++site_idx) {
            res.push_back(find_segment_cone_intersections(site_idx, segment));
        }
        return res;
    }

    static auto clip_all_segment_cone_intersections(const std::vector<Segment_cone_intersections> &isects, T tmin,
                                                    T tmax) {
        std::vector<Segment_cone_intersections> res;
        res.reserve(isects.size());
        std::transform(isects.begin(), isects.end(), std::back_inserter(res),
                       [=](const Segment_cone_intersections &isects) { return isects.clip(tmin, tmax); });
        return res;
    }

    void process_i_trace(const Internal_trace &tr) {
        auto intersect_f = traits.intersect_3_object();
        auto pline_f = traits.construct_parametric_line_3_object();

        auto bi_ray = tr.bisect_line;

        // Clip the bisector ray with the face on mesh
        bool has_edge_isect = false;
        mesh_halfedge_descriptor edge_hd;
        Pline_3 edge;
        for (auto hd : halfedges_around_face(tr.face_hd, mesh)) {
            if (tr.prev_hd.has_value() && hd == *tr.prev_hd) continue;
            edge = edge_segment(hd, mesh, vpm);
            // T t_edge, _;
            auto res = intersect_f(bi_ray, edge, true);
            // auto res = CGAL::isect(bi_ray, edge, t_edge, _, true);
            // if (res && !is_zero(t_edge)) {
            if (res) {
                auto [t_edge, _] = *res;
                if (is_zero(t_edge)) continue;
                has_edge_isect = true;
                edge_hd = hd;
                bi_ray = pline_f(bi_ray, -INF, t_edge);
                break;
            }
        }
        CGAL_assertion(has_edge_isect);

        // Clip the bisector ray with the cones
        T t_min, t_max;
        index_t k_min, k_max;
        metric_halfedge_descriptor_opt h_min, h_max;
        bool has_cone_isect = isect(tr.k0, tr.k1, bi_ray, t_min, t_max, k_min, k_max, h_min, h_max);
        bi_ray = pline_f(bi_ray, t_min, t_max);
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

            auto isects = find_segment_cone_intersections(site_idx, bi_ray);
            auto [c, m] = site(site_idx);
            for (auto [fd, isect_overlap] : isects) {
                Cone_descriptor k2{site_idx, fd};
                if (k2 == tr.k_prev) continue;

                // Pline_3 s_overlap;
                // if (!isect(k2, bi_ray, s_overlap)) continue;

                Plane_3 bi_plane = get_bisect_plane(tr.k0, k2);
                if (bi_plane.is_degenerate()) continue;

                auto s_overlap = construct_parametric_line(bi_ray, isect_overlap.t_min, isect_overlap.t_max);

                // T tb;
                auto tb = intersect_f(s_overlap, bi_plane);
                if (!tb) continue;

                auto dist = *tb - bi_ray.t_min();
                CGAL_assertion(dist >= 0);

                if (dist < dist_min) {
                    dist_min = dist;
                    tb_min = *tb;
                    bi_plane_min = bi_plane;
                    k2_min = k2;
                }
            }
        }

        // Find intersection with the boundary
        // TODO

        auto face_id = get(face_index_map, face(tr.face_hd, mesh));
        if (dist_min < INF) {
            // Found a 3-site bisector
            Internal_vertex_id vid(cone_index(tr.k0), cone_index(tr.k1), cone_index(k2_min), face_id);
            if (auto vd_it = vert_map.find(vid); vd_it != vert_map.cend()) {
                voronoi.connect(tr.v_vd, vd_it->second, voronoi.fd0, voronoi.fd0);
                return;
            }
            auto bisect_dir_02 = cross_product(bi_plane_min.orthogonal_vector(), tr.face_plane.orthogonal_vector());
            auto [c0, m0] = site(tr.k0.site_idx);
            auto [c1, m1] = site(tr.k1.site_idx);
            auto p0 = metric_face_plane(m0, tr.k0.face);
            auto p1 = metric_face_plane(m1, tr.k1.face);
            if (is_negative(scalar_product(
                    bisect_dir_02, p1.orthogonal_vector() / abs(p1.d()) - p0.orthogonal_vector() / abs(p0.d())))) {
                bisect_dir_02 = -bisect_dir_02;
            }
            auto pt_start = bi_ray(tb_min);
            auto bisect_line_02 = Pline_3::ray(pt_start, bisect_dir_02);
            auto v_vd = voronoi.add_vertex(pt_start, Three_site_bisector_info{tr.face_hd, tr.k0, tr.k1, k2_min},
                                           tr.face_plane.orthogonal_vector());
            voronoi.connect(tr.v_vd, v_vd, voronoi.fd0, voronoi.fd0);
            vert_map[vid] = v_vd;
            // vd.set_target(tr.v_hd, v_vd);

            // auto v_hd_02 = vd.add_loop(v_vd);
            i_traces.push_back({
                bisect_line_02,
                bi_plane_min,
                tr.face_plane,
                edge_hd,
                std::nullopt,
                tr.k0,
                k2_min,
                tr.k1,
                v_vd,
            });

            auto bisect_plane_12 = get_bisect_plane(tr.k1, k2_min);
            CGAL_assertion(!bisect_plane_12.is_degenerate());

            auto bisect_dir_12 = cross_product(bisect_plane_12.orthogonal_vector(), tr.face_plane.orthogonal_vector());
            if (is_negative(scalar_product(
                    bisect_dir_12, p0.orthogonal_vector() / abs(p0.d()) - p1.orthogonal_vector() / abs(p1.d())))) {
                bisect_dir_12 = -bisect_dir_12;
            }
            auto bisect_line_12 = Pline_3::ray(pt_start, bisect_dir_12);
            // auto v_hd_12 = vd.add_loop(v_vd);
            i_traces.push_back({
                bisect_line_12,
                bisect_plane_12,
                tr.face_plane,
                edge_hd,
                std::nullopt,
                tr.k1,
                k2_min,
                tr.k0,
                v_vd,
            });
        } else if (h_max) {
            // The bisector ray intersects the cone
            auto k0_next = k_max == 0 ? tr.k1 : tr.k0;
            auto k1_prev = k_max == 0 ? tr.k0 : tr.k1;
            auto [c1_prev, m1_prev] = site(k1_prev.site_idx);
            Cone_descriptor k1_next{k1_prev.site_idx, face(opposite(*h_max, m1_prev.graph), m1_prev.graph)};
            Internal_vertex_id vid(cone_index(k0_next), cone_index(k1_next), cone_index(k1_prev), face_id);
            if (auto vd_it = vert_map.find(vid); vd_it != vert_map.cend()) {
                voronoi.connect(tr.v_vd, vd_it->second, voronoi.fd0, voronoi.fd0);
                return;
            }

            Plane_3 bi_plane = get_bisect_plane(k0_next, k1_next);
            auto bisect_dir = cross_product(bi_plane.orthogonal_vector(), tr.face_plane.orthogonal_vector());
            auto n = m1_prev.cone_face_orthogonal_vector(*h_max);
            if (sign(scalar_product(n, bi_ray.d())) != sign(scalar_product(n, bisect_dir))) {
                bisect_dir = -bisect_dir;
            }
            auto bisect_line = Pline_3::ray(bi_ray.p_max(), bisect_dir);
            auto v_vd = voronoi.add_vertex(bi_ray.p_max(),
                                           Two_site_bisector_info{tr.face_hd, k0_next, k1_next.site_idx, *h_max});
            voronoi.connect(tr.v_vd, v_vd, voronoi.fd0, voronoi.fd0);
            vert_map[vid] = v_vd;
            // auto v_hd = vd.add_loop(v_vd);
            // vd.set_target(tr.v_hd, v_vd);
            // vd.set_next_double_sided(tr.v_hd, v_hd);
            i_traces.push_back({
                bisect_line,
                bi_plane,
                tr.face_plane,
                edge_hd,
                std::nullopt,
                k0_next,
                k1_next,
                k1_prev,
                v_vd,
            });
        } else {
            // The bisector leaves the face on mesh
            Boundary_vertex_id bvid(cone_index(tr.k0), cone_index(tr.k1),
                                    get(edge_index_map, CGAL::edge(edge_hd, mesh)));
            if (auto vd_it = b_vert_map.find(bvid); vd_it != b_vert_map.cend()) {
                voronoi.connect(tr.v_vd, vd_it->second, voronoi.fd0, voronoi.fd0);
                return;
            }
            auto edge_next_hd = opposite(edge_hd, mesh);
            auto b_plane = mesh_face_plane(edge_next_hd);
            auto bisect_dir = cross_product(tr.bisect_plane.orthogonal_vector(), b_plane.orthogonal_vector());
            auto orient = orientation(b_plane.orthogonal_vector(), bisect_dir, edge.d());
            auto bisect_line = Pline_3::ray(bi_ray.p_max(), orient == POSITIVE ? bisect_dir : -bisect_dir);
            auto v_vd = voronoi.add_vertex(bi_ray.p_max(), Boundary_bisector_info{tr.face_hd, tr.k0, tr.k1});
            voronoi.connect(tr.v_vd, v_vd, voronoi.fd0, voronoi.fd0);
            b_vert_map[bvid] = v_vd;
            // auto v_hd = vd.add_loop(v_vd);
            // vd.set_target(tr.v_hd, v_vd);
            // vd.set_next_double_sided(tr.v_hd, v_hd);

            if (is_border(edge_next_hd, mesh)) {
                return;
            }
            i_traces.push_back({
                bisect_line,
                tr.bisect_plane,
                b_plane,
                edge_next_hd,
                edge_next_hd,
                tr.k0,
                tr.k1,
                {},
                v_vd,
            });
        }
    }
};
}  // namespace SSM_restricted_voronoi_diagram
}  // namespace CGAL

#endif  // SSM_RESTRICTED_VORONOI_DIAGRAM_SSM_RESTRICTED_VORONOI_DIAGRAM_H
