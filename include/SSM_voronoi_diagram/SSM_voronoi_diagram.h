#ifndef SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H
#define SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H

#include <CGAL/Default.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Origin.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

#include <boost/container_hash/hash.hpp>
#include <boost/graph/graph_traits.hpp>
#include <deque>
#include <list>
#include <optional>
#include <unordered_map>
#include <variant>

#define TRAIT_FUNC(ret_type, name, functor)                  \
    template <class... T>                                    \
    ret_type name(T&&... args) const {                       \
        return m_traits.functor()(std::forward<T>(args)...); \
    }

#define PROPERTY(type, name)        \
    const type& name() const {      \
        return m_##name;            \
    }                               \
                                    \
    type& name() {                  \
        return m_##name;            \
    }                               \
                                    \
    void set_##name(type name) {    \
        m_##name = std::move(name); \
    }

namespace CGAL::SSM_voronoi_diagram {
template <class Traits, class VoronoiDiagramVertexPointPMap = Default, class VoronoiDiagramVertexIndexPMap = Default>
class SSM_voronoi_diagram {
#pragma region public types
   public:
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;
    using Polygon_2 = typename Traits::Polygon_2;
    using Parametric_line_2 = typename Traits::Parametric_line_2;

    using Metric = typename Traits::Metric;
    using Metric_traits = typename Traits::Metric_traits;
    using Metric_edge_iterator = typename Metric::Edge_const_iterator;

    class Metric_data {
       public:
        Metric_data(Metric m, Metric_traits traits) : m_polygon(std::move(m)), m_traits(std::move(traits)) {}

        auto& polygon() const { return m_polygon; }
        auto& traits() const { return m_traits; }

        auto any_intersection(const Vector_2& d) { return _any_intersection(m_polygon, d); }

        auto any_intersected_edge(const Vector_2& d) { return _any_intersected_edge(m_polygon, d); }

        auto index(Metric_edge_iterator ed) const { return std::distance(m_polygon.edges_begin(), ed); }

        Metric_edge_iterator next_edge(Metric_edge_iterator ed) const {
            return ++ed == m_polygon.edges_end() ? m_polygon.edges_begin() : ed;
        }

        Metric_edge_iterator prev_edge(Metric_edge_iterator ed) const {
            return ed == m_polygon.edges_begin() ? --m_polygon.edges_end() : --ed;
        }

       private:
        Metric m_polygon;
        Metric_traits m_traits;

        TRAIT_FUNC(auto, _any_intersection, metric_any_intersection_object)
        TRAIT_FUNC(auto, _any_intersected_edge, metric_any_intersected_edge_object)
    };

    using Metric_list = std::list<Metric_data>;
    using Metric_iterator = typename Metric_list::iterator;

    class Site {
       public:
        Site(Point_2 p, Metric_iterator m) : m_point(std::move(p)), m_metric(std::move(m)) {}

        PROPERTY(Point_2, point)
        PROPERTY(Metric_iterator, metric)

       private:
        Point_2 m_point;
        Metric_iterator m_metric;
    };

    using Site_list = std::list<Site>;
    using Site_const_iterator = typename Site_list::const_iterator;

    class Cone_descriptor {
       public:
        Cone_descriptor(Site_const_iterator s = {}, Metric_edge_iterator e = {})
            : m_site(std::move(s)), m_edge(std::move(e)) {}

        PROPERTY(Site_const_iterator, site)
        PROPERTY(Metric_edge_iterator, edge)

        bool operator==(const Cone_descriptor& other) const { return m_site == other.m_site && m_edge == other.m_edge; }

       private:
        Site_const_iterator m_site;
        Metric_edge_iterator m_edge;
    };

    using Boundary = Polygon_2;
    using Boundary_edge_iterator = typename Boundary::Edge_const_iterator;

    using Voronoi_diagram_graph = typename Traits::Voronoi_diagram;
    using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram_graph>;
    using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
    using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
    using vd_halfedge_descriptor = typename vd_graph_traits::halfedge_descriptor;
    using vd_face_descriptor = typename vd_graph_traits::face_descriptor;

    enum Vertex_type : std::size_t {
        BOUNDARY,
        BOUNDARY_CONE,
        BOUNDARY_BISECTOR,
        TWO_SITE_BISECTOR,
        THREE_SITE_BISECTOR,
    };

    struct Boundary_vertex_info {
        Boundary_edge_iterator ed;
        Cone_descriptor k;
    };

    struct Boundary_cone_info {
        Boundary_edge_iterator ed;
        Cone_descriptor k;
    };

    struct Boundary_bisector_info {
        Boundary_edge_iterator ed;
        Cone_descriptor k0, k1;
    };

    struct Two_site_bisector_info {
        Cone_descriptor k0;
        Site_const_iterator c1;
        Metric_edge_iterator e1;
    };

    struct Three_site_bisector_info {
        Cone_descriptor k0, k1, k2;
    };

    using Vertex_info = std::variant<Boundary_vertex_info, Boundary_cone_info, Boundary_bisector_info,
                                     Two_site_bisector_info, Three_site_bisector_info>;

    using Voronoi_diagram_vertex_point_pmap =
        typename Default::Get<VoronoiDiagramVertexPointPMap,
                              typename boost::property_map<Voronoi_diagram_graph, vertex_point_t>::const_type>::type;
    using Voronoi_diagram_vertex_index_pmap =
        typename Default::Get<VoronoiDiagramVertexIndexPMap,
                              typename boost::property_map<Voronoi_diagram_graph, vertex_index_t>::const_type>::type;

    struct Voronoi_diagram {
        using Vertex_info_property = CGAL::dynamic_vertex_property_t<Vertex_info>;
        using Vertex_info_map = typename boost::property_map<Voronoi_diagram_graph, Vertex_info_property>::type;

        Voronoi_diagram_graph graph;
        Voronoi_diagram_vertex_point_pmap vpm;
        Voronoi_diagram_vertex_index_pmap vertex_index_map;
        Vertex_info_map vertex_info_map;

        Voronoi_diagram(const Traits& traits)
            : vpm(get(vertex_point, graph)),
              vertex_index_map(get(vertex_index, graph)),
              vertex_info_map(get(Vertex_info_property{}, graph)),
              m_traits(traits) {}

        /**
         * @brief The copy constructor is deleted to avoid the property map ownership issue.
         */
        Voronoi_diagram(const Voronoi_diagram&) = delete;

        vd_vertex_descriptor add_vertex(const Point_2& p, const Vertex_info& info) {
            auto vd = CGAL::add_vertex(graph);
            put(vpm, vd, p);
            // put(vertex_index_map, vd, num_vertices(graph) - 1);
            put(vertex_info_map, vd, info);
            return vd;
        }

        Vector_2 normalized(const Vector_2& v) const { return v / approximate_sqrt(v.squared_length()); }

        void insert_halfedge_loop(vd_halfedge_descriptor hd) {
            auto vt = target(hd, graph), vs = source(hd, graph);
            auto hd_cur = halfedge(vt, graph);
            if (hd_cur == vd_graph_traits::null_halfedge()) {
                set_halfedge(vt, hd, graph);
                return;
            }

            if (hd_cur != opposite(next(hd_cur, graph), graph)) {
                // At least 2 halfedges around the target vertex
                auto pt = get(vpm, vt), ps = get(vpm, vs);
                auto v = ps - pt;

                auto v_cur = normalized(get(vpm, source(hd_cur, graph)) - pt);

                for (;;) {
                    auto hd_next = opposite(next(hd_cur, graph), graph);
                    auto v_next = normalized(get(vpm, source(hd_next, graph)) - pt);

                    // A monotonic angle function in [-2, 2]
                    auto angle = [this](const auto& v1, const auto& v2) {
                        auto cos_theta = scalar_product(v1, v2);
                        auto ori = orientation(v1, v2);
                        if (ori == ZERO) {
                            return cos_theta < 0 ? cos_theta + 1 : -cos_theta - 1;
                        } else {
                            // The halfedge loop around the vertex should be clockwise for the halfedge loop around the
                            // face to be counter-clockwise, hence the normal should point inward here
                            return ori == POSITIVE ? cos_theta + 1 : -cos_theta - 1;
                        }
                    };

                    if (angle(v_cur, v) < angle(v_cur, v_next)) break;

                    hd_cur = hd_next;
                    v_cur = v_next;
                }
            }

            set_next(hd, next(hd_cur, graph), graph);
            set_next(hd_cur, opposite(hd, graph), graph);
        }

        vd_halfedge_descriptor connect(vd_vertex_descriptor v0, vd_vertex_descriptor v1, vd_face_descriptor fd01,
                                       vd_face_descriptor fd10) {
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

       protected:
        TRAIT_FUNC(FT, scalar_product, compute_scalar_product_2_object)
        TRAIT_FUNC(Orientation, orientation, orientation_2_object)

        const Traits& m_traits;
    };

    using Voronoi_diagram_ptr = std::shared_ptr<Voronoi_diagram>;
    using Const_voronoi_diagram_ptr = std::shared_ptr<const Voronoi_diagram>;
#pragma endregion

#pragma region protected types
   protected:
    using index_t = std::ptrdiff_t;

    struct Cone_index {
        index_t site_idx, edge_idx;

        Cone_index(index_t site_idx, index_t edge_idx) : site_idx(site_idx), edge_idx(edge_idx) {}

        bool operator==(const Cone_index& other) const {
            return site_idx == other.site_idx && edge_idx == other.edge_idx;
        }

        bool operator>(const Cone_index& other) const {
            return site_idx > other.site_idx || (site_idx == other.site_idx && edge_idx > other.edge_idx);
        }
    };

    struct Cone_index_hash {
        std::size_t operator()(const Cone_index& k) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, k.site_idx);
            boost::hash_combine(seed, k.edge_idx);
            return seed;
        }
    };

    struct Internal_vertex_id {
        Cone_index k0, k1, k2;

        Internal_vertex_id(Cone_index ci0, Cone_index ci1, Cone_index ci2) : k0(ci0), k1(ci1), k2(ci2) {
            if (k0 > k1) std::swap(k0, k1);
            if (k1 > k2) std::swap(k1, k2);
            if (k0 > k1) std::swap(k0, k1);
        }

        bool operator==(const Internal_vertex_id& other) const {
            return k0 == other.k0 && k1 == other.k1 && k2 == other.k2;
        }
    };

    struct Internal_vertex_id_hash {
        std::size_t operator()(const Internal_vertex_id& v) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, Cone_index_hash{}(v.k0));
            boost::hash_combine(seed, Cone_index_hash{}(v.k1));
            boost::hash_combine(seed, Cone_index_hash{}(v.k2));
            return seed;
        }
    };

    struct Boundary_vertex_id {
        Cone_index k0, k1;

        Boundary_vertex_id(Cone_index ci0, Cone_index ci1) : k0(ci0), k1(ci1) {
            if (k0 > k1) std::swap(k0, k1);
        }

        bool operator==(const Boundary_vertex_id& other) const { return k0 == other.k0 && k1 == other.k1; }

        friend std::ostream& operator<<(std::ostream& os, const Boundary_vertex_id& v) {
            os << "Boundary(" << v.k0 << ", " << v.k1 << ", E" << v.mesh_halfedge_idx << ")";
            return os;
        }
    };

    struct Boundary_vertex_id_hash {
        std::size_t operator()(const Boundary_vertex_id& v) const {
            std::size_t seed = 0;
            boost::hash_combine(seed, Cone_index_hash{}(v.k0));
            boost::hash_combine(seed, Cone_index_hash{}(v.k1));
            return seed;
        }
    };

    struct Cone_line_intersection {
        FT tmin;
        std::optional<FT> tmax;
    };

    struct Intervals;

    class Intervals_iterator {
       public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = typename Intervals::value_type;
        using difference_type = std::ptrdiff_t;
        using iterator = Intervals_iterator;

        Intervals_iterator(Intervals* data, std::size_t idx) : data(data), idx(idx) {}

        value_type operator*() const { return (*data)[idx]; }

        iterator& operator++() {
            ++idx;
            return *this;
        }

        iterator operator++(int) {
            iterator tmp = *this;
            ++idx;
            return tmp;
        }

        bool operator==(const iterator& other) const { return data == other.data && idx == other.idx; }

       private:
        Intervals* data = nullptr;
        std::size_t idx;
    };

    struct Intervals {
        using value_type = std::pair<Metric_edge_iterator, Cone_line_intersection>;
        using iterator = Intervals_iterator;

        value_type operator[](std::size_t idx) const {
            CGAL_assertion(is_valid() && idx < num_cones());
            return std::make_pair(
                eds[idx],
                Cone_line_intersection{ts[idx], idx + 1 >= ts.size() ? std::nullopt : std::make_optional(ts[idx + 1])});
        }

        iterator begin() { return iterator(this, 0); }

        iterator end() { return iterator(this, eds.size()); }

        bool empty() const { return eds.empty(); }

        bool is_infinite() const { return ts.size() == eds.size(); }

        bool is_valid() const { return is_infinite() || ts.size() == eds.size() + 1; }

        std::size_t num_cones() const { return eds.size(); }

        std::size_t num_endpoints() const { return ts.size(); }

        auto& edges() const { return eds; }

        auto& endpoints() const { return ts; }

        FT t_min() const {
            CGAL_assertion(!empty());
            return ts.front();
        }

        std::optional<FT> t_max() const {
            CGAL_assertion(!empty());
            return is_infinite() ? std::nullopt : std::make_optional(ts.back());
        }

        iterator find_cone(Metric_edge_iterator ed) {
            // if (idx_map.empty()) {
            //     for (std::size_t i = 0; i < eds.size(); ++i) {
            //         idx_map[eds[i]] = i;
            //     }
            // }

            // auto it = idx_map.find(ed);
            // if (it == idx_map.end()) {
            //     return end();
            // }

            auto it = std::find(eds.begin(), eds.end(), ed);

            return iterator(this, std::distance(eds.begin(), it));
        }

        void clear() {
            ts.clear();
            eds.clear();
            // idx_map.clear();
        }

        void add_edge(Metric_edge_iterator ed) { eds.push_back(ed); }

        void add_intersection(FT t) {
            CGAL_assertion(empty() || t >= t_max());
            ts.push_back(t);
        }

        void clip(FT tmin, FT tmax, Intervals& res) const {
            res.clear();

            // Empty cases
            if (empty() || tmin > tmax || tmin > t_max() || tmax < t_min()) return;

            auto it0 = std::lower_bound(ts.begin(), ts.end(), tmin);
            auto it1 = std::upper_bound(ts.begin(), ts.end(), tmax);

            if (tmin < *it0) {
                // tmin is in range of it0 - 1 and it0
                res.add_intersection(tmin);
                auto i0 = std::distance(ts.begin(), it0);
                CGAL_assertion(i0 > 0);
                res.add_edge(eds[i0 - 1]);
            }

            for (auto it = it0; it != it1; ++it) {
                auto i = std::distance(ts.begin(), it);
                res.add_intersection(*it);
                if (*it < tmax) res.add_edge(eds[i]);
            }

            if (res.ts.empty() || tmax > res.ts.back()) {
                // tmax is in range of it1 and it1 + 1
                res.add_intersection(tmax);
            }
        }

       private:
        std::vector<FT> ts;                     // N+1
        std::vector<Metric_edge_iterator> eds;  // N

        // mutable std::unordered_map<Metric_edge_iterator, size_t> idx_map;
    };

    enum Endpoint_type { CCW = 1, SITE = 0, CW = -1, UNKNOWN = -2 };

    struct Internal_trace {
        Parametric_line_2 bisector;
        Cone_descriptor k0, k1;
        std::optional<Cone_descriptor> k_prev;
        Endpoint_type prev_type;
        vd_vertex_descriptor v_vd;

        Internal_trace(Parametric_line_2 bi, Cone_descriptor k0, Cone_descriptor k1,
                       std::optional<Cone_descriptor> k_prev, vd_vertex_descriptor v_prev)
            : bisector(std::move(bi)), k0(std::move(k0)), k1(std::move(k1)), k_prev(std::move(k_prev)), v_vd(v_prev) {}
    };

    struct Boundary_trace {
        vd_vertex_descriptor v_vd;
        Parametric_line_2 b_line;

        Boundary_trace(vd_vertex_descriptor v, Parametric_line_2 b) : v_vd(v), b_line(std::move(b)) {}
    };

#pragma endregion

#pragma region public methods
   public:
    SSM_voronoi_diagram(const Polygon_2& boundary, Traits traits = {})
        : m_boundary(boundary), m_traits(std::move(traits)) {}

    Metric_iterator add_metric(Metric m) { return m_metrics.emplace(m_metrics.end(), std::move(m), m_traits); }

    Site_const_iterator add_site(Point_2 p, Metric_iterator m) {
        return m_sites.emplace(m_sites.end(), std::move(p), m);
    }

    FT find_nearest_site(const Point_2& p, Cone_descriptor& cone) const {
        FT d_min;
        bool init = true;

        for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it) {
            auto res = site_it->metric()->any_intersection(construct_vector(site_it->point(), p));
            if (!res) continue;
            auto [pm, ed] = *res;
            FT weight = 1.0 / squared_distance(Point_2(ORIGIN), pm);

            auto d = approximate_sqrt(squared_distance(p, site_it->point()) * weight);
            if (d < d_min || init) {
                d_min = d;
                cone.set_site(site_it);
                cone.set_edge(ed);
            }
        }

        return d_min;
    }

    auto trace_boundary(Boundary_edge_iterator b_ed, Cone_descriptor k0, vd_vertex_descriptor prev_vd,
                        bool add_border_edges = false, bool add_cone_vertices = false,
                        vd_face_descriptor fd0 = vd_graph_traits::null_face(),
                        vd_face_descriptor fd1 = vd_graph_traits::null_face()) {
        auto b_line = construct_parametric_line(*b_ed);
        FT tmin = 0, tmax = 1;

        // Intervals of boundary segment clipped by cones of all sites
        find_all_intervals(b_line, m_intervals_vec_cache, tmin, tmax);
        auto& intervals_map = m_intervals_vec_cache;
        auto& intervals = intervals_map[site_index(k0.site())];
        auto interval_it = intervals.find_cone(k0.edge());
        CGAL_assertion_msg(interval_it != intervals.end(), "Initial cone k0 not found in intervals");

        Cone_descriptor k_prev;
        for (;;) {
            // Find the nearest two site bisector
            FT dist_min, tb_min;
            Line_2 bi_line_min;
            Cone_descriptor k1_min;
            bool found = false;
            for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it) {
                if (site_it == k0.site()) continue;

                auto interval = *interval_it;
                FT tmin_overlap = max(interval.tmin, tmin);
                FT tmax_overlap = min(interval.tmax, tmax);
                intervals_map[site_index(site_it)].clip(tmin_overlap, tmax_overlap, m_intervals_cache);
                auto& intervals_overlap = m_intervals_cache;

                for (auto [ed, interval_overlap] : intervals_overlap) {
                    Cone_descriptor k1(site_it, ed);
                    if (k1 == k_prev) continue;

                    Line_2 bi_line = get_bisector(k0, k1);
                    if (is_degenerate(bi_line)) continue;

                    auto isect = intersect(b_line, bi_line);
                    if (!isect || std::holds_alternative<Colinear>(*isect)) continue;

                    FT tb = std::get<FT>(*isect);
                    // TODO: marginal case: tb is on the terminal points
                    if (tb <= interval_overlap.tmin || tb > interval_overlap.tmax) continue;

                    FT dist = tb - tmin;
                    CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
                    if (dist < dist_min || !found) {
                        dist_min = dist;
                        tb_min = tb;
                        bi_line_min = bi_line;
                        k1_min = k1;
                        found = true;
                        break;  // goto next site
                    } else if (dist == dist_min) {
                        CGAL_assertion_msg(false, "TODO: handle multiple bisectors");
                    }
                }
            }

            if (found) {
                // Found a bisector
                auto& k1 = k1_min;
                Boundary_vertex_id v_id(cone_index(k0), cone_index(k1));

                Point_2 p = construct_point_on(b_line, tb_min);
                auto v_vd = m_voronoi->add_vertex(p, Boundary_bisector_info{b_ed, k0, k1});
                if (add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
                    m_voronoi->connect(prev_vd, v_vd, fd0, fd1);
                }
                prev_vd = v_vd;

                Vector_2 d = construct_vector(bi_line_min);
                if (orientation(construct_vector(b_line), d) == RIGHT_TURN) {
                    d = opposite_vector(d);
                }
                auto bisector = construct_parametric_line(p, d);
                m_i_traces.emplace_back(bisector, k0, k1, std::nullopt, v_vd);

                // Switch current cone to k1
                k_prev = k0;
                k0 = k1;
                tmin = tb_min;

                auto& intervals = intervals_map[site_index(k0.site())];
                interval_it = intervals.find_cone(k0.edge());
                CGAL_assertion(interval_it != intervals.end());
            } else {
                // No bisector found, terminate if this is the last interval
                interval_it++;
                auto& intervals = intervals_map[site_index(k0.site())];
                if (interval_it == intervals.end()) {
                    break;
                }

                // Otherwise there is a boundary-cone intersection
                // Otherwise switch to the next cone
                auto [m_ed, isect] = *interval_it;
                if (add_cone_vertices) {
                    Point_2 p = construct_point_on(b_line, isect.tmin);
                    auto v_vd = m_voronoi->add_vertex(p, Boundary_cone_info{b_ed, k0});
                    if (add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
                        m_voronoi->connect(prev_vd, v_vd, fd0, fd1);
                    }
                    prev_vd = v_vd;
                }

                k_prev = k0;
                k0.set_edge(m_ed);
                tmin = isect.tmin;
            }
        }

        return std::make_pair(k0, prev_vd);
    }

    void trace_all_boundaries(Boundary_edge_iterator b_ed0, Boundary_edge_iterator b_ed1, bool add_border_edges = true,
                              bool add_cone_vertices = false) {
        auto b_ed = b_ed0;
        Cone_descriptor k0;
        find_nearest_site(construct_source(*b_ed), k0);

        vd_vertex_descriptor prev_vd = vd_graph_traits::null_vertex(), vd0;

        for (; b_ed != b_ed1; ++b_ed) {
            if (add_border_edges) {
                auto v_vd = m_voronoi->add_vertex(construct_source(*b_ed), Boundary_vertex_info{b_ed, k0});
                if (prev_vd != vd_graph_traits::null_vertex()) {
                    m_voronoi->connect(prev_vd, v_vd, m_dummy_face, vd_graph_traits::null_face());
                } else {
                    vd0 = v_vd;
                }
                prev_vd = v_vd;
            }
            std::tie(k0, prev_vd) = trace_boundary(b_ed, k0, prev_vd, add_border_edges, add_cone_vertices, m_dummy_face,
                                                   vd_graph_traits::null_face());
        }
        if (add_border_edges) {
            m_voronoi->connect(prev_vd, vd0, m_dummy_face, vd_graph_traits::null_face());
        }
    }

    void trace_all_boundaries(bool add_border_edges = true, bool add_cone_vertices = false) {
        auto b_ed0 = m_boundary.edges_begin();
        auto b_ed1 = m_boundary.edges_end();
        trace_all_boundaries(b_ed0, b_ed1, add_border_edges, add_cone_vertices);
    }

    void reset() {
        m_voronoi = std::make_shared<Voronoi_diagram>(m_traits);
        m_dummy_face = add_face(m_voronoi->graph);

        m_i_traces.clear();
        m_i_vertices.clear();
        m_b_traces.clear();

        // if (b_trace_timer.is_running()) b_trace_timer.stop();
        // if (i_trace_timer.is_running()) i_trace_timer.stop();
        // b_trace_timer.reset();
        // i_trace_timer.reset();
    }

    bool step() {
        if (m_i_traces.empty()) {
            return false;
        }
        auto tr = std::move(m_i_traces.front());  // DFS
        m_i_traces.pop_front();
        // i_trace_timer.start();
        process_i_trace(tr);
        // i_trace_timer.stop();
        return true;
    }

    void build(bool add_border_edges = true, bool add_cone_vertices = false) {
        reset();
        trace_all_boundaries(add_border_edges, add_cone_vertices);
        bool stat;
        do {
            stat = step();
        } while (stat);
        trace_faces();
        CGAL_postcondition(is_valid_face_graph(m_voronoi->graph, true));

        // vout << IO::level(1) << "profiling: boundary trace: count = " << b_trace_timer.intervals()
        //      << ", time = " << b_trace_timer.time() << "s, speed = " << b_trace_timer.time() /
        //      b_trace_timer.intervals()
        //      << "s/trace" << std::endl;
        // vout << IO::level(1) << "profiling: internal trace: count = " << i_trace_timer.intervals()
        //      << ", time = " << i_trace_timer.time() << "s, speed = " << i_trace_timer.time() /
        //      i_trace_timer.intervals()
        //      << "s/trace" << std::endl;
    }
#pragma endregion

#pragma region protected methods
   protected:
    const Polygon_2& m_boundary;
    Metric_list m_metrics;
    Site_list m_sites;
    Traits m_traits;
    Voronoi_diagram_ptr m_voronoi;

    std::deque<Internal_trace> m_i_traces;
    std::unordered_map<Internal_vertex_id, vd_vertex_descriptor, Internal_vertex_id_hash> m_i_vertices;
    std::unordered_multimap<Boundary_vertex_id, Boundary_trace, Boundary_vertex_id_hash> m_b_traces;

    vd_face_descriptor m_dummy_face;

    Intervals m_intervals_cache;
    std::vector<Intervals> m_intervals_vec_cache;

    TRAIT_FUNC(Vector_2, construct_vector, construct_vector_2_object)
    TRAIT_FUNC(Point_2, construct_point, construct_point_2_object)
    TRAIT_FUNC(Point_2, construct_point_on, construct_point_on_2_object)
    TRAIT_FUNC(Point_2, construct_source, construct_source_2_object)
    TRAIT_FUNC(Point_2, construct_target, construct_target_2_object)
    TRAIT_FUNC(Parametric_line_2, construct_parametric_line, construct_parametric_line_2_object)
    TRAIT_FUNC(Line_2, construct_line, construct_line_2_object)
    TRAIT_FUNC(Vector_2, opposite_vector, construct_opposite_vector_2_object)
    TRAIT_FUNC(auto, intersect, intersect_2_object)

    TRAIT_FUNC(FT, squared_distance, compute_squared_distance_2_object)
    TRAIT_FUNC(FT, determinant, compute_determinant_2_object)
    TRAIT_FUNC(FT, ca, compute_a_2_object)
    TRAIT_FUNC(FT, cb, compute_b_2_object)
    TRAIT_FUNC(FT, cc, compute_c_2_object)
    TRAIT_FUNC(FT, cx, compute_x_2_object)
    TRAIT_FUNC(FT, cy, compute_y_2_object)
    TRAIT_FUNC(FT, scalar_product, compute_scalar_product_2_object)

    TRAIT_FUNC(bool, is_degenerate, is_degenerate_2_object)

    struct Colinear {
        // operator =
        bool operator==(const Colinear&) const { return true; }
    };  // namespace CGAL::SSM_voronoi_diagram

    using Parametric_line_intersection = std::variant<std::pair<FT, FT>, Colinear>;

    std::size_t site_index(Site_const_iterator site) const { return std::distance(m_sites.cbegin(), site); }

    Cone_index cone_index(Cone_descriptor k) const {
        return Cone_index(site_index(k.site()), k.site()->metric()->index(k.edge()));
    }

    /**
     * @brief Find the intersection of a parametric line with a ray with direction d and starting at the origin.
     *
     * @param l
     * @param d
     * @return std::optional<FT>
     */
    std::optional<Parametric_line_intersection> intersect_ray(const Vector_2& lp, const Vector_2& ld,
                                                              const Vector_2& d) {
        FT det = determinant(ld, d);
        if (is_zero(det)) {
            if (is_zero(determinant(lp, d))) {
                return Colinear{};
            }
            return std::nullopt;  // Parallel
        }

        FT denorm = 1.0 / det;
        FT tr = denorm * determinant(ld, lp);
        if (tr < 0) return std::nullopt;

        FT ts = denorm * determinant(d, lp);
        return std::make_pair(ts, tr);
    }

    std::optional<Parametric_line_intersection> intersect_ray(const Vector_2& lp, const Vector_2& ld, const Vector_2& d,
                                                              FT tmin, std::optional<FT> tmax) {
        auto isect = intersect_ray(lp, ld, d);
        if (!isect || std::holds_alternative<Colinear>(*isect)) return isect;

        auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect);
        if (ts < tmin || (tmax && ts >= *tmax)) return std::nullopt;  // TODO: check if t1 == tmax
        return isect;
    }

    struct Next_interval_info {
        Cone_descriptor k_next;
        FT ts;
        Endpoint_type type;
    };

    std::optional<Next_interval_info> find_next_interval(Cone_descriptor k_cur, const Parametric_line_2& segment,
                                                         Endpoint_type prev_type, FT tmin,
                                                         std::optional<FT> tmax = std::nullopt) {
        if (prev_type == SITE) {
            return std::nullopt;
        }

        auto p = construct_vector(k_cur.site()->point(), construct_point(segment));
        auto d = construct_vector(segment);

        auto process_isect = [&](const auto& isect, Endpoint_type type) -> std::optional<Next_interval_info> {
            if (!isect) return std::nullopt;
            if (std::holds_alternative<Colinear>(*isect)) {
                CGAL_assertion_msg(false, "TODO: case 2.b, 2.c");
            }
            auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect);
            if (is_zero(ts)) {
                CGAL_assertion_msg(false, "TODO: case 2.a");
            } else if (is_zero(tr)) {
                CGAL_assertion_msg(false, "TODO: case 1.b");
            } else {
                auto next_edge = type == CCW ? k_cur.site()->metric()->next_edge(k_cur.edge())
                                             : k_cur.site()->metric()->prev_edge(k_cur.edge());
                return Next_interval_info{Cone_descriptor(k_cur.site(), next_edge), ts, type};
            }
        };

        if (prev_type == CW || prev_type == UNKNOWN) {
            auto isect_next =
                intersect_ray(p, d, construct_vector(ORIGIN, construct_source(*k_cur.edge())), tmin, tmax);
            auto res = process_isect(isect_next, CW);
            if (res)
                return res;
            else if (prev_type == CW)
                return std::nullopt;
        }

        auto isect_next = intersect_ray(p, d, construct_vector(ORIGIN, construct_target(*k_cur.edge())), tmin, tmax);
        return process_isect(isect_next, CCW);
    }
    auto find_intervals(Site_const_iterator site, const Parametric_line_2& segment, Intervals& res, FT tmin,
                        std::optional<FT> tmax = std::nullopt) {
        res.clear();

        if (tmax && tmin > *tmax) return;

        auto pmin = construct_point_on(segment, tmin);
        auto p = construct_vector(site->point(), construct_point(segment));
        auto d = construct_vector(segment);

        // TODO: check case 3

        auto isect = site->metric()->any_intersected_edge(construct_vector(site->point(), pmin));
        // auto isect1 = metric_any_intersected_face(m.data, construct_vector(c.point, pmax));

        CGAL_assertion(!!isect);
        Metric_edge_iterator ed0 = *isect;

        // Find next cone of ed0
        Point_2 p0 = construct_source(*ed0);
        auto isect0 = intersect(p, d, construct_vector(ORIGIN, p0), tmin, tmax);
        int next_cone_idx = -1;
        if (isect0) {
            if (std::holds_alternative<Colinear>(*isect0)) {
                // TODO: case 2.b, 2.c
                CGAL_assertion(false);
            } else {
                auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect0);
                if (is_zero(ts)) {
                    // TODO: case 2.a

                    CGAL_assertion(false);
                } else if (is_zero(tr)) {
                    // TODO: case 1.b
                    CGAL_assertion(false);
                } else {
                    // case 1.a
                    res.add_intersection(tmin);
                    res.add_edge(ed0);

                    res.add_intersection(ts);

                    next_cone_idx = 0;
                    ed0 = site->metric()->prev_edge(ed0);
                    res.add_edge(ed0);
                }
            }
        } else {
            Point_2 p1 = construct_target(*ed0);
            auto isect1 = intersect(p, d, construct_vector(ORIGIN, p1), tmin, tmax);
            if (!isect1) {
                // Segment is in one cone
                res.add_intersection(tmin);
                res.add_edge(ed0);
                if (tmax) res.add_intersection(*tmax);
                return;
            }

            if (std::holds_alternative<Colinear>(*isect1)) {
                // TODO: case 2.b, 2.c
                CGAL_assertion(false);
            } else {
                auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect1);
                if (is_zero(ts)) {
                    // TODO: case 2.a
                    CGAL_assertion(false);
                } else if (is_zero(tr)) {
                    // TODO: case 1.b
                    CGAL_assertion(false);
                } else {
                    // case 1.a
                    res.add_intersection(tmin);
                    res.add_edge(ed0);

                    res.add_intersection(ts);

                    next_cone_idx = 1;
                    ed0 = site->metric()->next_edge(ed0);
                    res.add_edge(ed0);
                }
            }
        }

        while (true) {
            Point_2 p_next = next_cone_idx == 0 ? construct_source(*ed0) : construct_target(*ed0);
            auto isect_next = intersect(p, d, construct_vector(ORIGIN, p_next), tmin, tmax);
            if (!isect_next) break;
            auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect_next);
            CGAL_assertion(ts > tmin && tr > 0);
            res.add_intersection(ts);

            ed0 = next_cone_idx == 0 ? site->metric()->prev_edge(ed0) : site->metric()->next_edge(ed0);
            res.add_edge(ed0);
        }
        if (tmax) res.add_intersection(*tmax);
    }

    void find_all_intervals(const Parametric_line_2& segment, std::vector<Intervals>& res, FT tmin,
                            std::optional<FT> tmax = std::nullopt) {
        res.resize(m_sites.size());
        auto it = res.begin();
        for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it, ++it) {
            find_intervals(site_it, segment, *it, tmin, tmax);
        }
    }

    Vector_2 orthogonal_vector(const Line_2& l) const { return construct_vector(ca(l), cb(l)); }

    Line_2 get_bisector(const Cone_descriptor& k0, const Cone_descriptor& k1) {
        Line_2 l0 = construct_line(*(k0.edge())), l1 = construct_line(*(k1.edge()));
        Vector_2 nd0 = orthogonal_vector(l0) / abs(cc(l0)), nd1 = orthogonal_vector(l1) / abs(cc(l1));
        Vector_2 AB = nd0 - nd1;
        FT C = scalar_product(nd1, construct_vector(ORIGIN, k0.site()->point())) -
               scalar_product(nd0, construct_vector(ORIGIN, k1.site()->point()));
        return construct_line(cx(AB), cy(AB), C);
    }

    void process_i_trace(const Internal_trace& tr) {
        FT tmin = 0;

        // Clip the bisector with the cone k0 and k1
        int next_cone_idx = -1;
        std::optional<Next_interval_info> k_next_info;
        auto info0 = find_next_interval(tr.k0, tr.bisector,
                                        tr.k_prev == std::make_optional(tr.k0) ? tr.prev_type : UNKNOWN, tmin);
        auto info1 = find_next_interval(tr.k1, tr.bisector,
                                        tr.k_prev == std::make_optional(tr.k1) ? tr.prev_type : UNKNOWN, tmin);
        if (info0) {
            k_next_info = info0;
            next_cone_idx = 0;
        }
        if (info1 && (!k_next_info || info1->ts < k_next_info->ts)) {
            k_next_info = info1;
            next_cone_idx = 1;
        }

        std::optional<FT> tmax;
        if (k_next_info) tmax = k_next_info->ts;

        // Clip the bisector with the boundary
        Boundary_vertex_id v_id(cone_index(tr.k0), cone_index(tr.k1));
        vd_vertex_descriptor b_vd = vd_graph_traits::null_vertex();
        FT b_dist_min;
        auto [it_begin, it_end] = m_b_traces.equal_range(v_id);
        for (auto it = it_begin; it != it_end; ++it) {
            auto& b_trace = it->second;
            auto isect = intersect(tr.bisector, b_trace.b_line);
            if (!isect || std::holds_alternative<Colinear>(*isect)) continue;

            auto [ts, tb] = std::get<std::pair<FT, FT>>(*isect);
            if (ts < tmin || (tmax && ts >= *tmax)) continue;

            FT dist = ts - tmin;
            if (dist < b_dist_min || b_vd == vd_graph_traits::null_vertex()) {
                b_dist_min = dist;
                b_vd = b_trace.v_vd;
                tmax = ts;
            }
        }

        // Find the nearest three site bisector
        FT dist_min, tb_min;
        Line_2 bi_line_min;
        Cone_descriptor k2_min;
        bool found = false;
        for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it) {
            if (site_it == tr.k0.site() || site_it == tr.k1.site()) continue;

            find_intervals(site_it, tr.bisector, m_intervals_cache, tmin, tmax);
            auto& intervals = m_intervals_cache;

            for (auto [ed, interval_overlap] : intervals) {
                Cone_descriptor k2(site_it, ed);
                if (k2 == tr.k_prev) continue;

                Line_2 bi_line = get_bisector(tr.k0, k2);
                if (is_degenerate(bi_line)) continue;

                auto isect = intersect(tr.bisector, bi_line);
                if (!isect || std::holds_alternative<Colinear>(*isect)) continue;

                FT tb = std::get<FT>(*isect);
                if (tb <= interval_overlap.tmin || tb > interval_overlap.tmax) continue;

                FT dist = tb - tmin;
                CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
                if (dist < dist_min || !found) {
                    dist_min = dist;
                    tb_min = tb;
                    bi_line_min = bi_line;
                    k2_min = k2;
                    found = true;
                    break;
                } else if (dist == dist_min) {
                    CGAL_assertion_msg(false, "TODO: handle multiple bisectors");
                }
            }
        }

        if (found) {
            // Found a 3-site bisector
            Internal_vertex_id vid(cone_index(tr.k0), cone_index(tr.k1), cone_index(k2_min));
            if (auto vd_it = m_i_vertices.find(vid); vd_it != m_i_vertices.cend()) {
                m_voronoi->connect(tr.v_vd, vd_it->second, m_dummy_face, m_dummy_face);
                return;
            }

            Point_2 p = construct_point_on(tr.bisector, tb_min);
            auto v_vd = m_voronoi->add_vertex(p, Three_site_bisector_info{tr.k0, tr.k1, k2_min});
            m_voronoi->connect(tr.v_vd, v_vd, m_dummy_face, m_dummy_face);
            m_i_vertices.emplace(vid, v_vd);

            auto d_02 = construct_vector(bi_line_min);
            Line_2 l0 = construct_line(*(tr.k0.edge())), l1 = construct_line(*(tr.k1.edge()));
            if (is_negative(
                    scalar_product(d_02, orthogonal_vector(l1) / abs(cc(l1)) - orthogonal_vector(l0) / abs(cc(l0))))) {
                d_02 = opposite_vector(d_02);
            }
            auto bisector_02 = construct_parametric_line(p, d_02);
            m_i_traces.emplace_back(bisector_02, tr.k0, k2_min, tr.k1, v_vd);

            auto bisector_line_12 = get_bisector(tr.k1, k2_min);
            auto d_12 = construct_vector(bisector_line_12);
            if (is_negative(
                    scalar_product(d_12, orthogonal_vector(l0) / abs(cc(l0)) - orthogonal_vector(l1) / abs(cc(l1))))) {
                d_12 = opposite_vector(d_12);
            }
            auto bisector_12 = construct_parametric_line(p, d_12);
            m_i_traces.emplace_back(bisector_12, tr.k1, k2_min, tr.k0, v_vd);
        } else if (b_vd != vd_graph_traits::null_vertex()) {
            // Bisector intersects with the boundary
            m_voronoi->connect(tr.v_vd, b_vd, m_dummy_face, m_dummy_face);
        } else if (k_next_info) {
            // Found a 2-site bisector
            auto k0_next = next_cone_idx == 0 ? tr.k1 : tr.k0;
            auto k1_prev = next_cone_idx == 0 ? tr.k0 : tr.k1;
            Cone_descriptor k1_next = k_next_info->k_next;
            Internal_vertex_id v_id(cone_index(k0_next), cone_index(k1_prev), cone_index(k1_next));
            if (auto vd_it = m_i_vertices.find(v_id); vd_it != m_i_vertices.cend()) {
                m_voronoi->connect(tr.v_vd, vd_it->second, m_dummy_face, m_dummy_face);
                return;
            }

            Point_2 p = construct_point_on(tr.bisector, k_next_info->ts);
            auto v_vd = m_voronoi->add_vertex(p, Two_site_bisector_info{k0_next, k1_next.site(), k1_next.edge()});
            m_voronoi->connect(tr.v_vd, v_vd, m_dummy_face, m_dummy_face);
            m_i_vertices.emplace(v_id, v_vd);

            auto d = construct_vector(get_bisector(k0_next, k1_next));
            auto m = k_next_info->type == CCW ? construct_target(*k1_prev.edge()) : construct_source(*k1_prev.edge());
            auto md = construct_vector(ORIGIN, m);
            if (orientation(md, construct_vector(tr.bisector)) != orientation(md, d)) {
                d = opposite_vector(d);
            }
            auto bisector = construct_parametric_line(p, d);
            m_i_traces.emplace_back(bisector, k0_next, k1_next, k1_prev, v_vd);
        } else {
            CGAL_assertion_msg(false, "Unbounded interval trace");
        }
    }

    void trace_faces() {
        for (auto hd : halfedges(m_voronoi->graph)) {
            if (face(hd, m_voronoi->graph) != m_dummy_face) continue;
            auto fd = add_face(m_voronoi->graph);
            set_halfedge(fd, hd, m_voronoi->graph);
            for (auto hd_inner : halfedges_around_face(hd, m_voronoi->graph)) {
                set_face(hd_inner, fd, m_voronoi->graph);
            }
        }
        remove_face(m_dummy_face, m_voronoi->graph);
    }
};
#pragma endregion
}  // namespace CGAL::SSM_voronoi_diagram
#endif  // SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H
