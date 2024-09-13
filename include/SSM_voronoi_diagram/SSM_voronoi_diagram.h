#ifndef SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H
#define SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H

#include <CGAL/Default.h>
#include <CGAL/Origin.h>
#include <CGAL/Polygon_2.h>

#include <boost/graph/graph_traits.hpp>
#include <list>
#include <optional>
#include <variant>

#define TRAIT_FUNC(ret_type, name, functor)                  \
    template <class... T>                                    \
    ret_type name(T&&... args) {                             \
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
    using Site_iterator = typename Site_list::iterator;

    class Cone_descriptor {
       public:
        Cone_descriptor(Site_iterator s, Metric_edge_iterator e) : m_site(std::move(s)), m_edge(std::move(e)) {}

        PROPERTY(Site_iterator, site)
        PROPERTY(Metric_edge_iterator, edge)

       private:
        Site_iterator m_site;
        Metric_edge_iterator m_edge;
    };

    using Voronoi_diagram = typename Traits::Voronoi_diagram;
    using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram>;
    using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
    using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
    using vd_halfedge_descriptor = typename vd_graph_traits::halfedge_descriptor;
    using vd_face_descriptor = typename vd_graph_traits::face_descriptor;

   protected:
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
                res.add_intersection(*it, eds[i]);
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

    SSM_voronoi_diagram(const Polygon_2& boundary, Traits traits = {})
        : m_boundary(boundary), m_traits(std::move(traits)) {}

    Metric_iterator add_metric(Metric m) { return m_metrics.emplace(m_metrics.end(), std::move(m), m_traits); }

    Site_iterator add_site(Point_2 p, Metric_iterator m) { return m_sites.emplace(m_sites.end(), std::move(p), m); }

    FT find_nearest_site(const Point_2& p, Cone_descriptor& cone) const {
        FT d_min;
        bool init = true;

        for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it) {
            auto res = site_it->metric()->any_intersection(construct_vector(site_it->point(), p));
            if (!res) continue;
            auto [pm, ed] = *res;
            FT weight = 1.0 / squared_distance(ORIGIN, pm);

            auto d = approximate_sqrt(squared_distance(p, site_it->point()) * weight);
            if (d < d_min || init) {
                d_min = d;
                cone.set_site(site_it);
                cone.set_edge(ed);
            }
        }

        return d_min;
    }

    auto trace_boundary(Metric_edge_iterator ed, Cone_descriptor k0, vd_vertex_descriptor prev_vd,
                        bool add_border_edges = false, bool add_cone_vertices = false,
                        vd_face_descriptor fd0 = vd_graph_traits::null_face(),
                        vd_face_descriptor fd1 = vd_graph_traits::null_face()) {
        auto b_line = construct_parametric_line(*ed);
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
            std::vector<Cone_descriptor> k1_min;
            bool found = false;
            for (auto site_it = m_sites.begin(); site_it != m_sites.end(); ++site_it) {
                if (site_it == k0.site()) continue;

                auto interval = *interval_it;
                FT tmin_overlap = max(interval.t_min(), tmin);
                FT tmax_overlap = min(interval.t_max(), tmax);
                intervals_map[site_index(site_it)].clip(tmin_overlap, tmax_overlap, m_intervals_cache);
                auto& intervals_overlap = m_intervals_cache;

                for (auto& [ed, interval_overlap] : intervals_overlap) {
                    Cone_descriptor k1(site_it, ed);
                    if (k1 == k_prev) continue;

                    Line_2 bi_line = get_bisector(k0, k1);
                    if (is_degenerate(bi_line)) continue;

                    auto isect = intersect(b_line, bi_line);
                    if (!isect || std::holds_alternative<Colinear>(*isect)) continue;

                    FT tb = std::get<FT>(*isect);
                    // TODO: marginal case: tb is on the terminal points
                    if (tb <= interval_overlap.tmin() || tb > interval_overlap.tmax()) continue;

                    FT dist = tb - tmin;
                    CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
                    if (dist < dist_min || !found) {
                        dist_min = dist;
                        tb_min = tb;
                        bi_line_min = bi_line;
                        k1_min.clear();
                        k1_min.push_back(k1);
                        found = true;
                        break;  // goto next site
                    } else if (dist == dist_min) {
                        k1_min.push_back(k1);
                    }
                }
            }

            if (!found) {
                // No intersection found

            } else {
            }
        }
    }

   protected:
    const Polygon_2& m_boundary;
    Metric_list m_metrics;
    Site_list m_sites;
    Traits m_traits;

    Intervals m_intervals_cache;
    std::vector<Intervals> m_intervals_vec_cache;

    TRAIT_FUNC(Vector_2, construct_vector, construct_vector_2_object)
    TRAIT_FUNC(Point_2, construct_point, construct_point_2_object)
    TRAIT_FUNC(Point_2, construct_point_on, construct_point_on_2_object)
    TRAIT_FUNC(Point_2, construct_source, construct_source_2_object)
    TRAIT_FUNC(Point_2, construct_target, construct_target_2_object)
    TRAIT_FUNC(Parametric_line_2, construct_parametric_line, construct_parametric_line_2_object)
    TRAIT_FUNC(Line_2, construct_line, construct_line_2_object)

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

    std::size_t site_index(Site_iterator site) const { return std::distance(m_sites.begin(), site); }

    /**
     * @brief Find the intersection of a parametric line with a ray with direction d and starting at the origin.
     *
     * @param l
     * @param d
     * @return std::optional<FT>
     */
    std::optional<Parametric_line_intersection> intersect(const Vector_2& lp, const Vector_2& ld, const Vector_2& d) {
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

    std::optional<Parametric_line_intersection> intersect(const Vector_2& lp, const Vector_2& ld, const Vector_2& d,
                                                          FT tmin, std::optional<FT> tmax) {
        auto isect = intersect(lp, ld, d);
        if (!isect || std::holds_alternative<Colinear>(*isect)) return isect;

        auto [ts, tr] = std::get<std::pair<FT, FT>>(*isect);
        if (ts < tmin || (tmax && ts >= *tmax)) return std::nullopt;  // TODO: check if t1 == tmax
        return isect;
    }

    auto find_intervals(Site_iterator site, const Parametric_line_2& segment, Intervals& res, FT tmin,
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

    Vector_2 orthogonal_vertor(const Line_2& l) { return construct_vector(ca(l), cb(l)); }

    Line_2 get_bisector(const Cone_descriptor& k0, const Cone_descriptor& k1) {
        Line_2 l0 = construct_line(*(k0.edge())), l1 = construct_line(*(k1.edge()));
        Vector_2 nd0 = orthogonal_vertor(l0) / abs(cc(l0)), nd1 = orthogonal_vertor(l1) / abs(cc(l1));
        Vector_2 AB = nd0 - nd1;
        FT C = scalar_product(nd1, construct_vector(ORIGIN, k0.site()->point())) -
               scalar_product(nd0, construct_vector(ORIGIN, k1.site()->point()));
        return construct_line(cx(AB), cy(AB), C);
    }
};
}  // namespace CGAL::SSM_voronoi_diagram
#endif  // SSM_VORONOI_DIAGRAM_SSM_VORONOI_DIAGRAM_H
