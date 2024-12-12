#ifndef UTILS_GRAPH_HELPER_H
#define UTILS_GRAPH_HELPER_H
namespace CGAL {

/**
 * @brief A monotonic function mapping the angle between two vectors in range [0, 2pi) to [-2, 2)
 */
template <class GT>
inline typename GT::FT angle_func(const typename GT::Vector_2 &v1, const typename GT::Vector_2 &v2, const GT &traits) {
    auto cos_theta = traits.compute_scalar_product_2_object()(v1, v2);
    auto ori = traits.orientation_2_object()(v1, v2);
    if (ori == ZERO) {
        return cos_theta < 0 ? cos_theta + 1 : -cos_theta - 1;
    } else {
        return ori == NEGATIVE ? cos_theta + 1 : -cos_theta - 1;
    }
}

template <class Kernel>
class Projection_traits_with_scalar_product_3 : public Kernel {
   public:
    using typename Kernel::Vector_3;
    using Vector_2 = typename Kernel::Vector_3;

    auto compute_scalar_product_2_object() const {
        return [this](const Vector_2 &v1, const Vector_2 &v2) {
            auto dot = Kernel().compute_scalar_product_3_object();
            return dot(v1, v2) - dot(v1, m_n) * dot(v2, m_n) / dot(m_n, m_n);
        };
    }

    auto orientation_2_object() const {
        return [this](const Vector_2 &v1, const Vector_2 &v2) { return Kernel::orientation_3_object()(v1, v2, m_n); };
    }

    explicit Projection_traits_with_scalar_product_3(const typename Kernel::Vector_3 &n) : m_n(n) {}

   private:
    Vector_3 m_n;
};

template <class G, class VPM, class VNM, class GT>
void insert_halfedge_loop(typename boost::graph_traits<G>::halfedge_descriptor hd, G &graph, const VNM &vnm,
                          const VPM &vpm, const GT &traits) {
    using vd_graph_traits = boost::graph_traits<G>;
    using FT = typename GT::FT;

    auto vt = target(hd, graph), vs = source(hd, graph);
    auto hd_cur = halfedge(vt, graph);
    if (hd_cur == vd_graph_traits::null_halfedge()) {
        set_halfedge(vt, hd, graph);
        return;
    }

    auto hd_next = opposite(next(hd_cur, graph), graph);
    if (hd_cur != hd_next) {
        // At least 2 halfedges around the target vertex
        auto pt = get(vpm, vt), ps = get(vpm, vs);

        auto n = get(vnm, vt);
        // The halfedge loop around the vertex should be clockwise for the halfedge loop around the face
        // to be counter-clockwise, hence the normal should point inward here
        Projection_traits_with_scalar_product_3<GT> proj_traits{-n};

        auto hd0 = hd_cur;
        auto v0 = get(vpm, source(hd0, graph)) - pt;
        auto v = ps - pt;
        auto t = angle_func(v0, v, proj_traits);
        FT t_cur = -2;  // angle_func(v0, v, proj_traits);

        // std::clog << "***" << std::endl;
        // std::clog << "v0: " << v0 << std::endl;
        // std::clog << "v: " << v << std::endl;
        // std::clog << "t: " << t << std::endl;
        // std::clog << "---" << std::endl;
        do {
            auto v_next = get(vpm, source(hd_next, graph)) - pt;
            auto t_next = angle_func(v0, v_next, proj_traits);

            // std::clog << "t range: " << t_cur << " " << t_next << std::endl;
            // std::clog << "---" << std::endl;

            if (t >= t_cur && t < t_next) {
                break;
            }

            t_cur = t_next;
            hd_cur = hd_next;
            hd_next = opposite(next(hd_cur, graph), graph);
        } while (hd_next != hd0);
        // v is between v_{n-1} and v0
        CGAL_assertion(t >= t_cur && t < 2);
    }

    set_next(hd, next(hd_cur, graph), graph);
    set_next(hd_cur, opposite(hd, graph), graph);
}

template <class G, class GT, class VNM, class VPM>
typename boost::graph_traits<G>::halfedge_descriptor connect_vertices(
    typename boost::graph_traits<G>::vertex_descriptor v0, typename boost::graph_traits<G>::vertex_descriptor v1,
    G &graph, const VNM &vnm, const VPM &vpm, const GT &traits = {}) {
    using vd_graph_traits = boost::graph_traits<G>;

    CGAL_precondition(v0 != v1);
    CGAL_precondition(v0 != vd_graph_traits::null_vertex() && v1 != vd_graph_traits::null_vertex());

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
    insert_halfedge_loop(hd10, graph, vnm, vpm, traits);
    insert_halfedge_loop(hd01, graph, vnm, vpm, traits);

    return hd01;
}
}  // namespace CGAL
#endif  // UTILS_GRAPH_HELPER_H
