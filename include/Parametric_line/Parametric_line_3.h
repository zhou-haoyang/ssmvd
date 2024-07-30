#ifndef PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
#define PARAMETRIC_LINE_PARAMETRIC_LINE_3_H

#include <CGAL/Origin.h>

#include <algorithm>
#include <limits>
#include <utility>
#include <optional>

namespace CGAL {
template <class R_>
struct Parametric_line_3 {
   public:
    using R = R_;
    using FT = typename R::FT;
    using T = typename R::T;
    using Point_3 = typename R::Point_3;
    using Vector_3 = typename R::Vector_3;
    static constexpr T INF = R::INF;

    Parametric_line_3(Point_3 p = {}, Vector_3 d = {}, T t0 = -INF, T t1 = INF)
        : m_p(std::move(p)), m_d(std::move(d)), m_tmin(t0), m_tmax(t1) {
        CGAL_precondition(m_tmin <= m_tmax);
    }

    static Parametric_line_3 segment(const Point_3 &p0, const Point_3 &p1) {
        return Parametric_line_3(p0, p1 - p0, 0, 1);
    }

    static Parametric_line_3 ray(const Point_3 &p0, const Vector_3 &d) { return Parametric_line_3(p0, d, 0, INF); }

    Point_3 operator()(FT t) const { return point(t); }

    Point_3 point(FT t) const { return R().construct_point_on_3_object()(*this, t); }

    Parametric_line_3 opposite() const { return R().construct_opposite_parametric_line_3_object()(*this); }

    Parametric_line_3 clipped(T t_min = -INF, T t_max = INF) const {
        return R().construct_parametric_line_3_object()(*this, t_min, t_max);
    }

    Point_3 p_min() const { return R().construct_pmin_3_object()(*this); }

    Point_3 p_max() const { return R().construct_pmax_3_object()(*this); }

    T t_min() const { return m_tmin; }

    T t_max() const { return m_tmax; }

    Point_3 p() const { return m_p; }

    Vector_3 d() const { return m_d; }

    //    private:
    Point_3 m_p;
    Vector_3 m_d;
    T m_tmin = -INF, m_tmax = INF;
};

template <class K_, class T_ = double>
class Parametric_line_traits_3 {
   public:
    using K = K_;
    using T = T_;
    using FT = typename K::FT;
    using Point_3 = typename K::Point_3;
    using Vector_3 = typename K::Vector_3;
    using Plane_3 = typename K::Plane_3;
    using Self = Parametric_line_traits_3<K, T>;
    using Parametric_line_3 = Parametric_line_3<Self>;

    static constexpr auto INF = std::numeric_limits<T>::infinity();

    struct Construct_parametric_segment_3 {
        Parametric_line_3 operator()(const Point_3 &p0, const Point_3 &p1) const {
            return Parametric_line_3(p0, p1 - p0, 0, 1);
        }
    };

    struct Construct_parametric_ray_3 {
        Parametric_line_3 operator()(const Point_3 &p0, const Vector_3 &d) const {
            return Parametric_line_3(p0, d, 0, INF);
        }
    };

    struct Construct_parametric_line_3 {
        Parametric_line_3 operator()(const Point_3 &p, const Vector_3 &d, T tmin = -INF, T tmax = INF) const {
            return Parametric_line_3(p, d, tmin, tmax);
        }

        Parametric_line_3 operator()(const Point_3 &p0, const Point_3 &p1) const {
            return Parametric_line_3(p0, p1 - p0, -INF, INF);
        }

        Parametric_line_3 operator()(const Parametric_line_3 &pline, T tmin = -INF, T tmax = INF) const {
            return Parametric_line_3(pline.p(), pline.d(), std::max(pline.t_min(), tmin),
                                     std::min(pline.t_max(), tmax));
        }
    };

    struct Construct_opposite_parametric_line_3 {
        Parametric_line_3 operator()(const Parametric_line_3 &l) const {
            return Parametric_line_3(l.p(), -l.d(), -l.t_max(), -l.t_min());
        }
    };

    struct Construct_source_3 : public K::Construct_point_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.p(); }
    };

    struct Construct_vector_3 : public K::Construct_vector_3 {
        Vector_3 operator()(const Parametric_line_3 &l) const { return l.d(); }
    };

    struct Construct_point_on_3 : public K::Construct_point_on_3 {
        Point_3 operator()(const Parametric_line_3 &l, FT t) const { return l.p() + t * l.d(); }
    };

    struct Construct_pmin_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.p() + l.t_min() * l.d(); }
    };

    struct Construct_pmax_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.p() + l.t_max() * l.d(); }
    };

    struct Compute_tmin_3 {
        T operator()(const Parametric_line_3 &l) const { return l.t_min(); }
    };

    struct Compute_tmax_3 {
        T operator()(const Parametric_line_3 &l) const { return l.t_max(); }
    };

    struct Is_point_3 {
        bool operator()(const Parametric_line_3 &l) const { return iszero(l.t_min() - l.t_max()); }
    };

    struct Intersect_3 : public K::Intersect_3 {
        std::optional<T> operator()(const Parametric_line_3 &l, const Plane_3 &plane) const {
            auto scalar_product = K().compute_scalar_product_3_object();

            auto n = K().construct_orthogonal_vector_3_object()(plane);
            auto D = K().compute_d_3_object()(plane);
            auto nd = scalar_product(n, l.d());
            auto np = scalar_product(n, K().construct_vector_3_object()(ORIGIN, l.p())) + D;

            if (is_zero(nd)) {
                if (is_zero(np)) return INF;
                return std::nullopt;
            }
            T t = -np / nd;
            if (t >= l.t_min() && t <= l.t_max()) return t;
            return std::nullopt;
        }

        std::optional<std::pair<T, T>> operator()(const Parametric_line_3 &l0, const Parametric_line_3 &l1,
                                                  bool coplanar = false) {
            auto d0 = l0.d(), d1 = l1.d();
            auto p0 = l0.p(), p1 = l1.p();
            auto n = cross_product(d0, d1);
            if (is_zero(n.squared_length())) {
                // parallel
                return std::nullopt;
            }

            auto tn0 = cross_product((p1 - p0), d1);
            if (!coplanar && !is_zero(scalar_product(tn0, n))) {
                // not coplanar
                return std::nullopt;
            }
            T t0 = tn0.x() / n.x();
            auto tn1 = cross_product((p1 - p0), d0);
            T t1 = tn1.x() / n.x();

            if (t0 < l0.t_min() || t0 > l0.t_max() || t1 < l1.t_min() || t1 > l1.t_max()) {
                return std::nullopt;
            }
            return std::make_pair(t0, t1);
        }
    };

    auto construct_parametric_segment_3_object() const { return Construct_parametric_segment_3(); }

    auto construct_parametric_ray_3_object() const { return Construct_parametric_ray_3(); }

    auto construct_parametric_line_3_object() const { return Construct_parametric_line_3(); }

    auto construct_opposite_parametric_line_3_object() const { return Construct_opposite_parametric_line_3(); }

    auto construct_source_3_object() const { return Construct_source_3(); }

    auto construct_vector_3_object() const { return Construct_vector_3(); }

    auto construct_point_on_3_object() const { return Construct_point_on_3(); }

    auto construct_pmin_3_object() const { return Construct_pmin_3(); }

    auto construct_pmax_3_object() const { return Construct_pmax_3(); }

    auto compute_tmin_3_object() const { return Compute_tmin_3(); }

    auto compute_tmax_3_object() const { return Compute_tmax_3(); }

    auto is_point_3_object() const { return Is_point_3(); }

    auto intersect_3_object() const { return Intersect_3(); }
};

// template <class K, class T = double>
// class Parametric_line_3 : public Parametric_line_3_traits<K, T>::Parametric_line_3 {
//    public:
//     using R = Parametric_line_3_traits<K, T>;
//     using Rep = typename R::Parametric_line_3;
//     using Point_3 = typename R::Point_3;
//     using Vector_3 = typename R::Vector_3;
//     using FT = typename R::FT;
//     static constexpr T INF = R::INF;

//     Parametric_line_3(Point_3 p = {}, Vector_3 d = {}, T t0 = -INF, T t1 = INF) : Rep(p, d, t0, t1) {}

//     Point_3 point(FT t) const { return R().construct_point_on_3_object()(*this, t); }

//     Parametric_line_3 opposite() const { return R().construct_opposite_parametric_line_3_object()(*this); }

//     Point_3 pmin() const { return R().construct_pmin_3_object()(*this); }

//     Point_3 pmax() const { return R().construct_pmax_3_object()(*this); }
// };
}  // namespace CGAL

#endif  // PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
