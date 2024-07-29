#ifndef PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
#define PARAMETRIC_LINE_PARAMETRIC_LINE_3_H

#include <CGAL/Origin.h>
#include <algorithm>
#include <limits>
#include <utility>
#include <optional>

namespace CGAL {
template <class K, class T = double>
class Parametric_line_3 {
   public:
    static constexpr auto INF = std::numeric_limits<T>::infinity();

    using Kernel = K;
    using FT = typename K::FT;
    using Point_3 = typename K::Point_3;
    using Vector_3 = typename K::Vector_3;
    using Plane_3 = typename K::Plane_3;

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
        Parametric_line_3 operator()(const Point_3 &p0, const Point_3 &p1) const {
            return Parametric_line_3(p0, p1 - p0, -INF, INF);
        }
    };

    struct Construct_opposite_parametric_line_3 {
        Parametric_line_3 operator()(const Parametric_line_3 &l) const {
            return Parametric_line_3(l.m_p, -l.m_d, -l.m_tmax, -l.m_tmin);
        }
    };

    struct Construct_point_on_3 {
        Point_3 operator()(const Parametric_line_3 &l, FT t) const { return l.m_p + t * l.m_d; }
    };

    struct Construct_pmin_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.m_p + l.m_tmin * l.m_d; }
    };

    struct Construct_pmax_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.m_p + l.m_tmax * l.m_d; }
    };

    struct Is_point_3 {
        bool operator()(const Parametric_line_3 &l) const { return iszero(l.m_tmin - l.m_tmax); }
    };

    struct Intersect_3 {
        std::optional<T> operator()(const Parametric_line_3 &l, const Plane_3 &plane) const {
            auto scalar_product = K().scalar_product_3_object();

            auto n = K().construct_orthogonal_vector_3_object()(plane);
            auto D = K().compute_d_3_object()(plane);
            auto nd = scalar_product(n, l.m_d);
            auto np = scalar_product(n, K().construct_vector_3_object(l.m_p, ORIGIN)) + D;

            if (is_zero(nd)) {
                if (is_zero(np)) return INF;
                return std::nullopt;
            }
            T t = -np / nd;
            if (t >= l.m_tmin && t <= l.m_tmax) return t;
            return std::nullopt;
        }
    };

    Parametric_line_3(Point_3 p = {}, Vector_3 d = {}, T t0 = -INF, T t1 = INF)
        : m_p(std::move(p)), m_d(std::move(d)), m_tmin(std::min(t0, t1)), m_tmax(std::max(t0, t1)) {}

    auto construct_parametric_segment_3_object() const { return Construct_parametric_segment_3(); }

    auto construct_parametric_ray_3_object() const { return Construct_parametric_ray_3(); }

    auto construct_parametric_line_3_object() const { return Construct_parametric_line_3(); }

    auto construct_opposite_parametric_line_3_object() const { return Construct_opposite_parametric_line_3(); }

    auto construct_point_on_3_object() const { return Construct_point_on_3(); }

    auto construct_pmin_3_object() const { return Construct_pmin_3(); }

    auto construct_pmax_3_object() const { return Construct_pmax_3(); }

    auto is_point_3_object() const { return Is_point_3(); }

    Point_3 operator()(FT t) const { return m_p + t * m_d; }

    Parametric_line_3 opposite() const { return construct_opposite_parametric_line_3_object()(*this); }

    Point_3 pmin() const { return construct_pmin_3_object()(*this); }

    Point_3 pmax() const { return construct_pmax_3_object()(*this); }

   private:
    Point_3 m_p;
    Vector_3 m_d;
    T m_tmin = -INF, m_tmax = INF;
};
}  // namespace CGAL

#endif  // PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
