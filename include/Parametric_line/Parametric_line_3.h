#ifndef PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
#define PARAMETRIC_LINE_PARAMETRIC_LINE_3_H

#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Origin.h>

#include <algorithm>
#include <utility>
#include <optional>

namespace CGAL {
template <class R_>
struct Parametric_line_3 {
   public:
    using R = R_;
    using FT = typename R::FT;
    using Point_3 = typename R::Point_3;
    using Vector_3 = typename R::Vector_3;

    Parametric_line_3(Point_3 p = {}, Vector_3 d = {}) : m_p(std::move(p)), m_d(std::move(d)) {}

    Point_3 operator()(FT t) const { return point(t); }

    friend std::ostream &operator<<(std::ostream &os, const Parametric_line_3 &l) {
        return os << "{(" << l.p() << ") + t(" << l.d() << ")}";
    }

    Point_3 point(FT t) const { return R().construct_point_on_3_object()(*this, t); }

    Parametric_line_3 opposite() const { return R().construct_opposite_parametric_line_3_object()(*this); }

    FT parameter(const Point_3 &p) const { return R().compute_parameter_3_object()(*this, p); }

    Point_3 p() const { return m_p; }

    Vector_3 d() const { return m_d; }

    bool has_on(const Point_3 &p) const { return R().has_on_3_object()(*this, p); }

    friend FT squared_distance(const Parametric_line_3 &l, const Point_3 &p) {
        return R().compute_squared_distance_3_object()(l, p);
    }

   private:
    Point_3 m_p;
    Vector_3 m_d;
};

template <class K_, class T_ = double>
class Parametric_line_traits_3 {
   public:
    using K = K_;
    using FT = typename K::FT;
    using Point_3 = typename K::Point_3;
    using Vector_3 = typename K::Vector_3;
    using Plane_3 = typename K::Plane_3;
    using Self = Parametric_line_traits_3<K>;
    using Parametric_line_3 = Parametric_line_3<Self>;

    struct Construct_parametric_line_3 {
        Parametric_line_3 operator()(const Point_3 &p, const Vector_3 &d) const { return Parametric_line_3(p, d); }

        Parametric_line_3 operator()(const Point_3 &p0, const Point_3 &p1) const {
            return Parametric_line_3(p0, p1 - p0);
        }
    };

    struct Construct_opposite_parametric_line_3 {
        Parametric_line_3 operator()(const Parametric_line_3 &l) const { return Parametric_line_3(l.p(), -l.d()); }
    };

    struct Construct_point_3 : public K::Construct_point_3 {
        using K::Construct_point_3::operator();
        Point_3 operator()(const Parametric_line_3 &l) const { return l.p(); }
    };

    struct Construct_vector_3 : public K::Construct_vector_3 {
        using K::Construct_vector_3::operator();
        Vector_3 operator()(const Parametric_line_3 &l) const { return l.d(); }
    };

    struct Construct_point_on_3 : public K::Construct_point_on_3 {
        using K::Construct_point_on_3::operator();
        Point_3 operator()(const Parametric_line_3 &l, FT t) const { return l.p() + t * l.d(); }
    };

    struct Construct_target_3 {
        Point_3 operator()(const Parametric_line_3 &l) const { return l.p() + l.d(); }
    };

    struct Intersect_3 : public K::Intersect_3 {
        std::optional<FT> operator()(const Parametric_line_3 &l, const Plane_3 &plane) const {
            auto scalar_product = K().compute_scalar_product_3_object();

            auto n = K().construct_orthogonal_vector_3_object()(plane);
            auto D = K().compute_d_3_object()(plane);
            auto nd = scalar_product(n, l.d());
            auto np = scalar_product(n, K().construct_vector_3_object()(ORIGIN, l.p())) + D;

            if (is_zero(nd)) {
                // TODO: line on the plane
                // if (is_zero(np)) return INF;
                return std::nullopt;
            }
            FT t = -np / nd;
            return t;
        }

        std::optional<std::pair<FT, FT>> operator()(const Parametric_line_3 &l0, const Parametric_line_3 &l1,
                                                    bool coplanar = false) {
            auto d0 = l0.d(), d1 = l1.d();
            auto p0 = l0.p(), p1 = l1.p();
            auto n = cross_product(d0, d1);
            FT n_sqr = n.squared_length();

            if (is_zero(n_sqr)) {
                // parallel
                return std::nullopt;
            }

            auto tn0 = cross_product((p1 - p0), d1);
            auto tn0_dot_n = scalar_product(tn0, n);
            if (!coplanar && !is_zero(tn0_dot_n)) {
                // not coplanar
                return std::nullopt;
            }
            // Solve t0 . n = tn0
            FT t0 = tn0_dot_n / n_sqr;
            // FT t0 = tn0.x() / n.x();

            auto tn1 = cross_product((p1 - p0), d0);
            // Solve t1 . n = tn1
            FT t1 = scalar_product(tn1, n) / n_sqr;
            // FT t1 = tn1.x() / n.x();

            return std::make_pair(t0, t1);
        }
    };

    struct Compute_parameter_3 {
        FT operator()(const Parametric_line_3 &l, const Point_3 &p) const {
            // l.p() + t * l.d() = p
            auto v = K().construct_vector_3_object()(l.p(), p);
            // TODO: robustness, e.g. if l.d() is zero
            return scalar_product(v, l.d()) / l.d().squared_length();
        }
    };

    struct Compute_Squared_distance_3 {
        FT operator()(const Parametric_line_3 &l, const Point_3 &p) const {
            auto v = K().construct_vector_3_object()(l.p(), p);
            auto dp = K().compute_scalar_product_3_object()(l.d(), v);
            return v.squared_length() - dp * dp / l.d().squared_length();
        }
    };

    struct Has_on_3 : public K::Has_on_3 {
        using K::Has_on_3::operator();
        bool operator()(const Parametric_line_3 &l, const Point_3 &p) const {
            return K().collinear_3_object()(l.p(), l.p() + l.d(), p);
        }
    };

    auto construct_parametric_line_3_object() const { return Construct_parametric_line_3(); }

    auto construct_opposite_parametric_line_3_object() const { return Construct_opposite_parametric_line_3(); }

    auto construct_point_3_object() const { return Construct_point_3(); }

    auto construct_vector_3_object() const { return Construct_vector_3(); }

    auto construct_point_on_3_object() const { return Construct_point_on_3(); }

    auto construct_target_3_object() const { return Construct_target_3(); }

    auto intersect_3_object() const { return Intersect_3(); }

    auto compute_parameter_3_object() const { return Compute_parameter_3(); }

    auto compute_squared_distance_3_object() const { return Compute_Squared_distance_3(); }

    auto has_on_3_object() const { return Has_on_3(); }
};
}  // namespace CGAL

#endif  // PARAMETRIC_LINE_PARAMETRIC_LINE_3_H
