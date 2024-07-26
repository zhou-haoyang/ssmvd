#ifndef SSM_RESTRICTED_VORONOI_DIAGRAM_PARAMETRIC_LINE_3_H
#define SSM_RESTRICTED_VORONOI_DIAGRAM_PARAMETRIC_LINE_3_H

#include <CGAL/Origin.h>

#include <algorithm>
#include <limits>

namespace CGAL {
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

    bool is_point() const { return is_zero(t_min - t_max); }
};

template <class K, class T>
bool isect(const Parametric_line_3<K, T> &l, const typename K::Vector_3 &n, const typename K::FT &D, T &t) {
    auto nd = scalar_product(n, l.d), np = scalar_product(n, l.p - ORIGIN) + D;
    if (is_zero(nd)) {
        t = std::numeric_limits<T>::infinity();
        return is_zero(np);
    }
    t = -np / nd;
    return t >= l.t_min && t <= l.t_max;
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
}  // namespace CGAL
#endif  // SSM_RESTRICTED_VORONOI_DIAGRAM_PARAMETRIC_LINE_3_H
