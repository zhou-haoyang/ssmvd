#ifndef PARAMETRIC_LINE_PARAMETRIC_LINE_2_H
#define PARAMETRIC_LINE_PARAMETRIC_LINE_2_H

#include <CGAL/Origin.h>
#include <CGAL/determinant.h>

#include <algorithm>
#include <optional>
#include <utility>
#include <variant>

namespace CGAL {
template <class R_> struct Parametric_line_2
{
public:
  using R = R_;
  using FT = typename R::FT;
  using Point_2 = typename R::Point_2;
  using Vector_2 = typename R::Vector_2;

  Parametric_line_2(Point_2 p = {}, Vector_2 d = {})
      : m_p(std::move(p))
      , m_d(std::move(d)) {}

  Parametric_line_2(const Point_2& p0, const Point_2& p1)
      : m_p(p0)
      , m_d(p1 - p0) {}

  Point_2 operator()(FT t) const { return point(t); }

  friend std::ostream& operator<<(std::ostream& os, const Parametric_line_2& l) {
    return os << "{(" << l.p() << ") + t(" << l.d() << ")}";
  }

  Point_2 point(FT t) const { return R().construct_point_on_2_object()(*this, t); }

  Parametric_line_2 opposite() const { return R().construct_opposite_parametric_line_2_object()(*this); }

  Point_2 p() const { return m_p; }

  Vector_2 d() const { return m_d; }

private:
  Point_2 m_p;
  Vector_2 m_d;
};

template <class K_, class T_ = double> class Parametric_line_traits_2
{
public:
  using K = K_;
  using FT = typename K::FT;
  using Point_2 = typename K::Point_2;
  using Vector_2 = typename K::Vector_2;
  using Line_2 = typename K::Line_2;
  using Segment_2 = typename K::Segment_2;
  using Self = Parametric_line_traits_2<K>;
  using Parametric_line_2 = Parametric_line_2<Self>;

  struct Colinear
  {
    bool operator==(const Colinear&) const { return true; }
  };

  using Parameter_pair = std::pair<FT, FT>;

  struct Construct_parametric_line_2
  {
    Parametric_line_2 operator()(const Point_2& p, const Vector_2& d) const { return Parametric_line_2(p, d); }

    Parametric_line_2 operator()(const Point_2& p0, const Point_2& p1) const { return Parametric_line_2(p0, p1 - p0); }

    Parametric_line_2 operator()(const Segment_2& s) const { return Parametric_line_2(s.source(), s.target()); }
  };

  struct Construct_opposite_parametric_line_2
  {
    Parametric_line_2 operator()(const Parametric_line_2& l) const { return Parametric_line_2(l.p(), -l.d()); }
  };

  struct Construct_point_2 : public K::Construct_point_2
  {
    using K::Construct_point_2::operator();
    Point_2 operator()(const Parametric_line_2& l) const { return l.p(); }
  };

  struct Construct_vector_2 : public K::Construct_vector_2
  {
    using K::Construct_vector_2::operator();
    Vector_2 operator()(const Parametric_line_2& l) const { return l.d(); }
  };

  struct Construct_point_on_2 : public K::Construct_point_on_2
  {
    using K::Construct_point_on_2::operator();
    Point_2 operator()(const Parametric_line_2& l, FT t) const { return l.p() + t * l.d(); }
  };

  struct Intersect_2 : public K::Intersect_2
  {
    std::optional<std::variant<Parameter_pair, Colinear>> operator()(const Parametric_line_2& l0,
                                                                     const Parametric_line_2& l1) {
      Point_2 p0 = l0.p(), p1 = l1.p();
      Vector_2 d0 = l0.d(), d1 = l1.d();
      Vector_2 dp = p1 - p0;

      FT det = determinant(d0, d1);
      if(is_zero(det)) {
        if(is_zero(determinant(dp, d0)))
          return Colinear();
        return std::nullopt; // parallel
      }

      FT denorm = 1.0 / det;
      FT t0 = denorm * determinant(dp, d1);
      FT t1 = denorm * determinant(dp, d0);

      return std::make_pair(t0, t1);
    }

    std::optional<std::variant<FT, Colinear>> operator()(const Parametric_line_2& l0, const Line_2& l1) {
      Point_2 p = l0.p();
      Vector_2 d = l0.d();
      Vector_2 n(l1.a(), l1.b());
      FT c = l1.c();

      FT nd = scalar_product(n, d);
      FT num = scalar_product(n, p - ORIGIN) + c;
      if(is_zero(nd)) {
        if(is_zero(num))
          return Colinear();
        return std::nullopt; // parallel
      }

      return -num / nd;
    }
  };

  auto construct_parametric_line_2_object() const { return Construct_parametric_line_2(); }

  auto construct_opposite_parametric_line_2_object() const { return Construct_opposite_parametric_line_2(); }

  auto construct_point_2_object() const { return Construct_point_2(); }

  auto construct_vector_2_object() const { return Construct_vector_2(); }

  auto construct_point_on_2_object() const { return Construct_point_on_2(); }

  auto intersect_2_object() const { return Intersect_2(); }
};
} // namespace CGAL

#endif // PARAMETRIC_LINE_PARAMETRIC_LINE_2_H
