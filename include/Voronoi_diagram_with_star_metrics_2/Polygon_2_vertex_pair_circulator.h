#ifndef VORONOI_DIAGRAM_WITH_STAR_METRICS_2_POLYGON_2_VERTEX_PAIR_CIRCULATOR_H
#define VORONOI_DIAGRAM_WITH_STAR_METRICS_2_POLYGON_2_VERTEX_PAIR_CIRCULATOR_H

#include <CGAL/circulator.h>

namespace CGAL {
template <class Metric_vertex_circulator> class Polygon_2_vertex_pair_circulator
{
public:
  using value_type = std::pair<Metric_vertex_circulator, Metric_vertex_circulator>;
  using difference_type = typename Metric_vertex_circulator::difference_type;
  using size_type = typename Metric_vertex_circulator::size_type;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using reference = value_type&;
  using const_reference = const value_type&;
  using iterator_category = Bidirectional_circulator_tag;

  explicit Polygon_2_vertex_pair_circulator(Metric_vertex_circulator vertex = Metric_vertex_circulator{})
      : m_vertex(std::move(vertex)) {}

  bool operator==(const Polygon_2_vertex_pair_circulator& other) const { return m_vertex == other.m_vertex; }

  bool operator!=(const Polygon_2_vertex_pair_circulator& other) const { return m_vertex != other.m_vertex; }

  const value_type& operator*() const {
    auto v_next = m_vertex;
    ++v_next;
    m_value = std::make_pair(m_vertex, v_next);
    return m_value;
  }

  const value_type* operator->() const { return &(**this); }

  Polygon_2_vertex_pair_circulator& operator++() {
    ++m_vertex;
    return *this;
  }

  Polygon_2_vertex_pair_circulator operator++(int) {
    Polygon_2_vertex_pair_circulator tmp = *this;
    ++*this;
    return tmp;
  }

  Polygon_2_vertex_pair_circulator& operator--() {
    --m_vertex;
    return *this;
  }

  Polygon_2_vertex_pair_circulator operator--(int) {
    Polygon_2_vertex_pair_circulator tmp = *this;
    --*this;
    return tmp;
  }

  difference_type operator-(const Polygon_2_vertex_pair_circulator& other) const { return m_vertex - other.m_vertex; }

  Polygon_2_vertex_pair_circulator operator+(difference_type n) const {
    return Polygon_2_vertex_pair_circulator(m_vertex + n);
  }

  Polygon_2_vertex_pair_circulator operator-(difference_type n) const {
    return Polygon_2_vertex_pair_circulator(m_vertex - n);
  }

protected:
  Metric_vertex_circulator m_vertex;
  mutable value_type m_value;
};
} // namespace CGAL

#endif // VORONOI_DIAGRAM_WITH_STAR_METRICS_2_POLYGON_2_VERTEX_PAIR_CIRCULATOR_H
