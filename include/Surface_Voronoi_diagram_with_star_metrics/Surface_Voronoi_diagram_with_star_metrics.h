#ifndef SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_H
#define SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_H

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/tss.h>
#include <IO/Verbosity_level_ostream.h>
#include <Parametric_line/Parametric_line_3.h>
#include <Surface_Voronoi_diagram_with_star_metrics/Surface_Voronoi_diagram_with_star_metrics_traits.h>
#include <Utils/Graph_helper.h>

#include <CGAL/Default.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Origin.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Real_timer.h>
#include <CGAL/basic.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <CGAL/boost/graph/graph_traits_HalfedgeDS_default.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/config.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_utils_classes.h>

#include <BS_thread_pool.hpp>

#include <algorithm>
#include <deque>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <vector>

#define TRAIT_FUNC(ret_type, name, functor)                                                                            \
  template <class... T> static ret_type name(T&&... args) { return Traits().functor()(std::forward<T>(args)...); }

#define SOURCE_LOC __func__ << " (" << __LINE__ << ")"

namespace CGAL {
namespace Surface_Voronoi_diagram_with_star_metrics {
namespace PMP = Polygon_mesh_processing;

/*!
 * \brief Compute Voronoi diagrams on a surface mesh using star-metrics.
 *
 * This class generalizes the planar star-metric Voronoi diagram to surface meshes. Each site is
 * associated with a metric polyhedron (a discretized star-shaped metric). The class provides
 * operations to add sites/metrics, build the surface Voronoi diagram and inspect the resulting
 * face/edge/vertex information stored in a halfedge-based Voronoi graph.
 *
 * Template parameters mirror the Traits and property map types required to bridge the surface mesh,
 * metric polyhedron and the Voronoi graph.
 */
template <class Traits,
          class MeshVertexPointPMap = Default,
          class MeshFaceIndexPMap = Default,
          class MeshEdgeIndexPMap = Default,
          class MetricVertexPointPMap = Default,
          class MetricFaceIndexPMap = Default,
          class VoronoiDiagramVertexPointPMap = Default,
          class VoronoiDiagramVertexIndexPMap = Default>
class Surface_Voronoi_diagram_with_star_metrics
{
public:
#pragma region PublicTypes

  using Kernel = typename Traits::Kernel;
  using T = typename Traits::T;
  using FT = typename Traits::FT;
  using Point_3 = typename Traits::Point_3;
  using Vector_3 = typename Traits::Vector_3;
  using Ray_3 = typename Traits::Ray_3;
  using Plane_3 = typename Traits::Plane_3;
  using Line_3 = typename Traits::Line_3;
  using Parametric_line_3 = typename Traits::Parametric_line_3;

  using Surface_mesh = typename Traits::Surface_mesh;
  using Metric_polyhedron = typename Traits::Metric_polyhedron;
  using Voronoi_diagram_graph = typename Traits::Voronoi_diagram;
  using Metric_traits_data = typename Traits::Metric_traits_data;

  using index_t = std::ptrdiff_t;

  using mesh_graph_traits = typename boost::graph_traits<Surface_mesh>;
  using mesh_vertex_descriptor = typename mesh_graph_traits::vertex_descriptor;
  using mesh_face_descriptor = typename mesh_graph_traits::face_descriptor;
  using mesh_halfedge_descriptor = typename mesh_graph_traits::halfedge_descriptor;

  using metric_graph_traits = typename boost::graph_traits<Metric_polyhedron>;
  using metric_face_descriptor = typename metric_graph_traits::face_descriptor;
  using metric_face_iterator = typename metric_graph_traits::face_iterator;
  using metric_halfedge_descriptor = typename metric_graph_traits::halfedge_descriptor;
  using metric_edge_descriptor = typename metric_graph_traits::edge_descriptor;

  using vd_graph_traits = typename boost::graph_traits<Voronoi_diagram_graph>;
  using vd_vertex_descriptor = typename vd_graph_traits::vertex_descriptor;
  using vd_edge_descriptor = typename vd_graph_traits::edge_descriptor;
  using vd_halfedge_descriptor = typename vd_graph_traits::halfedge_descriptor;
  using vd_face_descriptor = typename vd_graph_traits::face_descriptor;

  using Face_component_property = CGAL::dynamic_face_property_t<index_t>;
  using Face_component_map = typename boost::property_map<Voronoi_diagram_graph, Face_component_property>::type;

  using Mesh_vertex_point_pmap =
      typename Default::Get<MeshVertexPointPMap,
                            typename boost::property_map<Surface_mesh, vertex_point_t>::const_type>::type;
  using Mesh_face_index_pmap =
      typename Default::Get<MeshFaceIndexPMap,
                            typename boost::property_map<Surface_mesh, face_index_t>::const_type>::type;
  using Mesh_edge_index_pmap =
      typename Default::Get<MeshEdgeIndexPMap,
                            typename boost::property_map<Surface_mesh, edge_index_t>::const_type>::type;

  using Metric_vertex_point_pmap =
      typename Default::Get<MetricVertexPointPMap,
                            typename boost::property_map<Metric_polyhedron, vertex_point_t>::const_type>::type;
  using Metric_face_index_pmap =
      typename Default::Get<MetricFaceIndexPMap,
                            typename boost::property_map<Metric_polyhedron, face_index_t>::const_type>::type;

  using Voronoi_diagram_vertex_point_pmap =
      typename Default::Get<VoronoiDiagramVertexPointPMap,
                            typename boost::property_map<Voronoi_diagram_graph, vertex_point_t>::const_type>::type;
  using Voronoi_diagram_vertex_index_pmap =
      typename Default::Get<VoronoiDiagramVertexIndexPMap,
                            typename boost::property_map<Voronoi_diagram_graph, vertex_index_t>::const_type>::type;

  static constexpr T INF = std::numeric_limits<T>::infinity();

  struct Site
  {
    Point_3 point;
    index_t metric_idx;
    bool enabled = true;
  };

  struct Cone_index
  {
    index_t site_idx;
    size_t face_idx;

    bool operator==(const Cone_index& other) const { return site_idx == other.site_idx && face_idx == other.face_idx; }
    bool operator>(const Cone_index& other) const {
      return site_idx > other.site_idx || (site_idx == other.site_idx && face_idx > other.face_idx);
    }

    // stream operator
    friend std::ostream& operator<<(std::ostream& os, const Cone_index& ci) {
      os << "Cone(S" << ci.site_idx << ", M" << ci.face_idx << ")";
      return os;
    }
  };

  struct Cone_index_hash
  {
    std::size_t operator()(const Cone_index& k) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, k.site_idx);
      boost::hash_combine(seed, k.face_idx);
      return seed;
    }
  };

  struct Internal_vertex_id
  {
    Cone_index k0, k1, k2;
    size_t mesh_face_idx;

    Internal_vertex_id(Cone_index ci0, Cone_index ci1, Cone_index ci2, index_t mesh_face_idx)
        : k0(ci0)
        , k1(ci1)
        , k2(ci2)
        , mesh_face_idx(mesh_face_idx) {
      if(k0 > k1)
        std::swap(k0, k1);
      if(k1 > k2)
        std::swap(k1, k2);
      if(k0 > k1)
        std::swap(k0, k1);
    }

    bool operator==(const Internal_vertex_id& other) const {
      return k0 == other.k0 && k1 == other.k1 && k2 == other.k2 && mesh_face_idx == other.mesh_face_idx;
    }
  };

  struct Internal_vertex_id_hash
  {
    std::size_t operator()(const Internal_vertex_id& v) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, Cone_index_hash{}(v.k0));
      boost::hash_combine(seed, Cone_index_hash{}(v.k1));
      boost::hash_combine(seed, Cone_index_hash{}(v.k2));
      boost::hash_combine(seed, v.mesh_face_idx);
      return seed;
    }
  };

  struct Boundary_vertex_id
  {
    Cone_index k0, k1;
    size_t mesh_halfedge_idx;

    Boundary_vertex_id(Cone_index ci0, Cone_index ci1, index_t mesh_halfedge_idx)
        : k0(ci0)
        , k1(ci1)
        , mesh_halfedge_idx(mesh_halfedge_idx) {
      if(k0 > k1)
        std::swap(k0, k1);
    }

    bool operator==(const Boundary_vertex_id& other) const {
      return k0 == other.k0 && k1 == other.k1 && mesh_halfedge_idx == other.mesh_halfedge_idx;
    }

    friend std::ostream& operator<<(std::ostream& os, const Boundary_vertex_id& v) {
      os << "Boundary(" << v.k0 << ", " << v.k1 << ", E" << v.mesh_halfedge_idx << ")";
      return os;
    }
  };

  struct Boundary_vertex_id_hash
  {
    std::size_t operator()(const Boundary_vertex_id& v) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, Cone_index_hash{}(v.k0));
      boost::hash_combine(seed, Cone_index_hash{}(v.k1));
      boost::hash_combine(seed, v.mesh_halfedge_idx);
      return seed;
    }
  };

  struct Bisector_segment_id
  {
    Cone_index k0, k1;
    size_t mesh_face_index;

    Bisector_segment_id(Cone_index ci0, Cone_index ci1, index_t mesh_face_index)
        : k0(ci0)
        , k1(ci1)
        , mesh_face_index(mesh_face_index) {
      if(k0 > k1)
        std::swap(k0, k1);
    }

    bool operator==(const Bisector_segment_id& other) const {
      return k0 == other.k0 && k1 == other.k1 && mesh_face_index == other.mesh_face_index;
    }

    friend std::ostream& operator<<(std::ostream& os, const Bisector_segment_id& v) {
      os << "Bisector(" << v.k0 << ", " << v.k1 << ", F" << v.mesh_face_index << ")";
      return os;
    }
  };

  struct Bisector_segment_id_hash
  {
    std::size_t operator()(const Bisector_segment_id& v) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, Cone_index_hash{}(v.k0));
      boost::hash_combine(seed, Cone_index_hash{}(v.k1));
      boost::hash_combine(seed, v.mesh_face_index);
      return seed;
    }
  };

  struct Cone_descriptor
  {
    index_t site_idx = -1;
    metric_face_descriptor face;

    bool operator==(const Cone_descriptor& other) const { return site_idx == other.site_idx && face == other.face; }

    bool is_valid() const { return site_idx >= 0 && face != metric_graph_traits::null_face(); }
  };

  enum Vertex_type : std::size_t {
    BOUNDARY,
    BOUNDARY_CONE,
    BOUNDARY_BISECTOR,
    TWO_SITE_BISECTOR,
    THREE_SITE_BISECTOR,
  };

  struct Boundary_vertex_info
  {
    mesh_vertex_descriptor vd;
    Cone_descriptor k;
  };

  struct Boundary_cone_info
  {
    mesh_halfedge_descriptor hd;
    Cone_descriptor k;
  };

  struct Boundary_bisector_info
  {
    mesh_halfedge_descriptor hd;
    Cone_descriptor k0, k1;
  };

  struct Two_site_bisector_info
  {
    mesh_halfedge_descriptor hd;
    Cone_descriptor k0;
    index_t c1;
    metric_halfedge_descriptor hd1;
  };

  struct Three_site_bisector_info
  {
    mesh_halfedge_descriptor hd;
    Cone_descriptor k0, k1, k2;
  };

  using Vertex_info = std::variant<Boundary_vertex_info,
                                   Boundary_cone_info,
                                   Boundary_bisector_info,
                                   Two_site_bisector_info,
                                   Three_site_bisector_info>;

  struct Bisector_edge_info
  {};

  struct Boundary_edge_info
  {
    mesh_halfedge_descriptor hd;
  };

  using Edge_info = std::variant<Bisector_edge_info, Boundary_edge_info>;

  struct Halfedge_info
  {
    Cone_descriptor k;

    /*!
     * \brief The face of surface mesh containing the halfedge.
     *
     * If the halfedge is on the boundary, the face corresponds to the mesh boundary halfedge.
     */
    mesh_face_descriptor mesh_fd;
  };

  struct Voronoi_diagram_data
  {
    using Vertex_info_property = CGAL::dynamic_vertex_property_t<Vertex_info>;
    using Vertex_info_map = typename boost::property_map<Voronoi_diagram_graph, Vertex_info_property>::type;
    using Vertex_normal_property = CGAL::dynamic_vertex_property_t<Vector_3>;
    using Vertex_normal_map = typename boost::property_map<Voronoi_diagram_graph, Vertex_normal_property>::type;
    using Edge_info_property = CGAL::dynamic_edge_property_t<Edge_info>;
    using Edge_info_map = typename boost::property_map<Voronoi_diagram_graph, Edge_info_property>::type;
    using Halfedge_info_property = CGAL::dynamic_halfedge_property_t<Halfedge_info>;
    using Halfedge_info_map = typename boost::property_map<Voronoi_diagram_graph, Halfedge_info_property>::type;

    Voronoi_diagram_graph graph;
    Voronoi_diagram_vertex_point_pmap vpm;
    Voronoi_diagram_vertex_index_pmap vertex_index_map;
    Vertex_info_map vertex_info_map;
    Vertex_normal_map vertex_normal_map;
    Edge_info_map edge_info_map;
    Halfedge_info_map halfedge_info_map;

    Voronoi_diagram_data()
        : vpm(get(vertex_point, graph))
        , vertex_index_map(get(vertex_index, graph))
        , vertex_info_map(get(Vertex_info_property{}, graph))
        , vertex_normal_map(get(Vertex_normal_property{}, graph))
        , edge_info_map(get(Edge_info_property{}, graph))
        , halfedge_info_map(get(Halfedge_info_property{}, graph)) {}

    /*!
     * \brief The copy constructor is deleted to avoid the property map ownership issue.
     */
    Voronoi_diagram_data(const Voronoi_diagram_data&) = delete;

    bool is_boundary_edge(vd_edge_descriptor ed) const {
      return std::holds_alternative<Boundary_edge_info>(get(edge_info_map, ed));
    }

    bool is_boundary_edge(vd_halfedge_descriptor hd) const { return is_boundary_edge(edge(hd, graph)); }

    bool is_bisector_edge(vd_edge_descriptor ed) const {
      return std::holds_alternative<Bisector_edge_info>(get(edge_info_map, ed));
    }

    bool is_bisector_edge(vd_halfedge_descriptor hd) const { return is_bisector_edge(edge(hd, graph)); }

    vd_vertex_descriptor add_vertex(const Point_3& p, const Vertex_info& info, const Vector_3& n = {}) {
      auto vd = CGAL::add_vertex(graph);
      put(vpm, vd, p);
      // put(vertex_index_map, vd, num_vertices(graph) - 1);
      put(vertex_info_map, vd, info);
      put(vertex_normal_map, vd, n);
      return vd;
    }

    void print_halfedge_loop(vd_vertex_descriptor vd) const {
      for(auto hd : halfedges_around_target(vd, graph)) {
        std::cout << source(hd, graph) << " ";
      }
    }

    vd_halfedge_descriptor connect(vd_vertex_descriptor v0,
                                   vd_vertex_descriptor v1,
                                   vd_face_descriptor fd01,
                                   vd_face_descriptor fd10,
                                   Edge_info info,
                                   Halfedge_info hinfo01,
                                   Halfedge_info hinfo10) {
      auto n0 = get(vertex_normal_map, v0);
      auto n1 = get(vertex_normal_map, v1);
      auto hd01 = connect_vertices_3(v0, v1, graph, n0, n1, vpm, Kernel{});
      auto hd10 = opposite(hd01, graph);
      auto ed = edge(hd01, graph);

      set_face(hd01, fd01, graph);
      set_face(hd10, fd10, graph);
      put(edge_info_map, ed, std::move(info));
      put(halfedge_info_map, hd01, std::move(hinfo01));
      put(halfedge_info_map, hd10, std::move(hinfo10));
      return hd01;
    }

    /*!
     * \brief Remove halfedge `hd` from the halfedge loop around the target vertex `target(hd, graph)`.
     * Return true if the target vertex is not orphaned after the removal.
     * \param hd
     * \return true
     * \return false
     */
    bool remove_halfedge(vd_halfedge_descriptor hd) {
      // Find the previous halfedge in the target vertex loop
      auto hd_cur = hd;
      for(;;) {
        auto hd_next = opposite(next(hd_cur, graph), graph);
        if(hd_next == hd)
          break;
        hd_cur = hd_next;
      }
      auto vt = target(hd, graph);
      if(hd_cur == hd) {
        // hd is the only halfedge around the target vertex
        set_halfedge(vt, vd_graph_traits::null_halfedge(), graph);
        return false;
        // if (remove_orphaned_vertex)
        //     remove_vertex(vt, graph);
        // else
        //     set_halfedge(vt, vd_graph_traits::null_halfedge(), graph);
      } else {
        // Skip hd in the halfedge around target loop
        set_next(hd_cur, next(hd, graph), graph);
        if(halfedge(vt, graph) == hd)
          set_halfedge(vt, hd_cur, graph);
        return true;
      }
    }

    /*!
     * \brief Remove `halfedge(ed, graph)` and `opposite(halfedge(ed, graph), graph)` from the loop of halfedge
     * around target. This does not remove the edge but isolate the edge from the graph.
     * \param ed
     */
    void detach_edge(vd_edge_descriptor ed) {
      auto hd = halfedge(ed, graph), ohd = opposite(hd, graph);
      remove_halfedge(hd);
      remove_halfedge(ohd);
      set_next(hd, ohd, graph);
      set_next(ohd, hd, graph);
    }

    /*!
     * \brief Detach the edge from the graph then remove the edge.
     *
     * \param ed
     * \param remove_orphaned_vertex
     */
    void remove_edge(vd_edge_descriptor ed) {
      detach_edge(ed);
      CGAL::remove_edge(ed, graph);
    }

    void flood_fill_faces(auto& map,
                          const auto& default_value,
                          const auto& seed_faces,
                          bool flood_through_bisector = false,
                          bool flood_through_boundary = true) {
      std::queue<vd_face_descriptor> queue;

      for(auto fd : seed_faces) {
        assert(get(map, fd) != default_value);
        queue.push(fd);
      }

      while(!queue.empty()) {
        auto fd = queue.front();
        queue.pop();
        for(auto hd : CGAL::halfedges_around_face(halfedge(fd, graph), graph)) {
          if(!flood_through_bisector && is_bisector_edge(edge(hd, graph)))
            continue;
          if(!flood_through_boundary && is_boundary_edge(edge(hd, graph)))
            continue;

          auto nfd = face(opposite(hd, graph), graph);
          if(get(map, nfd) == default_value) {
            put(map, nfd, get(map, fd));
            queue.push(nfd);
          }
        }
      }
    }

    template <class ComponentMap>
    typename boost::property_traits<ComponentMap>::value_type trace_components(ComponentMap c) const {
      using comp_t = typename boost::property_traits<ComponentMap>::value_type;
      const comp_t UNVISITED = std::numeric_limits<comp_t>::max();

      // Initialize all faces to UNVISITED
      for(auto fd : faces(graph)) {
        put(c, fd, UNVISITED);
      }

      comp_t num_components = 0;
      std::deque<vd_face_descriptor> queue;

      for(auto fd : faces(graph)) {
        if(get(c, fd) != UNVISITED)
          continue;
        put(c, fd, num_components);
        flood_fill_faces(c, UNVISITED, std::vector<vd_face_descriptor>{fd}, false, true);
        ++num_components;
      }

      return num_components;
    }

    auto trace_components() {
      Face_component_map map{get(Face_component_property{}, graph)};
      std::size_t num_components = trace_components(map);
      return std::make_pair(num_components, map);
    }

    auto face_bounded_side(vd_face_descriptor fd, const Point_3& p) const {
      auto n = PMP::compute_face_normal(fd, graph);
      std::vector<Point_3> vertices;
      for(auto vd : vertices_around_face(halfedge(fd, graph), graph)) {
        vertices.push_back(get(vpm, vd));
      }
      return CGAL::bounded_side_2(vertices.begin(), vertices.end(), p, CGAL::Projection_traits_3<Kernel>(n));
    }

    bool check_halfedge_info(bool has_boundary_edges = true) const {
      for(auto fd : faces(graph)) {
        auto cell = halfedge(fd, graph);
        auto info = get(halfedge_info_map, cell);
        auto site = info.k.site_idx;
        auto mesh_fd = info.mesh_fd;
        for(auto hd : halfedges_around_face(halfedge(fd, graph), graph)) {
          auto hd_info = get(halfedge_info_map, hd);
          if(info.k.site_idx != site)
            return false;
          if(has_boundary_edges && info.mesh_fd != mesh_fd)
            return false;
        }
      }
      return true;
    }

    bool is_valid_polygon_mesh() const {
      for(auto fd : faces(graph)) {
        auto n = PMP::compute_face_normal(fd, graph);
        auto dt = Projection_traits_3<Kernel>(n);
        std::vector<Point_3> vertices;
        for(auto vd : vertices_around_face(halfedge(fd, graph), graph)) {
          vertices.push_back(get(vpm, vd));
        }
        if(!CGAL::is_simple_2(vertices.begin(), vertices.end(), dt)) {
          return false;
        }
        if(CGAL::orientation_2(vertices.begin(), vertices.end(), dt) != COUNTERCLOCKWISE) {
          return false;
        }
      }
      return true;
    }

    template <class... Ts> struct overloaded : Ts...
    {
      using Ts::operator()...;
    };

    // explicit deduction guide (not needed as of C++20)
    template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

    bool check_topology() const {
      for(auto vd : vertices(graph)) {
        std::visit(overloaded{
                       [&](const Boundary_vertex_info& info) {
                         for(auto hd : halfedges_around_target(vd, graph)) {
                           CGAL_assertion(is_boundary_edge(hd));
                         }
                       },
                       [&](const Boundary_cone_info& info) {
                         CGAL_assertion(degree(vd, graph) == 2);
                         auto hd0 = halfedge(vd, graph), hd1 = next(hd0, graph);
                         CGAL_assertion(is_boundary_edge(hd0) && is_boundary_edge(hd1));
                       },
                       [&](const Boundary_bisector_info& info) {
                         CGAL_assertion(degree(vd, graph) == 4 || degree(vd, graph) == 3);
                         int n_boundary = 0, n_bisector = 0;
                         for(auto hd : halfedges_around_target(vd, graph)) {
                           if(is_bisector_edge(hd))
                             ++n_bisector;
                           if(is_boundary_edge(hd))
                             ++n_boundary;
                         }
                         CGAL_assertion(n_boundary == 2 && (n_bisector == 2 || n_bisector == 1));
                       },
                       [&](const Two_site_bisector_info& info) {
                         CGAL_assertion(degree(vd, graph) == 2);
                         auto hd0 = halfedge(vd, graph), hd1 = next(hd0, graph);
                         CGAL_assertion(is_bisector_edge(hd0) && is_bisector_edge(hd1));
                       },
                       [&](const Three_site_bisector_info& info) {
                         CGAL_assertion(degree(vd, graph) == 3);

                         for(auto hd : halfedges_around_target(vd, graph)) {
                           CGAL_assertion(is_bisector_edge(hd));
                         }
                       },
                   },
                   get(vertex_info_map, vd));
      }
      return true;
    }

    index_t cell_site(vd_face_descriptor fd) const { return get(halfedge_info_map, halfedge(fd, graph)).k.site_idx; }

    void prune_orphan_vertices() {
      std::vector<vd_vertex_descriptor> orphans;
      for(auto vd : vertices(graph)) {
        if(degree(vd, graph) == 0) {
          orphans.push_back(vd);
        }
      }
      for(auto vd : orphans) {
        CGAL::remove_vertex(vd, graph);
      }
    }
  };

  using const_voronoi_diagram_ptr = std::shared_ptr<const Voronoi_diagram_data>;
  using Voronoi_diagram_ptr = std::shared_ptr<Voronoi_diagram_data>;

  struct Internal_trace
  {
    Parametric_line_3 bisect_line;
    Plane_3 bisect_plane, face_plane;
    mesh_halfedge_descriptor face_hd;
    mesh_halfedge_descriptor prev_hd;
    Cone_descriptor k0, k1, k_prev;
    /*!
     * \brief When the previous itrace goes through a cone and enter the neighboring cone,
     * the metric_prev_hd is the halfedge of current metric that the previous itrace enter the cone.
     */
    metric_halfedge_descriptor metric_prev_hd;
    vd_vertex_descriptor v_vd;
  };

  struct Metric_data
  {
    Metric_polyhedron graph;
    Metric_vertex_point_pmap vpm;
    Metric_face_index_pmap face_index_map;
    Metric_traits_data data;

    Metric_data(Metric_polyhedron&& graph, Metric_vertex_point_pmap&& vpm, Metric_face_index_pmap&& fim)
        : graph(std::move(graph))
        , vpm(std::move(vpm))
        , face_index_map(std::move(fim))
        , data(std::move(construct_metric_traits_data(this->graph))) {}

    Metric_data(Metric_polyhedron&& graph)
        : graph(std::move(graph))
        , vpm(std::move(get(vertex_point, this->graph)))
        , face_index_map(std::move(get(face_index, this->graph)))
        , data(std::move(construct_metric_traits_data(this->graph))) {}

    auto cone_face_bases(metric_halfedge_descriptor hd) const {
      auto v0 = get(vpm, source(hd, graph)) - ORIGIN;
      auto v1 = get(vpm, target(hd, graph)) - ORIGIN;
      return std::make_pair(v0, v1);
    }

    /*!
     * \brief Return the normal of the 2D cone face, pointing inward
     *
     * \param hd
     * \return Vector_3
     */
    Vector_3 cone_face_orthogonal_vector(metric_halfedge_descriptor hd) const {
      // TODO: cache the result
      auto [v0, v1] = cone_face_bases(hd);
      Vector_3 n = cross_product(v0, v1);
      return n;
    }
  };
#pragma endregion

#pragma region PrivateTypes
protected:
  struct Cone_line_intersection
  {
    T t_min, t_max;
    metric_halfedge_descriptor ed_min, ed_max;
  };

  struct Segment_cone_intersections;

  class Segment_cone_intersections_iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = typename Segment_cone_intersections::value_type;
    using difference_type = std::ptrdiff_t;
    using iterator = Segment_cone_intersections_iterator;

    Segment_cone_intersections_iterator(Segment_cone_intersections* data, std::size_t idx)
        : data(data)
        , idx(idx) {}

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
    Segment_cone_intersections* data = nullptr;
    std::size_t idx;
  };

  struct Segment_cone_intersections
  {
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
      if(idx_map.empty()) {
        for(std::size_t i = 0; i < fds.size(); ++i) {
          idx_map[fds[i]] = i;
        }
      }

      auto it = idx_map.find(fd);
      if(it == idx_map.end()) {
        return end();
      }

      return iterator(this, it->second);
    }

    void clear() {
      ts.clear();
      eds.clear();
      fds.clear();
      idx_map.clear();
    }

    void add_face(metric_face_descriptor fd) { fds.push_back(fd); }

    void add_intersection(T t, metric_halfedge_descriptor hd) {
      CGAL_assertion(empty() || t >= t_max());
      ts.push_back(t);
      eds.push_back(hd);
    }

    void clip(T tmin, T tmax, Segment_cone_intersections& res) const {
      res.clear();

      // Empty cases
      if(empty() || tmin > tmax || tmin > t_max() || tmax < t_min())
        return;

      auto it0 = std::lower_bound(ts.begin(), ts.end(), tmin);
      auto it1 = std::upper_bound(ts.begin(), ts.end(), tmax);

      if(tmin < *it0) {
        // tmin is in range of it0 - 1 and it0
        res.add_intersection(tmin, metric_graph_traits::null_halfedge());
        auto i0 = std::distance(ts.begin(), it0);
        CGAL_assertion(i0 > 0);
        res.add_face(fds[i0 - 1]);
      }

      for(auto it = it0; it != it1; ++it) {
        auto i = std::distance(ts.begin(), it);
        res.add_intersection(*it, eds[i]);
        if(*it < tmax)
          res.add_face(fds[i]);
      }

      if(res.ts.empty() || tmax > res.ts.back()) {
        // tmax is in range of it1 and it1 + 1
        res.add_intersection(tmax, metric_graph_traits::null_halfedge());
      }
    }

  private:
    std::vector<T> ts;                           // N+1
    std::vector<metric_halfedge_descriptor> eds; // N+1
    std::vector<metric_face_descriptor> fds;     // N

    mutable std::unordered_map<metric_face_descriptor, size_t> idx_map;
  };

  struct Disconnected_component
  {
    std::unordered_set<vd_face_descriptor> faces;
    std::unordered_set<index_t> neighbor_sites;
    std::unordered_set<vd_halfedge_descriptor> trace_inward_edges;
    std::unordered_set<vd_edge_descriptor> boundary_edges;
    std::unordered_set<vd_edge_descriptor> bisector_edges;
  };
#pragma endregion

#pragma region PublicInterface
public:
  Surface_Voronoi_diagram_with_star_metrics(const Surface_mesh& mesh,
                                            Mesh_vertex_point_pmap vpm,
                                            Mesh_face_index_pmap face_index_map,
                                            Mesh_edge_index_pmap edge_index_map,
                                            Traits traits = Traits())
      : mesh(mesh)
      , vpm(std::move(vpm))
      , face_index_map(std::move(face_index_map))
      , edge_index_map(std::move(edge_index_map))
      , traits(traits) {
    reset();
  }

  Surface_Voronoi_diagram_with_star_metrics(const Surface_mesh& mesh,
                                            Traits traits = Traits())
      : Surface_Voronoi_diagram_with_star_metrics(
            mesh, get(vertex_point, mesh), get(face_index, mesh), get(edge_index, mesh), traits) {}

  void add_site(const Point_3& p, index_t metric_idx) { m_sites.push_back({p, metric_idx}); }

  /*!
   * \brief Add a metric attached to sites added later.
   * To avoid copy, use add_metric(std::move(...), ...) to move the metric.
   *
   * \param m
   * \param vpm
   * \param fim
   * \return index_t
   */
  index_t add_metric(Metric_polyhedron m, Metric_vertex_point_pmap vpm, Metric_face_index_pmap fim) {
    auto idx = m_metrics.size();
    m_metrics.emplace_back(std::move(m), std::move(vpm), std::move(fim));
    return idx;
  }

  index_t add_metric(Metric_polyhedron m) {
    auto idx = m_metrics.size();
    m_metrics.emplace_back(std::move(m));
    return idx;
  }

  void reserve_sites(std::size_t n) { m_sites.reserve(n); }

  void reserve_metrics(std::size_t n) { m_metrics.reserve(n); }

  void clear_sites() { m_sites.clear(); }

  void clear_metrics() { m_metrics.clear(); }

  std::size_t num_sites() const { return m_sites.size(); }

  std::size_t num_metrics() const { return m_metrics.size(); }

  std::size_t num_i_traces() const { return i_traces.size(); }

  auto site_cbegin() const { return m_sites.begin(); }

  auto site_cend() const { return m_sites.end(); }

  auto metric_cbegin() const { return m_metrics.begin(); }

  auto metric_cend() const { return m_metrics.end(); }

  auto& metrics() { return m_metrics; }

  auto& metrics() const { return m_metrics; }

  auto i_trace_cbegin() const { return i_traces.begin(); }

  auto i_trace_cend() const { return i_traces.end(); }

  const_voronoi_diagram_ptr voronoi_diagram_ptr() const { return voronoi; }

  Voronoi_diagram_ptr voronoi_diagram_ptr() { return voronoi; }

  const Voronoi_diagram_data& voronoi_diagram() const { return *voronoi; }

  Voronoi_diagram_data& voronoi_diagram() { return *voronoi; }

  const Real_timer& i_timer() const { return i_trace_timer; }

  const Real_timer& b_timer() const { return b_trace_timer; }

  bool read_sites(std::istream& is) {
    std::size_t n_metrics;
    is >> n_metrics;

    clear_metrics();
    reserve_metrics(n_metrics);
    for(size_t i = 0; i < n_metrics; ++i) {
      Metric_polyhedron P;
      if(!IO::read_OFF(is, P))
        return false;
      add_metric(P);
    }

    std::size_t n_sites;
    is >> n_sites;
    clear_sites();
    reserve_sites(n_sites);
    for(size_t i = 0; i < n_sites; ++i) {
      double x, y, z;
      index_t site_idx;
      is >> x >> y >> z >> site_idx;
      add_site(Point_3(x, y, z), site_idx);
    }

    return true;
  }

  bool write_sites(std::ostream& os) const {
    os << num_metrics() << std::endl;
    for(const auto& m : m_metrics) {
      if(!IO::write_OFF(os, m.graph))
        return false;
    }

    os << num_sites() << std::endl;
    for(const auto& s : m_sites) {
      os << s.point.x() << " " << s.point.y() << " " << s.point.z() << " " << s.metric_idx << std::endl;
    }

    return true;
  }

  std::pair<Site&, Metric_data&> site(index_t idx) {
    auto& c = m_sites[idx];
    return {c, m_metrics[c.metric_idx]};
  }

  std::pair<const Site&, const Metric_data&> site(index_t idx) const {
    auto& c = m_sites[idx];
    return {c, m_metrics[c.metric_idx]};
  }

  void disable_site(index_t idx) { m_sites[idx].enabled = false; }

  void enable_site(index_t idx) { m_sites[idx].enabled = true; }

  void disable_all_sites() {
    for(auto& s : m_sites) {
      s.enabled = false;
    }
  }

  void enable_all_sites() {
    for(auto& s : m_sites) {
      s.enabled = true;
    }
  }

  void disable_sites(auto it_begin, auto it_end) {
    for(auto it = it_begin; it != it_end; ++it) {
      disable_site(*it);
    }
  }

  void enable_sites(auto it_begin, auto it_end) {
    for(auto it = it_begin; it != it_end; ++it) {
      enable_site(*it);
    }
  }

  T find_nearest_site(const Point_3& p, Cone_descriptor& m_cone) const {
    T d_min = INF;

    for(index_t i = 0; i < m_sites.size(); ++i) {
      auto [c, m] = site(i);
      if(!c.enabled)
        continue;

      auto res = metric_any_intersection(m.data, p - c.point);
      if(!res)
        continue;
      auto [pm, fd] = *res;

      FT weight = 1.0 / (pm - ORIGIN).squared_length();

      // auto it = m_weights.find(&metrics[m_idx].graph);
      // if (it == m_weights.end()) {
      //     continue;
      // }
      auto d = approximate_sqrt((p - c.point).squared_length() * weight);
      if(d < d_min) {
        d_min = d;
        m_cone.site_idx = i;
        m_cone.face = fd;
      }
    }
    CGAL_assertion(d_min < INF);
    return d_min;
  }

  auto trace_boundary(mesh_halfedge_descriptor bh,
                      Cone_descriptor k0,
                      vd_vertex_descriptor prev_vd,
                      bool same_side = true,
                      bool opposite_side = true,
                      bool add_border_edges = false,
                      bool add_cone_vertices = false,
                      vd_face_descriptor fd0 = vd_graph_traits::null_face(),
                      vd_face_descriptor fd1 = vd_graph_traits::null_face(),
                      FT t_min = 0,
                      FT t_max = 1) {
    CGAL_precondition(k0.is_valid());
    CGAL_precondition(t_max >= t_min);

    vout << IO::level(1) << SOURCE_LOC << ": tracing boundary " << bh << " with cone " << cone_index(k0) << std::endl;

    auto b_line = mesh_edge_segment(bh);
    // FT t_min = 0, t_max = 1;

    auto bh_opposite = opposite(bh, mesh);
    if(is_border(bh, mesh)) {
      same_side = false;
    }
    if(is_border(bh_opposite, mesh)) {
      opposite_side = false;
    }

    Vector_3 b_normal(NULL_VECTOR);
    Plane_3 b_plane, b_plane_opposite;
    if(same_side) {
      b_plane = mesh_face_plane(bh);
      b_normal = b_plane.orthogonal_vector();
    }
    if(opposite_side) {
      b_plane_opposite = mesh_face_plane(bh_opposite);
      b_normal += b_plane_opposite.orthogonal_vector();
    }

    find_all_segment_cone_intersections(b_line, t_min, t_max, isects_vec_cache);
    auto& isects = isects_vec_cache;

    // Sub-interval of the boundary segment that intersects k0
    auto isect_iter = isects[k0.site_idx].find_face(k0.face);
    CGAL_assertion_msg(isect_iter != isects[k0.site_idx].end(), "Boundary must intersect the cone");

    // Find all boundary-cone / boundary-bisector intersections
    Cone_descriptor k_prev;
    for(;;) {
      // Find nearest intersection of 2-site bisector plane with the boundary segment
      T dist_min = INF;
      T tb_min;
      Plane_3 bi_plane_min;
      Cone_descriptor k1_min;

      for(index_t k1_site = 0; k1_site < m_sites.size(); ++k1_site) {
        if(k1_site == k0.site_idx) {
          continue;
        }

        if(!m_sites[k1_site].enabled) {
          continue;
        }

        // Iterate over all intervals of k1_site that intersect current interval of k0
        auto isect = *isect_iter;
        isects[k1_site].clip(max(isect.second.t_min, t_min), min(isect.second.t_max, t_max), isects_cache);
        auto& isects_overlap = isects_cache;
        for(auto [fd, isect_overlap] : isects_overlap) {
          Cone_descriptor k1{k1_site, fd};
          if(k1 == k_prev)
            continue;

          Plane_3 bi_plane = get_bisect_plane(k0, k1);
          if(bi_plane.is_degenerate())
            continue;

          auto tb = intersect(b_line, bi_plane);
          if(!tb)
            continue;

          // TODO: marginal case: tb is on the terminal points
          if(*tb < isect_overlap.t_min || *tb >= isect_overlap.t_max)
            continue;
          T dist = *tb - t_min;
          CGAL_assertion_msg(dist >= 0, "Intersection must be on the segment");
          if(dist < dist_min) {
            dist_min = dist;
            tb_min = *tb;
            bi_plane_min = bi_plane;
            k1_min = k1;
          }
        }
      }

      if(dist_min < INF) {
        // A 2-site bisector intersects the boundary segment
        Boundary_vertex_id bvid(cone_index(k0), cone_index(k1_min), get(edge_index_map, edge(bh, mesh)));
        auto pt_start = b_line(tb_min);
        vd_vertex_descriptor v_vd;

        {
#ifndef CGAL_HAS_NO_THREADS
          std::lock_guard lock(vd_mutex);
#endif

          if(auto vd_it = b_vert_map.find(bvid); vd_it != b_vert_map.cend()) {
            v_vd = vd_it->second;
          } else {
            v_vd = voronoi->add_vertex(pt_start, Boundary_bisector_info{bh, k0, k1_min}, b_normal);
            b_vert_map[bvid] = v_vd;
          }
        }

        vout << IO::level(1) << SOURCE_LOC << ": boundary 2-site bisector " << bvid << " intersects at point "
             << pt_start << std::endl;

        if(add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
#ifndef CGAL_HAS_NO_THREADS
          std::lock_guard lock(vd_mutex);
#endif
          voronoi->connect(prev_vd, v_vd, fd0, fd1, Boundary_edge_info{bh}, Halfedge_info{k0, face(bh, mesh)},
                           Halfedge_info{k0, face(bh_opposite, mesh)});
        }
        prev_vd = v_vd;

        if(same_side) {
          auto bisect_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane.orthogonal_vector());
          auto orient = orientation(b_plane.orthogonal_vector(), b_line.d(), bisect_dir);
          auto bisect_line = construct_parametric_line(pt_start, orient == POSITIVE ? bisect_dir : -bisect_dir);

#ifndef CGAL_HAS_NO_THREADS
          std::lock_guard lock(i_trace_mutex);
#endif
          i_traces.push_back({
              bisect_line,
              bi_plane_min,
              b_plane,
              bh,
              bh,
              k0,
              k1_min,
              {},
              metric_graph_traits::null_halfedge(),
              v_vd,
          });
        }
        if(opposite_side) {
          auto bisect_dir = cross_product(bi_plane_min.orthogonal_vector(), b_plane_opposite.orthogonal_vector());
          auto orient = orientation(b_plane_opposite.orthogonal_vector(), -b_line.d(), bisect_dir);
          auto bisect_line = construct_parametric_line(pt_start, orient == POSITIVE ? bisect_dir : -bisect_dir);

#ifndef CGAL_HAS_NO_THREADS
          std::lock_guard lock(i_trace_mutex);
#endif
          i_traces.push_back({
              bisect_line,
              bi_plane_min,
              b_plane_opposite,
              bh_opposite,
              bh_opposite,
              k1_min,
              k0,
              {},
              metric_graph_traits::null_halfedge(),
              v_vd,
          });
        }

        // Switch current cone to new cone k1
        t_min = tb_min;
        k_prev = k0;
        k0 = k1_min;

        isect_iter = isects[k0.site_idx].find_face(k0.face);
        CGAL_assertion_msg(isect_iter != isects[k0.site_idx].end(), "Boundary must intersect the cone");
      } else {
        // The bondary segment move to the neighboring cone
        isect_iter++;
        // If this is the last interval
        if(isect_iter == isects[k0.site_idx].end())
          break;

        // Otherwise switch to the next cone
        auto isect = *isect_iter;
        if(add_cone_vertices) {
#ifndef CGAL_HAS_NO_THREADS
          std::lock_guard lock(vd_mutex);
#endif
          auto v_vd = voronoi->add_vertex(b_line(isect.second.t_min), Boundary_cone_info{bh, k0}, b_normal);
          if(add_border_edges && prev_vd != vd_graph_traits::null_vertex()) {
            voronoi->connect(prev_vd, v_vd, fd0, fd1, Boundary_edge_info{bh}, Halfedge_info{k0, face(bh, mesh)},
                             Halfedge_info{k0, face(bh_opposite, mesh)});
          }
          prev_vd = v_vd;
        }

        k_prev = k0;
        k0.face = isect.first;
        t_min = isect.second.t_min;
      }
    }
    vout << IO::level(1) << SOURCE_LOC << ": boundary trace leaves from cone " << cone_index(k0) << std::endl;
    return std::make_pair(k0, prev_vd);
  }

  // void trace_all_boundaries(mesh_vertex_descriptor vd, bool add_border_edges = true, bool add_cone_vertices =
  // false) {
  //     vout << IO::level(1) << SOURCE_LOC << ": start" << std::endl;

  //     Cone_descriptor k0;
  //     find_nearest_site(get(vpm, vd), k0);

  //     using edge_bool_t = CGAL::dynamic_edge_property_t<bool>;
  //     using edge_visited_map = typename boost::property_map<Surface_mesh, edge_bool_t>::type;

  //     auto edge_visited = get(edge_bool_t{}, mesh);

  //     std::vector<std::pair<mesh_halfedge_descriptor, Cone_descriptor>> queue;
  //     for (auto hd : halfedges_around_source(vd, mesh)) {
  //         queue.emplace_back(hd, k0);
  //         // put(edge_visited, edge(hd, mesh), true);
  //     }

  //     while (!queue.empty()) {
  //         auto [hd, k0] = queue.back();
  //         queue.pop_back();

  //         if (get(edge_visited, edge(hd, mesh))) continue;

  //         if (is_border(hd, mesh)) {
  //             // hd is a border halfedge, trace the boundary loop
  //             auto bhd = hd;
  //             auto kb = k0;  // Current cone
  //             vd_vertex_descriptor prev_vd = vd_graph_traits::null_vertex(), vd0;
  //             do {
  //                 if (add_border_edges) {
  //                     auto b_plane = mesh_face_plane(opposite(bhd, mesh));
  //                     auto v_vd = voronoi->add_vertex(get(vpm, source(bhd, mesh)), Boundary_vertex_info{hd, k0},
  //                                                     b_plane.orthogonal_vector());
  //                     if (prev_vd != vd_graph_traits::null_vertex()) {
  //                         voronoi->connect(prev_vd, v_vd, vd_graph_traits::null_face(), dummy_face);
  //                     } else {
  //                         vd0 = v_vd;
  //                     }
  //                     prev_vd = v_vd;
  //                 }

  //                 std::tie(kb, prev_vd) = trace_boundary(bhd, kb, prev_vd, false, true, add_border_edges,
  //                                                        add_cone_vertices, vd_graph_traits::null_face(),
  //                                                        dummy_face);
  //                 put(edge_visited, edge(bhd, mesh), true);
  //                 bhd = next(bhd, mesh);
  //             } while (bhd != hd);
  //             if (add_border_edges) voronoi->connect(prev_vd, vd0, vd_graph_traits::null_face(), dummy_face);
  //         } else if (is_border(opposite(hd, mesh), mesh)) {
  //             continue;  // We handle the border halfedge in the previous case
  //         } else {
  //             auto [k, _] = trace_boundary(hd, k0, vd_graph_traits::null_vertex(), true, true, false, false);
  //             put(edge_visited, edge(hd, mesh), true);
  //             for (auto hd_inner : halfedges_around_source(target(hd, mesh), mesh)) {
  //                 if (get(edge_visited, edge(hd_inner, mesh))) continue;
  //                 queue.emplace_back(hd_inner, k);
  //                 // put(edge_visited, edge(hd_inner, mesh), true);
  //             }
  //         }
  //     }
  // }

  void trace_all_boundaries(mesh_vertex_descriptor vd,
                            bool add_mesh_vertices = true,
                            bool add_boundary_edges = true,
                            bool add_cone_vertices = false) {
    b_trace_timer.start();

    using vd_vertex_t = CGAL::dynamic_vertex_property_t<vd_vertex_descriptor>;
    using vd_vertex_map = typename boost::property_map<Surface_mesh, vd_vertex_t>::type;
    auto vd_map = get(vd_vertex_t{}, mesh);

    if(add_mesh_vertices) {
      for(auto vd : vertices(mesh)) {
        auto v_vd = voronoi->add_vertex(get(vpm, vd), Boundary_vertex_info{vd, Cone_descriptor{}},
                                        PMP::compute_vertex_normal(vd, mesh, parameters::vertex_point_map(vpm)));
        put(vd_map, vd, v_vd);
      }
    }

    Cone_descriptor k0;
    find_nearest_site(get(vpm, vd), k0);

    using edge_bool_t = CGAL::dynamic_edge_property_t<bool>;
    using edge_visited_map = typename boost::property_map<Surface_mesh, edge_bool_t>::type;

    auto edge_visited = get(edge_bool_t{}, mesh);

    std::function<void(mesh_halfedge_descriptor, Cone_descriptor)> trace;
    std::mutex edge_visited_mutex;

    // Non-recursive trace implementation using an explicit stack.
    trace = [&](mesh_halfedge_descriptor start_hd, Cone_descriptor start_k0) {
      try {
        // Local stack for DFS-like traversal when running single-threaded.
        std::vector<std::pair<mesh_halfedge_descriptor, Cone_descriptor>> stack;
        stack.emplace_back(start_hd, start_k0);

        while(!stack.empty()) {
          {
#ifndef CGAL_HAS_NO_THREADS
            std::lock_guard lock(exception_mutex);
#endif
            if(exception)
              return;
          }

          auto [hd, k0] = stack.back();
          stack.pop_back();

          // Check and mark visited atomically
          {
            std::lock_guard lock(edge_visited_mutex);
            if(get(edge_visited, edge(hd, mesh)))
              continue;
            put(edge_visited, edge(hd, mesh), true);
          }

          auto prev_vd = vd_graph_traits::null_vertex();
          if(add_mesh_vertices) {
            prev_vd = get(vd_map, source(hd, mesh));
          }

          bool same_side = !is_border(hd, mesh);
          bool opposite_side = !is_border(opposite(hd, mesh), mesh);
          auto fd01 = same_side ? dummy_face : vd_graph_traits::null_face();
          auto fd10 = opposite_side ? dummy_face : vd_graph_traits::null_face();

          auto [k, v_vd] = trace_boundary(hd, k0, prev_vd, same_side, opposite_side, add_boundary_edges,
                                          add_cone_vertices, fd01, fd10);

          if(add_boundary_edges && add_mesh_vertices) {
#ifndef CGAL_HAS_NO_THREADS
            std::lock_guard lock(vd_mutex);
#endif
            voronoi->connect(v_vd, get(vd_map, target(hd, mesh)), fd01, fd10, Boundary_edge_info{hd},
                             Halfedge_info{k, face(hd, mesh)}, Halfedge_info{k, face(opposite(hd, mesh), mesh)});
          }

          // Iterate neighbors. If multi-threaded, spawn a detached task per neighbor so each runs its own stack;
          // otherwise push neighbor onto local stack to continue iteratively.
          for(auto hd_inner : halfedges_around_source(target(hd, mesh), mesh)) {
#ifndef CGAL_HAS_NO_THREADS
            detach_task(std::bind(trace, hd_inner, k));
#else
            stack.emplace_back(hd_inner, k);
#endif
          }
        }
      } catch(...) {
#ifndef CGAL_HAS_NO_THREADS
        std::lock_guard lock(exception_mutex);
#endif
        exception = std::current_exception();
      }
    };

    for(auto hd : halfedges_around_source(vd, mesh)) {
#ifndef CGAL_HAS_NO_THREADS
      detach_task(std::bind(trace, hd, k0));
#else
      trace(hd, k0);
#endif
    }

// wait for all tasks to finish
#ifndef CGAL_HAS_NO_THREADS
    wait_for_tasks();
#endif

    if(exception)
      std::rethrow_exception(exception);

    b_trace_timer.stop();
  }

  void
  trace_all_boundaries(bool add_mesh_vertices = true, bool add_boundary_edges = true, bool add_cone_vertices = false) {
    if(is_empty(mesh))
      return;
    trace_all_boundaries(*vertices(mesh).first, add_mesh_vertices, add_boundary_edges, add_cone_vertices);
  }

  void reload() { vpm = get(vertex_point, mesh); }

  void reset() {
    voronoi = std::make_shared<Voronoi_diagram_data>();
    dummy_face = add_face(voronoi->graph);

    i_traces.clear();
    vert_map.clear();
    b_vert_map.clear();
    processed_b_verts.clear();
    bounding_vertices.clear();

    if(b_trace_timer.is_running())
      b_trace_timer.stop();
    if(i_trace_timer.is_running())
      i_trace_timer.stop();
    b_trace_timer.reset();
    i_trace_timer.reset();

    exception = nullptr;
  }

  std::optional<Internal_trace> step() {
    if(i_traces.empty()) {
      return std::nullopt;
    }
    auto tr = std::move(i_traces.back());
    i_traces.pop_back();
    i_trace_timer.start();

    // pool.detach_task([=, this]() { process_i_trace(tr, false); });
    // pool.wait();
    process_i_trace(tr, false);

    i_trace_timer.stop();
    return tr;
  }

  void process_i_traces() {
    i_trace_timer.start();
    for(auto& trace : i_traces) {
      detach_task([=, this]() { process_i_trace(trace); });
    }
    wait_for_tasks();
    i_traces.clear();
    i_trace_timer.stop();
  }

  void process_i_traces_single_thread() {
    while(step())
      ;
  }

  void build(bool add_mesh_vertices = true, bool add_boundary_edges = true, bool add_cone_vertices = false) {
    reset();
    trace_all_boundaries(add_mesh_vertices, add_boundary_edges, add_cone_vertices);

    process_i_traces();

    trace_faces();
    CGAL_postcondition(is_valid_polygon_mesh(voronoi->graph, verbosity() > 0));
    CGAL_postcondition(voronoi->check_halfedge_info(add_boundary_edges));
    CGAL_postcondition(voronoi->check_topology());
    CGAL_postcondition(voronoi->is_valid_polygon_mesh());

    vout << IO::level(0) << "profiling: boundary trace: count = " << b_trace_timer.intervals()
         << ", time = " << b_trace_timer.time() << "s, speed = " << b_trace_timer.time() / b_trace_timer.intervals()
         << "s/trace" << std::endl;
    vout << IO::level(0) << "profiling: internal trace: count = " << i_trace_timer.intervals()
         << ", time = " << i_trace_timer.time() << "s, speed = " << i_trace_timer.time() / i_trace_timer.intervals()
         << "s/trace" << std::endl;
  }

#pragma region DisconnectedComponentRemoval
  void trace_principal_components(const auto& seed_sites,
                                  const auto& seed_points,
                                  const auto& seed_mesh_faces,
                                  auto& component_map) {
    std::unordered_multimap<mesh_face_descriptor, index_t> seed_face_map;
    for(index_t s = 0; s < seed_mesh_faces.size(); ++s) {
      seed_face_map.emplace(seed_mesh_faces[s], s);
    }

    std::vector<vd_face_descriptor> seed_faces;
    auto& G = voronoi->graph;
    for(auto fd : faces(G)) {
      auto hd = halfedge(fd, G);
      auto mesh_fd = get(voronoi->halfedge_info_map, hd).mesh_fd;
      put(component_map, fd, -1);

      auto range = seed_face_map.equal_range(mesh_fd);

      for(auto it = range.first; it != range.second; ++it) {
        auto p = seed_points[it->second];
        // if (bounded_side_3(G, fd, m_sites[site].point, voronoi->vpm) != CGAL::ON_UNBOUNDED_SIDE) {
        //     put(component_map, fd, site);
        //     seed_faces.push_back(fd);
        //     break;
        // }

        if(voronoi->face_bounded_side(fd, p) != CGAL::ON_UNBOUNDED_SIDE) {
          put(component_map, fd, seed_sites[it->second]);
          seed_faces.push_back(fd);
          break;
        }
      }
    }

    voronoi->flood_fill_faces(component_map, -1, seed_faces);
  }

  Face_component_map
  trace_principal_components(const auto& seed_sites, const auto& seed_points, const auto& seed_mesh_faces) {
    Face_component_map component_map{get(Face_component_property{}, voronoi->graph)};
    trace_principal_components(seed_sites, seed_points, seed_mesh_faces, component_map);
    return component_map;
  }

  void remove_all_disconnected_components(const auto&... args) {
    auto cmap = get(Face_component_property{}, voronoi->graph);
    for(;;) {
      trace_principal_components(args..., cmap);
      // Check if there is any disconnected component
      bool has_disconnected = false;
      for(auto fd : faces(voronoi->graph)) {
        if(get(cmap, fd) < 0) {
          has_disconnected = true;
          break;
        }
      }
      if(!has_disconnected) {
        break;
      }
      remove_disconnected_components(cmap);
    }
  }

  void remove_disconnected_components(const auto& cmap, bool single_thread = false) {
    scan_disconnected_components(cmap, m_disconnected_components_cache);

    dummy_face = add_face(voronoi->graph);
    for(auto& comp : m_disconnected_components_cache) {
      add_traces(comp);
      if(single_thread) {
        process_i_traces_single_thread();
      } else {
        process_i_traces();
      }
      // process_i_traces();
    }

    for(auto& comp : m_disconnected_components_cache) {
      prune_bounding_vertices(comp);
    }

    retrace_faces();
    voronoi->prune_orphan_vertices();

    CGAL_postcondition(CGAL::is_valid_face_graph(voronoi->graph, verbosity() > 0));
    CGAL_postcondition(voronoi->check_halfedge_info(true));
  }

  void retrace_faces() {
    auto& G = voronoi->graph;
    for(auto hd : halfedges(G)) {
      if(face(hd, G) != vd_graph_traits::null_face()) {
        CGAL::set_face(hd, dummy_face, G);
      }
    }

    std::vector<vd_face_descriptor> faces_to_remove;
    for(auto fd : faces(G)) {
      if(fd != dummy_face) {
        faces_to_remove.push_back(fd);
      }
    }

    for(auto fd : faces_to_remove) {
      CGAL::remove_face(fd, G);
    }

    trace_faces();
  }
#pragma endregion

  static int verbosity() { return vout.output_level(); }

  static void set_verbosity(int level) { vout.output_level(level); }
#pragma endregion

protected:
#pragma region PrivateVariables
  Traits traits;
  std::vector<Site> m_sites;
  std::vector<Metric_data> m_metrics;

  const Surface_mesh& mesh;
  Mesh_vertex_point_pmap vpm;
  Mesh_face_index_pmap face_index_map;
  Mesh_edge_index_pmap edge_index_map;

  std::shared_ptr<Voronoi_diagram_data> voronoi;
  vd_face_descriptor dummy_face;

  std::deque<Internal_trace> i_traces;
  std::unordered_map<Internal_vertex_id, vd_vertex_descriptor, Internal_vertex_id_hash> vert_map;

  std::unordered_map<Boundary_vertex_id, vd_vertex_descriptor, Boundary_vertex_id_hash> b_vert_map;
  std::unordered_multimap<vd_vertex_descriptor, mesh_halfedge_descriptor> processed_b_verts;

  std::exception_ptr exception;

  /// Disconnected component removal
  std::vector<Disconnected_component> m_disconnected_components_cache;
  std::unordered_multimap<Bisector_segment_id, vd_vertex_descriptor, Bisector_segment_id_hash> bounding_vertices;

/// Multithreading
#ifndef CGAL_HAS_NO_THREADS
  static inline BS::thread_pool pool;
  std::mutex vd_mutex, i_trace_mutex;
  std::mutex exception_mutex;
#endif

  inline CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Segment_cone_intersections, isects_cache);
  inline CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(std::vector<Segment_cone_intersections>, isects_vec_cache);

  inline CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(IO::Verbosity_level_ostream, vout);

  Real_timer b_trace_timer;
  Real_timer i_trace_timer;

#pragma endregion

#pragma region PrivateMethods

#ifndef CGAL_HAS_NO_THREADS
  void detach_task(auto&& func) { pool.detach_task(std::forward<decltype(func)>(func)); }
  void wait_for_tasks() { pool.wait(); }
#else
  void detach_task(auto&& func) { func(); }
  void wait_for_tasks() {
    // No-op
  }
#endif

  Cone_index cone_index(const Cone_descriptor& k) const {
    auto [c, m] = site(k.site_idx);
    return {k.site_idx, m.face_index_map[k.face]};
  }

  template <class FaceGraph, class VPMap>
  static Plane_3 supporting_plane(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                                  const FaceGraph& g,
                                  const VPMap& vpm) {
    auto i0 = source(hd, g), i1 = target(hd, g), i2 = target(next(hd, g), g);
    auto v0 = get(vpm, i0), v1 = get(vpm, i1), v2 = get(vpm, i2);
    return {v0, v1, v2};
  }

  template <class FaceGraph, class VPMap>
  Parametric_line_3 edge_segment(typename boost::graph_traits<FaceGraph>::halfedge_descriptor hd,
                                 const FaceGraph& g,
                                 const VPMap& vpm) const {
    auto i0 = source(hd, g), i1 = target(hd, g);
    return construct_parametric_line(get(vpm, i0), get(vpm, i1));
  }

  Parametric_line_3 mesh_edge_segment(mesh_halfedge_descriptor hd) const { return edge_segment(hd, mesh, vpm); }

  Plane_3 metric_face_plane(const Metric_data& m, metric_face_descriptor fd) const {
    return supporting_plane(halfedge(fd, m.graph), m.graph, m.vpm);
  }

  Plane_3 mesh_face_plane(mesh_halfedge_descriptor hd) const { return supporting_plane(hd, mesh, vpm); }

  Plane_3 get_bisect_plane(const Cone_descriptor& k0, const Cone_descriptor& k1) const {
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

  TRAIT_FUNC(Ray_3, construct_ray, construct_ray_3_object)
  TRAIT_FUNC(Point_3, construct_point, construct_point_3_object)
  TRAIT_FUNC(Point_3, construct_point_on, construct_point_on_3_object)
  TRAIT_FUNC(Vector_3, construct_vector, construct_vector_3_object)
  TRAIT_FUNC(Parametric_line_3, construct_parametric_line, construct_parametric_line_3_object)
  TRAIT_FUNC(FT, scalar_product, compute_scalar_product_3_object)
  TRAIT_FUNC(Orientation, orientation, orientation_3_object)
  TRAIT_FUNC(auto, intersect, intersect_3_object)
  TRAIT_FUNC(Metric_traits_data, construct_metric_traits_data, construct_metric_data_object)
  TRAIT_FUNC(auto, metric_any_intersection, metric_any_intersection_object)
  TRAIT_FUNC(auto, metric_any_intersected_face, metric_any_intersected_face_object)

  void trace_faces() {
    for(auto hd : halfedges(voronoi->graph)) {
      if(face(hd, voronoi->graph) != dummy_face)
        continue;
      auto fd = add_face(voronoi->graph);
      set_halfedge(fd, hd, voronoi->graph);
      for(auto hd_inner : halfedges_around_face(hd, voronoi->graph)) {
        set_face(hd_inner, fd, voronoi->graph);
      }
    }
    remove_face(dummy_face, voronoi->graph);
  }

  struct Cone_intersection
  {
    metric_face_descriptor fd;
    metric_halfedge_descriptor hd;
    FT t;
  };

  std::optional<Cone_intersection> find_next_cone(const Metric_data& m,
                                                  metric_face_descriptor fd,
                                                  metric_halfedge_descriptor hd_prev,
                                                  const Vector_3& p,
                                                  const Vector_3& d,
                                                  FT tmin,
                                                  std::optional<FT> tmax = std::nullopt) const {
    vout << IO::level(3) << SOURCE_LOC << ": start: line=" << p << " + " << d << " * t, range=" << tmin << " - "
         << (tmax ? std::to_string(*tmax) : "inf") << std::endl;
    auto hd0 = halfedge(fd, m.graph);
    for(auto hd : halfedges_around_face(hd0, m.graph)) {
      if(hd == hd_prev)
        continue;
      auto [v0, v1] = m.cone_face_bases(hd);
      auto n = m.cone_face_orthogonal_vector(hd);

      vout << IO::level(4) << SOURCE_LOC << ": test segment-cone face intersection: v0=" << v0 << ", v1=" << v1
           << ", n=" << n << std::endl;

      FT nd = scalar_product(n, d);
      FT np = scalar_product(n, p);
      if(is_zero(nd)) {
        // TODO
        CGAL_assertion_msg(!is_zero(np), "TODO: handle the case that the line lies on the face");
        vout << IO::level(4) << SOURCE_LOC << ": segment is parallel to cone facet" << std::endl;
        continue;
      }

      FT ti = -np / nd;
      if(ti < tmin || (tmax && ti > *tmax)) {
        vout << IO::level(4) << SOURCE_LOC << ": intersection ti = " << ti << " is out of range" << std::endl;
        continue;
      }

      Vector_3 vi = p + ti * d;

      auto ori0 = orientation(n, v0, vi);
      auto ori1 = orientation(n, vi, v1);
      if(!(ori0 == POSITIVE && ori1 == POSITIVE)) {
        // TODO
        CGAL_assertion_msg(ori0 != ZERO && ori1 != ZERO,
                           "TODO: marginal case: segment leaves from a lateral edge of the cone");
        // p_i is outside the 2D cone
        vout << IO::level(4) << SOURCE_LOC << ": intersection is outside the cone facet" << std::endl;
        continue;
      }

      // The segment leaves the cone from edge hd
      hd_prev = opposite(hd, m.graph);
      fd = face(hd_prev, m.graph);

      vout << IO::level(3) << SOURCE_LOC << ": end: found: t=" << ti << std::endl;

      return Cone_intersection{fd, hd_prev, ti};
    }
    vout << IO::level(3) << SOURCE_LOC << ": end: not found" << std::endl;
    return std::nullopt;
  }

  auto find_next_cone(const Cone_descriptor& k,
                      metric_halfedge_descriptor hd_prev,
                      const Parametric_line_3& l,
                      FT tmin,
                      std::optional<FT> tmax = std::nullopt) const {
    auto [c, m] = site(k.site_idx);
    auto p = construct_vector(c.point, construct_point(l));
    auto d = construct_vector(l);
    return find_next_cone(m, k.face, hd_prev, p, d, tmin, tmax);
  }

  /*!
   * \brief Find all intervals within a segment or a ray divided by cones of a site
   *
   * \param site_idx
   * \param segment
   * \param tmin
   * \param tmax
   * \param res
   */
  void find_segment_cone_intersections(
      index_t site_idx, const Parametric_line_3& segment, FT tmin, FT tmax, Segment_cone_intersections& res) const {
    CGAL_precondition(tmax > tmin);

    auto [c, m] = site(site_idx);
    auto pmin = construct_point_on(segment, tmin);
    auto pmax = construct_point_on(segment, tmax);
    auto p = construct_vector(c.point, construct_point(segment));
    auto d = construct_vector(segment);

    vout << IO::level(2) << SOURCE_LOC << ": start: site=" << c.point << ", segment=" << segment << std::endl;

    auto isect0 = metric_any_intersected_face(m.data, construct_vector(c.point, pmin));
    auto isect1 = metric_any_intersected_face(m.data, construct_vector(c.point, pmax));

    CGAL_assertion(isect0 && isect1);
    auto fd = *isect0;

    // TODO: fix fd0 if the query ray intersects the edge of the metric polyhedron

    res.clear();
    res.add_intersection(tmin, metric_graph_traits::null_halfedge());
    res.add_face(fd);

    metric_halfedge_descriptor hd_prev;
    while(fd != *isect1) {
      auto isect_next = find_next_cone(m, fd, hd_prev, p, d, tmin, tmax);
      CGAL_assertion_msg(!!isect_next, "Can not find the next cone");
      auto [fd_next, hd, t] = *isect_next;
      res.add_intersection(t, hd);
      res.add_face(fd_next);
      // TODO: is `tmin = isect_next->t;` needed?
      fd = fd_next;
      hd_prev = hd;
    }
    res.add_intersection(tmax, metric_graph_traits::null_halfedge());
  }

  void find_all_segment_cone_intersections(const Parametric_line_3& segment,
                                           FT tmin,
                                           FT tmax,
                                           std::vector<Segment_cone_intersections>& res) const {
    res.resize(m_sites.size());
    for(index_t site_idx = 0; site_idx < m_sites.size(); ++site_idx) {
      if(!m_sites[site_idx].enabled)
        continue;
      find_segment_cone_intersections(site_idx, segment, tmin, tmax, res[site_idx]);
    }
  }

  void _process_i_trace(const Internal_trace& tr, bool immediate = true) {
    // Check if the bisector has already been traced (from another direction)
    {
#ifndef CGAL_HAS_NO_THREADS
      std::lock_guard lock(vd_mutex);
#endif
      auto v_hd = halfedge(tr.v_vd, voronoi->graph);
      auto range = processed_b_verts.equal_range(tr.v_vd);
      for(auto it = range.first; it != range.second; ++it) {
        if(it->second == v_hd)
          return;
      }
    }
    // Range of the bisector segment
    FT tmin = 0, tmax;

    // Clip the bisector ray with the face on mesh
    mesh_halfedge_descriptor edge_hd;
    Parametric_line_3 edge;
    {
      bool has_edge_isect = false;
      for(auto hd : halfedges_around_face(tr.face_hd, mesh)) {
        if(hd == tr.prev_hd)
          continue;
        edge = mesh_edge_segment(hd);
        auto res = intersect(tr.bisect_line, edge, true);
        if(res) {
          auto [t_edge, t_mesh] = *res; // TODO: leaving from vertex
          if(t_edge <= 0 || t_mesh < 0 || t_mesh >= 1)
            continue;

          has_edge_isect = true;
          edge_hd = hd;
          tmax = t_edge;
          break;
        }
      }
      CGAL_assertion(has_edge_isect);
    }

    // Clip the bisector ray with bounding vertices
    vd_vertex_descriptor bounding_vd = vd_graph_traits::null_vertex();
    {
      Bisector_segment_id bisect_id{cone_index(tr.k0), cone_index(tr.k1),
                                    static_cast<index_t>(get(face_index_map, face(tr.face_hd, mesh)))};
      auto [vd_begin, vd_end] = bounding_vertices.equal_range(bisect_id);
      for(auto it = vd_begin; it != vd_end; ++it) {
        auto vd = it->second;
        if(vd == tr.v_vd)
          continue;
        // if (vd == tr.v_vd || !tr.bisect_line.has_on(pt)) continue;
        auto t = tr.bisect_line.parameter(get(voronoi->vpm, vd));
        if(t < tmin || t > tmax)
          continue;
        tmax = t;
        bounding_vd = vd;
      }
    }

    // Check if the bisector leaves the cone k0 or k1
    std::optional<Cone_intersection> cone_isect = std::nullopt;
    int cone_idx_next = -1;
    {
      auto isect0 = find_next_cone(tr.k0, tr.metric_prev_hd, tr.bisect_line, tmin);
      auto isect1 = find_next_cone(tr.k1, tr.metric_prev_hd, tr.bisect_line, tmin);

      if(isect0 && isect0->t < tmax) {
        cone_isect = *isect0;
        cone_idx_next = 0;
        tmax = isect0->t;
      }
      if(isect1 && isect1->t < tmax) {
        cone_isect = *isect1;
        cone_idx_next = 1;
        tmax = isect1->t;
      }
    }

    // Find the nearest 3-site bisector intersection
    T dist_min = INF;
    T tb_min;
    Plane_3 bi_plane_min;
    Cone_descriptor k2_min;

    for(index_t k2_site = 0; k2_site < m_sites.size(); ++k2_site) {
      if(k2_site == tr.k0.site_idx || k2_site == tr.k1.site_idx) {
        continue;
      }

      if(!m_sites[k2_site].enabled) {
        continue;
      }

      find_segment_cone_intersections(k2_site, tr.bisect_line, tmin, tmax, isects_cache);
      auto& isects = isects_cache;
      auto [c, m] = site(k2_site);
      for(auto [fd, isect_overlap] : isects) {
        Cone_descriptor k2{k2_site, fd};
        if(k2 == tr.k_prev)
          continue;

        Plane_3 bi_plane = get_bisect_plane(tr.k0, k2);
        if(bi_plane.is_degenerate())
          continue;

        auto tb = intersect(tr.bisect_line, bi_plane);
        if(!tb)
          continue;

        // TODO: marginal case: tb is on the terminal points
        if(*tb < isect_overlap.t_min || *tb >= isect_overlap.t_max)
          continue;

        auto dist = *tb - tmin;
        CGAL_assertion(dist >= 0);

        if(dist < dist_min) {
          dist_min = dist;
          tb_min = *tb;
          bi_plane_min = bi_plane;
          k2_min = k2;
        }
      }
    }

    auto mesh_fd = face(tr.face_hd, mesh);
    auto face_id = get(face_index_map, mesh_fd);
    if(dist_min < INF) {
      // Found a 3-site bisector
      Internal_vertex_id vid(cone_index(tr.k0), cone_index(tr.k1), cone_index(k2_min), face_id);
      Point_3 pt_start;
      vd_vertex_descriptor v_vd;
      {
#ifndef CGAL_HAS_NO_THREADS
        std::lock_guard lock(vd_mutex);
#endif

        if(auto vd_it = vert_map.find(vid); vd_it != vert_map.cend()) {
          voronoi->connect(tr.v_vd, vd_it->second, dummy_face, dummy_face, Bisector_edge_info{},
                           Halfedge_info{tr.k0, mesh_fd}, Halfedge_info{tr.k1, mesh_fd});
          return;
        }

        pt_start = construct_point_on(tr.bisect_line, tb_min);

        v_vd = voronoi->add_vertex(pt_start, Three_site_bisector_info{tr.face_hd, tr.k0, tr.k1, k2_min},
                                   tr.face_plane.orthogonal_vector());
        voronoi->connect(tr.v_vd, v_vd, dummy_face, dummy_face, Bisector_edge_info{}, Halfedge_info{tr.k0, mesh_fd},
                         Halfedge_info{tr.k1, mesh_fd});
        vert_map[vid] = v_vd;
      }

      // Add branch bisector of site 0 and 2
      auto bisect_dir_02 = cross_product(bi_plane_min.orthogonal_vector(), tr.face_plane.orthogonal_vector());
      auto [c0, m0] = site(tr.k0.site_idx);
      auto [c1, m1] = site(tr.k1.site_idx);
      auto p0 = metric_face_plane(m0, tr.k0.face);
      auto p1 = metric_face_plane(m1, tr.k1.face);
      if(is_negative(scalar_product(bisect_dir_02,
                                    p1.orthogonal_vector() / abs(p1.d()) - p0.orthogonal_vector() / abs(p0.d()))))
      {
        bisect_dir_02 = -bisect_dir_02;
      }
      auto bisect_line_02 = construct_parametric_line(pt_start, bisect_dir_02);

      Internal_trace tr02{
          bisect_line_02,
          bi_plane_min,
          tr.face_plane,
          edge_hd,
          mesh_graph_traits::null_halfedge(),
          tr.k0,
          k2_min,
          tr.k1,
          metric_graph_traits::null_halfedge(),
          v_vd,
      };

      if(immediate)
        detach_task([=, this]() { process_i_trace(tr02); });
      else
        i_traces.push_back(tr02);

      // Add branch bisector of site 1 and 2
      auto bisect_plane_12 = get_bisect_plane(tr.k1, k2_min);
      CGAL_assertion(!bisect_plane_12.is_degenerate());

      auto bisect_dir_12 = cross_product(bisect_plane_12.orthogonal_vector(), tr.face_plane.orthogonal_vector());
      if(is_negative(scalar_product(bisect_dir_12,
                                    p0.orthogonal_vector() / abs(p0.d()) - p1.orthogonal_vector() / abs(p1.d()))))
      {
        bisect_dir_12 = -bisect_dir_12;
      }
      auto bisect_line_12 = construct_parametric_line(pt_start, bisect_dir_12);
      Internal_trace tr12{
          bisect_line_12,
          bisect_plane_12,
          tr.face_plane,
          edge_hd,
          mesh_graph_traits::null_halfedge(),
          k2_min,
          tr.k1,
          tr.k0,
          metric_graph_traits::null_halfedge(),
          v_vd,
      };

      if(immediate)
        detach_task([=, this]() { process_i_trace(tr12); });
      else
        i_traces.push_back(tr12);
    } else if(cone_isect) {
      // The bisector ray intersects the cone
      auto k0_next = cone_idx_next == 0 ? tr.k1 : tr.k0;
      auto k1_prev = cone_idx_next == 0 ? tr.k0 : tr.k1;
      auto [c1_prev, m1_prev] = site(k1_prev.site_idx);
      Cone_descriptor k1_next{k1_prev.site_idx, cone_isect->fd};
      if(cone_idx_next == 0) {
        std::swap(k0_next, k1_next);
      }

      Internal_vertex_id vid(cone_index(k0_next), cone_index(k1_next), cone_index(k1_prev), face_id);
      Point_3 p;
      vd_vertex_descriptor v_vd;
      {
#ifndef CGAL_HAS_NO_THREADS
        std::lock_guard lock(vd_mutex);
#endif

        if(auto vd_it = vert_map.find(vid); vd_it != vert_map.cend()) {
          voronoi->connect(tr.v_vd, vd_it->second, dummy_face, dummy_face, Bisector_edge_info{},
                           Halfedge_info{tr.k0, mesh_fd}, Halfedge_info{tr.k1, mesh_fd});
          return;
        }

        p = construct_point_on(tr.bisect_line, tmax);
        v_vd = voronoi->add_vertex(p, Two_site_bisector_info{tr.face_hd, k0_next, k1_next.site_idx, cone_isect->hd},
                                   tr.face_plane.orthogonal_vector());
        voronoi->connect(tr.v_vd, v_vd, dummy_face, dummy_face, Bisector_edge_info{}, Halfedge_info{tr.k0, mesh_fd},
                         Halfedge_info{tr.k1, mesh_fd});
        vert_map[vid] = v_vd;
      }

      Plane_3 bi_plane = get_bisect_plane(k0_next, k1_next);
      auto bisect_dir = cross_product(bi_plane.orthogonal_vector(), tr.face_plane.orthogonal_vector());
      auto n = m1_prev.cone_face_orthogonal_vector(cone_isect->hd);
      if(sign(scalar_product(n, tr.bisect_line.d())) != sign(scalar_product(n, bisect_dir))) {
        bisect_dir = -bisect_dir;
      }
      auto bisect_line = construct_parametric_line(p, bisect_dir);

      Internal_trace tr1{
          bisect_line, bi_plane, tr.face_plane, edge_hd,        mesh_graph_traits::null_halfedge(),
          k0_next,     k1_next,  k1_prev,       cone_isect->hd, v_vd,
      };

      if(immediate)
        detach_task([=, this]() { process_i_trace(tr1); });
      else
        i_traces.push_back(tr1);
    } else if(bounding_vd != vd_graph_traits::null_vertex()) {
#ifndef CGAL_HAS_NO_THREADS
      std::lock_guard lock(vd_mutex);
#endif
      voronoi->connect(tr.v_vd, bounding_vd, dummy_face, dummy_face, Bisector_edge_info{},
                       Halfedge_info{tr.k0, mesh_fd}, Halfedge_info{tr.k1, mesh_fd});
      CGAL_assertion(degree(bounding_vd, voronoi->graph) == 2);

    } else {
      // The bisector leaves the face on mesh
      Boundary_vertex_id bvid(cone_index(tr.k0), cone_index(tr.k1), get(edge_index_map, CGAL::edge(edge_hd, mesh)));
#ifndef CGAL_HAS_NO_THREADS
      std::lock_guard lock(vd_mutex);
#endif
      auto vd_it = b_vert_map.find(bvid);
      CGAL_assertion_msg(vd_it != b_vert_map.cend(), "Bisector leaves the face from an unknown vertex");
      voronoi->connect(tr.v_vd, vd_it->second, dummy_face, dummy_face, Bisector_edge_info{},
                       Halfedge_info{tr.k0, mesh_fd}, Halfedge_info{tr.k1, mesh_fd});
      processed_b_verts.emplace(vd_it->second, edge_hd);
    }
  }

  void process_i_trace(const Internal_trace& tr, bool immediate = true) {
    {
#ifndef CGAL_HAS_NO_THREADS
      std::lock_guard lock(exception_mutex);
#endif
      if(exception)
        return;
    }
    try {
      _process_i_trace(tr, immediate);
    } catch(...) {
#ifndef CGAL_HAS_NO_THREADS
      std::lock_guard lock(exception_mutex);
#endif
      exception = std::current_exception();
    }
  }

  void _scan_disconnected_components(const auto& cmap,
                                     std::vector<Disconnected_component>& components,
                                     const auto& faces) const {
    std::unordered_set<vd_face_descriptor> visited_faces;
    auto& G = voronoi->graph;
    for(auto fd : faces) {
      if(get(cmap, fd) != -1)
        continue;

      if(visited_faces.contains(fd))
        continue;

      // Found a disconnected component
      Disconnected_component comp;

      std::queue<vd_face_descriptor> queue;
      queue.push(fd);

      // Search for all faces in the disconnected component
      while(!queue.empty()) {
        auto fd = queue.front();
        queue.pop();

        if(visited_faces.contains(fd))
          continue;
        visited_faces.insert(fd);

        comp.faces.insert(fd);

        for(auto hd : CGAL::halfedges_around_face(halfedge(fd, G), G)) {
          auto nfd = face(opposite(hd, G), G);
          auto site = get(cmap, nfd);
          if(site == -1) {
            // nfd is part of the disconnected component
            queue.push(nfd);
          } else {
            // nfd is part of a principal component around current disconnected component
            comp.neighbor_sites.insert(site);
          }

          auto ed = edge(hd, G);
          if(voronoi->is_boundary_edge(ed)) {
            comp.boundary_edges.insert(ed);
          } else if(voronoi->is_bisector_edge(ed)) {
            comp.bisector_edges.insert(ed);
          }
        }
      }

      // Check if the disconnected component group contains a site from a neighboring principal component
      auto cmap_copy = cmap;
      bool split = false;
      for(auto fd : comp.faces) {
        auto site = voronoi->cell_site(fd);
        if(comp.neighbor_sites.contains(site)) {
          // If so, mark it as a principal component so that it will not be removed
          split = true;
          put(cmap_copy, fd, site);
        }
      }

      // The rest of the disconnected component group will be removed
      if(split) {
        _scan_disconnected_components(cmap_copy, components, comp.faces);
        continue;
      }

      // Otherwise, add the disconnected component group to the list
      // Search for all inward trace edges
      for(auto fd : comp.faces) {
        for(auto vd : CGAL::vertices_around_face(halfedge(fd, G), G)) {
          for(auto hd : CGAL::halfedges_around_target(vd, G)) {
            if(voronoi->is_boundary_edge(hd) || comp.faces.contains(face(hd, G)) ||
               comp.faces.contains(face(opposite(hd, G), G)))
            {
              continue;
            }
            auto vt = target(hd, G);
            CGAL_assertion(std::holds_alternative<Three_site_bisector_info>(get(voronoi->vertex_info_map, vt)));
            comp.trace_inward_edges.insert(hd);
          }
        }
      }

      components.push_back(std::move(comp));
    }
  }

  void scan_disconnected_components(const auto& cmap, std::vector<Disconnected_component>& components) const {
    components.clear();
    _scan_disconnected_components(cmap, components, faces(voronoi->graph));
  }

  void scan_disconnected_components(const auto& cmap) const {
    m_disconnected_components_cache.clear();
    _scan_disconnected_components(cmap, m_disconnected_components_cache, faces(voronoi->graph));
  }

  void add_traces(const Disconnected_component& comp) {
    b_trace_timer.start();
    disable_all_sites();
    enable_sites(comp.neighbor_sites.cbegin(), comp.neighbor_sites.cend());
    bounding_vertices.clear();
    processed_b_verts.clear();

    auto& G = voronoi->graph;

    for(auto ed : comp.bisector_edges) {
      voronoi->remove_edge(ed);
    }

    for(auto ed : comp.boundary_edges) {
      voronoi->detach_edge(ed);
    }

    for(auto ed : comp.boundary_edges) {
      auto hd = halfedge(ed, voronoi->graph);
      auto vs = source(hd, voronoi->graph), vt = target(hd, voronoi->graph);
      auto mesh_hd = std::get<Boundary_edge_info>(get(voronoi->edge_info_map, ed)).hd;
      auto b_line = mesh_edge_segment(mesh_hd);
      auto ps = get(voronoi->vpm, vs), pt = get(voronoi->vpm, vt);
      FT ts = b_line.parameter(ps), tt = b_line.parameter(pt);
      if(ts > tt) {
        std::swap(vs, vt);
        std::swap(ps, pt);
        std::swap(ts, tt);
      }

      Cone_descriptor k0;
      find_nearest_site(ps, k0);
      vout << CGAL::IO::level(2) << SOURCE_LOC << ": tracing boundary " << mesh_hd << " with cone " << cone_index(k0)
           << " between point " << ps << " and " << pt << std::endl;
      auto [k1, vd] = trace_boundary(mesh_hd, k0, vs, true, true, true, false, dummy_face, dummy_face, ts, tt);

      voronoi->connect(vd, vt, dummy_face, dummy_face, Boundary_edge_info{mesh_hd}, Halfedge_info{k1, face(mesh_hd, G)},
                       Halfedge_info{k1, face(opposite(mesh_hd, G), G)});
    }
    for(auto hd : comp.trace_inward_edges) {
      add_itrace_from_voronoi_edge(hd);
    }
    b_trace_timer.stop();
  }

  void add_itrace_from_voronoi_edge(vd_halfedge_descriptor hd) {
    auto& vg = voronoi->graph;

    auto [k0, mesh_fd] = get(voronoi->halfedge_info_map, hd);
    auto [k1, _] = get(voronoi->halfedge_info_map, opposite(hd, vg));
    auto mesh_hd = halfedge(mesh_fd, mesh);
    auto vs = source(hd, vg), vt = target(hd, vg);
    auto ps = get(voronoi->vpm, vs), pt = get(voronoi->vpm, vt);

    Bisector_segment_id bisect_id{cone_index(k0), cone_index(k1), get(face_index_map, mesh_fd)};
    i_traces.push_back({
        construct_parametric_line(pt, construct_vector(ps, pt)),
        {},
        mesh_face_plane(mesh_hd),
        mesh_hd,
        mesh_graph_traits::null_halfedge(), // vt was a 3-site bisector inside the cone k0 and k1
        k0,
        k1,
        {},
        metric_graph_traits::null_halfedge(),
        vt,
    });
    bounding_vertices.emplace(bisect_id, vt);
  }

  auto join_edges(vd_vertex_descriptor vd) {
    auto& graph = voronoi->graph;
    auto hd0 = halfedge(vd, graph);
    if(hd0 == vd_graph_traits::null_halfedge()) {
      return hd0;
    }

    CGAL_precondition(degree(vd, voronoi->graph) == 2);

    auto hd1 = next(hd0, graph);
    auto hd0_opposite = opposite(hd0, graph);
    CGAL_precondition(CGAL::next(CGAL::opposite(hd1, graph), graph) == CGAL::opposite(hd0, graph));
    auto ed0 = edge(hd0, graph), ed1 = edge(hd1, graph);
    voronoi->detach_edge(ed0);
    voronoi->detach_edge(ed1);
    auto vs = source(hd0, graph), vt = target(hd1, graph);
    auto hd_new = voronoi->connect(vs, vt, CGAL::face(hd0, graph), CGAL::face(hd0_opposite, graph),
                                   get(voronoi->edge_info_map, ed0), get(voronoi->halfedge_info_map, hd0),
                                   get(voronoi->halfedge_info_map, hd0_opposite));
    CGAL::remove_edge(ed0, graph);
    CGAL::remove_edge(ed1, graph);
    return hd_new;
  }

  void prune_bounding_vertices(const Disconnected_component& comp) {
    for(auto ed : comp.boundary_edges) {
      auto hd = halfedge(ed, voronoi->graph);
      auto vs = source(hd, voronoi->graph), vt = target(hd, voronoi->graph);
      if(std::holds_alternative<Boundary_bisector_info>(get(voronoi->vertex_info_map, vs))) {
        join_edges(vs);
      }
      if(std::holds_alternative<Boundary_bisector_info>(get(voronoi->vertex_info_map, vt))) {
        join_edges(vt);
      }
      CGAL::remove_edge(ed, voronoi->graph);
    }

    for(auto hd : comp.trace_inward_edges) {
      auto vs = source(hd, voronoi->graph), vt = target(hd, voronoi->graph);
      join_edges(vt);
    }
  }

#pragma endregion
};
} // namespace Surface_Voronoi_diagram_with_star_metrics
} // namespace CGAL

#endif // SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_SURFACE_VORONOI_DIAGRAM_WITH_STAR_METRICS_H
