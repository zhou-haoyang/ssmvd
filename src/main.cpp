
// #include <CGAL/Kernel/global_functions_3.h>
// #include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh.h>
// #include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
// #include <CGAL/Surface_mesh_parameterization/parameterize.h>
// #include <CGAL/Polygon_mesh_processing/measure.h>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Random.h>
// #include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Barycentric_coordinates_2/triangle_coordinates_2.h>

#include <CGAL/Simple_cartesian.h>
// #include <CGAL/AABB_tree.h>
// #include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
// #include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/convex_hull_3.h>

#include <ssmrvd.hpp>

#include <Eigen/Dense>

#include <igl/stb/read_image.h>
#include <igl/stb/write_image.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/writePLY.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>

// #include <boost/lexical_cast.hpp>

#include <cstdlib>
#include <iostream>
#include <unordered_map>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Triangle_2 Triangle_2;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef CGAL::Polyhedron_3<Kernel> Metric_polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Voronoi_diagram;
namespace RVD = CGAL::SSM_restricted_voronoi_diagram;
typedef RVD::SSM_restricted_voronoi_diagram<Kernel, Metric_polyhedron, Triangle_mesh, Voronoi_diagram>
    Restricted_voronoi_diagram;
typedef Restricted_voronoi_diagram::Cone_descriptor Cone_descriptor;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::vertex_descriptor vertex_descriptor;
typedef Graph_traits::face_iterator face_iterator;
typedef Graph_traits::face_descriptor face_descriptor;
typedef Graph_traits::halfedge_descriptor halfedge_descriptor;
using Index = std::ptrdiff_t;

using Eigen::Matrix2X;
using Eigen::Matrix2Xi;
using Eigen::Vector2i;
using Vector2F = Eigen::Vector2<FT>;
using Eigen::Vector2;

template <typename T>
auto normalize(T const &V) {
    auto const slen = V.squared_length();
    auto const d = CGAL::approximate_sqrt(slen);
    return V / d;
}

class App {
   public:
    App() {
        viewer.plugins.push_back(&imgui);
        imgui.widgets.push_back(&menu);
    }

    virtual ~App() {}

    void launch() {
        viewer.callback_init = [&](auto &viewer) { return init(); };

        viewer.callback_pre_draw = [&](auto &viewer) { return draw(); };

        viewer.callback_key_pressed = [&](auto &viewer, unsigned int key, int mod) { return key_pressed(key, mod); };

        menu.callback_draw_viewer_menu = [&]() { return draw_viewer_menu(); };

        viewer.launch();
    }

   protected:
    igl::opengl::glfw::Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    igl::opengl::glfw::imgui::ImGuiPlugin imgui;

    virtual bool init() { return false; }
    virtual bool draw() { return false; }
    virtual bool draw_viewer_menu() { return false; }
    virtual bool key_pressed(unsigned int key, int mod) { return false; }
};

class GVDApp : public App {
   public:
    GVDApp() : App(), gvd(sm) {}

   protected:
    // TextureGenerator gvd;
    Triangle_mesh sm;
    Metric_polyhedron metric;
    Restricted_voronoi_diagram gvd;

    std::string filename;
    Index n_sites = 10;
    // Index n_metrics = 0;
    // Index dir_path_idx = -2;
    // Eigen::VectorXd metrics_r, metrics_theta;

    // Index n_metrics_3d = 4;
    // Eigen::MatrixX3d metrics_3d;

    FT isoline_max = 1;
    Index isoline_steps = 10;

    Index tex_w = 256, tex_h = 256;
    Eigen::MatrixX<Index> tex_site;
    Eigen::MatrixX<FT> tex_dist;
    Eigen::MatrixX<uint8_t> tex_r, tex_g, tex_b, tex_a;
    Eigen::MatrixX3<uint8_t> site_colors;

    int mesh_view, vd_view, site_view;

    bool load_mesh(const std::string &filename) {
        Eigen::MatrixXd V, TC, N;   // #V by 3, #TC by 2, #N by 3
        Eigen::MatrixXi F, FT, FN;  // #F by 3
        if (!igl::readOBJ(filename, V, TC, N, F, FT, FN)) {
            std::cerr << "Invalid input file." << std::endl;
            return false;
        }

        sm = Triangle_mesh();
        auto [uv_pmap, uv_new] = sm.add_property_map<halfedge_descriptor, Point_2>("h:uv");
        auto [n_pmap, n_new] = sm.add_property_map<halfedge_descriptor, Vector_3>("h:n");

        std::vector<vertex_descriptor> v_map(V.rows());
        for (Index i = 0; i < V.rows(); ++i) {
            v_map[i] = sm.add_vertex(Point_3(V(i, 0), V(i, 1), V(i, 2)));
        }

        for (Index i = 0; i < F.rows(); ++i) {
            auto v0 = v_map[F(i, 0)], v1 = v_map[F(i, 1)], v2 = v_map[F(i, 2)];
            auto fi = sm.add_face(v0, v1, v2);
            auto h0 = sm.halfedge(fi), h1 = sm.next(h0), h2 = sm.next(h1);

            if (TC.rows() > 0 && FT.rows() > 0) {
                Eigen::RowVector2d uv0 = TC.row(FT(i, 0)), uv1 = TC.row(FT(i, 1)), uv2 = TC.row(FT(i, 2));
                std::unordered_map<vertex_descriptor, Point_2> uv_map{
                    {v0, Point_2(uv0.x(), uv0.y())}, {v1, Point_2(uv1.x(), uv1.y())}, {v2, Point_2(uv2.x(), uv2.y())}};
                uv_pmap[h0] = uv_map[sm.source(h0)];
                uv_pmap[h1] = uv_map[sm.target(h0)];
                uv_pmap[h2] = uv_map[sm.target(h1)];
            }

            if (N.rows() > 0 && FN.rows() > 0) {
                Eigen::RowVector3d n0 = N.row(FN(i, 0)), n1 = N.row(FN(i, 1)), n2 = N.row(FN(i, 2));
                n_pmap[h2] = Vector_3(n0.x(), n0.y(), n0.z());
                n_pmap[h0] = Vector_3(n1.x(), n1.y(), n1.z());
                n_pmap[h1] = Vector_3(n2.x(), n2.y(), n2.z());
            }
        }

        // duplicate vertices, uv, and norms for each face
        Eigen::MatrixXd VF, UVF, NF;
        Eigen::MatrixXi FF;

        VF.resize(F.rows() * 3, 3);
        FF.resize(F.rows(), 3);
        if (TC.rows() > 0) {
            UVF.resize(F.rows() * 3, 2);
        }
        if (N.rows() > 0) {
            NF.resize(F.rows() * 3, 3);
        }
        for (Index i = 0; i < F.rows(); ++i) {
            for (Index j = 0; j < 3; ++j) {
                VF.row(i * 3 + j) = V.row(F(i, j));
                if (TC.rows() > 0) UVF.row(i * 3 + j) = TC.row(FT(i, j));
                if (N.rows() > 0) NF.row(i * 3 + j) = N.row(FN(i, j));
            }
            FF.row(i) << i * 3, i * 3 + 1, i * 3 + 2;
        }

        auto &view = viewer.data(mesh_view);
        view.clear();
        view.set_mesh(VF, FF);
        view.set_colors(Eigen::RowVector3d(1, 1, 1));
        view.set_uv(UVF);
        view.set_normals(NF);
        igl::stb::read_image("texture.png", tex_r, tex_g, tex_b, tex_a);
        view.set_texture(tex_r, tex_g, tex_b, tex_a);
        view.show_texture = true;

        gvd.reload();
        return true;
    }

    auto triangle_vert_idx(face_descriptor f) const {
        auto h = sm.halfedge(f);
        return std::make_tuple(sm.source(h), sm.target(h), sm.target(sm.next(h)));
    }

    void random_sites() {
        auto [fbegin, fend] = sm.faces();
        std::vector<face_descriptor> faces(fbegin, fend);

        Eigen::VectorXi face_idx = Eigen::ArrayXi::Random(n_sites);
        Eigen::MatrixX2d bary = (Eigen::ArrayX2d::Random(n_sites, 2) + 1) / 2;
        Eigen::MatrixX3d site_pos(n_sites, 3);

        gvd.clear_sites();
        for (Index i = 0; i < n_sites; ++i) {
            auto f = faces[face_idx(i) % faces.size()];
            if (bary(i, 0) + bary(i, 1) >= 1) {
                bary(i, 0) = 1 - bary(i, 0);
                bary(i, 1) = 1 - bary(i, 1);
            }
            auto [i0, i1, i2] = triangle_vert_idx(f);
            auto pos = CGAL::barycenter(sm.point(i0), bary(i, 0), sm.point(i1), bary(i, 1), sm.point(i2));
            gvd.add_site(pos, 0);
            site_pos.row(i) << pos.x(), pos.y(), pos.z();
        }
        site_colors.setRandom(n_sites, 3);

        auto &view = viewer.data(site_view);
        view.set_points(site_pos, Eigen::RowVector3d(0.8, 0, 0));
        view.point_size = 10;
        std::vector<std::string> labels;
        for (Index i = 0; i < n_sites; ++i) {
            labels.push_back(std::to_string(i));
        }
        auto VL = site_pos.rowwise() + Eigen::RowVector3d(0, 0, 0.01);
        view.set_labels(VL, labels);
        view.show_custom_labels = true;
        view.label_size = 4;
        // view.label_color = Eigen::RowVector4d(0.8, 0, 0, 1);
    }

    void render_texture() {
        tex_a = Eigen::MatrixX<uint8_t>::Zero(tex_w, tex_h);
        tex_r = Eigen::MatrixX<uint8_t>::Zero(tex_w, tex_h);
        tex_g = Eigen::MatrixX<uint8_t>::Zero(tex_w, tex_h);
        tex_b = Eigen::MatrixX<uint8_t>::Zero(tex_w, tex_h);

        for (Index i = 0; i < tex_a.rows(); ++i) {
            for (Index j = 0; j < tex_a.cols(); ++j) {
                Index site = tex_site(i, j);
                FT dist = tex_dist(i, j);
                if (site < 0) continue;
                tex_a(i, j) = (uint8_t)std::round(std::clamp(
                    (1. - (std::floor(dist / isoline_max * isoline_steps) / isoline_steps)) * 255., 0., 255.));
                tex_r(i, j) = site_colors(site, 0);
                tex_g(i, j) = site_colors(site, 1);
                tex_b(i, j) = site_colors(site, 2);
            }
        }

        auto &view = viewer.data(mesh_view);
        view.set_texture(tex_r, tex_g, tex_b, tex_a);
        view.show_texture = true;
    }

    Point_2 pixel_idx_to_uv(Index i, Index j) const { return Point_2((i + 0.5) / tex_w, (j + 0.5) / tex_h); }

    auto pixel_uv_to_idx(FT x, FT y) const {
        return std::make_pair(std::clamp((Index)(x * tex_w), 0L, tex_w - 1),
                              std::clamp((Index)(y * tex_h), 0L, tex_h - 1));
    }

    auto pixel_uv_to_idx(const Point_2 &uv) const { return pixel_uv_to_idx(uv.x(), uv.y()); }

    void gen_tex_rvd() {
        tex_site = Eigen::MatrixX<Index>::Constant(tex_w, tex_h, -1);
        tex_dist = Eigen::MatrixX<FT>::Zero(tex_w, tex_h);

        auto [uv_map, uv_exist] = sm.property_map<halfedge_descriptor, Point_2>("h:uv");
        assert(uv_exist);

        for (auto face : sm.faces()) {
            auto h0 = sm.halfedge(face), h1 = sm.next(h0), h2 = sm.next(h1);
            vertex_descriptor i0 = sm.source(h0), i1 = sm.target(h0), i2 = sm.target(h1);
            Point_3 p0 = sm.point(i0), p1 = sm.point(i1), p2 = sm.point(i2);
            Point_2 uv0 = uv_map[h0], uv1 = uv_map[h1], uv2 = uv_map[h2];

            auto bb = Triangle_2(uv0, uv1, uv2).bbox();
            auto [i_min, j_min] = pixel_uv_to_idx(bb.xmin(), bb.ymin());
            auto [i_max, j_max] = pixel_uv_to_idx(bb.xmax(), bb.ymax());

            for (Index i = i_min; i <= i_max; ++i) {
                for (Index j = j_min; j <= j_max; ++j) {
                    Point_2 uv = pixel_idx_to_uv(i, j);
                    std::array<FT, 3> bary;

                    CGAL::Barycentric_coordinates::triangle_coordinates_2(uv0, uv1, uv2, uv, bary.begin());
                    if (bary[0] < 0 || bary[1] < 0 || bary[2] < 0) {
                        continue;
                    }
                    Point_3 p = CGAL::barycenter(p0, bary[0], p1, bary[1], p2, bary[2]);

                    Cone_descriptor m_cone;
                    auto dist = gvd.find_nearest_site(p, m_cone);

                    tex_site(i, j) = m_cone.site_idx;
                    tex_dist(i, j) = dist;
                }
            }
        }

        render_texture();
    }

    bool init() override {
        viewer.core().lighting_factor = 0.0;

        mesh_view = viewer.append_mesh();
        vd_view = viewer.append_mesh();
        site_view = viewer.append_mesh();

        filename = "isohemisphere.obj";

        Eigen::MatrixX3d metrics_3d(4, 3);
        metrics_3d << -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1;
        metrics_3d *= std::sqrt(6.) / 4.;
        std::vector<Point_3> m_points;
        for (Index i = 0; i < metrics_3d.rows(); ++i) {
            m_points.emplace_back(metrics_3d(i, 0), metrics_3d(i, 1), metrics_3d(i, 2));
        }
        CGAL::convex_hull_3(m_points.begin(), m_points.end(), metric);
        gvd.add_metric(metric);

        if (load_mesh(filename)) {
            random_sites();
        }

        return false;
    }

    bool draw() override { return false; }

    bool draw_viewer_menu() override {
        menu.draw_viewer_menu();

        // if (ImGui::CollapsingHeader("Polygon Metrics", ImGuiTreeNodeFlags_DefaultOpen)) {
        //     if (ImGui::InputScalar("Number of Metrics", ImGuiDataType_S64, &n_metrics)) {
        //         reset_metrics_2d();
        //     }

        //     for (Index i = 0; i < metrics_theta.size(); ++i) {
        //         ImGui::InputScalar(("Theta " + boost::lexical_cast<std::string>(i)).c_str(), ImGuiDataType_Double,
        //                            &metrics_theta(i));
        //         ImGui::InputScalar(("R " + boost::lexical_cast<std::string>(i)).c_str(), ImGuiDataType_Double,
        //                            &metrics_r(i));
        //     }
        // }

        // if (ImGui::CollapsingHeader("Polyhedron Metrics", ImGuiTreeNodeFlags_DefaultOpen)) {
        //     if (ImGui::InputScalar("Number of Metric Points", ImGuiDataType_S64, &n_metrics_3d)) {
        //         metrics_3d.resize(n_metrics_3d, 3);
        //         metrics_3d.setRandom();
        //     }

        //     for (Index i = 0; i < metrics_3d.rows(); ++i) {
        //         ImGui::InputScalar(("X " + boost::lexical_cast<std::string>(i)).c_str(), ImGuiDataType_Double,
        //                            &metrics_3d(i, 0));
        //         ImGui::InputScalar(("Y " + boost::lexical_cast<std::string>(i)).c_str(), ImGuiDataType_Double,
        //                            &metrics_3d(i, 1));
        //         ImGui::InputScalar(("Z " + boost::lexical_cast<std::string>(i)).c_str(), ImGuiDataType_Double,
        //                            &metrics_3d(i, 2));
        //     }
        // }

        if (ImGui::CollapsingHeader("VD", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputText("filename", filename);
            if (ImGui::Button("Open Mesh")) {
                filename = igl::file_dialog_open();
                load_mesh(filename);
            }
            if (ImGui::Button("Load Mesh")) {
                load_mesh(filename);
            }

            ImGui::InputScalar("Number of Sites", ImGuiDataType_S64, &n_sites);

            if (ImGui::Button("Random Sites")) {
                random_sites();
            }

            ImGui::InputScalar("Texture Width", ImGuiDataType_S64, &tex_w);
            ImGui::InputScalar("Texture Height", ImGuiDataType_S64, &tex_h);

            if (ImGui::Button("Generate RVD")) {
                gen_tex_rvd();
            }

            auto iso_max_changed = ImGui::InputScalar("Isoline Max", ImGuiDataType_Double, &isoline_max);
            auto iso_steps_changed = ImGui::InputScalar("Isoline Steps", ImGuiDataType_S64, &isoline_steps);

            if (iso_max_changed || iso_steps_changed) {
                render_texture();
            }

            if (ImGui::Button("Save Texture")) {
                auto path = igl::file_dialog_save();
                igl::stb::write_image(path, tex_r, tex_g, tex_b, tex_a);
            }

            // if (ImGui::Button("Save Mesh")) {
            //     auto path = igl::file_dialog_save();
            //     igl::writePLY(path, viewer.data().V, viewer.data().F, viewer.data().V_normals, viewer.data().V_uv);
            // }
        }
        return false;
    }

    void update_vd() {
        std::unordered_map<vertex_descriptor, Index> v_map;
        auto &vd = gvd.vd.graph;
        Eigen::MatrixX3d V(vd.number_of_vertices(), 3);
        Index i = 0;
        std::vector<std::string> labels;
        for (auto v : vd.vertices()) {
            v_map[v] = i;
            auto p = vd.point(v);
            V.row(i) << p.x(), p.y(), p.z();
            auto info = get(gvd.vd.vertex_info_map, v);
            if (auto d = std::get_if<Restricted_voronoi_diagram::Boundary_vertex_info>(&info)) {
                labels.push_back(std::format("V{}: {}", v.idx(), d->k.site_idx));
            } else if (auto d = std::get_if<Restricted_voronoi_diagram::Boundary_bisector_info>(&info)) {
                labels.push_back(std::format("B{}: {} {}", v.idx(), d->k0.site_idx, d->k1.site_idx));
            } else if (auto d = std::get_if<Restricted_voronoi_diagram::Boundary_cone_info>(&info)) {
                labels.push_back(std::format("C{}: {}", v.idx(), d->k.site_idx));
            } else {
                labels.push_back("");
            }
            i++;
        }

        Eigen::MatrixX2i E(vd.number_of_edges(), 2);
        i = 0;
        for (auto e : vd.edges()) {
            auto he = vd.halfedge(e);
            E.row(i) << v_map[vd.source(he)], v_map[vd.target(he)];
            i++;
        }

        auto &view = viewer.data(vd_view);
        view.set_points(V, Eigen::RowVector3d(0, 0, 1));
        view.point_size = 10;
        view.set_edges(V, E, Eigen::RowVector3d(0, 1, 0));
        view.line_width = 4;
        auto VL = V.rowwise() + Eigen::RowVector3d(0, 0, 0.01);
        view.set_labels(VL, labels);
        view.show_custom_labels = true;
        view.label_size = 2;
    }

    bool key_pressed(unsigned int key, int mod) override {
        switch (key) {
            case 'b':
                gvd.trace_boundary(CGAL::Polygon_mesh_processing::longest_border(sm).first);
                update_vd();
                break;
            case 'g':
                gvd.step();
                update_vd();
                break;
            case 'r':
                gvd.reset();
                update_vd();
                break;
        }
        return false;
    }
};

int main(int argc, char **argv) {
    GVDApp app;
    app.launch();
    return EXIT_SUCCESS;
}