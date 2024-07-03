
#include <ssmrvd.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polyhedron_3.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
namespace RVD = CGAL::SSM_restricted_voronoi_diagram;

// typedef RVD::SSM_restricted_voronoi_diagram<Kernel, Polyhedron, SurfaceMesh> SSM_restricted_voronoi_diagram;

int main(int argc, char** argv) {
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/nefertiti.off");
    SurfaceMesh sm;
    if (!CGAL::IO::read_polygon_mesh(filename, sm)) {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    // a halfedge on the border
    auto bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
    // SMP::parameterize(sm, bhd, uv_map);

    RVD::SSM_restricted_voronoi_diagram rvd(sm);
}
