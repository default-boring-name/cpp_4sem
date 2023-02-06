#include <gmsh.h>
#include <vector>
#include <set>
double lc = 1e-1;
struct Circle
{
    int boundary_points[4] = {0};
    int centre = 0;
    int boundary_curves[4] = {0};
    void make_boundary(){
        for (int i = 0; i < 4; ++i){
            boundary_curves[i] = gmsh::model::geo::addCircleArc(boundary_points[i], centre, boundary_points[(i + 1) % 4]);            
        }
    }
};
struct Torus
{
    Circle boundary_circles[4];
    int boundary_curves[16];
    int boundary_surfaces[16];
    int top_centre;
    int down_centre;
    int centre;
    Torus(double x, double y, double z, double R, double r):centre(gmsh::model::geo::addPoint(x, y, z, lc)), top_centre(gmsh::model::geo::addPoint(x, y, z + r, lc)), down_centre(gmsh::model::geo::addPoint(x, y, z - r, lc)){
        boundary_circles[0].centre = gmsh::model::geo::addPoint(x + R, y, z, lc);
        boundary_circles[0].boundary_points[0] = gmsh::model::geo::addPoint(R + x, y, z + r, lc);
        boundary_circles[0].boundary_points[1] = gmsh::model::geo::addPoint(R + r + x, y, z, lc);
        boundary_circles[0].boundary_points[2] = gmsh::model::geo::addPoint(R + x, y, z - r, lc);
        boundary_circles[0].boundary_points[3] = gmsh::model::geo::addPoint(R - r + x, y, z, lc);
    
        boundary_circles[1].centre = gmsh::model::geo::addPoint(x, R, y, lc);
        boundary_circles[1].boundary_points[0] = gmsh::model::geo::addPoint(x, y + R, z + r, lc);
        boundary_circles[1].boundary_points[1] = gmsh::model::geo::addPoint(x, y + R + r, z, lc);
        boundary_circles[1].boundary_points[2] = gmsh::model::geo::addPoint(x, y + R, z - r, lc);
        boundary_circles[1].boundary_points[3] = gmsh::model::geo::addPoint(x, y + R - r, z, lc);
    
        boundary_circles[2].centre = gmsh::model::geo::addPoint(x - R, y, 0, lc);
        boundary_circles[2].boundary_points[0] = gmsh::model::geo::addPoint(x - R, y, z + r, lc);
        boundary_circles[2].boundary_points[1] = gmsh::model::geo::addPoint(x - R - r, y, z, lc);
        boundary_circles[2].boundary_points[2] = gmsh::model::geo::addPoint(x - R, y, z - r, lc);
        boundary_circles[2].boundary_points[3] = gmsh::model::geo::addPoint(x - R + r, y, z, lc);
    
        boundary_circles[3].centre = gmsh::model::geo::addPoint(x, -R, y, lc);
        boundary_circles[3].boundary_points[0] = gmsh::model::geo::addPoint(x, y - R, z + r, lc);
        boundary_circles[3].boundary_points[1] = gmsh::model::geo::addPoint(x, y - R - r, z, lc);
        boundary_circles[3].boundary_points[2] = gmsh::model::geo::addPoint(x, y - R, z - r, lc);
        boundary_circles[3].boundary_points[3] = gmsh::model::geo::addPoint(x, y - R + r, z, lc);

    }
    void make_wire(){
        for (int i = 0; i < 4; ++i){
            boundary_circles[i].make_boundary();
        }
        for (int i = 0; i < 4; ++i){
             boundary_curves[i * 4] = gmsh::model::geo::addCircleArc(boundary_circles[i].boundary_points[0], top_centre, boundary_circles[(i + 1) % 4].boundary_points[0]);
             boundary_curves[i * 4 + 2] = gmsh::model::geo::addCircleArc(boundary_circles[i].boundary_points[2], down_centre, boundary_circles[(i + 1) % 4].boundary_points[2]);
             boundary_curves[i * 4 + 1] = gmsh::model::geo::addCircleArc(boundary_circles[i].boundary_points[1], centre, boundary_circles[(i + 1) % 4].boundary_points[1]);
             boundary_curves[i * 4 + 3] = gmsh::model::geo::addCircleArc(boundary_circles[i].boundary_points[3], centre, boundary_circles[(i + 1) % 4].boundary_points[3]);
        }
    }
    void make_surface(){
        for (int i = 0; i < 4; ++i){
            for (int j = 0; j < 4; ++j){
                int curve_loop = gmsh::model::geo::addCurveLoop({boundary_circles[i].boundary_curves[j], boundary_curves[i * 4 + (j + 1) % 4], 
                                                                   -boundary_circles[(i + 1) % 4].boundary_curves[j], -boundary_curves[i * 4 + j]});
                boundary_surfaces[i * 4 + j] = gmsh::model::geo::addSurfaceFilling({curve_loop});
            }
        }
    }
};

struct ThickTorus{
    Torus inner;
    Torus outer;
    int volume = 0;
    ThickTorus(double x, double y, double z, double R, double r_in, double r_out):inner(x, y, z, R, r_in),outer(x, y, z, R, r_out){
    }
    void make_wire(){
        inner.make_wire();
        outer.make_wire();
    }
    void make_surface(){
        inner.make_surface();
        outer.make_surface();
    }
    void make_volume(){
        std::vector<int> surfaces;
        for (int i = 0; i < 16; ++i){
            surfaces.push_back(-inner.boundary_surfaces[i]);
            surfaces.push_back(outer.boundary_surfaces[i]);
        }
        int surface_loop = gmsh::model::geo::addSurfaceLoop(surfaces);
        volume = gmsh::model::geo::addVolume({surface_loop});
    }
};

 int main(int argc, char **argv){
    gmsh::initialize(); 
    gmsh::model::add("torus");

    double r_in = 0.6;
    double r_out = 0.4;
    double R = 1;

    ThickTorus torus(0, 0, 0, R, r_in, r_out);

    torus.make_wire();
    torus.make_surface();
    torus.make_volume();
    
    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("torus.msh");

    std::set<std::string> args(argv, argv + argc);
      if(!args.count("-nopopup")) gmsh::fltk::run();
    gmsh::finalize();
    return 0;
 }
