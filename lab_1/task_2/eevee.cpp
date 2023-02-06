#include <gmsh.h>
#include <set>
#include <cmath>
#include <vector>
#include <iostream>

int main(int argc, char **argv){
    gmsh::initialize();
    gmsh::model::add("eevee");
    try{
        gmsh::merge("eevee.stl");
    }catch (...){
        gmsh::finalize();
        return 0;
    }
    
    double angle = 20;
    bool forceParametrizablePatches = true;
    bool includeBoundary = true;

    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary,
                                        forceParametrizablePatches);
    gmsh::model::mesh::createGeometry();
    
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for (auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "2");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);
    gmsh::write("eevee.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

    gmsh::finalize();
    return 0;
}
