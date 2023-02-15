#include  <iostream>
#include  <cmath>
#include  <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

class CalcNode{
friend class CalcMesh;
friend class Tail;
friend class Ears;
protected:
    double x;
    double y;
    double z;

    double vx;
    double vy;
    double vz;

    double smth;
public:
    CalcNode():CalcNode(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0){}
    CalcNode(double x, double y, double z, double vx, double vy, double vz, double smth): x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), smth(smth){}
};

class Element{
friend class CalcMesh;
protected:
    unsigned long nodesIds[4];
};
class BodyPart{
protected:
    std::vector<size_t> nodes;
    std::vector<CalcNode>& global_nodes;
    double phi;
    double omega;
    double root_x;
    double root_y;
    double root_z;
public:
    BodyPart(std::vector<CalcNode>& global_nodes, double root_x, double root_y, double root_z, double omega):
        global_nodes(global_nodes), phi(0), root_x(root_x), root_y(root_y), root_z(root_z), omega(omega){}
    void addNode(size_t node_id){
        nodes.push_back(node_id);
    }
    virtual void move(double tau){
    }
};
class Tail: public BodyPart{
public:
    Tail(std::vector<CalcNode>& global_nodes, double root_x, double root_y, double root_z, double omega):
        BodyPart(global_nodes, root_x, root_y, root_z, omega){}
    void move(double tau){
        for(unsigned int i = 0; i < nodes.size(); i++) {
            double len_x = (global_nodes[nodes[i]].x - root_x); 
            double len_y = (global_nodes[nodes[i]].y - root_y); 
            double len = sqrt((len_x * len_x + len_y * len_y )) / 10;
            global_nodes[nodes[i]].x += -len_y * len * omega * tau;
            global_nodes[nodes[i]].y += len_x * len *  omega * tau;

            global_nodes[nodes[i]].vx = -len_y * len * 10 * omega * tau;
            global_nodes[nodes[i]].vy = len_x * len * 10 *  omega * tau;
            global_nodes[nodes[i]].smth = 2 * cos(len / 10 - 5 * phi);
        }
        phi += omega * tau;
        if (abs(phi) > M_PI / 4){
            omega *= -1;
        }
    }
};
class Ears: public BodyPart{
public:
    Ears(std::vector<CalcNode>& global_nodes, double root_x, double root_y, double root_z, double omega):
        BodyPart(global_nodes, root_x, root_y, root_z, omega){}
    void move(double tau){
        for(unsigned int i = 0; i < nodes.size(); i++) {
            double len_y = (global_nodes[nodes[i]].y - root_y); 
            double len_z = (global_nodes[nodes[i]].z - root_z); 
            double len_x = (global_nodes[nodes[i]].x - root_x); 

            global_nodes[nodes[i]].z += len_y * len_z /10 * omega * tau;
            global_nodes[nodes[i]].y += (-len_z - copysign(1.0, len_x) * len_x)* len_z / 10 *  omega * tau ;
            global_nodes[nodes[i]].x += copysign(1.0, len_x) * len_y * len_z / 10 *  omega * tau ;

            global_nodes[nodes[i]].vz = len_y * len_z * omega * tau * 10;
            global_nodes[nodes[i]].vy = (-len_z - copysign(1.0, len_x) * len_x)* len_z *  omega * tau * 10;
            global_nodes[nodes[i]].vx = copysign(1.0, len_x) * len_y * len_z  *  omega * tau * 10;
        }
        phi += omega * tau;
        if (phi > M_PI / 8 || phi < 0){
            omega *= -1;
        }
    }
};

class CalcMesh{
protected:
    std::vector<CalcNode> nodes;
    std::vector<Element> elements;
    Tail tail;
    Ears ears;
    double timer;
    double nose[3];
    double left_foot[3];
    double rigth_foot[3];
public:
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints):tail(nodes, -47, -192, 0.0, 0.5), ears(nodes, -48, -222, 46.4, 0.2) {
        timer = 0;

        nose[0] = -47.65;
        nose[1] = -229.22;
        nose[2] = 30.92;

        left_foot[0]= -43.37;
        left_foot[1]= -223.98;
        left_foot[2]= 0;

        rigth_foot[0]= -51.93;
        rigth_foot[1]= -223.98;
        rigth_foot[2]= 0;

        nodes.resize(nodesCoords.size() / 3);
        
        for (unsigned int i = 0; i < nodesCoords.size() / 3;  ++i){
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            nodes[i] = CalcNode(pointX, pointY, pointZ, 0.0, 0.0, 0.0, 0.0);
            if (pointY > -192)
                tail.addNode(i);
            if (pointZ > 47)
                ears.addNode(i);
        }

        elements.resize(tetrsPoints.size() / 4);
        for (unsigned int i = 0; i < tetrsPoints.size() / 4; ++i){
            for (int j = 0; j < 4; ++j){
                elements[i].nodesIds[j] = tetrsPoints[i * 4 + j] - 1;
            }
        }
    }
    void doTimeStep(double tau){
        for(unsigned int i = 0; i < nodes.size(); i++) {
         nodes[i].smth = 1 + cos((abs(nodes[i].x - nose[0]) + abs(nodes[i].y - nose[1]) + abs(nodes[i].z - nose[2])) / 10 - 2 * timer) 
                                + 0.5 * exp(-abs(nodes[i].y - left_foot[1]) / 10) * 
                                        cos((abs(nodes[i].x - left_foot[0]) + abs(nodes[i].y - left_foot[1]) + abs(nodes[i].z - left_foot[2])) / 10 - 3 * timer) 
                                + 0.5 * exp(-abs(nodes[i].y - rigth_foot[1]) / 10)
                                      * cos((abs(nodes[i].x - rigth_foot[0]) + abs(nodes[i].y - rigth_foot[1]) + abs(nodes[i].z - rigth_foot[2])) / 10 - 3 * timer);
        }
        tail.move(tau);
        ears.move(tau);
        timer += tau;

    }
    void snapshot(unsigned int snap_number){
        auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        auto dumpPoints = vtkSmartPointer<vtkPoints>::New();

        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);
        for(unsigned int i = 0; i < nodes.size(); i++) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
            double _vel[3] = {nodes[i].vx * 10, nodes[i].vy * 10, nodes[i].vz * 10};
            vel->InsertNextTuple(_vel);
            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        std::string fileName = "eevee-step-" + std::to_string(snap_number) + ".vtu"; 
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};
int main()
{
    double tau = 0.05;

    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("eve");
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
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodesTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodesTags, nodesCoord, parametricCoord);


    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    for (unsigned i = 0; i < elementTypes.size(); ++i){
        if (elementTypes[i] == GMSH_TETR_CODE)
            tetrsNodesTags = &elementNodeTags[i];
    }
    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }
    CalcMesh mesh(nodesCoord, *tetrsNodesTags);
    gmsh::finalize();
   
    for (unsigned step = 0; step < 6 * 20; step++){
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    return 0;
}
