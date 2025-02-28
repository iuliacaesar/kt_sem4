#include <gmsh.h>
#include <set>
int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("torus");

    double R=15.0, r_out=5.0, r_in=3.0;
    int outTorus = gmsh::model::occ::addTorus(0, 0, 0, R, r_out);
    int inTorus = gmsh::model::occ::addTorus(0, 0, 0, R, r_in);
    std::vector<std::pair<int, int> > ov;
    std::vector<std::vector<std::pair<int, int> > > ovv;
    gmsh::model::occ::cut({{3, outTorus}}, {{3, inTorus}}, ov, ovv);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.5);
    

    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("torus.msh");
        
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();


  
    gmsh::finalize();
  
    return 0;
  }
