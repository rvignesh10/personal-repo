#include "../Headers/mesh.hpp"
#include "../Headers/geometries.hpp"

int main(){
    Mesh<1> mesh1;
    Mesh<2> mesh2;
    Mesh<1> mesh3;
    mesh1.Make1DCartesian(4,0,0.,2.);
    mesh2.Make2DCartesian(2,2,0., 1., 0., 1.,triangle);
    mesh3.Make2DCartesian(3,3,0.,1.,3.,4.,quadrilateral);

    Element<2> *el;
    el = mesh2.GetElement(8);
    el->printElementNodes();
    return 0;
}