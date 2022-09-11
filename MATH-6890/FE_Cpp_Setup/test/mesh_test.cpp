#include "../Headers/mesh.hpp"
#include "../Headers/geometries.hpp"

int main(){
    Mesh<2,0,1> mesh1;
    Mesh<2,2,2> mesh2;
    Mesh<2,2,3> mesh3;
    mesh1.Make1DCartesian(0.,2.);
    mesh2.Make2DCartesian(0., 1., 0., 1.,triangle);
    mesh3.Make2DCartesian(0.,1.,3.,4.,quadrilateral);
}