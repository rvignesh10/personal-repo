#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"

int main(){
    Mesh<3,0,1> mesh;
    mesh.Make1DCartesian(0.,2.);
    //mesh.Make2DCartesian(-3.,3.,-3.,3.,quadrilateral);
    // mesh.Make2DCartesian(-1.,1.,-1.,1.,triangle);
    H1_FiniteElementSpace<1> fespace(&mesh);

    return 0;
}