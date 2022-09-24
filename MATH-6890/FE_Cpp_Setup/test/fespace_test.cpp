#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"

int main(){
    Mesh<2> mesh;
    //mesh.Make1DCartesian(3,0,0.,2.);
    //mesh.Make2DCartesian(2,2,-3.,3.,-3.,3.,quadrilateral);
    mesh.Make2DCartesian(2,2,-1.,1.,-1.,1.,triangle);
    H1_FiniteElementSpace<2> fespace(&mesh);

    return 0;
}