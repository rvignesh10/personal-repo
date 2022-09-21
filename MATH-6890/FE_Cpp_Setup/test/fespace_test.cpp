#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"

int main(){
    Mesh<4,0,3> mesh;
    mesh.Make1DCartesian(0.,2.);
    H1_FiniteElementSpace<3> fespace(&mesh);

    return 0;
}