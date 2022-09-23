#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/Integrator.hpp"
#include "../Headers/ADSolver.hpp"

int main(){
    Mesh<2> mesh;
    mesh.Make1DCartesian(50,0,-2.,2.);
    //mesh.Make2DCartesian(2,2,-1.,1.,-1.,1.,triangle);
    H1_FiniteElementSpace<2> fespace(&mesh);
    ADSolver<2> solver;
    solver.AddDomainIntegrator(new DiffusionIntegrator<2>(.1, &fespace));
    //solver.AddDomainIntegrator(new MassIntegrator<2>(&fespace));
    solver.writeSparse();
    //solver.AddDomainIntegrator(new DiffusionIntegrator<2>(1., &fespace));
    return 0;
}