#include "../Headers/MatrixAlgebra.hpp"
#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/Integrator.hpp"
#include "../Headers/ADSolver.hpp"

void calc_a(double x, double y, Vector<double> &a);

int main(){
    Vector<double> a(2);
    Mesh<2> mesh;
    //mesh.Make1DCartesian(50,0,-2.,2.);
    mesh.Make2DCartesian(30,30,-3.,3.,-3.,3.,quadrilateral);
    H1_FiniteElementSpace<2> fespace(&mesh);
    ADSolver<2> solver;
    solver.AddDomainIntegrator(new AdvectionIntegrator<2>(&calc_a, &fespace));
    solver.AddDomainIntegrator(new DiffusionIntegrator<2>(.1, &fespace));
    solver.AddDomainIntegrator(new MassIntegrator<2>(&fespace));
    solver.writeSparse();
    return 0;
}

void calc_a(double x, double y, Vector<double> &a){
    int l = a.getLength_();
    if (l == 2){
        a.setValue(0, 0.2*x);
        a.setValue(1, 0.01*y);
    }
    else{
        a.setValue(0, 0.);
        a.setValue(1, 0.);
        std::cerr << "incorrect dimension for a \n";
    }
}

