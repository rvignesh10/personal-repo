#include "../Headers/MatrixAlgebra.hpp"
#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/BiLinearForm.hpp"
#include "../Headers/LinearForm.hpp"
#include "../Headers/ADSolver.hpp"

void calc_a(double x, double y, Vector<double> &a);

int main(){
    
    const int degree = 1;
    Mesh<degree> mesh;
    mesh.Make1DCartesian(50,0,0.,1.);

    H1_FiniteElementSpace<degree> fespace(&mesh);


    ADSolver<degree> solver;
    solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(&calc_a, &fespace));
    solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(.1, &fespace));
    solver.AddDomainIntegrator(new MassIntegrator<degree>(&fespace));
    solver.AddLinearForm(new LinearForm<degree>(1., &fespace));
    solver.writeSparse();

    return 0;
}

void calc_a(double x, double y, Vector<double> &a){
    int l = a.getLength_();
    if (l == 2){
        a.setValue(0, 1+x);
        a.setValue(1, 0.01*y);
    }
    else{
        a.setValue(0, 0.);
        a.setValue(1, 0.);
        std::cerr << "incorrect dimension for a \n";
    }
}

