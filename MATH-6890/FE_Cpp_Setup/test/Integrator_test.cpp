#include "../Headers/MatrixAlgebra.hpp"
#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/BiLinearForm.hpp"
#include "../Headers/LinearForm.hpp"
#include "../Headers/ADSolver.hpp"

void calc_a(double x, double y, Vector<double> &a);
void set_Bc(double x, double y, double &q);

int main(){
    
    const int degree = 1;
    Mesh<degree> mesh;
    mesh.Make1DCartesian(100,0,0.,1.);

    H1_FiniteElementSpace<degree> fespace(&mesh);
    fespace.setBdrCondition(&set_Bc);

    ADSolver<degree> solver;
    solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(&calc_a, &fespace));
    solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(.1, &fespace));
    solver.AddLinearForm(new LinearForm<degree>(1., &fespace));
    solver.writeSparse();

    return 0;
}

void calc_a(double x, double y, Vector<double> &a){
    int l = a.getLength_();
    if (l == 2){
        a.setValue(0, 1.+x);
        a.setValue(1, 0.0);
    }
    else{
        a.setValue(0, 1.+x);
        a.setValue(1, 0.);
        std::cerr << "incorrect dimension for a \n";
    }
}

void set_Bc(double x, double y, double &q){
    if(x < 1e-14){
        q = 0.;
    }
    else if ((x-1.)<1e-14){
        q = 1.;
    }
}
