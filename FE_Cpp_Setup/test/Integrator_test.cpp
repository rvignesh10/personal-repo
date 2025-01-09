#include "../Headers/MatrixAlgebra.hpp"
#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/element.hpp"
#include "../Headers/BiLinearForm.hpp"
#include "../Headers/LinearForm.hpp"
#include "../Headers/FlowSolver.hpp"

void calc_a(double x, double y, Vector<double> &a);
void calc_k(double x, double y, double &kappa);
void calc_tau1D(double h, double ax, double Pee, double &tau);
void set_Bc(double x, double y, double &q);

int main(){
    
    const int degree = 1;
    int num_elem = 4;
    Mesh<degree> mesh;
    mesh.Make1DCartesian(100,0,0.,1.);
    // mesh.Init(1, segment, num_elem, num_elem, 0);
    // for (int i=0; i<num_elem; i++){
    //     POINT p[2];
    //     if(i==0){
    //         p[0].setCoordinates(0.,0.);
    //         p[0].setAttribute(1);
    //         p[0].setIdx(1);
    //         p[1].setCoordinates(0.5,0.);
    //         p[1].setAttribute(0);
    //         p[1].setIdx(2);
    //         mesh.AddElement(i, p, 2, 1);
    //     }
    //     else if(i==1){
    //         p[0].setCoordinates(0.5,0.);
    //         p[0].setAttribute(0);
    //         p[0].setIdx(2);
    //         p[1].setCoordinates(0.75,0.);
    //         p[1].setAttribute(0);
    //         p[1].setIdx(3);
    //         mesh.AddElement(i, p, 2, 0);
    //     }
    //     else if(i==2){
    //         p[0].setCoordinates(0.75,0.);
    //         p[0].setAttribute(0);
    //         p[0].setIdx(3);
    //         p[1].setCoordinates(0.875,0.);
    //         p[1].setAttribute(0);
    //         p[1].setIdx(4);
    //         mesh.AddElement(i, p, 2, 0);
    //     }
    //     else if(i==3){
    //         p[0].setCoordinates(0.875,0.);
    //         p[0].setAttribute(0);
    //         p[0].setIdx(4);
    //         p[1].setCoordinates(1.,0.);
    //         p[1].setAttribute(1);
    //         p[1].setIdx(5);
    //         mesh.AddElement(i, p, 2, 1);
    //     }
    // }

    // Element<degree> *el;
    // el = mesh.GetElement(4);
    // el->printElementNodes();

    

    H1_FiniteElementSpace<degree> fespace(&mesh);
    fespace.setBdrCondition(&set_Bc);
    Vector<double> a(2);
    a.setValue(0,1e0);
    a.setValue(1,0.);
    double kappa  = 1e-5;
    FlowSolver<degree> solver(steady, AD);
    solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(a, &fespace));
    solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(kappa, &fespace));
    solver.AddLinearForm(new LinearForm<degree>(0., &fespace));
    solver.AddSUPGStabilization(new SUPGStabilization<degree>(AD,a,kappa,&calc_tau1D,0.,&fespace));
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

void calc_k(double x, double y, double &kappa){
    kappa = 0.1*(1.+x);
}

void calc_tau1D(double h, double ax, double Pee, double &tau){
    tau = (h/abs(ax))*(((1+exp(-2.0*Pee))/(1-exp(-2.0*Pee)))-(1/Pee));
}

void set_Bc(double x, double y, double &q){
    if(x < 1e-14){
        q = 0.;
    }
    else if ((x-1.)<1e-14){
        q = 1.;
    }
}
