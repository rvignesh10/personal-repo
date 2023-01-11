#include "../Headers/MatrixAlgebra.hpp"
#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/element.hpp"
#include "../Headers/BiLinearForm.hpp"
#include "../Headers/LinearForm.hpp"
#include "../Headers/FlowSolver.hpp"
#include <string.h>
#include <math.h>

void calc_tau1D_AD(double h, double ax, double Pee, double &tau);
void calc_tau1D_ADR(double h, double ax, double kappa, double c, double &tau);
void set_Bc_1(double x, double y, double &q);
void set_Bc_2(double x, double y, double &q);

int main(int argc, char *argv[]){
    
    const int degree = 1;
    int num_elem;
    Mesh<degree> mesh;

    if((strcmp(argv[1],"1")==0 || strcmp(argv[1],"2")==0 || strcmp(argv[1],"3")==0) && strcmp(argv[2],"a")==0){
        std::cout << "choice is: " << argv[1] << " " << argv[2] << "\n";
        mesh.Make1DCartesian(8,0,0.,1.);
    }
    else if ((strcmp(argv[1],"1")==0 || strcmp(argv[1],"2")==0) && strcmp(argv[2],"b")==0){
        std::cout << "choice is: " << argv[1] << " " << argv[2] << "\n";
        num_elem = 4;
        mesh.Init(1, segment, num_elem, num_elem, 0);
        for (int i=0; i<num_elem; i++){
            POINT p[2];
            if(i==0){
                p[0].setCoordinates(0.,0.);
                p[0].setAttribute(1);
                p[0].setIdx(1);
                p[1].setCoordinates(0.5,0.);
                p[1].setAttribute(0);
                p[1].setIdx(2);
                mesh.AddElement(i, p, 2, 1);
            }
            else if(i==1){
                p[0].setCoordinates(0.5,0.);
                p[0].setAttribute(0);
                p[0].setIdx(2);
                p[1].setCoordinates(0.75,0.);
                p[1].setAttribute(0);
                p[1].setIdx(3);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==2){
                p[0].setCoordinates(0.75,0.);
                p[0].setAttribute(0);
                p[0].setIdx(3);
                p[1].setCoordinates(0.875,0.);
                p[1].setAttribute(0);
                p[1].setIdx(4);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==3){
                p[0].setCoordinates(0.875,0.);
                p[0].setAttribute(0);
                p[0].setIdx(4);
                p[1].setCoordinates(1.,0.);
                p[1].setAttribute(1);
                p[1].setIdx(5);
                mesh.AddElement(i, p, 2, 1);
            }
        }
        mesh.FinalizeMesh();
    }
    else if((strcmp(argv[1],"1")==0 || strcmp(argv[1],"2")==0) && strcmp(argv[2],"c")==0){
        std::cout << "choice is: " << argv[1] << " " << argv[2] << "\n";
        num_elem = 2;
        mesh.Init(1, segment, num_elem, num_elem, 0);
        for (int i=0; i<num_elem; i++){
            POINT p[2];
            if(i==0){
                p[0].setCoordinates(0.,0.);
                p[0].setAttribute(1);
                p[0].setIdx(1);
                p[1].setCoordinates(0.875,0.);
                p[1].setAttribute(0);
                p[1].setIdx(2);
                mesh.AddElement(i, p, 2, 1);
            }
            else if(i==1){
                p[0].setCoordinates(0.875,0.);
                p[0].setAttribute(0);
                p[0].setIdx(2);
                p[1].setCoordinates(1.,0.);
                p[1].setAttribute(1);
                p[1].setIdx(3);
                mesh.AddElement(i, p, 2, 1);
            }
        }
        mesh.FinalizeMesh();         
    }
    else if ((strcmp(argv[1],"3")==0) && strcmp(argv[2],"b")==0){
        std::cout << "choice is: " << argv[1] << " " << argv[2] << "\n";
        num_elem = 6;
        mesh.Init(1, segment, num_elem, num_elem, 0);
        for (int i=0; i<num_elem; i++){
            POINT p[2];
            if(i==0){
                p[0].setCoordinates(0.,0.);
                p[0].setAttribute(1);
                p[0].setIdx(1);
                p[1].setCoordinates(0.125,0.);
                p[1].setAttribute(0);
                p[1].setIdx(2);
                mesh.AddElement(i, p, 2, 1);
            }
            else if(i==1){
                p[0].setCoordinates(0.125,0.);
                p[0].setAttribute(0);
                p[0].setIdx(2);
                p[1].setCoordinates(0.25,0.);
                p[1].setAttribute(0);
                p[1].setIdx(3);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==2){
                p[0].setCoordinates(0.25,0.);
                p[0].setAttribute(0);
                p[0].setIdx(3);
                p[1].setCoordinates(0.5,0.);
                p[1].setAttribute(0);
                p[1].setIdx(4);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==3){
                p[0].setCoordinates(0.5,0.);
                p[0].setAttribute(0);
                p[0].setIdx(4);
                p[1].setCoordinates(0.75,0.);
                p[1].setAttribute(0);
                p[1].setIdx(5);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==4){
                p[0].setCoordinates(0.75,0.);
                p[0].setAttribute(0);
                p[0].setIdx(5);
                p[1].setCoordinates(0.875,0.);
                p[1].setAttribute(0);
                p[1].setIdx(6);
                mesh.AddElement(i, p, 2, 0);
            }
            else if(i==5){
                p[0].setCoordinates(0.875,0.);
                p[0].setAttribute(0);
                p[0].setIdx(6);
                p[1].setCoordinates(1.,0.);
                p[1].setAttribute(1);
                p[1].setIdx(7);
                mesh.AddElement(i, p, 2, 1);
            }
        }
        mesh.FinalizeMesh();
    }
    else if ((strcmp(argv[2],"d"))==0){
        std::cout << "choice is: " << argv[1] << "\n";
        mesh.Make1DCartesian(100,0,0.,1.);
        // std::cerr << "Wrong options chosen \n";
    }

    if (strcmp(argv[1],"1")==0){
        H1_FiniteElementSpace<degree> fespace(&mesh);
        fespace.setBdrCondition(&set_Bc_1);
        Vector<double> a(2);
        a.setValue(0,1e-0);
        a.setValue(1,0.0);
        double kappa  = 1e-4;
        double source = 0.;
        FlowSolver<degree> solver(steady, AD);
        solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(a, &fespace));
        solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(kappa, &fespace));
        solver.AddLinearForm(new LinearForm<degree>(source, &fespace));
        solver.AddStabilization(new SUPGStabilization<degree>(a,kappa,&calc_tau1D_AD,source,&fespace));
        solver.writeSparse();
    }
    else if (strcmp(argv[1],"2")==0){
        H1_FiniteElementSpace<degree> fespace(&mesh);
        fespace.setBdrCondition(&set_Bc_2);
        Vector<double> a(2);
        a.setValue(0,1e-0);
        a.setValue(1,0.);
        double kappa  = 1e-4;
        double source = 1.;
        FlowSolver<degree> solver(steady, AD);
        solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(a, &fespace));
        solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(kappa, &fespace));
        solver.AddLinearForm(new LinearForm<degree>(source, &fespace));
        solver.AddStabilization(new SUPGStabilization<degree>(a,kappa,&calc_tau1D_AD,source,&fespace));
        solver.writeSparse();        
    }
    else if (strcmp(argv[1],"3")==0){
        H1_FiniteElementSpace<degree> fespace(&mesh);
        fespace.setBdrCondition(&set_Bc_1);
        Vector<double> a(2);
        a.setValue(0,1e0);
        a.setValue(1,0.);
        double kappa  = 1e-1;
        double c = 1e2;
        double source = 10.;
        FlowSolver<degree> solver(steady, ADR);
        solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(a, &fespace));
        solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(kappa, &fespace));
        solver.AddDomainIntegrator(new ReactionIntegrator<degree>(c, &fespace));
        solver.AddLinearForm(new LinearForm<degree>(source, &fespace));
        solver.AddStabilization(new VMSStabilization<degree>(a,kappa,c,&calc_tau1D_ADR, source, &fespace));
        solver.writeSparse();   
    }
    else{
        H1_FiniteElementSpace<degree> fespace(&mesh);
        fespace.setBdrCondition(&set_Bc_1);
        Vector<double> a(2);
        a.setValue(0,1.e0);
        a.setValue(1,0.);
        double kappa  = 1.e-1;
        double c = 1.e2;
        double source = 10.;
        FlowSolver<degree> solver(steady, ADR);
        solver.AddDomainIntegrator(new AdvectionIntegrator<degree>(a, &fespace));
        solver.AddDomainIntegrator(new DiffusionIntegrator<degree>(kappa, &fespace));
        solver.AddDomainIntegrator(new ReactionIntegrator<degree>(c, &fespace));
        solver.AddLinearForm(new LinearForm<degree>(source, &fespace));
        solver.AddStabilization(new VMSStabilization<degree>(a,kappa,c,&calc_tau1D_ADR, source, &fespace));
        solver.writeSparse();         
    }

    return 0;
}


void calc_tau1D_AD(double h, double ax, double Pee, double &tau){
    tau =  (h/abs(ax))*(((1+exp(-2.0*Pee))/(1-exp(-2.0*Pee)))-(1/Pee));
}

void calc_tau1D_ADR(double h, double ax, double kappa, double c, double &tau){
    tau = 1./sqrt( pow((2.*ax/h),2.) + 9.*pow((4.*kappa/(h*h)),2.) + c*c );
}

void set_Bc_1(double x, double y, double &q){
    if(x < 1e-14){
        q = 0.;
    }
    else if ((x-1.)<1e-14){
        q = 1.;
    }
}

void set_Bc_2(double x, double y, double &q){
    if(x < 1e-14){
        q = 0.;
    }
    else if ((x-1.)<1e-14){
        q = 0.;
    }
}
