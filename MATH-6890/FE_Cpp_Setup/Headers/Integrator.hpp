#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"

enum IntegratorType {invalidIntegrator=-1, Advection=0, Diffusion=1, Mass=2};

template<int degree>
class Integrator {
protected:
    H1_FiniteElementSpace<degree> *fes;
    AppendList *K; AppendList *head;
public :
    IntegratorType IntegType;
    Integrator();
    Integrator(const Integrator<degree> &source);
    Integrator(H1_FiniteElementSpace<degree> *fespace_ptr);
    AppendList* returnHead(){return head;}
    void AddDiffusionIntegrator(double kappa);
    double NumericalIntegration1D(Vector<double> &q_wt, Vector<double> &funcVal){return q_wt.dotProduct(funcVal);}
};

template<int degree>
Integrator<degree>::Integrator(){
    fes = nullptr;
    K = nullptr; head = nullptr;
    IntegType = invalidIntegrator;
}

template<int degree>
Integrator<degree>::Integrator(const Integrator<degree> &source){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList;
    *fes = *source.fes;
    *K = *source.K;
    *head = *source.head;
    IntegType = source.IntegType;
}

template<int degree>
Integrator<degree>::Integrator(H1_FiniteElementSpace<degree> *fespace_ptr){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList;
    head = new AppendList;
    head = K;
    fes = fespace_ptr;
}

template<int degree>
void Integrator<degree>::AddDiffusionIntegrator(double kappa){
    
    std::cout << "Forming Diffusion matrix -} w,x kappa u,x d_omega \n";
    int NE = fes->mesh_file->GetNE();

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;

    for(int i=0; i< NE; i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Matrix<double> D(e->sizeof_p, e->sizeof_p);
        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv_detJac(num_q);
            for(int j=0; j<num_q; j++){
                q_inv_detJac.setValue(j, 1./q_detJac.getValue(j));
            }
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_p; m++){
                    Vector<double> dNldxi(e->sizeof_q), dNmdxi(e->sizeof_q);
                    dNdxi_copy.getRow(l, dNldxi);
                    dNdxi_copy.getRow(m, dNmdxi);
                    Vector<double> func(e->sizeof_q);
                    q_inv_detJac.ElementMultiplication(dNldxi, dNldxi);
                    dNldxi.ElementMultiplication(dNmdxi, func);
                    func.Scale(kappa, func);
                    D.setValue(l,m,NumericalIntegration1D(q_wt, func));

                    int row = e->node_idx.getValue(l);  
                    int col = e->node_idx.getValue(m);
                    if (gl_ien.getValue(l) == -1){
                        K->i = row; 
                        K->j = col;
                        K->value = 0.;
                        K->next = new AppendList;
                        K = K->next;
                    }
                    else{
                        K->i = row; 
                        K->j = col;
                        K->value = D.getValue(l,m);
                        K->next = new AppendList;
                        K = K->next;
                    }
                }
            }
        }     
    }

    K = nullptr;
}


/* ------------------------------------------------------------------------------------------- */

template<int degree>
class DiffusionIntegrator : public Integrator<degree> {
public :
    DiffusionIntegrator(double kappa, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
    Integrator<degree>::AddDiffusionIntegrator(kappa);
    this->Integrator<degree>::IntegType = Diffusion;
    }
};


/* ------------------------------------------------------------------------------------------- */

template<int degree>
class MassIntegrator : public Integrator<degree> {
public :
    MassIntegrator(H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr) {
        this->Integrator<degree>::IntegType = Mass;
    }
};


#endif 