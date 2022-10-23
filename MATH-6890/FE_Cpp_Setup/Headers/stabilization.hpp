#ifndef STAB_INTEG
#define STAB_INTEG

#include "fe_space.hpp"
#include "BiLinearForm.hpp"
#include "LinearForm.hpp"
#include "MatrixAlgebra.hpp"
#include <functional>

enum Residual { noRes=-1, A=0, D=1, R=2, AD=3, ADR=4};
enum StabilizationType {noStabilization=-1, SUPG=0, GLS=1, VMS=2};

template<int degree>
class Stabilization : public BiLinearForm<degree>, public LinearForm<degree>{
protected:
    H1_FiniteElementSpace<degree> *fes;
    AppendList *K; AppendList *head;
    AppendList1D *RHS; AppendList1D *r_head;
public:
    StabilizationType sType;
    Stabilization();
    Stabilization(H1_FiniteElementSpace<degree> *fespace_ptr);
    AppendList* returnHead(){return head;}
    AppendList1D* returnHead_r(){return r_head;}
    void AddStabilizationAdvection_W_Advection_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &));
    void AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &));
    void AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, void(*func)(double, double, double, double &), double source);
    void tauInvoke(double h, double ax, double Pee, double &tauEx, void(*func)(double, double, double, double &)){func(h,ax,Pee,tauEx);}
    void Assemble(Vector<int> e_node_idx, Matrix<double> M);
    void Assemble(Vector<int> e_node_idx, Vector<double> lf);
}; 

template<int degree>
Stabilization<degree>::Stabilization(){
    fes = nullptr;
    K = nullptr;
    head = nullptr;
    RHS = nullptr;
    r_head = nullptr;
}

template<int degree>
Stabilization<degree>::Stabilization(H1_FiniteElementSpace<degree> *fespace_ptr){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList;
    head = new AppendList; 
    head = K;
    RHS = new AppendList1D; 
    r_head = new AppendList1D;
    r_head = RHS;
    fes = fespace_ptr;
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &)){
    std::cout << "Adding stabilization term (-a . (grad)(w), -tau*(a . (grad)(u))) ... \n";

    Matrix<double> N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        e->printElementNodes();
        Vector<double> xq, yq, q_wt, q_detJac;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);

        // populating the a_vector at each of the integration points
        Matrix<double> a_(2, e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            a_.setColumn(l, a);
        }

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.ElementMultiplication(kappa_num, kappa_num);
            kappa_num.Scale(tauEx, kappa_num);

            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, dNdx, dNdx, D);
            Stabilization<degree>::Assemble(e->node_idx, D);
        }
    }
    // K = nullptr;
    // RHS = nullptr;
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, 
                                                                    void(*func)(double, double, double, double &)){
    std::cout << " Adding stabilization term (-a . (grad)(w), tau*kappa*(nabla).(grad)(u)) ... \n";

    Matrix<double> N_copy, dNdxi_copy, dNdeta_copy, d2Ndxi2_copy, d2Ndeta2_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);

        // populating the a_vector at each of the integration points
        Matrix<double> a_(2, e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            a_.setColumn(l, a);
        }

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);

            dNdx.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            tauEx = -1.*tauEx;
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*kappa, kappa_num);

            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, dNdx, d2Ndx2, D);
            Stabilization<degree>::Assemble(e->node_idx, D);
        }
    }
    // K = nullptr;
    // RHS = nullptr;
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, 
                                                                 void(*func)(double, double, double, double &), double source){
    std::cout << "Adding Stabilization to RHS through the term (-a . (grad)(w), tau*source) ... \n";

    Matrix<double> N_copy, dNdxi_copy, dNdeta_copy, d2Ndxi2_copy, d2Ndeta2_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);

        // populating the a_vector at each of the integration points
        Matrix<double> a_(2, e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            a_.setColumn(l, a);
        }

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            std::cout << "Element ID = " << i << "\t h = " << h << "\t Pee = " << Pee << "\t tauEx = " << tauEx << "\n";
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            tauEx = -1.*tauEx;

            Vector<double> stable_lf(e->sizeof_p);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*source, kappa_num);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, N, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, stable_lf);
        }
    }
}

template<int degree>
void Stabilization<degree>::Assemble(Vector<int> e_node_idx, Matrix<double> M){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        for(int j=0; j<e_node_idx.getLength_(); j++){
            int row = e_node_idx.getValue(i);
            int col = e_node_idx.getValue(j);
            K->i = row;
            K->j = col;
            K->value = M.getValue(i,j);
            K->next = new AppendList;
            K = K->next; 
        }
    }
}

template<int degree>
void Stabilization<degree>::Assemble(Vector<int> e_node_idx, Vector<double> lf){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        int row = e_node_idx.getValue(i);
        double v = lf.getValue(i);
        RHS->i = row;
        RHS->value = v;
        RHS->next = new AppendList1D;
        RHS = RHS->next;
    }
}

/* ---------------------------------------------------------------------------------------------------------- */
template<int degree>
class SUPGStabilization : public Stabilization<degree> {
public:
    SUPGStabilization(Residual rType, Vector<double> a, double kappa, void(*func)(double, double, double, double &),
                      double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr){
        this->Stabilization<degree>::sType = SUPG;
        std::cout << "Adding SUPG stabilization (-a . (grad)(w), -tau * Residual(u)) ... \n";
        if (rType == AD){
            Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(a,kappa,func);
            Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(a,kappa,func);
            Stabilization<degree>::AddStabilizationAdvection_W_source_S(a,kappa,func,source);            
        }
        else if (rType == ADR){

        }
    }
};

/* ---------------------------------------------------------------------------------------------------------- */
template<int degree>
class VMSStabilization : public Stabilization<degree> {
public:
    VMSStabilization(Residual rType, Vector<double> a, double kappa, void(*func)(double, double, double, double &),
                     double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr) {
        this->Stabilization<degree>::sType = VMS;
        if (rType == AD){
            std::cout << "Adding VMS Stabilization (-a . (grad)(w) - kappa*(nabla)^2(u) , -tau * Residual(u)) ... \n";

        }
        else if (rType == ADR){
            std::cout << "Adding VMS Stabilization (-a . (grad)(w) - kappa*(nabla)^2(u) + cu , -tau * Residual(u)) ... \n";
        }
    }
};

#endif