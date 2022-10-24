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
    void AddStabilizationAdvection_W_Advection_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &));
    void AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddStabilizationAdvection_W_Reaction_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, void(*func)(double, double, double, double &), double source);
    void AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &), double source);

    void AddVMSStabilizationDiffusion_W_Advection_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &));
    void AddVMSStabilizationDiffusion_W_Advection_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddVMSStabilizationDiffusion_W_Diffusion_U(Vector<double> a, double kappa, void(*func)(double, double, double, double &));
    void AddVMSStabilizationDiffusion_W_Diffusion_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddVMSStabilizationDiffusion_W_Reaction_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));

    void AddVMSStabilizationDiffusion_W_source_S(Vector<double> a, double kappa, void(*func)(double, double, double, double &), double source);
    void AddVMSStabilizationDiffusion_W_source_S(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &), double source);

    void AddStabilizationReaction_W_Advection_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));
    void AddStabilizationReaction_W_Diffusion_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));
    void AddStabilizationReaction_W_Reaction_U(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &));
    void AddStabilizationReaction_W_source_S(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &), double source);

    void tauInvoke(double h, double ax, double Pee, double &tauEx, void(*func)(double, double, double, double &)){func(h,ax,Pee,tauEx);}
    void tauInvoke(double h, double ax, double Pee, double c, double &tauEx, void(*func)(double, double, double, double, double &)){func(h,ax,Pee,c,tauEx);}
    void Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> qbc_e, Matrix<double> M);
    void Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> lf);
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
    std::cout << " \t Adding stabilization term (-a . (grad)(w), -tau*(a . (grad)(u))) ... \n";

    Matrix<double> N_copy, dNdxi_copy, dNdeta_copy;
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
        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

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
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(Vector<double> a, double kappa, double c, 
                                                                    void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding stabilization term (-a . (grad)(w), -tau*(a . (grad)(u))) ... \n";

    Matrix<double> N_copy, dNdxi_copy, dNdeta_copy;
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
        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

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
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.ElementMultiplication(kappa_num, kappa_num);
            kappa_num.Scale(tauEx, kappa_num);

            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, dNdx, dNdx, D);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, 
                                                                    void(*func)(double, double, double, double &)){
    std::cout << " \t Adding stabilization term (-a . (grad)(w), tau*kappa*(nabla).(grad)(u)) ... \n";

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
        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

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
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
    }
    // K = nullptr;
    // RHS = nullptr;
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(Vector<double> a, double kappa, double c, 
                                                                    void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding stabilization term (-a . (grad)(w), tau*kappa*(nabla).(grad)(u)) ... \n";

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
        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

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
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            tauEx = -1.*tauEx;
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*kappa, kappa_num);

            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, dNdx, d2Ndx2, D);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_Reaction_U(Vector<double> a, double kappa, double c, 
                                                                   void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding Stabilization through the term (-a . (grad)(w) , -tau*c*(u) ) \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);

            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*c, kappa_num);

            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, dNdx, N, D);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, 
                                                                 void(*func)(double, double, double, double &), double source){
    std::cout << " \t Adding Stabilization to RHS through the term (-a . (grad)(w), tau*source) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            tauEx = -1.*tauEx;

            Vector<double> stable_lf(e->sizeof_p);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*source, kappa_num);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, N, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, stable_lf);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationAdvection_W_source_S(Vector<double> a, double kappa, double c, 
                                                                 void(*func)(double, double, double, double, double &), double source){
    std::cout << " \t Adding Stabilization to RHS through the term (-a . (grad)(w), tau*source) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            tauEx = -1.*tauEx;

            Vector<double> stable_lf(e->sizeof_p);
            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(tauEx*source, kappa_num);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, N, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, stable_lf);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_Advection_U(Vector<double> a, double kappa, 
                                                                       void(*func)(double, double, double, double &)){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), -tau* a . (grad)(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
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

            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(kappa*tauEx, kappa_num);

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, d2Ndx2, dNdx, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }    
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_Advection_U(Vector<double> a, double kappa, double c, 
                                                                       void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), -tau* a . (grad)(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
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
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);

            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(kappa*tauEx, kappa_num);

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, d2Ndx2, dNdx, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }    
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_Diffusion_U(Vector<double> a, double kappa, 
                                                                       void(*func)(double, double, double, double &)){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), tau* kappa(nabla)^2(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            double kappa_num = -1.*kappa*tauEx*kappa;

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, d2Ndx2, d2Ndx2, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_Diffusion_U(Vector<double> a, double kappa, double c, 
                                                                       void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), tau* kappa(nabla)^2(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            double kappa_num = -1.*kappa*tauEx*kappa;

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, d2Ndx2, d2Ndx2, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_Reaction_U(Vector<double> a, double kappa, double c, 
                                                                      void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), -tau* c*(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, tauEx, c, func);
            double kappa_num = kappa*tauEx*c;

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, d2Ndx2, N, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_source_S(Vector<double> a, double kappa, 
                                                                    void(*func)(double, double, double, double &), double source){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), tau* source ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double Pee = abs(a.getValue(0))*h/(2.*kappa);
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), Pee, tauEx, func);
            double kappa_num = -1.*kappa*tauEx*source;

            Vector<double> stable_lf(e->sizeof_p);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, d2Ndx2, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, stable_lf);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddVMSStabilizationDiffusion_W_source_S(Vector<double> a, double kappa, double c, 
                                                                    void(*func)(double, double, double, double, double &), double source){
    std::cout << " \t Adding VMS Stabilization (-kappa(nabla)^2(w), tau* source ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            double kappa_num = -1.*kappa*tauEx*source;

            Vector<double> stable_lf(e->sizeof_p);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, d2Ndx2, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, stable_lf);
        }
    }
}

template<int degree>
void Stabilization<degree>::AddStabilizationReaction_W_Advection_U(Vector<double> a, double kappa, double c, 
                                                                   void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding Stabilization ( c(w) , -tau * a . (grad)(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> N; N_copy.copy(N);
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);

            Vector<double> kappa_num(e->sizeof_q);
            a_.getRow(0, kappa_num);
            kappa_num.Scale(-1.*c*tauEx, kappa_num);

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, N, dNdx, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }    
}

template<int degree>
void Stabilization<degree>::AddStabilizationReaction_W_Diffusion_U(Vector<double> a, double kappa, double c, 
                                                                   void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding Stabilization ( c(w) , tau*kappa*(nabla)^2(u) ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> N; N_copy.copy(N);
            Matrix<double> d2Ndx2; d2Ndxi2_copy.copy(d2Ndx2);
            Vector<double> q_inv_sq(num_q);
            q_inv.ElementMultiplication(q_inv, q_inv_sq);
            d2Ndx2.ElementMultiplication(q_inv_sq, d2Ndx2, 1);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            double kappa_num = c*tauEx*kappa;

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, N, d2Ndx2, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }    
}

template<int degree>
void Stabilization<degree>::AddStabilizationReaction_W_Reaction_U(Vector<double> a, double kappa, double c, 
                                                                  void(*func)(double, double, double, double, double &)){
    std::cout << " \t Adding Stabilization ( c(w) , -tau*c(u) ) ... \n ";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);

            double kappa_num = -1.*c*tauEx*c;

            Matrix<double> M(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_num, q_wt, q_detJac, N, N, M);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, ql_bc, M);
        }
    }    
}

template<int degree>
void Stabilization<degree>::AddStabilizationReaction_W_source_S(Vector<double> a, double kappa, double c, 
                                                                void(*func)(double, double, double, double, double &),double source){
    std::cout << " \t Adding Stabilization ( c(w) , tau*source ) ... \n";

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

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if (fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns2Der(d2Ndxi2_copy, d2Ndeta2_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }

            Matrix<double> N; N_copy.copy(N);

            double h = e->get_h();
            double tauEx;
            Stabilization<degree>::tauInvoke(h, a.getValue(0), kappa, c, tauEx, func);
            double kappa_num = c*tauEx*source;

            Vector<double> stable_lf(e->sizeof_p);
            LinearForm<degree>::NumericalIntegration(kappa_num, q_wt, q_detJac, N, stable_lf);
            Stabilization<degree>::Assemble(e->node_idx, gl_ien, stable_lf);
        }
    }
}


template<int degree>
void Stabilization<degree>::Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, 
                                     Vector<double>  qbc_e, Matrix<double> M){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        for(int j=0; j<e_node_idx.getLength_(); j++){
            int row = e_node_idx.getValue(i);
            int col = e_node_idx.getValue(j);

            if (g_node_idx.getValue(i) == -1 && g_node_idx.getValue(j) == -1){
                K->i = row;
                K->j = col;
                K->value = 0.;
                K->next = new AppendList;
                K = K->next;  
            }
            else if(g_node_idx.getValue(i) != -1 && g_node_idx.getValue(j)==-1){
                RHS->i = row;
                RHS->value = -1.*qbc_e.getValue(j)*M.getValue(i,j); // you subtract this on RHS
                RHS->next = new AppendList1D;
                RHS = RHS->next;
            }
            else{
                K->i = row;
                K->j = col;
                K->value = M.getValue(i,j);
                K->next = new AppendList;
                K = K->next; 
            }
        }
    }
}

template<int degree>
void Stabilization<degree>::Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> lf){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        int row = e_node_idx.getValue(i);
        double v = lf.getValue(i);
        if(g_node_idx.getValue(i) == -1){
            RHS->i = row;
            RHS->value = 0.;
            RHS->next = new AppendList1D;
            RHS = RHS->next;
        }
        else{
            RHS->i = row;
            RHS->value = -1*v; // you subtract this on the RHS
            RHS->next = new AppendList1D;
            RHS = RHS->next;
        }
    }
}

/* ---------------------------------------------------------------------------------------------------------- */
template<int degree>
class SUPGStabilization : public Stabilization<degree> {
public:
    SUPGStabilization(Vector<double> a, double kappa, void(*func)(double, double, double, double &),
                      double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr){
        this->Stabilization<degree>::sType = SUPG;
        std::cout << "Adding SUPG stabilization for Advection-Diffusion equation --> (-a . (grad)(w), -tau * Residual(u)) ... \n";

        Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(a,kappa,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(a,kappa,func);
        Stabilization<degree>::AddStabilizationAdvection_W_source_S(a,kappa,func,source);   
    }

    SUPGStabilization(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &),
                      double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr){
        this->Stabilization<degree>::sType = SUPG;
        std::cout << "Adding SUPG stabilization for Advection-Diffusion-Reaction equation --> (-a . (grad)(w), -tau * Residual(u)) ... \n";

        Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Reaction_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_source_S(a,kappa,c,func,source);  
    }
};

/* ---------------------------------------------------------------------------------------------------------- */
template<int degree>
class VMSStabilization : public Stabilization<degree> {
public:
    VMSStabilization(Vector<double> a, double kappa, void(*func)(double, double, double, double &),
                     double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr) {
        this->Stabilization<degree>::sType = VMS;
        std::cout << "Adding VMS Stabilization for Advection-Diffusion equation --> (-a . (grad)(w) - kappa*(nabla)^2(u) , -tau * Residual(u)) ... \n";
        Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(a,kappa,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(a,kappa,func);
        Stabilization<degree>::AddStabilizationAdvection_W_source_S(a,kappa,func,source);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_Advection_U(a,kappa,func);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_Diffusion_U(a,kappa,func);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_source_S(a,kappa,func,source);
    }
    VMSStabilization(Vector<double> a, double kappa, double c, void(*func)(double, double, double, double, double &),
                     double source, H1_FiniteElementSpace<degree> *fespace_ptr) : Stabilization<degree>(fespace_ptr){
        std::cout << "Adding VMS Stabilization for Advection-Diffusion-Reaction equation --> (-a . (grad)(w) - kappa*(nabla)^2(u) + cu , -tau * Residual(u)) ... \n";
        Stabilization<degree>::AddStabilizationAdvection_W_Advection_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Diffusion_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_Reaction_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationAdvection_W_source_S(a,kappa,c,func,source);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_Advection_U(a,kappa,c,func);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_Diffusion_U(a,kappa,c,func);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_Reaction_U(a,kappa,c,func);
        Stabilization<degree>::AddVMSStabilizationDiffusion_W_source_S(a,kappa,c,func,source);
        Stabilization<degree>::AddStabilizationReaction_W_Advection_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationReaction_W_Diffusion_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationReaction_W_Reaction_U(a,kappa,c,func);
        Stabilization<degree>::AddStabilizationReaction_W_source_S(a,kappa,c,func,source);
    }
};

#endif