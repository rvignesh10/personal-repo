#ifndef BILINEARFORM_HPP
#define BILINEARFORM_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"
#include <functional>

enum IntegratorType {invalidIntegrator=-1, Advection=0, Diffusion=1, Mass=2};

template<int degree>
class BiLinearForm {
protected:
    H1_FiniteElementSpace<degree> *fes;
    AppendList *K; AppendList *head;
    AppendList1D *RHS; AppendList1D *r_head;
public :
    IntegratorType IntegType;
    BiLinearForm();
    BiLinearForm(const BiLinearForm<degree> &source);
    BiLinearForm(H1_FiniteElementSpace<degree> *fespace_ptr);
    AppendList* returnHead(){return head;}
    AppendList1D* returnHead_r(){return r_head;}
    void AddDiffusionIntegrator(double kappa);
    void AddDiffusionIntegrator(Matrix<double> &kappa);
    void AddDiffusionIntegrator(void(*func)(double , double, double & ));
    void AddMassIntegrator();
    void AddAdvectionIntegrator(double a);
    void AddAdvectionIntegrator(Vector<double> &a);
    void AddAdvectionIntegrator(void (*func)(double, double, Vector<double> &));
    void invoke(double x, double y, Vector<double> &a, 
                void (*func)(double, double, Vector<double> &));
    void invoke(Vector<double> xq, Vector<double> yq, Vector<double> &k, void(*func)(double, double, double &));
    void NumericalIntegration2D(double kappa, Vector<double> q_wt, Vector<double> jac_det,
                                Matrix<double> N1, Matrix<double> N2, Matrix<double> &M);
    void NumericalIntegration2D(Vector<double> a, Vector<double> q_wt, Vector<double> jac_det,
                                Matrix<double> N1, Matrix<double> N2, Matrix<double> &M);
    void Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> qbc_e, Matrix<double> &M);
};

template<int degree>
BiLinearForm<degree>::BiLinearForm(){
    fes = nullptr;
    K = nullptr; head = nullptr;
    RHS = nullptr; r_head = nullptr;
    IntegType = invalidIntegrator;
}

template<int degree>
BiLinearForm<degree>::BiLinearForm(const BiLinearForm<degree> &source){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList;
    RHS = new AppendList1D;
    *fes = *source.fes;
    *K = *source.K;
    *head = *source.head;
    *RHS = *source.RHS;
    *r_head = *source.r_head;
    IntegType = source.IntegType;
}

template<int degree>
BiLinearForm<degree>::BiLinearForm(H1_FiniteElementSpace<degree> *fespace_ptr){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList;
    head = new AppendList;
    RHS = new AppendList1D;
    r_head = new AppendList1D;
    head = K;
    r_head = RHS;
    fes = fespace_ptr;
}

template<int degree>
void BiLinearForm<degree>::AddDiffusionIntegrator(double kappa){
    std::cout << "Forming Diffusion matrix = -} w,x kappa u,x d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
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
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdx, dNdx, D);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
        else{
            if (e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            Matrix<double> D2(e->sizeof_p, e->sizeof_p); 
            // integrate Ni,x kappa Nj,x
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdx, dNdx, D);
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdy, dNdy, D2);
            D.Add(D2, D);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }     
    }

    K = nullptr;
    RHS = nullptr;
}

template<int degree>
void BiLinearForm<degree>::AddDiffusionIntegrator(Matrix<double> &kappa){
    std::cout << "Forming Diffusion matrix =  -} w,x kappa u,x d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if(fes->mesh_file->GetDim() == 1){
            std::cerr << "Not supported for 1D problems - use other overload \n";
        }
        else{
            if (e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            Matrix<double> D2(e->sizeof_p, e->sizeof_p);
            Matrix<double> D3(e->sizeof_p, e->sizeof_p);
            Matrix<double> D4(e->sizeof_p, e->sizeof_p); 
            // integrate Ni,x kappa Nj,x
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,0), q_wt, q_detJac, dNdx, dNdx, D);
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,1), q_wt, q_detJac, dNdx, dNdy, D2);
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,1), q_wt, q_detJac, dNdy, dNdx, D3);
            BiLinearForm<degree>::NumericalIntegration2D(-1.*kappa.getValue(1,1), q_wt, q_detJac, dNdy, dNdy, D4);
            D.Add(D2, D);
            D.Add(D3, D);
            D.Add(D4, D);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }     
    }

    K = nullptr;
    RHS = nullptr;
}

template<int degree>
void BiLinearForm<degree>::AddDiffusionIntegrator(void(*func)(double, double, double &)){
    std::cout << "Forming Diffusion matrix = -} w,x kappa u,x d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        Vector<double> xq, yq, kappa_q;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);

        this->invoke(xq, yq, kappa_q, func);

        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
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
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_q, q_wt, q_detJac, dNdx, dNdx, D);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }
        else{
            if (e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            Matrix<double> D2(e->sizeof_p, e->sizeof_p); 
            // integrate Ni,x kappa Nj,x
            BiLinearForm<degree>::NumericalIntegration2D(kappa_q, q_wt, q_detJac, dNdx, dNdx, D);
            BiLinearForm<degree>::NumericalIntegration2D(kappa_q, q_wt, q_detJac, dNdy, dNdy, D2);
            D.Add(D2, D);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, D);
        }     
    }

    K = nullptr;
    RHS = nullptr;
}

template<int degree>
void BiLinearForm<degree>::AddMassIntegrator(){
    std::cout << "Forming Mass matrix  =  } w u d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(fes->mesh_file->GetDim() == 1){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy); 
        }
        else{
            if (e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
        }
        Matrix<double> N1; N_copy.copy(N1);
        Matrix<double> N2; N_copy.copy(N2);
        Matrix<double> M(e->sizeof_p, e->sizeof_p);
        BiLinearForm<degree>::NumericalIntegration2D(1., q_wt, q_detJac, N1, N2, M);
        BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, M);
    }
    K = nullptr;
}

template<int degree>
void BiLinearForm<degree>::AddAdvectionIntegrator(double a){
    std::cout << "Forming Advection matrix = -} w,x (au) d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac, aq;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        aq.setSize(2); // hard coded because max 2

        // populating the a_vector at each of the integration points
        Matrix<double> a_(2, e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            aq.setValue(0, -a); aq.setValue(1, -a);
            a_.setColumn(l, aq);
        }

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> dNdx;
            dNdxi_copy.copy(dNdx);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a_.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
        else{
            if(e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a_.getRow(0,ax);
            a_.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N, A2);
            A.Add(A2, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
    }
    K = nullptr;
    RHS = nullptr;
}

template<int degree>
void BiLinearForm<degree>::AddAdvectionIntegrator(Vector<double> &a){
    std::cout << "Forming Advection matrix = -} w,x (au) d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);

        // populating the a_vector at each of the integration points
        Matrix<double> a_(fes->mesh_file->GetDim(), e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            a.Scale(-1.,a);
            a_.setColumn(l, a);
        }

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            Matrix<double> N; N_copy.copy(N);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a_.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
        else{
            if(e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a_.getRow(0,ax);
            a_.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N, A2);
            A.Add(A2, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
    }
    K = nullptr;
    RHS = nullptr;
}


template<int degree>
void BiLinearForm<degree>::AddAdvectionIntegrator(void (*func)(double, double, Vector<double> &)){
    std::cout << "Forming Advection matrix = -} w,x (au) d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac, aq;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        aq.setSize(2); // hard coded to max dimension

        // populating the a_vector at each of the integration points
        Matrix<double> a(2, e->sizeof_q); // hard coded to max dimension
        for(int l=0; l<e->sizeof_q; l++){
            BiLinearForm<degree>::invoke(xq.getValue(l), yq.getValue(l), aq, func);
            aq.Scale(-1,aq);
            a.setColumn(l, aq);
        }

        Vector<int> gl_ien(e->sizeof_p);
        Vector<double> ql_bc(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        fes->qbc.getRow(i, ql_bc);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> dNdx; dNdxi_copy.copy(dNdx);
            dNdx.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
        else{
            if(e->geometry == triangle){
                fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            else{
                fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            }
            // need to find Ni,x and Ni,y first
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            Matrix<double> dNdy(e->sizeof_p, e->sizeof_q);
            for(int l=0; l<e->sizeof_p; l++){
                for(int m=0; m<e->sizeof_q; m++){
                    Vector<double> dN_natural(fes->mesh_file->GetDim());
                    dN_natural.setValue(0,dNdxi_copy.getValue(l,m));
                    dN_natural.setValue(1,dNdeta_copy.getValue(l,m));
                    e->q->Jinv.Multiply(dN_natural, dN_natural);
                    dNdx.setValue(l,m,dN_natural.getValue(0));
                    dNdy.setValue(l,m,dN_natural.getValue(1));
                }
            }
            Matrix<double> N; N_copy.copy(N);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a.getRow(0,ax);
            a.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            BiLinearForm<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N, A);
            BiLinearForm<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N, A2);
            A.Add(A2, A);
            BiLinearForm<degree>::Assemble(e->node_idx, gl_ien, ql_bc, A);
        }
    }
    K = nullptr;
    RHS = nullptr;
}

template<int degree>
void BiLinearForm<degree>::NumericalIntegration2D(double kappa, Vector<double> q_wt, Vector<double> jac_det, 
                                                Matrix<double> N1, Matrix<double> N2, Matrix<double> &M){
    q_wt.ElementMultiplication(jac_det, q_wt);
    int m1, n1; N2.getSize_(m1, n1);
    Matrix<double> N2_T(n1, m1);
    N2.Transpose(N2_T);
    N2_T.Multiply(kappa, N2_T);
    N1.ElementMultiplication(q_wt, N1, 1);
    if(M.checkSquare()){
        int md, nd; M.getSize_(md, nd);
        if (md == m1){
            N1.Multiply(N2_T, M);
        }
        else{
            std::cerr << "Dimension mismatch between dNdxi and D, numerical integration error \n";
        }
    }
    else{
        std::cerr << "error in dimensions of D matrix \n";
    }
}

template<int degree>
void BiLinearForm<degree>::NumericalIntegration2D(Vector<double> a, Vector<double> q_wt, Vector<double> jac_det,
                                Matrix<double> N1, Matrix<double> N2, Matrix<double> &M){
    q_wt.ElementMultiplication(jac_det, q_wt);
    int m1, n1; N2.getSize_(m1, n1);
    Matrix<double> N2_T(n1, m1);
    N2.Transpose(N2_T);
    N2_T.ElementMultiplication(a, N2_T, 2);
    N1.ElementMultiplication(q_wt, N1, 1);
    if(M.checkSquare()){
        int md, nd; M.getSize_(md, nd);
        if (md == m1){
            N1.Multiply(N2_T, M);
        }
        else{
            std::cerr << "Dimension mismatch between dNdxi and D, numerical integration error \n";
        }
    }
    else{
        std::cerr << "error in dimensions of D matrix \n";
    }
}

template<int degree>
void BiLinearForm<degree>::Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, 
                                    Vector<double> qbc_e, Matrix<double> &M){
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
                std::cout << "we do get here, " << row << "\n";
                RHS->value = -1.*qbc_e.getValue(j)*M.getValue(i,j);
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
void BiLinearForm<degree>::invoke(double x, double y, Vector<double> &a, 
                void (*func)(double, double, Vector<double> &)){
    func(x,y,a);
}

template<int degree>
void BiLinearForm<degree>::invoke(Vector<double> xq, Vector<double> yq, Vector<double> &kq,
                                  void(*func)(double, double, double &)){
    kq.setSize(xq.getLength_());
    for(int i=0; i<xq.getLength_(); i++){
        double val;
        func(xq.getValue(i), yq.getValue(i), val);
        if (val > 0.){
            kq.setValue(i, -1.*val);
        }
        else{
            std::cerr << "kappa values are becoming negative \n";
        }
    }
}

/* ------------------------------------------------------------------------------------------- */

template<int degree>
class DiffusionIntegrator : public BiLinearForm<degree> {
public :
    DiffusionIntegrator(double kappa, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->BiLinearForm<degree>::IntegType = Diffusion;
        BiLinearForm<degree>::AddDiffusionIntegrator(kappa);
    }
    DiffusionIntegrator(Matrix<double> kappa, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->BiLinearForm<degree>::IntegType = Diffusion;
        BiLinearForm<degree>::AddDiffusionIntegrator(kappa);
    }
    DiffusionIntegrator(void(*func)(double, double, double &), 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->BiLinearForm<degree>::IntegType = Diffusion;
        BiLinearForm<degree>::AddDiffusionIntegrator(func);
    }
};


/* ------------------------------------------------------------------------------------------- */

template<int degree>
class AdvectionIntegrator : public BiLinearForm<degree> {
public :
    AdvectionIntegrator(double a, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->IntegType = Advection;
        BiLinearForm<degree>::AddAdvectionIntegrator(a);
    }
    AdvectionIntegrator(Vector<double> &a,
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->IntegType = Advection;
        BiLinearForm<degree>::AddAdvectionIntegrator(a);
    }
    AdvectionIntegrator(void (*func)(double, double, Vector<double> &), 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr){
        this->IntegType = Advection;
        BiLinearForm<degree>::AddAdvectionIntegrator(func);
    }
};

/* ------------------------------------------------------------------------------------------- */

template<int degree>
class MassIntegrator : public BiLinearForm<degree> {
public :
    MassIntegrator(H1_FiniteElementSpace<degree> *fespace_ptr) : BiLinearForm<degree>(fespace_ptr) {
        this->BiLinearForm<degree>::IntegType = Mass;
        this->AddMassIntegrator();
    }
};


#endif 