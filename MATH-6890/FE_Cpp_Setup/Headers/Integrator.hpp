#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"
#include <functional>

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
    void AddDiffusionIntegrator(Matrix<double> &kappa);
    void AddMassIntegrator();
    void AddAdvectionIntegrator(double a);
    void AddAdvectionIntegrator(Vector<double> &a);
    void AddAdvectionIntegrator(void (*func)(double, double, Vector<double> &));
    void invoke(double x, double y, Vector<double> &a, 
                void (*func)(double, double, Vector<double> &));
    void NumericalIntegration2D(double kappa, Vector<double> &q_wt, Vector<double> &jac_det,
                                Matrix<double> &N1, Matrix<double> &N2, Matrix<double> &M);
    void NumericalIntegration2D(Vector<double> &a, Vector<double> &q_wt, Vector<double> &jac_det,
                                Matrix<double> &N1, Matrix<double> &N2, Matrix<double> &M);
    void Assemble(Vector<int> &e_node_idx, Vector<int> &g_node_idx, Matrix<double> &M);
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
    std::cout << "Forming Diffusion matrix = -} w,x kappa u,x d_omega \n";

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
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            dNdxi_copy.ElementMultiplication(q_inv, dNdx, 1);
            Matrix<double> D(e->sizeof_p, e->sizeof_p);
            Integrator<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdx, dNdx, D);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, D);
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
            Integrator<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdx, dNdx, D);
            Integrator<degree>::NumericalIntegration2D(-1.*kappa, q_wt, q_detJac, dNdy, dNdy, D2);
            D.Add(D2, D);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, D);
        }     
    }

    K = nullptr;
}

template<int degree>
void Integrator<degree>::AddDiffusionIntegrator(Matrix<double> &kappa){
    std::cout << "Forming Diffusion matrix =  -} w,x kappa u,x d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

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
            Integrator<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,0), q_wt, q_detJac, dNdx, D);
            Integrator<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,1), q_wt, q_detJac, dNdx, dNdy, D2);
            Integrator<degree>::NumericalIntegration2D(-1.*kappa.getValue(0,1), q_wt, q_detJac, dNdy, dNdx, D3);
            Integrator<degree>::NumericalIntegration2D(-1.*kappa.getValue(1,1), q_wt, q_detJac, dNdy, D4);
            D.Add(D2, D);
            D.Add(D3, D);
            D.Add(D4, D);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, D);
        }     
    }

    K = nullptr;
}

template<int degree>
void Integrator<degree>::AddMassIntegrator(){
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
        Matrix<double> M(e->sizeof_p, e->sizeof_p);
        Integrator<degree>::NumericalIntegration2D(1., q_wt, q_detJac, N_copy, N_copy, M);
        Integrator<degree>::Assemble(e->node_idx, gl_ien, M);
    }
    K = nullptr;
}

template<int degree>
void Integrator<degree>::AddAdvectionIntegrator(double a){
    std::cout << "Forming Advection matrix = -} w,x (au) d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac, aq;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        aq.setSize(fes->mesh_file->GetDim());

        // populating the a_vector at each of the integration points
        Matrix<double> a_(fes->mesh_file->GetDim(), e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            aq.setValue(0, a); aq.setValue(1, a);
            a_.setColumn(l, aq);
        }

        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            dNdxi_copy.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a_.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, N_copy, dNdx, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
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
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a_.getRow(0,ax);
            a_.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N_copy, A);
            Integrator<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N_copy, A2);
            A.Add(A2, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
        }
    }
    K = nullptr;
}

template<int degree>
void Integrator<degree>::AddAdvectionIntegrator(Vector<double> &a){
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
            a_.setColumn(l, a);
        }

        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            dNdxi_copy.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a_.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, N_copy, dNdx, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
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
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a_.getRow(0,ax);
            a_.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N_copy, A);
            Integrator<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N_copy, A2);
            A.Add(A2, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
        }
    }
    K = nullptr;
}


template<int degree>
void Integrator<degree>::AddAdvectionIntegrator(void (*func)(double, double, Vector<double> &)){
    std::cout << "Forming Advection matrix = -} w,x (au) d_omega \n";

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i< fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac, aq;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        aq.setSize(fes->mesh_file->GetDim());

        // populating the a_vector at each of the integration points
        Matrix<double> a(fes->mesh_file->GetDim(), e->sizeof_q);
        for(int l=0; l<e->sizeof_q; l++){
            Integrator<degree>::invoke(xq.getValue(l), yq.getValue(l), aq, func);
            a.setColumn(l, aq);
        }

        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        if(e->geometry == segment){
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
            int num_q = q_detJac.getLength_();
            Vector<double> q_inv(num_q);
            for(int j=0; j<num_q; j++){
                q_inv.setValue(j, 1./q_detJac.getValue(j));
            }
            Matrix<double> dNdx(e->sizeof_p, e->sizeof_q);
            dNdxi_copy.ElementMultiplication(q_inv, dNdx, 1);
            Vector<double> ax(e->sizeof_q);
            a.getRow(0, ax);
            ax.Scale(-1., ax);
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, N_copy, dNdx, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
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
            Matrix<double> A(e->sizeof_p, e->sizeof_p);
            Matrix<double> A2(e->sizeof_p, e->sizeof_p);
            Vector<double> ax(e->sizeof_q), ay(e->sizeof_q);
            a.getRow(0,ax);
            a.getRow(1,ay);
            ax.Scale(-1., ax);
            ay.Scale(-1., ay);
            Integrator<degree>::NumericalIntegration2D(ax, q_wt, q_detJac, dNdx, N_copy, A);
            Integrator<degree>::NumericalIntegration2D(ay, q_wt, q_detJac, dNdy, N_copy, A2);
            A.Add(A2, A);
            Integrator<degree>::Assemble(e->node_idx, gl_ien, A);
        }
    }
    K = nullptr;
}

template<int degree>
void Integrator<degree>::NumericalIntegration2D(double kappa, Vector<double> &q_wt, Vector<double> &jac_det, 
                                                Matrix<double> &N1, Matrix<double> &N2, Matrix<double> &M){
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
void Integrator<degree>::NumericalIntegration2D(Vector<double> &a, Vector<double> &q_wt, Vector<double> &jac_det,
                                Matrix<double> &N1, Matrix<double> &N2, Matrix<double> &M){
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
void Integrator<degree>::Assemble(Vector<int> &e_node_idx, Vector<int> &g_node_idx, Matrix<double> &M){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        for(int j=0; j<e_node_idx.getLength_(); j++){
            int row = e_node_idx.getValue(i);
            int col = e_node_idx.getValue(j);
            if (g_node_idx.getValue(i) == -1){
                K->i = row;
                K->j = col;
                K->value = 0.;
                K->next = new AppendList;
                K = K->next;    
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
void Integrator<degree>::invoke(double x, double y, Vector<double> &a, 
                void (*func)(double, double, Vector<double> &)){
    func(x,y,a);
}

/* ------------------------------------------------------------------------------------------- */

template<int degree>
class DiffusionIntegrator : public Integrator<degree> {
public :
    DiffusionIntegrator(double kappa, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
        this->Integrator<degree>::IntegType = Diffusion;
        Integrator<degree>::AddDiffusionIntegrator(kappa);
    }
    DiffusionIntegrator(Matrix<double> kappa, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
        this->Integrator<degree>::IntegType = Diffusion;
        Integrator<degree>::AddDiffusionIntegrator(kappa);
    }
};


/* ------------------------------------------------------------------------------------------- */

template<int degree>
class AdvectionIntegrator : public Integrator<degree> {
public :
    AdvectionIntegrator(double a, 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
        this->IntegType = Advection;
        Integrator<degree>::AddAdvectionIntegrator(a);
    }
    AdvectionIntegrator(Vector<double> &a,
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
        this->IntegType = Advection;
        Integrator<degree>::AddAdvectionIntegrator(a);
    }
    AdvectionIntegrator(void (*func)(double, double, Vector<double> &), 
                        H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr){
        this->IntegType = Advection;
        Integrator<degree>::AddAdvectionIntegrator(func);
    }
};

/* ------------------------------------------------------------------------------------------- */

template<int degree>
class MassIntegrator : public Integrator<degree> {
public :
    MassIntegrator(H1_FiniteElementSpace<degree> *fespace_ptr) : Integrator<degree>(fespace_ptr) {
        this->Integrator<degree>::IntegType = Mass;
        this->AddMassIntegrator();
    }
};


#endif 