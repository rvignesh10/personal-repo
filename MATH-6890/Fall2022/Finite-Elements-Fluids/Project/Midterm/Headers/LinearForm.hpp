#ifndef LINEARFORM_HPP
#define LINEARFORM_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"
#include <functional>

template<int degree>
class LinearForm {
protected:
    H1_FiniteElementSpace<degree> *fes;
    AppendList1D *K; AppendList1D *head;
public:
    LinearForm();
    LinearForm(const LinearForm &source);
    LinearForm(double source, H1_FiniteElementSpace<degree> *fespace_ptr);
    LinearForm(void(*func)(double, double, double &), H1_FiniteElementSpace<degree> *fespace_ptr);
    void invoke(double x, double y, double &source, void(*func)(double, double, double &));
    void AddLinearForm(double source);
    void AddLinearForm(void(*func)(double, double, double &));
    void NumericalIntegration(double source, Vector<double> q_wt, Vector<double> q_Jac, Matrix<double> Ni, Vector<double> &lf);
    void NumericalIntegration(Vector<double> source, Vector<double> q_wt, Vector<double> q_Jac, Matrix<double> Ni, Vector<double> &lf);
    void Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> lf);
    AppendList1D* returnHead(){return head;}
};

template<int degree>
LinearForm<degree>::LinearForm(){
    fes = nullptr;
    K = nullptr;
    head = nullptr;
}

template<int degree>
LinearForm<degree>::LinearForm(const LinearForm &source){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList1D;
    head = new AppendList1D;
    *fes = *source.fes;
    *K = *source.K;
    *head = *source.head;
}

template<int degree>
LinearForm<degree>::LinearForm(double source, H1_FiniteElementSpace<degree> *fespace_ptr){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList1D;
    head = new AppendList1D;
    head = K;
    fes = fespace_ptr;
    LinearForm<degree>::AddLinearForm(source);
}

template<int degree>
LinearForm<degree>::LinearForm(void(*func)(double, double, double &), H1_FiniteElementSpace<degree> *fespace_ptr){
    fes = new H1_FiniteElementSpace<degree>;
    K = new AppendList1D;
    head = new AppendList1D;
    head = K;
    fes = fespace_ptr;
    LinearForm<degree>::AddLinearForm(func);
}

template<int degree>
void LinearForm<degree>::AddLinearForm(double source){
    std::cout << "Adding Linear form ... } w source d_omega \n"; 

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    
    for(int i=0; i<fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> q_wt, q_detJac;
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);
        if(e->geometry == segment){ 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        else if(e->geometry == triangle){
            fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        else{
            fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        Vector<double> lf(e->sizeof_p);
        Matrix<double> N; N_copy.copy(N);
        LinearForm<degree>::NumericalIntegration(source, q_wt, q_detJac, N, lf);
        LinearForm<degree>::Assemble(e->node_idx, gl_ien, lf);
    }
    K = nullptr;
}

template<int degree>
void LinearForm<degree>::AddLinearForm(void(*func)(double, double, double &)){
    std::cout << "Adding Linear form ... } w source d_omega \n"; 

    Matrix<double>N_copy, dNdxi_copy, dNdeta_copy;
    for(int i=0; i<fes->mesh_file->GetNE(); i++){
        Element<degree>* e = fes->mesh_file->GetElement(i+1);
        Vector<double> xq, yq, q_wt, q_detJac;
        e->getQuadrature(1, xq);
        e->getQuadrature(2, yq);
        e->getQuadrature(3, q_wt);
        e->getQuadrature(4, q_detJac);
        Vector<int> gl_ien(e->sizeof_p);
        fes->Global_IEN.getRow(i, gl_ien);

        Vector<double> source(e->sizeof_q);
        for(int j=0; j<e->sizeof_q; j++){
            double s;
            LinearForm<degree>::invoke(xq.getValue(j), yq.getValue(j), s);
            source.setValue(j, s);
        }

        if(e->geometry == segment){ 
            fes->h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        else if(e->geometry == triangle){
            fes->h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        else{
            fes->h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(N_copy, dNdxi_copy, dNdeta_copy);
        }
        Vector<double> lf(e->sizeof_p);
        Matrix<double> N; N_copy.copy(N);
        LinearForm<degree>::NumericalIntegration(source, q_wt, q_detJac, N, lf);
        LinearForm<degree>::Assemble(e->node_idx, gl_ien, lf);
    }
    K = nullptr;
}

template<int degree>
void LinearForm<degree>::invoke(double x, double y, double &source, void(*func)(double, double, double &)){
    func(x,y,source);
}

template<int degree>
void LinearForm<degree>::NumericalIntegration(double source, Vector<double> q_wt, 
                                              Vector<double> q_Jac, Matrix<double> N, Vector<double> &lf){
    Matrix<double> N_prime; N.copy(N_prime);
    N_prime.ElementMultiplication(q_Jac, N_prime, 1);
    N_prime.Multiply(q_wt, lf);
    lf.Scale(source, lf);
}

template<int degree>
void LinearForm<degree>::NumericalIntegration(Vector<double> source, Vector<double> q_wt, 
                                              Vector<double> q_Jac, Matrix<double> N, Vector<double> &lf){
    Matrix<double> N_prime; N.copy(N_prime);
    N_prime.ElementMultiplication(q_Jac, N_prime, 1);
    N_prime.ElementMultiplication(source, N_prime, 1);
    N_prime.Multiply(q_wt, lf);
}

template<int degree>
void LinearForm<degree>::Assemble(Vector<int> e_node_idx, Vector<int> g_node_idx, Vector<double> lf){
    for(int i=0; i<e_node_idx.getLength_(); i++){
        int row = e_node_idx.getValue(i);
        double v = lf.getValue(i);
        if(g_node_idx.getValue(i) == -1){
            K->i = row;
            K->value = 0.;
            K->next = new AppendList1D;
            K = K->next;
        }
        else{
            K->i = row;
            K->value = v;
            K->next = new AppendList1D;
            K = K->next;
        }
    }
}

#endif