#ifndef H1_ELEMENTS_HPP
#define H1_ELEMENTS_HPP

#include "element.hpp"
#include <math.h>

template<int degree>
class H1_FiniteElement_BiUnitSegment : public Element<degree> {
protected:
    Matrix<double> N;
    Matrix<double> dNdxi;
    Matrix<double> dNdeta;
public:
    H1_FiniteElement_BiUnitSegment();
    void Init_UnitSegment();
    void calcShapeFunction();
    void getShapeFns(Matrix<double> &m_in1, Matrix<double> &m_in2, Matrix<double> &m_in3){m_in1 = N; m_in2 = dNdxi; m_in3 = dNdeta;}
};

template<int degree>
H1_FiniteElement_BiUnitSegment<degree>::H1_FiniteElement_BiUnitSegment(){
//    Init_UnitSegment();
}

template<int degree>
void H1_FiniteElement_BiUnitSegment<degree>::Init_UnitSegment(){
    this->Element<degree>::Init_(segment);
    for (int i=0; i<degree+1; i++){
        POINT t(-1+(double)(i)*(2./(double)(degree)), 0.);
        this->Element<degree>::AddVertex(&t, i);
    }
    this->Element<degree>::setQuadrature(0, (-1./3.)*sqrt(5. + 2.*sqrt(10./7.)), 0., (322. - 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(1, (-1./3.)*sqrt(5. - 2.*sqrt(10./7.)), 0., (322. + 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(2, 0., 0., 128./225.);
    this->Element<degree>::setQuadrature(3, (1./3.)*sqrt(5. - 2.*sqrt(10./7.)), 0., (322. + 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(4, (1./3.)*sqrt(5. + 2.*sqrt(10./7.)), 0., (322. - 13.*sqrt(70.))/900.); 
    int n = this->sizeof_q;
    N.setSize(degree+1, n);  
    dNdxi.setSize(degree+1, n);
    dNdeta.setSize(degree+1,n);
    H1_FiniteElement_BiUnitSegment<degree>::calcShapeFunction(); 
}

template<int degree>
void H1_FiniteElement_BiUnitSegment<degree>::calcShapeFunction(){
    Vector<double> xi;
    Vector<double> quad;
    this->Element<degree>::getVertexLoc(1, xi);
    this->Element<degree>::getQuadrature(1, quad);
    for(int k=0; k<xi.getLength_(); k++){
        double center = xi.getValue(k);
        for (int j=0; j<quad.getLength_(); j++){
            int count = 0;
            double p = 1;
            Vector<double> Num(xi.getLength_()-1);
            Vector<double> Den(xi.getLength_()-1);
            for (int i=0; i<xi.getLength_(); i++){
                if (abs(center-xi.getValue(i))<=1e-12){
                    continue;
                }
                else{
                    Num.setValue(count, (quad.getValue(j)-xi.getValue(i)));
                    Den.setValue(count, (center - xi.getValue(i)));
                    p *= (quad.getValue(j)-xi.getValue(i))/(center - xi.getValue(i));
                    ++count;
                }
            }
            this->N.setValue(k,j,p);
            double sum = 0.;
            for(int l=0; l<xi.getLength_()-1; l++){
                p = 1;
                for(int u=0; u<xi.getLength_()-1; u++){
                    if(l==u){
                        p *= 1./Den.getValue(u);
                    }
                    else{
                        p *= Num.getValue(u)/Den.getValue(u);
                    }
                }
                sum += p;
            }
            this->dNdxi.setValue(k,j,sum);
        }
    }
}

/* ---------------------------------------------------------------------------------------------------------------------- */
template<int degree>
class H1_FiniteElement_UnitTriangle : public Element<degree> {
protected:
    Matrix<double> N;
    Matrix<double> dNdxi;
    Matrix<double> dNdeta;
public:
    H1_FiniteElement_UnitTriangle();
    void Init_UnitTriangle();
    void calcShapeFunction();
    void getShapeFns(Matrix<double> &m_in1, Matrix<double> &m_in2, Matrix<double> &m_in3){m_in1 = N; m_in2 = dNdxi; m_in3 = dNdeta;}
};

template<int degree>
H1_FiniteElement_UnitTriangle<degree>::H1_FiniteElement_UnitTriangle(){
//    Init_UnitTriangle();
}

template<int degree>
void H1_FiniteElement_UnitTriangle<degree>::Init_UnitTriangle(){
    this->Element<degree>::Init_(triangle);
    int count = 0;
    for (int l=0; l<degree+1; l++){
        for (int m=0; m<=l; m++){
            POINT t(0+(double)(m)*(1./degree), 1-(double)(l)*(1./degree));
            this->Element<degree>::AddVertex(&t,count);
            ++count;
        }
    }
    this->Element<degree>::setQuadrature(0, 2./3., 1./6., 1./3.);
    this->Element<degree>::setQuadrature(1, 1./6., 1./6., 1./3.);
    this->Element<degree>::setQuadrature(2, 1./6., 2./3., 1./3.);  
    int n = this->sizeof_q;
    switch (degree)
    {
    case 1:{
        this->N.setSize(3,n);
        this->dNdxi.setSize(3,n);
        this->dNdeta.setSize(3,n);
        break;
    }
    case 2:{
        this->N.setSize(6,n);
        this->dNdxi.setSize(6,n);
        this->dNdeta.setSize(6,n);
        break;
    }
    case 3:{
        this->N.setSize(10,n);
        this->dNdxi.setSize(10,n);
        this->dNdeta.setSize(10,n);
        break;
    }
    default:{
        std::cerr << "Polynomial order > 3 not supported \n";
        break;
    }    
    }
    H1_FiniteElement_UnitTriangle<degree>::calcShapeFunction();
}

template<int degree>
void H1_FiniteElement_UnitTriangle<degree>::calcShapeFunction(){

    Vector<double> quad_xi, quad_eta;
    this->Element<degree>::getQuadrature(1, quad_xi);
    this->Element<degree>::getQuadrature(2, quad_eta);

    switch (degree)
    {
    case 1:{
        for(int k=0; k<quad_xi.getLength_(); k++){
            double r = quad_xi.getValue(k);
            double s = quad_eta.getValue(k);
            this->N.setValue(0, k, s);
            this->N.setValue(1, k, 1.-s-r);
            this->N.setValue(2, k, r);

            this->dNdxi.setValue(0, k, 0.);
            this->dNdxi.setValue(1, k, -1.);
            this->dNdxi.setValue(2, k, 1.);

            this->dNdeta.setValue(0, k, 1.);
            this->dNdeta.setValue(1, k, -1.);
            this->dNdeta.setValue(2, k, 0.);
        }
        break;
    }
    case 2:{
        for(int k=0; k<quad_xi.getLength_(); k++){
            double r = quad_xi.getValue(k);
            double s = quad_eta.getValue(k);
            double t = 1. - r - s;
            this->N.setValue(0, k, s*(2.*s - 1.));
            this->N.setValue(1, k, 4.*s*t);
            this->N.setValue(2, k, 4.*r*s);
            this->N.setValue(3, k, t*(2.*t - 1.));
            this->N.setValue(4, k, 4.*r*t);
            this->N.setValue(5, k, r*(2.*r - 1.));

            this->dNdxi.setValue(0, k, 0.); this->dNdeta.setValue(0, k, 4.*s - 1.);
            this->dNdxi.setValue(1, k, -4.*s); this->dNdeta.setValue(1, k, 4.- 8.*s -4.*r);
            this->dNdxi.setValue(2, k, 4.*s); this->dNdeta.setValue(2, k, 4.*r);
            this->dNdxi.setValue(3, k, -3.+ 4.*s + 4.*r); this->dNdeta.setValue(3, k, -3.+ 4.*s + 4.*r);
            this->dNdxi.setValue(4, k, 4. - 4.*s - 8.*r); this->dNdeta.setValue(4, k, -4.*r);
            this->dNdxi.setValue(5, k, 4.*r - 1.); this->dNdeta.setValue(5, k, 0.);
        }
        break;
    }
    case 3:{
        for (int k=0; k<quad_xi.getLength_(); k++){
            double r = quad_xi.getValue(k);
            double s = quad_eta.getValue(k);
            double t = 1. - r - s;
        }
        break;
    }
    default:{
        std::cerr << "Polynomial order > 3 not supported \n";
        break;
    }
    }

}

/* ---------------------------------------------------------------------------------------------------------------------- */

template<int degree>
class H1_FiniteElement_BiUnitSquare : public Element<degree> {
protected:
    Matrix<double> N;
    Matrix<double> dNdxi;
    Matrix<double> dNdeta;
public:
    H1_FiniteElement_BiUnitSquare();
    void Init_BiUnitSquare();
    void calcShapeFunction();
    void getShapeFns(Matrix<double> &m_in1, Matrix<double> &m_in2, Matrix<double> &m_in3){m_in1 = N; m_in2 = dNdxi; m_in3 = dNdeta;}
};

template<int degree>
H1_FiniteElement_BiUnitSquare<degree>::H1_FiniteElement_BiUnitSquare(){
//    H1_FiniteElement_BiUnitSquare<degree>::Init_BiUnitSquare();
}
template<int degree>
void H1_FiniteElement_BiUnitSquare<degree>::Init_BiUnitSquare(){
    this->Element<degree>::Init_(quadrilateral);
    int count = 0;
    for(int j=0; j<degree+1; j++){
        for (int i=0; i<degree+1; i++){
            POINT t(-1+(double)(i)*(2./(double)degree), -1+(double)(j)*(2./(double)degree));
            this->Element<degree>::AddVertex(&t, count);
            ++count;
        }
    }
    this->Element<degree>::setQuadrature(0, -1./sqrt(3.), -1./sqrt(3.), 1.);
    this->Element<degree>::setQuadrature(1, 1./sqrt(3.), -1./sqrt(3.), 1.);
    this->Element<degree>::setQuadrature(2, -1./sqrt(3.), 1./sqrt(3.), 1.);
    this->Element<degree>::setQuadrature(3, 1./sqrt(3.), 1./sqrt(3.), 1.);
    int n = this->sizeof_q;
    N.setSize((degree+1)*(degree+1),n);
    dNdxi.setSize((degree+1)*(degree+1),n);
    dNdeta.setSize((degree+1)*(degree+1),n);
    H1_FiniteElement_BiUnitSquare<degree>::calcShapeFunction();
}

template<int degree>
void H1_FiniteElement_BiUnitSquare<degree>::calcShapeFunction(){
    Vector<double> xi(degree+1); Vector<double> eta(degree+1);
    for (int i=0; i<degree+1; i++){
        xi.setValue(i,-1+(double)(i)*(2./(double)(degree)));
        eta.setValue(i,-1+(double)(i)*(2./(double)(degree)));
    }
    Vector<double> quad_xi, quad_eta;
    this->Element<degree>::getQuadrature(1, quad_xi);
    this->Element<degree>::getQuadrature(2, quad_eta);
    Matrix<double> l(degree+1,quad_xi.getLength_());
    Matrix<double> dldxi(degree+1,quad_xi.getLength_());
    Matrix<double> n(degree+1,quad_eta.getLength_());
    Matrix<double> dndeta(degree+1,quad_eta.getLength_());

    for(int k=0; k<xi.getLength_(); k++){
        double center = xi.getValue(k);
        for (int j=0; j<quad_xi.getLength_(); j++){
            int count = 0;
            double p = 1;
            Vector<double> Num(xi.getLength_()-1);
            Vector<double> Den(xi.getLength_()-1);
            for (int i=0; i<xi.getLength_(); i++){
                if (abs(center-xi.getValue(i))<=1e-12){
                    continue;
                }
                else{
                    Num.setValue(count, (quad_xi.getValue(j)-xi.getValue(i)));
                    Den.setValue(count, (center - xi.getValue(i)));
                    p *= (quad_xi.getValue(j)-xi.getValue(i))/(center - xi.getValue(i));
                    ++count;
                }
            }
            l.setValue(k,j,p);
            double sum = 0.;
            for(int l=0; l<xi.getLength_()-1; l++){
                p = 1;
                for(int u=0; u<xi.getLength_()-1; u++){
                    if(l==u){
                        p *= 1./Den.getValue(u);
                    }
                    else{
                        p *= Num.getValue(u)/Den.getValue(u);
                    }
                }
                sum += p;
            }
            dldxi.setValue(k,j,sum);
        }
    }
    
    for(int k=0; k<eta.getLength_(); k++){
        double center = eta.getValue(k);
        for (int j=0; j<quad_eta.getLength_(); j++){
            int count = 0;
            double p = 1;
            Vector<double> Num(eta.getLength_()-1);
            Vector<double> Den(eta.getLength_()-1);
            for (int i=0; i<eta.getLength_(); i++){
                if (abs(center-eta.getValue(i))<=1e-12){
                    continue;
                }
                else{
                    Num.setValue(count, (quad_eta.getValue(j)-eta.getValue(i)));
                    Den.setValue(count, (center - eta.getValue(i)));
                    p *= (quad_eta.getValue(j)-eta.getValue(i))/(center - eta.getValue(i));
                    ++count;
                }
            }
            n.setValue(k,j,p);
            double sum = 0.;
            for(int l=0; l<eta.getLength_()-1; l++){
                p = 1;
                for(int u=0; u<eta.getLength_()-1; u++){
                    if(l==u){
                        p *= 1./Den.getValue(u);
                    }
                    else{
                        p *= Num.getValue(u)/Den.getValue(u);
                    }
                }
                sum += p;
            }
            dndeta.setValue(k,j,sum);
        }
    }

    int idx=0;
    for(int i=0; i<degree+1; i++){
        for(int j=0; j<degree+1; j++){
            for(int k=0; k<quad_xi.getLength_(); k++){
                this->N.setValue(idx,k,l.getValue(j,k)*n.getValue(i,k));
                this->dNdxi.setValue(idx,k,dldxi.getValue(j,k)*n.getValue(i,k));
                this->dNdeta.setValue(idx,k,l.getValue(j,k)*dndeta.getValue(i,k));
            }
            ++idx;
        }
    }

}

/* ---------------------------------------------------------------------------------------------------------------------- */

template<int degree>
class H1_FiniteElement : public H1_FiniteElement_BiUnitSquare<degree>, 
                         public H1_FiniteElement_UnitTriangle<degree>,
                         public H1_FiniteElement_BiUnitSegment<degree>{ };


#endif