#ifndef H1_ELEMENTS_HPP
#define H1_ELEMENTS_HPP

#include "element.hpp"
#include <math.h>

template<int degree>
class H1_FiniteElement_UnitSegment : public Element<degree> {
protected:
    Matrix<double> N;
    Matrix<double> dNdxi;
    Matrix<double> dNdeta;
public:
    H1_FiniteElement_UnitSegment();
    void Init_UnitSegment();
    void calcShapeFunction();
    void getShapeFns(Matrix<double> &m_in1, Matrix<double> &m_in2, Matrix<double> &m_in3){m_in1 = N; m_in2 = dNdxi; m_in3 = dNdeta;}
    double **getN(int &row, int &col){
        N.getSize_(row,col);
        return N.returnMatrix();}
    double **getdNdxi(){return dNdxi.returnMatrix();}
    double **getdNdeta(){return dNdeta.returnMatrix();}
};

template<int degree>
H1_FiniteElement_UnitSegment<degree>::H1_FiniteElement_UnitSegment(){
//    Init_UnitSegment();
}

template<int degree>
void H1_FiniteElement_UnitSegment<degree>::Init_UnitSegment(){
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
    calcShapeFunction();
}

template<int degree>
void H1_FiniteElement_UnitSegment<degree>::calcShapeFunction(){
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
    N.displayMatrix();
}

/* ---------------------------------------------------------------------------------------------------------------------- */
template<int degree>
class H1_FiniteElement_UnitTriangle : public Element<degree> {
public:
    H1_FiniteElement_UnitTriangle();
    void Init_UnitTriangle();
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
}

/* ---------------------------------------------------------------------------------------------------------------------- */

template<int degree>
class H1_FiniteElement_BiUnitSquare : public Element<degree> {
public:
    H1_FiniteElement_BiUnitSquare();
    void Init_BiUnitSquare();
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
}

/* ---------------------------------------------------------------------------------------------------------------------- */

template<int degree>
class H1_FiniteElement : public H1_FiniteElement_BiUnitSquare<degree>, 
                         public H1_FiniteElement_UnitTriangle<degree>,
                         public H1_FiniteElement_UnitSegment<degree>{ };


#endif