#ifndef H1_ELEMENTS_HPP
#define H1_ELEMENTS_HPP

#include "element.hpp"
#include <math.h>

template<const int degree>
class H1_FiniteElement_1D : public Element {
public:
    H1_FiniteElement_1D();
    void calcShapeFn();
    void calcDShapeFn();
};

template<const int degree>
H1_FiniteElement_1D<degree>::H1_FiniteElement_1D(){
    this->Element<degree>::Init_(segment);
    for (int i=0; i<degree+1; i++){
        POINT t(0, -1+(double)(i)*(2./(double)(degree)));
        this->Element<degree>::AddVertex(t, i);
    }
    this->Element<degree>::setQuadrature(0, (-1./3.)*sqrt(5. + 2.*sqrt(10./7.)), 0., (322. - 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(1, (-1./3.)*sqrt(5. - 2.*sqrt(10./7.)), 0., (322. + 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(2, 0., 0., 128./225.);
    this->Element<degree>::setQuadrature(3, (1./3.)*sqrt(5. - 2.*sqrt(10./7.)), 0., (322. + 13.*sqrt(70.))/900.);
    this->Element<degree>::setQuadrature(4, (1./3.)*sqrt(5. + 2.*sqrt(10./7.)), 0., (322. - 13.*sqrt(70.))/900.);
}

#endif