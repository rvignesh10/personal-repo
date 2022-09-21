#ifndef FESPACE_HPP
#define FESPACE_HPP

#include "H1_FiniteElements.hpp"
#include "mesh.hpp"


template<int degree>
class H1_FiniteElementSpace{
private:  
    H1_FiniteElement<degree> h1_fe;
public:
    template<int nx, int ny>
    H1_FiniteElementSpace(Mesh<nx,ny,degree> *mesh);
};

template<int degree>
template<int nx, int ny>
H1_FiniteElementSpace<degree>::H1_FiniteElementSpace(Mesh<nx,ny,degree> *mesh){

    Vector<double> wt;
    Matrix<double> Ncopy, dNdxicopy, dNdetacopy;

    switch (mesh->GetDim())
    {
    case 1:{
        h1_fe.H1_FiniteElement_UnitSegment<degree>::Init_UnitSegment(); 
        h1_fe.H1_FiniteElement_UnitSegment<degree>::getQuadrature(3, wt);
        h1_fe.getShapeFns(Ncopy, dNdxicopy, dNdetacopy);

        break;
    }
    case 2:{
        if (mesh->GetGeometry() == triangle){
            h1_fe.H1_FiniteElement_UnitTriangle<degree>::Init_UnitTriangle(); 
        }
        else{
            h1_fe.H1_FiniteElement_BiUnitSquare<degree>::Init_BiUnitSquare(); 
        }
        break;
    }
    default:{
        std::cerr << "Dimensions > 2 are not supported \n";
        break;
    }
    } 

    int ne = mesh->GetNE();
    for (int i=0; i<ne; i++){  
        mesh->Mesh<nx,ny,degree>::ElemTransformation(i,Ncopy,dNdxicopy,dNdetacopy,wt);
    }

    std::cout << "deleted Ncopy, and other shit \n";
}


#endif