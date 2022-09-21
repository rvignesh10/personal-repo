#ifndef FESPACE_HPP
#define FESPACE_HPP

#include "H1_FiniteElements.hpp"
#include "mesh.hpp"


template<int degree>
class H1_FiniteElementSpace{
private:  
    H1_FiniteElement<degree> h1_fe;
    Matrix<int> Global_IEN;
public:
    template<int nx, int ny>
    H1_FiniteElementSpace(Mesh<nx,ny,degree> *mesh);
};

template<int degree>
template<int nx, int ny>
H1_FiniteElementSpace<degree>::H1_FiniteElementSpace(Mesh<nx,ny,degree> *mesh){

    Vector<double> wt;
    Matrix<double> Ncopy, dNdxicopy, dNdetacopy;

    int num_pts = 0;

    switch (mesh->GetDim())
    {
    case 1:{
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::Init_UnitSegment(); 
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getQuadrature(3, wt);
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(Ncopy, dNdxicopy, dNdetacopy);
        num_pts = h1_fe.H1_FiniteElement_BiUnitSegment<degree>::sizeof_p;
        break;
    }
    case 2:{
        if (mesh->GetGeometry() == triangle){
            h1_fe.H1_FiniteElement_UnitTriangle<degree>::Init_UnitTriangle();
            h1_fe.H1_FiniteElement_UnitTriangle<degree>::getQuadrature(3,wt);
            h1_fe.H1_FiniteElement_UnitTriangle<degree>::getShapeFns(Ncopy, dNdxicopy, dNdetacopy);
            num_pts = h1_fe.H1_FiniteElement_UnitTriangle<degree>::sizeof_p;
        }
        else{
            h1_fe.H1_FiniteElement_BiUnitSquare<degree>::Init_BiUnitSquare();
            h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getQuadrature(3,wt);
            h1_fe.H1_FiniteElement_BiUnitSquare<degree>::getShapeFns(Ncopy, dNdxicopy, dNdetacopy); 
            num_pts = h1_fe.H1_FiniteElement_BiUnitSquare<degree>::sizeof_p;
        }
        break;
    }
    default:{
        std::cerr << "Dimensions > 2 are not supported \n";
        break;
    }
    } 

    int ne = mesh->GetNE();
    Global_IEN.setSize(ne, num_pts);

    for (int i=0; i<ne; i++){  
        mesh->Mesh<nx,ny,degree>::ElemTransformation(i,Ncopy,dNdxicopy,dNdetacopy,wt);
        Vector<int> l_ien;
        mesh->Mesh<nx,ny,degree>::GetElemIEN(i, l_ien);
        Global_IEN.setRow(i, l_ien);
    }

    // std::cout << "Global_IEN is:  \n";
    // Global_IEN.displayMatrix();
}


#endif