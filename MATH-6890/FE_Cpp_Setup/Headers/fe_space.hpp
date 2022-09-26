#ifndef FESPACE_HPP
#define FESPACE_HPP

#include "H1_FiniteElements.hpp"
#include "mesh.hpp"
#include <fstream>
#include <sstream>

template<int degree>
class H1_FiniteElementSpace{
protected:
public:
    H1_FiniteElement<degree> h1_fe;
    Matrix<int> Global_IEN;
    Mesh<degree> *mesh_file;
    H1_FiniteElementSpace(){ mesh_file = new Mesh<degree>; }
    H1_FiniteElementSpace(Mesh<degree> *mesh);

};

template<int degree>
H1_FiniteElementSpace<degree>::H1_FiniteElementSpace(Mesh<degree> *mesh){

    mesh_file = mesh;
    Vector<double> wt;
    Matrix<double> Ncopy, dNdxicopy, dNdetacopy;

    int num_pts = 0;

    switch (mesh_file->GetDim())
    {
    case 1:{
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::Init_UnitSegment(); 
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getQuadrature(3, wt);
        h1_fe.H1_FiniteElement_BiUnitSegment<degree>::getShapeFns(Ncopy, dNdxicopy, dNdetacopy);
        num_pts = h1_fe.H1_FiniteElement_BiUnitSegment<degree>::sizeof_p;
        break;
    }
    case 2:{
        if (mesh_file->GetGeometry() == triangle){
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

    int ne = mesh_file->GetNE();
    Global_IEN.setSize(ne, num_pts);

    for (int i=0; i<ne; i++){  
        mesh_file->Mesh<degree>::ElemTransformation(i,Ncopy,dNdxicopy,dNdetacopy,wt);
        Vector<int> l_ien;
        mesh_file->Mesh<degree>::GetElemIEN(i, l_ien);
        Global_IEN.setRow(i, l_ien);
    }

    std::stringstream fileName;
    fileName << "../solveFiles/FiniteElementSpace.txt" << std::flush;
    std::cout << "Finite Element Space created : " << fileName.str() << "\n";

    std::fstream fileCreate;
    fileCreate.open(fileName.str(), std::fstream::out);
    fileCreate << "degree: " << degree << std::endl;
    fileCreate << "#Elements: " << ne << std::endl;
    fileCreate << "#Nodes: " << mesh_file->GetNNodes() << std::endl;
    fileCreate.close();
}


#endif