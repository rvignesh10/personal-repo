#ifndef MESH_HPP
#define MESH_HPP

#include "MatrixAlgebra.hpp"
#include "element.hpp"

/* only generating maximum 2D mesh */

template<int degree>
class Mesh{
private:
    Element<degree> *e;
    POINT *nodes;
    int dim, nx, ny;
    double x1, x2, y1, y2;
    int num_elem;
    int num_nodes;
    Geom geometry;
public:
    Mesh(){e = nullptr; nodes = nullptr; dim=0; nx=0; ny=0;x1 = 0.; x2 = 0.; y1=0.; y2 = 0.; num_elem = 0; num_nodes = 0; geometry = invalid;}
    void Make1DCartesian(int n1, int n2, double x_left, double x_right);
    void Make2DCartesian(int n1, int n2,double x_left, double x_right, double y_bottom, double y_top, Geom g);
    void MakeRectMesh();
    void MakeTriMesh();
    int GetNNodes(){return num_nodes;}
    int GetNE(){return num_elem;}
    int GetDim(){return dim;}
    Geom GetGeometry() {return geometry;}
    void ElemTransformation(int idx,Matrix<double> &N, Matrix<double> &dNdxi, Matrix<double> &dNdeta, Vector<double> &w);
    void GetElemIEN(int idx, Vector<int> &local_ien);
    Element<degree>* GetElement(int i);
    ~Mesh(){
        delete[] e;
        delete[] nodes;
    }
};

template<int degree>
void Mesh<degree>::Make1DCartesian(int n1, int n2,double x_left, double x_right){
    nx = n1; ny = n2;
    e = new Element<degree>[nx];
    num_elem = nx;
    x1 = x_left;
    x2 = x_right;
    y1 = 0.; y2 = 0.;
    dim = 1;
    geometry = segment;

    double dx = (x2-x1)/(double)(nx);

    nodes = new POINT[degree*num_elem + 1];
    num_nodes = degree*num_elem + 1;
    for (int i=0; i<nx*degree+1; i++){
        (nodes+i)->setCoordinates(x1+(double)(i)*(dx)/(double)(degree),0.);
        (nodes+i)->setIdx(i+1);
        
        if (i==0 || i==nx*degree){
            (nodes+i)->setAttribute(1);
        }
        else{
            (nodes+i)->setAttribute(0);
        }

    }
    

    int k = 0;
    for (int i=0; i<nx; i++){
        if (i==0 || i==nx-1){
            (e+i)->Element<degree>::Init(i+1, geometry, 1);
        }
        else{
            (e+i)->Element<degree>::Init(i+1, geometry, 0);
        }
        for (int j=0; j<=degree; j++){
            (e+i)->Element<degree>::AddVertex((nodes+k), j);
            ++k;
        }
        k-=1;
        
    }
}
template<int degree>
void Mesh<degree>::Make2DCartesian(int n1, int n2,double x_left, 
                           double x_right, double y_bottom, double y_top, Geom g){
    
    nx = n1; ny = n2;
    x1 = x_left; x2 = x_right;
    y1 = y_bottom; y2 = y_top;

    dim = 2;
    geometry = g;

    if (g == triangle){
        MakeTriMesh();
    }
    else{
        MakeRectMesh();
    }

}
template<int degree>
void Mesh<degree>::MakeRectMesh(){
    e = new Element<degree>[nx*ny];
    num_elem = nx*ny;
    double dx = (x2-x1)/((double)(nx));
    double dy = (y2-y1)/((double)(ny));

    nodes = new POINT[(nx*degree+1)*(ny*degree+1)];
    num_nodes = (nx*degree+1)*(ny*degree+1);
    int count = 0;
    for (int j=0; j<ny*degree+1; j++){
        for (int i=0; i<nx*degree+1; i++){
            (nodes+count)->setCoordinates(x1+(double)(i)*(dx)/((double)(degree)),
                                        y1+(double)(j)*(dy)/((double)(degree)));
            (nodes+count)->setIdx(count+1);
            if (i==0 || i==nx*degree || j==0 || j==ny*degree){
                (nodes+count)->setAttribute(1);
            }
            else{
                (nodes+count)->setAttribute(0);
            }
            ++count;
        }
    }

    int el_count = 0;
    for (int J=0; J<ny; J++){
        int j = J*degree;
        for (int I=0; I<nx; I++){
            int i = I*degree;
            count = 0;
            if (I==0 || I==nx-1 || J==0 || J==ny-1){
                // boundary elements marked 1
                (e+el_count)->Element<degree>::Init(el_count+1, geometry, 1);
            }
            else{
                (e+el_count)->Element<degree>::Init(el_count+1, geometry, 0);
            }
            for (int l=0; l<degree+1; l++){
                for (int m=0; m<degree+1; m++){
                    (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
                    ++count;
                }
            }
            ++el_count;
        }
    }
}
template<int degree>
void Mesh<degree>::MakeTriMesh(){
    e = new Element<degree>[2*nx*ny];
    num_elem = 2*nx*ny;
    double dx = (x2-x1)/(double)(nx);
    double dy = (y2-y1)/(double)(ny);

    nodes = new POINT[(nx*degree+1)*(ny*degree+1)];
    num_nodes = (nx*degree+1)*(ny*degree+1);
    int count = 0;
    for (int j=0; j<ny*degree+1; j++){
        for (int i=0; i<nx*degree+1; i++){
            (nodes+count)->setCoordinates(x1+(double)(i)*dx/(double)(degree),
                                        y1+(double)(j)*dy/(double)(degree));
            (nodes+count)->setIdx(count+1);
            if (i==0 || i==nx*degree || j==0 || j==ny*degree){
                (nodes+count)->setAttribute(1);
            }
            else{
                (nodes+count)->setAttribute(0);
            }
            ++count;
        }
    }

    int el_count = 0;
    for (int J=0; J<ny; J++){
        int j = J*degree;
        for (int I=0; I<nx; I++){
            int i = I*degree;
            // upper triangle
            (e+el_count)->Element<degree>::Init(el_count+1, geometry, 0); 
            count = 0;
            // for (int l=0; l<degree+1; l++){
            //     for (int m=0; m<=l; m++){
            //         (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
            //         ++count;
            //     }
            // }
            for (int m=degree; m>=0; m--){
                for (int l=degree; l>=m; l--){
                    (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
                    ++count;
                }
            }
            ++el_count;
            //lower triangle
            (e+el_count)->Element<degree>::Init(el_count+1, geometry, 0); 
            count = 0;
            for(int l=degree; l>=0; l--){
                for(int m=l; m<degree+1; m++){
                    (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
                    ++count;
                }
            }
            // for(int l=degree; l>=0; l--){
            //     for (int m=degree; m>=l; m--){
            //         (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
            //         ++count;
            //     }
            // }
            // for(int l=0; l<degree+1; l++){
            //     for(int m=0; m<=l; m++){
            //         (e+el_count)->Element<degree>::AddVertex(nodes+(j+l)*(nx*degree+1)+(i+m),count);
            //         ++count;
            //     }
            // }
            ++el_count;
        }
    }
}

template<int degree>
void Mesh<degree>::ElemTransformation(int idx, Matrix<double> &N, Matrix<double> &dNdxi, Matrix<double> &dNdeta,Vector<double> &w){
    // std::cout << "elem idx :" << idx << "\n";
    (e+idx)->Element<degree>::ElemTransformation(N,dNdxi,dNdeta,w);
}

template<int degree>
void Mesh<degree>::GetElemIEN(int idx, Vector<int> &local_ien){
    (e+idx)->Element<degree>::LocalIEN(local_ien);
}

/// @brief get element details, make sure to send i from 1 to NE and not from 0 to NE-1
/// @tparam degree 
/// @param i 
/// @return 
template<int degree>
Element<degree>* Mesh<degree>::GetElement(int i){
    if (i >=1 && i <= num_elem){
        return e+i-1;
    }
    else{
        std::cerr << "inaccessbile element ID = " << i << "\n";
        return nullptr;
    }
}

#endif