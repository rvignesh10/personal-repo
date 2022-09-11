#ifndef MESH_HPP
#define MESH_HPP

#include "MatrixAlgebra.hpp"
#include "element.hpp"

/* only generating maximum 2D mesh */

template<const int nx, const int ny, const int degree>
class Mesh{
private:
    Element<degree> *e;
    POINT *nodes;
    int dim;
    double x1, x2, y1, y2;
    int num_elem;
public:
    void Make1DCartesian(double x_left, double x_right);
    void Make2DCartesian(double x_left, double x_right, double y_bottom, double y_top, Geom g);
    void MakeRectMesh();
    void MakeTriMesh();
    int GetNE(){return num_elem;}
    ~Mesh(){
        delete[] e;
        delete[] nodes;
    }
};

template<const int nx, const int ny, const int degree>
void Mesh<nx,ny,degree>::Make1DCartesian(double x_left, double x_right){
    e = new Element<degree>[nx];
    num_elem = nx;
    x1 = x_left;
    x2 = x_right;
    y1 = 0.; y2 = 0.;
    dim = 1;
    double dx = (x2-x1)/(double)(nx);

    nodes = new POINT[degree*num_elem + 1];
    for (int i=0; i<nx*degree+1; i++){
        nodes[i].setCoordinates(x1+(double)(i)*(dx)/(double)(degree),0.);
        nodes[i].setIdx(i+1);
        if (i==0 || i==nx*degree){
            nodes[i].setAttribute(1);
        }
        else{
            nodes[i].setAttribute(0);
        }
    }
    

    int k = 0;
    for (int i=0; i<nx; i++){
        if (i==0 || i==nx-1){
            e[i].Init(i+1, segment, 1);
        }
        else{
            e[i].Init(i+1, segment, 0);
        }
        for (int j=0; j<=degree; j++){
            e[i].AddVertex(nodes[k], j);
            ++k;
        }
        k-=1;
        
    }
}
template<const int nx, const int ny, const int degree>
void Mesh<nx,ny,degree>::Make2DCartesian(double x_left, 
                           double x_right, double y_bottom, double y_top, Geom g){

    x1 = x_left; x2 = x_right;
    y1 = y_bottom; y2 = y_top;

    dim = 2;

    if (g == triangle){
        MakeTriMesh();
    }
    else{
        MakeRectMesh();
    }

}
template<const int nx, const int ny, const int degree>
void Mesh<nx,ny,degree>::MakeRectMesh(){
    e = new Element<degree>[nx*ny];
    num_elem = nx*ny;
    double dx = (x2-x1)/(double)(nx);
    double dy = (y2-y1)/(double)(ny);

    nodes = new POINT[(nx*degree+1)*(ny*degree+1)];
    int count = 0;
    for (int j=0; j<ny*degree+1; j++){
        for (int i=0; i<nx*degree+1; i++){
            nodes[count].setCoordinates(x1+(double)(i)*dx/(double)(degree),
                                        x2+(double)(j)*dy/(double)(degree));
            nodes[count].setIdx(count+1);
            if (i==0 || i==nx*degree || j==0 || j==ny*degree){
                nodes[count].setAttribute(1);
            }
            else{
                nodes[count].setAttribute(0);
            }
            ++count;
        }
    }

    int el_count = 0;
    for (int j=0; j<ny; j++){
        for (int i=0; i<nx; i++){
            count = 0;
            if (i==0 || i==nx-1 || j==0 || j==ny-1){
                // boundary elements marked 1
                e[el_count].Init(el_count+1, quadrilateral, 1);
            }
            else{
                e[el_count].Init(el_count+1, quadrilateral, 0);
            }
            for (int l=0; l<degree+1; l++){
                for (int m=0; m<degree+1; m++){
                    e[el_count].AddVertex(nodes[(j+l)*(nx*degree+1)+(i+m)],count);
                    ++count;
                }
            }
            ++el_count;
        }
    }
}
template<const int nx, const int ny, const int degree>
void Mesh<nx,ny,degree>::MakeTriMesh(){
    e = new Element<degree>[2*nx*ny];
    num_elem = 2*nx*ny;
    double dx = (x2-x1)/(double)(nx);
    double dy = (y2-y1)/(double)(ny);

    nodes = new POINT[(nx*degree+1)*(ny*degree+1)];
    int count = 0;
    for (int j=0; j<ny*degree+1; j++){
        for (int i=0; i<nx*degree+1; i++){
            nodes[count].setCoordinates(x1+(double)(i)*dx/(double)(degree),
                                        x2+(double)(j)*dy/(double)(degree));
            nodes[count].setIdx(count+1);
            if (i==0 || i==nx*degree || j==0 || j==ny*degree){
                nodes[count].setAttribute(1);
            }
            else{
                nodes[count].setAttribute(0);
            }
            ++count;
        }
    }

    int el_count = 0;
    for (int j=0; j<ny; j++){
        for (int i=0; i<nx; i++){
            // lower triangle
            e[el_count].Init(el_count+1, triangle, 0); 
            count = 0;
            for (int m=0; m<degree+1; m++){
                for (int l=m; l<degree+1; l++){
                    e[el_count].AddVertex(nodes[(j+l)*(nx*degree+1)+(i+m)],count);
                    ++count;
                }
            }
            ++el_count;
            //upper triangle
            e[el_count].Init(el_count+1, triangle, 0); 
            count = 0;
            for(int l=0; l<degree+1; l++){
                for (int m=0; m<l; m++){
                    e[el_count].AddVertex(nodes[(j+l)*(nx*degree+1)+(i+m)],count);
                    ++count;
                }
            }
            ++el_count;
        }
    }
}

#endif