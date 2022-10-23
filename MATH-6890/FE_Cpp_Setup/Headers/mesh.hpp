#ifndef MESH_HPP
#define MESH_HPP

#include "MatrixAlgebra.hpp"
#include "element.hpp"
#include <fstream>
#include <sstream>

/* only generating maximum 2D mesh */
/* only generating rectangular meshes */

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
    void Init(int dimension, Geom g, int n_elem, int n1, int n2);
    void Make1DCartesian(int n1, int n2, double x_left, double x_right);
    void Make2DCartesian(int n1, int n2,double x_left, double x_right, double y_bottom, double y_top, Geom g);
    void MakeRectMesh();
    void MakeTriMesh();
    void setBdrCondition_onElement(int e_idx, Vector<int> l_idx, Vector<double> &l_bdr, void (*func)(double, double, double &));
    void invoke(double x, double y, double &qbc, void (*func)(double, double, double &)){func(x,y,qbc);}
    int GetNNodes(){return num_nodes;}
    int GetNE(){return num_elem;}
    int GetDim(){return dim;}
    Geom GetGeometry() {return geometry;}
    void ElemTransformation(int idx,Matrix<double> N, Matrix<double> dNdxi, Matrix<double> dNdeta, Vector<double> w);
    void GetElemIEN(int idx, Vector<int> &local_ien);
    Element<degree>* GetElement(int i);
    void AddElement(int e_id, POINT *p, int num_pts, int e_attribute);
    void FinalizeMesh();
    ~Mesh(){
        delete[] e;
        delete[] nodes;
    }
};

template<int degree>
void Mesh<degree>::Init(int dimension, Geom g, int n_elem, int n1, int n2){
    dim = dimension;
    geometry = g;
    num_elem = n_elem;
    nx = n1;
    ny = n2;

    if (g == segment){
        if (nx != n_elem){
            std::cerr << "incompatible mesh formation \n";
        }
        else {
            num_nodes = degree*num_elem+1;
            e = new Element<degree>[num_elem];
            nodes = new POINT[num_nodes];
        }
    }
    else if(g == quadrilateral){
        if (nx*ny != n_elem){
            std::cerr << "incompatible mesh formation \n";
        }
        else{
            num_nodes = (nx*degree+1)*(ny*degree+1);
            e = new Element<degree>[num_elem];
            nodes = new POINT[num_nodes];
        }
    }
    else if(g == triangle){
        if (2*nx*ny != n_elem){
            std::cerr << "incompatible mesh formation \n";
        }
        else{
            num_nodes = (nx*degree+1)*(ny*degree+1);
            e = new Element<degree>[num_elem];
            nodes = new POINT[num_nodes];
        }
    }
}

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
    this->FinalizeMesh();
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
                if(i==0){
                    (nodes+count)->setAttribute(1);
                }
                else if(i==nx*degree){
                    (nodes+count)->setAttribute(2);
                }
                else if(j==0){
                    (nodes+count)->setAttribute(3);
                }
                else{
                    (nodes+count)->setAttribute(4);
                }
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
                if(i==0){
                    (nodes+count)->setAttribute(1);
                }
                else if(i==nx*degree){
                    (nodes+count)->setAttribute(2);
                }
                else if(j==0){
                    (nodes+count)->setAttribute(3);
                }
                else{
                    (nodes+count)->setAttribute(4);
                }
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
            ++el_count;
        }
    }
}

template<int degree>
void Mesh<degree>::ElemTransformation(int idx, Matrix<double> N, Matrix<double> dNdxi, Matrix<double> dNdeta,Vector<double> w){
    // std::cout << "elem idx :" << idx << "\n";
    (e+idx)->Element<degree>::ElemTransformation(N,dNdxi,dNdeta,w);
}

template<int degree>
void Mesh<degree>::GetElemIEN(int idx, Vector<int> &local_ien){
    (e+idx)->Element<degree>::LocalIEN(local_ien);
}

template<int degree>
void Mesh<degree>::setBdrCondition_onElement(int e_idx, Vector<int> l_idx, Vector<double> &l_bdr,void (*func)(double, double, double &)){
    Vector<double> x((e+e_idx)->Element<degree>::sizeof_p);
    Vector<double> y((e+e_idx)->Element<degree>::sizeof_p);
    (e+e_idx)->Element<degree>::getVertexLoc(1, x);
    (e+e_idx)->Element<degree>::getVertexLoc(2, y);
    for(int i=0; i<(e+e_idx)->Element<degree>::sizeof_p; i++){
        if(l_idx.getValue(i) == -1){
            double qbc;
            Mesh<degree>::invoke(x.getValue(i), y.getValue(i), qbc, func);
            l_bdr.setValue(i, qbc);
        }
        else{
            l_bdr.setValue(i, 0.);
        }
    }
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

// set element id from 1 to num_elem. 
template<int degree>
void Mesh<degree>::AddElement(int e_id, POINT *pt, int num_pt, int e_attribute){
    if (e_id < 0 || e_id >= num_elem){
        std::cerr << "wrong element idx trying to be accessed \n";
    }
    else{
        (e+e_id)->Element<degree>::Init(e_id, geometry, e_attribute);
        if (geometry == segment){
            if (num_pt!=2){
                std::cerr << "wrong number of points specified \n";
            }
            else{
                double xp[2];
                for (int i=0; i<num_pt; i++){
                    double xc, yc;
                    (pt+i)->getCoordinates(xc,yc);
                    xp[i] = xc; 
                }
                double h = abs(xp[1]-xp[0])/degree;
                double xmin, xmax;
                if (h <= 1e-16){
                    std::cerr << "mesh is too small \n";
                }
                else{
                    if(xp[0]<xp[1]){
                        xmin = xp[0];
                        xmax = xp[1];

                    }
                    else{
                        xmin = xp[1];
                        xmax = xp[0];
                    }
                }
                for(int i=0; i<= degree; i++){
                    (nodes+degree*(e_id)+i)->setCoordinates(xmin+(double)(i)*h,0.);
                    if(degree*(e_id)+i==0 || degree*(e_id)+i==num_nodes){
                        (nodes+degree*(e_id)+i)->setAttribute(1);
                    }
                    else{
                        (nodes+degree*(e_id)+i)->setAttribute(0);
                    }
                    (nodes+degree*(e_id)+i)->setIdx(degree*(e_id)+i+1);
                    (e+e_id)->Element<degree>::AddVertex(nodes+degree*(e_id)+i, i);
                }
            }
        }
        else if (geometry == quadrilateral){
            if(num_pt!=4){
                std::cerr<< "wrong number of points specified \n";
            }
            else {
                double xp[4][2];
            }

        }
        else if (geometry == triangle){

        }
    }
}

template<int degree>
void Mesh<degree>::FinalizeMesh(){
    if (geometry==segment){
        std::stringstream fileName;
        fileName << "../solveFiles/Mesh.txt" << std::flush;
        std::cout << "Mesh File created and stored ... \n";

        std::fstream fileCreate;
        fileCreate.open(fileName.str(), std::fstream::out);
        for(int i=0; i<num_nodes; i++){
            double x,y;
            (nodes+i)->getCoordinates(x,y);
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << x;
            std::string s = double2str.str();
            fileCreate << s << std::endl;
        }
        fileCreate.close();
    }
}

#endif