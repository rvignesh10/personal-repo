#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "MatrixAlgebra.hpp"
#include "geometries.hpp"

template<int degree>
class Element{
protected:
    int elem_id;
    int attribute; // 1 - boundary element, 0 - nonboundary element
public:
    Geom geometry;
    POINT *p;
    QuadratureNode *q;
    int sizeof_p, sizeof_q;
    Vector<int> node_idx;
    Element(){
        elem_id = 0;
        attribute = 0;
        geometry = invalid;
        p = nullptr;
        q = nullptr;
        sizeof_p = 0;
        sizeof_q = 0;
    }
    Element(Geom g); 
    Element(int id, Geom g, int att);
    void Init(int id, Geom g, int att);
    void Init_(Geom g);
    void setAttribute(int att){attribute = att;}
    void AddVertex(POINT *node, int i);
    void getVertexLoc(int dir, Vector<double> &loc);
    void setQuadrature(int i, double x, double y, double weight);
    void getQuadrature(int dir, Vector<double> &q_dir);
    void ElemTransformation(Matrix<double> &N, Matrix<double> &dNdxi, Matrix<double> &dNdeta, Vector<double> &w);
    void LocalIEN(Vector<int> &local_ien);
    void printElementNodes();
    ~Element(){
        delete[] p;
        delete[] q;
    }
};

/// @brief used to initialize elements with no id and attribute
/// @tparam degree 
/// @param g - geometry of this element
template<int degree>
Element<degree>::Element(Geom g){
    Init_(g);
}

/// @brief Element class default constructor
/// @tparam degree 
/// @param id - Assign this number as Element ID
/// @param g - Geometry of the element (point, segment, triangle or quadrilateral)
/// @param att - Attribute of the element (boundary or domain element)
template<int degree>
Element<degree>::Element(int id, Geom g, int att){
    Init(id, g, att);
}

/// @brief A user-defined initialization of each element class
/// @tparam degree 
/// @param id - Assign this number as Element ID
/// @param g - Geometry of the element (point, segment, triangle or quadrilateral)
/// @param att - Attribute of the element (boundary or domain element)
template<int degree>
void Element<degree>::Init(int id, Geom g, int att){
    elem_id = id;
    attribute = att;
    Init_(g);
}

/// @brief initialization without id and attribute
/// @tparam degree 
/// @param g - required geometry of the element
template<int degree>
void Element<degree>::Init_(Geom g){
    geometry = g;
    switch (geometry)
    {
        case point : {
            p = new POINT;
            q = new QuadratureNode;
            sizeof_p = 1;
            sizeof_q = 1;
            break;
        }
        case segment : {
            p = new POINT[degree+1];
            q = new QuadratureNode[5]; // integrates upto polynomial order of 9 accurately
            sizeof_p = degree + 1;
            sizeof_q = 5;
            break;
        }
        case triangle : {
            if (degree == 1){
                p = new POINT[3];
                q = new QuadratureNode[3];
                sizeof_p = 3;
                sizeof_q = 3;
            }
            else if (degree == 2){
                p = new POINT[6];
                q = new QuadratureNode[3];
                sizeof_p = 6;
                sizeof_q = 3;
            }
            else if (degree == 3){
                p = new POINT[10];
                q = new QuadratureNode[3];
                sizeof_p = 10;
                sizeof_q = 3;
            }
            else {
                std::cerr << "degree>3 not supported for triangular elements \n";
            }
            break;
        }
        case quadrilateral : {
            p = new POINT[(degree+1)*(degree+1)];
            q = new QuadratureNode[4];
            sizeof_p = (degree+1)*(degree+1);
            sizeof_q = 4; // integrates upto polynomial order of 7 accurately
            break;
        }
        default :{
            p = nullptr;
            q = nullptr;
            std::cerr << "Element type not supported \n";
            break;
        }
    }
}

/// @brief Adds the vertices associated with this element
/// @tparam degree 
/// @param node - node location is a point with idx, attribute and (x,y) physical location
/// @param i - the index of the p array
template<int degree>
void Element<degree>::AddVertex(POINT *node, int i){
    if (i < sizeof_p){
        (p+i)->setIdx(node->getIdx());
        (p+i)->setAttribute(node->getAttribute());
        double x,y;
        node->getCoordinates(x,y);
        (p+i)->setCoordinates(x,y);
    }
    else{
        std::cerr << "invalid node addition to geometry \n";
    }
}

/// @brief get the location of the node locations in this element
/// @tparam degree 
/// @param dir direction of the location needed (x or y)
/// @param loc Vector to assign these values to
template<int degree>
void Element<degree>::getVertexLoc(int dir, Vector<double> &loc){
    loc.setSize(sizeof_p);
    for (int i=0; i<sizeof_p; i++){
        double x,y;
        (p+i)->getCoordinates(x,y);
        if (dir==1){
            loc.setValue(i,x);
        }
        else{
            loc.setValue(i,y);
        }
    }
}

/// @brief This function sets the quadrature location and its weight for each element
/// @tparam degree 
/// @param i - the ith index of the quadrature point array
/// @param x - the x-component of quadrature point
/// @param y - the y-component of quadrature point
/// @param weight - integration weight associated with that location
template<int degree>
void Element<degree>::setQuadrature(int i, double x, double y, double weight){
    if (i < sizeof_q){
        (q+i)->int_x = x;
        (q+i)->int_y = y;
        (q+i)->int_wt= weight;
    } 
    else{
        std::cerr << "invalid Quadrature point access \n";
    }   
}

/// @brief  get the quadrature location (x or y) associated with this element
/// @tparam degree 
/// @param dir direction of quadrature location needed
/// @param q_dir Vector to assign those values to
template<int degree>
void Element<degree>::getQuadrature(int dir, Vector<double> &q_dir){
    q_dir.setSize(sizeof_q);
    for(int i=0; i<sizeof_q; i++){
        if (dir ==1){
            q_dir.setValue(i, (q+i)->int_x);
        }
        else if(dir==2){
            q_dir.setValue(i, (q+i)->int_y);
        } 
        else if (dir==3){
            q_dir.setValue(i, (q+i)->int_wt);
        }
        else if(dir ==4){
            q_dir.setValue(i, (q+i)->det_jacobian);
        } 
    }
}


template<int degree>
void Element<degree>::ElemTransformation(Matrix<double> &N, Matrix<double> &dNdxi, Matrix<double> &dNdeta, Vector<double> &w){
    Vector<double> x; Vector<double> y;
    this->getVertexLoc(1, x);
    this->getVertexLoc(2, y);
    Vector<double> N_qi(this->sizeof_p);
    Vector<double> dNdxi_qi(this->sizeof_p);
    Vector<double> dNdeta_qi(this->sizeof_p);

    for(int i=0; i<this->sizeof_q; i++){
        
        for(int j=0; j<this->sizeof_p; j++){
            N_qi.setValue(j, N.getValue(j,i));
            dNdxi_qi.setValue(j,dNdxi.getValue(j,i));
            dNdeta_qi.setValue(j,dNdeta.getValue(j,i));
        }
        double xq, yq;
        xq = N_qi.dotProduct(x);
        yq = N_qi.dotProduct(y);
        double dxdxi = dNdxi_qi.dotProduct(x);
        double dydxi = dNdxi_qi.dotProduct(y);
        double dxdeta = dNdeta_qi.dotProduct(x);
        double dydeta = dNdeta_qi.dotProduct(y);
        (q+i)->J.setSize(2,2);
        (q+i)->J.setValue(0,0,dxdxi);
        (q+i)->J.setValue(0,1,dxdeta);
        (q+i)->J.setValue(1,0,dydxi);
        (q+i)->J.setValue(1,1,dydeta);
        if (this->geometry == segment){
            if(abs(dxdxi)<1e-4){
                (q+i)->det_jacobian = 1e-4;
            } 
            else{
                (q+i)->det_jacobian = dxdxi;
            }     
        }
        else{
            if(abs(dxdxi*dydeta - dxdeta*dydxi) < 1e-4){
                (q+i)->det_jacobian = 1e-4;
            }
            else{
                (q+i)->det_jacobian = dxdxi*dydeta - dxdeta*dydxi;
            }
        }
        (q+i)->Jinv.setSize(2,2);
        (q+i)->Jinv.setValue(0,0,(1./(q+i)->det_jacobian)*dydeta);
        (q+i)->Jinv.setValue(0,1,(-1./(q+i)->det_jacobian)*dxdeta);
        (q+i)->Jinv.setValue(1,0,(-1./(q+i)->det_jacobian)*dydxi);
        (q+i)->Jinv.setValue(1,1,(1./(q+i)->det_jacobian)*dxdxi);
        this->setQuadrature(i, xq, yq, w.getValue(i));
    }

}

template<int degree>
void Element<degree>::LocalIEN(Vector<int> &local_ien){
    local_ien.setSize(sizeof_p);
    node_idx.setSize(sizeof_p);
    for(int i=0; i<sizeof_p; i++){
        node_idx.setValue(i, (p+i)->getIdx());
        if((p+i)->getAttribute() == 0){
            local_ien.setValue(i, (p+i)->getIdx());
        }
        else{
            local_ien.setValue(i, -1);
        }
        
    }
}

template<int degree>
void Element<degree>::printElementNodes(){
    std::cout << "ID" << "  " << "x" << "  " << "y \n"; 
    for (int i=0; i<sizeof_p; i++){
        double x,y;
        (p+i)->getCoordinates(x,y);
        int n = (p+i)->getIdx();
        std::cout << n << "  " << x << "  " << y << "\n";
    }
}
#endif