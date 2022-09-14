#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "geometries.hpp"

template<const int degree>
class Element{
protected:
    int elem_id;
    int attribute; // 1 - boundary element, 0 - nonboundary element
public:
    Geom geometry;
    POINT *p;
    QuadratureNode *q;
    int sizeof_p, sizeof_q;
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
    void setQuadrature(int i, double x, double y, double weight);
    void printElementNodes();
    ~Element(){
        delete[] p;
        delete[] q;
    }
};

/// @brief used to initialize elements with no id and attribute
/// @tparam degree 
/// @param g - geometry of this element
template<const int degree>
Element<degree>::Element(Geom g){
    Init_(g);
}

/// @brief Element class default constructor
/// @tparam degree 
/// @param id - Assign this number as Element ID
/// @param g - Geometry of the element (point, segment, triangle or quadrilateral)
/// @param att - Attribute of the element (boundary or domain element)
template<const int degree>
Element<degree>::Element(int id, Geom g, int att){
    Init(id, g, att);
}

/// @brief A user-defined initialization of each element class
/// @tparam degree 
/// @param id - Assign this number as Element ID
/// @param g - Geometry of the element (point, segment, triangle or quadrilateral)
/// @param att - Attribute of the element (boundary or domain element)
template<const int degree>
void Element<degree>::Init(int id, Geom g, int att){
    elem_id = id;
    attribute = att;
    Init_(g);
}

/// @brief initialization without id and attribute
/// @tparam degree 
/// @param g - required geometry of the element
template<const int degree>
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
                q = new QuadratureNode;
                sizeof_p = 3;
                sizeof_q = 1;
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
            sizeof_q = 4;
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
template<const int degree>
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

/// @brief This function sets the quadrature location and its weight for each element
/// @tparam degree 
/// @param i - the ith index of the quadrature point array
/// @param x - the x-component of quadrature point
/// @param y - the y-component of quadrature point
/// @param weight - integration weight associated with that location
template<const int degree>
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

template<const int degree>
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