#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "geometries.hpp"

template<const int degree>
class Element{
private:
    int elem_id;
    Geom geometry;
    int attribute; // 1 - boundary element, 0 - nonboundary element
protected:
    POINT *p;
public:
    Element(){
        elem_id = 0;
        attribute = 0;
        geometry = invalid;
        p = nullptr;
    }
    Element(int id, Geom g, int att);
    void Init(int id, Geom g, int att);
    void setAttribute(int att){attribute = att;}
    void AddVertex(POINT node, int i);
    void AddVertex(double x, double y, int i);
    ~Element(){
        delete[] p;
    }
};

template<const int degree>
Element<degree>::Element(int id, Geom g, int att){
    Init(id, g, att);
}
template<const int degree>
void Element<degree>::Init(int id, Geom g, int att){
    elem_id = id;
    attribute = att;
    geometry = g;
    switch (geometry)
    {
        case point : {
            p = new POINT;
            break;
        }
        case segment : {
            p = new POINT[degree+1];
            break;
        }
        case triangle : {
            if (degree == 1){
                p = new POINT[3];
            }
            else if (degree == 2){
                p = new POINT[6];
            }
            else if (degree == 3){
                p = new POINT[10];
            }
            else {
                std::cerr << "p>3 not supported for triangular elements \n";
            }
            break;
        }
        case quadrilateral : {
            p = new POINT[(degree+1)*(degree+1)];
            break;
        }
        default :{
            p = nullptr;
            std::cerr << "Element type not supported \n";
            break;
        }
    }
}

template<const int degree>
void Element<degree>::AddVertex(POINT node, int i){
    switch(geometry)
    {
        case point : {
            if (i==0){
                *p = node;
            }
            else{
                std::cerr << "invalid node addition to geometry \n";
            }
            break;
        }
        case segment : {
            if (i < degree+1){
                p[i] = node;
            }
            else{
                std::cerr << "invalid node addition to geometry \n";
            }
            break;
        }
        case triangle : {
            if (degree == 1){
                if (i<3){
                    p[i] = node;
                }
                else{
                    std::cerr << "invalid node addition to geometry \n";
                }
            }
            else if(degree == 2){
                if (i<6){
                    p[i] = node;
                }
                else{
                    std::cerr << "invalid node addition to geometry \n";
                }
            }
            else if (degree == 3){
                if (i<10){
                    p[i] = node;
                }
                else{
                    std::cerr << "invalid node addition to geometry \n";
                }
            }
            break;
        }
        case quadrilateral : {
            if (i < (degree+1)*(degree+1))
            {
                p[i] = node;
            }
            else{
                std::cerr << "invalid node addition to geometry \n";
            }
            break;
        }
        default : {
            std::cerr << "invalid node addition to geometry \n";
            break;
        }
    }
}

#endif