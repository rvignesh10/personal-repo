#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


enum Geom {invalid=-1, point=0, segment=1, triangle=2, quadrilateral=3};

/* assume all problems are at most 2D */

/* defining a point element with (x,y) values */
class POINT {
private:
    int idx, attribute;
    double x1,x2;
public:
    // POINT(){x1 = 0.; x2 = 0.;}
    POINT(int id, double x = 0., double y = 0., int att = 0){
        idx = id; x1 = x; x2 = y; attribute = att;
    }
    POINT(int id, double x = 0., double y = 0.){
        idx = id; x1 = x; x2 = y; attribute = 0;
    }
    POINT(double x = 0., double y = 0.){
        x1 = x; x2 = y;
        idx = 0;
    }
    void setCoordinates(double x, double y){
        x1 = x; x2 = y;
    }
    void getCoordinates(double &x, double &y){
        x = x1; y = x2;
    }
    void setIdx(int id){
        idx = id;
    }
    int getIdx() {return idx;}
    void setAttribute(int att){
        attribute = att;
    }
    int getAttribute(){return attribute;}
    
};

// class SEGMENT {
// private:
//     POINT p[2];
//     int bdr_seg;
// public:
//     void AddPoints(Point a, Point b){
//         p[0] = a; p[1] = b;
//     }
// };


#endif