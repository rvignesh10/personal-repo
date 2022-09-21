#include "../Headers/MatrixAlgebra.hpp"

int main(){
    Matrix<double> S(2.,3);
    Matrix<double> M1,M2;
    M1.setSize(3,2);
    double k = 1.;
    for(int i=0; i<3; i++){
        for (int j=0; j<2; j++){
            M1.setValue(i,j,k);
            ++k;
        }
    }
    M2.setSize(2,3);
    M1.Transpose(M2);
    M2.displayMatrix();
    
    return 0;
}