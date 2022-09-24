#ifndef ADSOLVER_HPP
#define ADSOLVER_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"
#include "Integrator.hpp"

#include <string.h>
#include <fstream>
#include <sstream>




template<int degree>
class ADSolver {
protected:
    Integrator<degree> *DomainInt; 
    int num;   
public :
    ADSolver();
    void AddDomainIntegrator(Integrator<degree> *i = nullptr);
    void readSparse();
    void writeSparse();
    void solve();
};


template<int degree>
ADSolver<degree>::ADSolver(){
    DomainInt = new Integrator<degree>[10];
    num = 0; 
}


template<int degree>
void ADSolver<degree>::AddDomainIntegrator(Integrator<degree> *i){

    if (num < 10){
        *(DomainInt+num) = *i;
        ++num;
        std::cout << "You have added a domain integrator and number of domain integrators present are: " << num << "\n";
    }
    else{
        std::cout << "limit reached for number of integrators to add \n";
    }

}

template<int degree>
void ADSolver<degree>::readSparse(){
    for (int it=0; it < num; it++){
        AppendList* h = DomainInt[it].Integrator<degree>::returnHead();
        std::cout << " --------- " << DomainInt[it].Integrator<degree>::IntegType << " -------------- \n";
        std::cout << "I" << "  " << "J" << "  " << "value \n";
        for (; h!=nullptr; h=h->next){
            std::cout << h->i << "  " << h->j << "  " << h->value << "\n";
        }

    }
}

template<int degree>
void ADSolver<degree>::writeSparse(){

    for (int it=0; it < num; it++){
        AppendList* h = DomainInt[it].Integrator<degree>::returnHead();
        std::stringstream fileName;
        fileName << "../solveFiles/Integrator_" << (it+1) << "_" <<DomainInt[it].Integrator<degree>::IntegType << ".txt" << std::flush;
        std::cout << fileName.str() << "\n";

        std::fstream fileCreate;
        fileCreate.open(fileName.str(), std::fstream::out);

        for (; h!=nullptr; h=h->next){
            int I = h->i; int J = h->j; double Val = h->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate << I << "," << J << "," << s << std::endl;
        }

        fileCreate.close();

    }
}

template<int degree>
void ADSolver<degree>::solve(){
    
}

#endif