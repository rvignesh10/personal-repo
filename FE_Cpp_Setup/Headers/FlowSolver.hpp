#ifndef FlowSolver_HPP
#define FlowSolver_HPP

#include "fe_space.hpp"
#include "MatrixAlgebra.hpp"
#include "BiLinearForm.hpp"
#include "LinearForm.hpp"
#include "stabilization.hpp"

#include <string.h>
#include <fstream>
#include <sstream>

enum FlowType { steady=0, unsteady=1};

template<int degree>
class FlowSolver {
protected:
    FlowType fType;
    Residual rType;
    BiLinearForm<degree> *DomainInt; 
    int num;   
    LinearForm<degree> *DomainLF;
    Stabilization<degree> *DomainHatInt;
public :
    FlowSolver();
    FlowSolver(FlowType f, Residual rType);
    void AddDomainIntegrator(BiLinearForm<degree> *i = nullptr);
    void AddLinearForm(LinearForm<degree> *i = nullptr);
    void AddStabilization(Stabilization<degree> *i = nullptr);
    void readSparse();
    void writeSparse();
};


template<int degree>
FlowSolver<degree>::FlowSolver(){
    DomainInt = new BiLinearForm<degree>[10];
    num = 0; 
    DomainLF = nullptr;
    DomainHatInt = nullptr;
    fType = steady;
    rType = noRes;
}

template<int degree>
FlowSolver<degree>::FlowSolver(FlowType f, Residual r){
    DomainInt = new BiLinearForm<degree>[10];
    num = 0; 
    DomainLF = nullptr;
    DomainHatInt = nullptr;
    fType = f;
    rType = r;
}

template<int degree>
void FlowSolver<degree>::AddDomainIntegrator(BiLinearForm<degree> *i){

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
void FlowSolver<degree>::AddLinearForm(LinearForm<degree> *i){
    DomainLF = new LinearForm<degree>; 
    DomainLF = i;
    std::cout << "You have added 1 Linear Form for the RHS \n";
}

template<int degree>
void FlowSolver<degree>::AddStabilization(Stabilization<degree> *i){
    DomainHatInt = new Stabilization<degree>;
    DomainHatInt = i;
    std::cout << "You have added Stabilization \n";
}

template<int degree>
void FlowSolver<degree>::readSparse(){
    // for (int it=0; it < num; it++){
    //     AppendList* h = DomainInt[it].BiLinearForm<degree>::returnHead();
    //     std::cout << " --------- " << DomainInt[it].BiLinearForm<degree>::IntegType << " -------------- \n";
    //     std::cout << "I" << "  " << "J" << "  " << "value \n";
    //     for (; h!=nullptr; h=h->next){
    //         std::cout << h->i << "  " << h->j << "  " << h->value << "\n";
    //     }
    // }
    AppendList1D *h = DomainLF->LinearForm<degree>::returnHead();
    std::cout << " -------- LinearForm ----------- \n";
    std::cout << "I" << "  " << "value \n";
    for(; h!=nullptr; h = h->next){
        std::cout << h->i << "  " << h->value << "\n";
    }
}

template<int degree>
void FlowSolver<degree>::writeSparse(){

    for (int it=0; it < num; it++){
        AppendList* h = DomainInt[it].BiLinearForm<degree>::returnHead();
        AppendList1D *hr = DomainInt[it].BiLinearForm<degree>::returnHead_r();

        std::stringstream fileName, fileName2;
        fileName << "../solveFiles/Integrator_" << (it+1) << "_" <<DomainInt[it].BiLinearForm<degree>::IntegType << ".txt" << std::flush;
        fileName2 << "../solveFiles/Integrator_" << (it+1) << "_" << DomainInt[it].BiLinearForm<degree>::IntegType << "_" << "RHS.txt" << std::flush;
        std::cout << fileName.str() << " and "<< fileName2.str() <<"\n";

        std::fstream fileCreate, fileCreate2;
        fileCreate.open(fileName.str(), std::fstream::out);
        fileCreate2.open(fileName2.str(), std::fstream::out);

        for (; h!=nullptr; h=h->next){
            int I = h->i; int J = h->j; double Val = h->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate << I << "," << J << "," << s << std::endl;
        }

        for(; hr!=nullptr; hr=hr->next){
            int I = hr->i; double Val = hr->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate2 << I << "," << s << std::endl;
        }

        fileCreate.close();
        fileCreate2.close();
    }

    if (fType == steady && DomainLF){
        std::stringstream fileName;
        fileName << "../solveFiles/LinearForm.txt" << std::flush;
        std::cout << fileName.str() << "\n";

        std::fstream fileCreate;
        fileCreate.open(fileName.str(), std::fstream::out);

        AppendList1D *h = this->DomainLF->LinearForm<degree>::returnHead();
        for(; h!=nullptr; h=h->next){
            int I = h->i; double Val = h->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate << I << "," << s << std::endl;
        }
        fileCreate.close();
    }

    if (DomainHatInt){
        AppendList *h = DomainHatInt->Stabilization<degree>::returnHead();
        AppendList1D *hr = DomainHatInt->Stabilization<degree>::returnHead_r();

        std::stringstream fileName1, fileName2;
        fileName1 << "../solveFiles/Stabilization_" << DomainHatInt->sType << ".txt" << std::flush;
        fileName2 << "../solveFiles/Stabilization_" << DomainHatInt->sType << "_RHS.txt" << std::flush;
        std::cout << fileName1.str() << " and "<< fileName2.str() <<"\n";

        std::fstream fileCreate, fileCreate2;
        fileCreate.open(fileName1.str(), std::fstream::out);
        fileCreate2.open(fileName2.str(), std::fstream::out);

        for (; h!=nullptr; h=h->next){
            int I = h->i; int J = h->j; double Val = h->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate << I << "," << J << "," << s << std::endl;
        }

        for(; hr!=nullptr; hr=hr->next){
            int I = hr->i; double Val = hr->value;
            std::ostringstream double2str;
            double2str << std::fixed;
            double2str << std::setprecision(16);
            double2str << Val;
            std::string s = double2str.str();
            fileCreate2 << I << "," << s << std::endl;
        }

        fileCreate.close();
        fileCreate2.close();        
    }

}

#endif