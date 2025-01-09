#define PY_SSIZE_T_CLEAN

#include "../Headers/H1_FiniteElements.hpp"
#include "../Headers/fe_space.hpp"
#include "../Headers/mesh.hpp"
#include "../Headers/Integrator.hpp"
#include "../Headers/ADSolver.hpp"
#include </Library/Frameworks/Python.framework/Versions/3.8/include/python3.8/Python.h>

int main(int argc, char const *argv[]){
    Mesh<2> mesh;
    mesh.Make1DCartesian(50,0,-2.,2.);
    H1_FiniteElementSpace<2> fespace(&mesh);
    ADSolver<2> solver;
    solver.AddDomainIntegrator(new DiffusionIntegrator<2>(.1, &fespace));
    solver.writeSparse();

    wchar_t *progam = Py_DecodeLocale(argv[0], NULL); 

    if (progam == NULL){
        std::cerr << "Fatal error .. Cannot Decode argv[0] \n";
        exit(1);
    }

    Py_SetProgramName(progam); 
    Py_Initialize();

    PyRun_SimpleString("exec(open('solve.py').read())");

    if (Py_FinalizeEx() < 0.){
        std::cerr << "ummfh \n";
        exit(120);
    }

    PyMem_RawFree(progam);

    return 0;
}