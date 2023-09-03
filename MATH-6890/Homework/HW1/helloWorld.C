#include <stdio.h>
#include <iostream>

int main(int argc, char *argv[]){
    printf("Hello World! \n");

    std::cout << "argument 1 is: " << argv[0] << std::endl;
    if (argc > 1){
        for (int i=1; i < argc; i++){
            std::cout << "argument " << i+1 << " is: " << argv[i] << std::endl;
        }
    }

    return 0;
}