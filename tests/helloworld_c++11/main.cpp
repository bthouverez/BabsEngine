#include <iostream>
#include "Utilities/shell_colors.h"

class LaClasse {
    int nombre;

public:
    LaClasse(int n) : nombre(n) {}
    LaClasse() : LaClasse(42) {} // Calling another constructor, c++11 feature
};



int main(int argc, char ** argv){

	ShellColors sc;
	sc.text_blue();
	std::cout << "############## TEST -- C++11 ##############" << std::endl;
	sc.reset();

   	LaClasse lc;

	sc.text_green();
    std::cout<< "############## TEST C++11 OK ##############" << std::endl;
    sc.reset();

    return 0;
}