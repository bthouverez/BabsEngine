#include <iostream>
#include "Utilities/shell_colors.h"

int main(int argc, char ** argv){

	ShellColors sc;
	sc.text_blue();
	std::cout << "############## TEST -- C++14 ##############" << std::endl;
	sc.reset();

    int big = 1'000'000'000;
    std::cout << "Big (C++14) : " << big << std::endl;

    int quatre = 2**2;
    std::cout << "Quatre : " << quatre << std::endl;

	sc.text_green();
    std::cout<< "############## TEST C++14 OK ##############" << std::endl;
    sc.reset();

    return 0;
}
