#include <iostream>
#include "Utilities/shell_colors.h"

int main(int argc, char ** argv){

	ShellColors sc;
	sc.text_blue();
	std::cout << "############## TEST -- C++14 ##############" << std::endl;
	sc.reset();

    int big = 1'000'000'000; // ' to make int readable, C++14 feature
    std::cout << "Big (C++14) : " << big << std::endl;

	sc.text_green();
    std::cout<< "############## TEST C++14 OK ##############" << std::endl;
    sc.reset();

    return 0;
}
