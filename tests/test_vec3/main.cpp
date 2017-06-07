#include <iostream>
#include "Utilities/shell_colors.h"
#include "Algebra/vec3.h"

int main(int argc, char ** argv){

	ShellColors sc;
	sc.text_blue();
	std::cout << "############## TEST -- Vec3 ##############" << std::endl;
	sc.reset();

    // Default construction
    Vec3 vector1;
    std::cout << "vector 1 " << vector1 << std::endl;

    // Construction with values
    Vec3 vector2(0.0, 5.0, 0.0);
    std::cout << "vector 2 " << vector2 << std::endl;

    // Copy construction
    Vec3 vector3(vector2);
    std::cout << "vector 3 " << vector3 << std::endl;

    // Construction from 2 points
    Vec3 vector4(vector3, vector2);
    std::cout << "vector 4 " << vector4 << std::endl;

    // Addition
    Vec3 vector5 = vector2 + vector3 + Vec3(-1.0, -1.0, -1.0);
    std::cout << "vector 5 " << vector5 << std::endl;


	sc.text_green();
    std::cout<< "############## TEST Vec3 OK ##############" << std::endl;
    sc.reset();

    return 0;
}