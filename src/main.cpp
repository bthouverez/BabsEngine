#include <iostream>

#include "Algebra/vec3.h"
#include "Viewer/glwidget.h"

#include "Viewer/window.h"

int main(int argc, char ** argv){
    std::cout<< "Hello Poulpe" << std::endl;
    Vec3 v(1,2, 3.6);
    std::cout<< v << std::endl;


    int big = 1'000'000'000;
    std::cout<< "Big (C++14) : " << big << std::endl;

    b_Window window(800, 600);
    window.init();
    window.run();
    window.release();


    std::cout<< "Bye Poulpe" << std::endl;
    return 0;
}
