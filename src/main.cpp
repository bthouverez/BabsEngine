#include <iostream>

#include "Algebra/vec3.h"
#include "Viewer/glwidget.h"


/*
void displayMe(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_TRIANGLES);
        glVertex3f(-0.5, -0.5, 0.0);
        glVertex3f(0.5, -0.5, 0.0);
        glVertex3f(0.0, 0.5, 0.0);
    glEnd();
    glFlush();
}
*/


int main(int argc, char ** argv){
	std::cout<< "Hello Poulpe" << std::endl;
	Vec3 v(1,2, 3.6);
	std::cout<< v << std::endl;

    glutInit(&argc, argv);
    glutCreateWindow("Hello world :D");

    /*
    glutInitDisplayMode(GLUT_SINGLE);
    glutInitWindowSize(300, 300);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Hello world :D");
    glutDisplayFunc(displayMe);
    glutMainLoop();
    */
    GLenum err = glewInit();
    if(err!=GLEW_OK) {
      // Problem: glewInit failed, something is seriously wrong.
      std::cout << "glewInit failed: " << glewGetErrorString(err) << std::endl;
      exit(1);
    }


	std::cout<< "Bye Poulpe" << std::endl;
	return 0;
}
