#include "glwidget.h"
#include <GL/glut.h>

#include <iostream>

GLWidget::GLWidget(QWidget *parent) : QOpenGLWidget(parent)
{

}

GLWidget::~GLWidget()
{

}


void GLWidget::initializeGL() {
    glClearColor(0.2, 0.2, 0.2, 1);
    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_LIGHT0);
    //glEnable(GL_LIGHTING);
}

void GLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0.0, 5.0,  0.0, 0.0, 0.0,  0.0, 1.0, 0.0 );

    glBegin(GL_TRIANGLES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(-0.5, -0.5, 0.0);
        glColor3f(1.0, 1.0, 0.0);
        glVertex3f(0.5, -0.5, 0);
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.5, 0.0);
    glEnd();
}

void GLWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (float)w/h, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, 5.0,  0.0, 0.0, 0.0,  0.0, 1.0, 0.0 );
}

void GLWidget::rotateX(const double &a){
    glRotatef(a, 1.0, 0.0, 0.0);
    update();
}

void GLWidget::rotateY(const double &a){
    glRotatef(a, 0.0, 1.0, 0.0);
   // update();
}

void GLWidget::rotateZ(const double &a){
    glRotatef(a, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    update();
    glMatrixMode(GL_MODELVIEW);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
      float dx = (event->x() - lastPos.x()) / 10.0f;
      float dy = (event->y() - lastPos.y())/ 10.0f;
      float dz = (event->y() - lastPos.y())/ 10.0f;

     if (event->buttons() & Qt::LeftButton)
     {
       glRotatef(dy*0.1, 1.0f, 0.0f, 0.0f);
       glRotatef(dx*0.1, 0.0f, 1.0f, 0.0f);
       glRotatef(dz*0.1, 0.0f, 0.0f, 1.0f);
     }

    lastPos = event->pos();
    update();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
    update();
}

