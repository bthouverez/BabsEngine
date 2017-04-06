#ifndef GL_WIDGET_HEADER
#define GL_WIDGET_HEADER


#include <QOpenGLWidget>
#include <QWheelEvent>
#include "shader.h"
#include "gl_version.h"
#include <GL/gl.h>
#include "../Algebra/vec3.h"
#include "../Algebra/mat4.h"
#include "camera.h"

class GLWidget : public QOpenGLWidget, protected OpenGLFunctionCore
{
	Q_OBJECT

public:
	GLWidget(QWidget *parent = 0)
		: QOpenGLWidget(parent)		
	{
		setFocusPolicy(Qt::StrongFocus);//for key events
	}
	~GLWidget();

    //QSize minimumSizeHint() const Q_DECL_OVERRIDE;
    //QSize sizeHint() const Q_DECL_OVERRIDE;


protected:
	void initializeGL() Q_DECL_OVERRIDE;
    void paintGL() Q_DECL_OVERRIDE;


    void wheelEvent(QWheelEvent *event);
    //void resizeGL(int width, int height) Q_DECL_OVERRIDE;

	void cleanup();

    GLuint m_vao;
    GLuint m_buffer;
    GLuint m_instance_buffer;
    GLuint m_program;

    b_Camera m_camera;

    float _scale;
    QMatrix4x4 _scalem;

    /*
	void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void wheelEvent(QWheelEvent * event) Q_DECL_OVERRIDE;
	void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
	void keyReleaseEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    */

};

#endif //GL_WIDGET_HEADER
