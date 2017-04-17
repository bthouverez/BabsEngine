#ifndef GL_WIDGET_HEADER
#define GL_WIDGET_HEADER

#include <vector>

#include "shader.h"
#include "../Algebra/vec3.h"
#include "../Algebra/mat4.h"
#include "camera.h"


typedef Yann_Shader Shader;

class GLWidget {

public:
	GLWidget();
	~GLWidget();

    //QSize minimumSizeHint() const Q_DECL_OVERRIDE;
    //QSize sizeHint() const Q_DECL_OVERRIDE;


	void initializeGL();
    void paintGL();



	void cleanup();

protected:
    GLuint m_vao;
    GLuint m_buffer;
    GLuint m_instance_buffer;
    GLuint m_program;

    b_Camera m_camera;

    float _scale;
    Mat4 _scalem;

    /*
	void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void wheelEvent(QWheelEvent * event) Q_DECL_OVERRIDE;
	void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
	void keyReleaseEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    */

};

#endif //GL_WIDGET_HEADER
