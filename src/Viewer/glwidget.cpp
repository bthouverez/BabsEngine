#include "glwidget.h"


GLWidget::GLWidget()
{   
}

GLWidget::~GLWidget()
{   
}

void GLWidget::cleanup()
{
    // liberer les ressources gl la.
}

void GLWidget::initializeGL() {

    // Important : pour appeler automatiquement ta fonction de cleanUp lors de la destruction du widget.
    //connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);

	// Important : toute appli opengl doit "charger" les fonctions opengl...
    //initializeOpenGLFunctions();

    _scale = 1.f;
    _scalem = Mat4();

    // shaders
    Shader shader("/home/babs/Documents/e-roma/Hybride/src/Viewer/vertex.glsl", "/home/babs/Documents/e-roma/Hybride/src/Viewer/fragment.glsl", std::cerr);
    shader.init();
    m_program = shader.getProgramID();

    // camera
    m_camera = b_Camera();
    m_camera.lookAt(Vec3(0.0, 0.0, 5.0),  Vec3(0.0, 0.0, 0.0),Vec3(0.0, 1.0, 0.0));
    m_camera.projection(800, 600, 45);

    // data points
    std::vector<Vec3> points;
    points.push_back(Vec3(-0.5, -0.5, 0));
    points.push_back(Vec3(0.5, -0.5, 0));
    points.push_back(Vec3(0.0, 0.5, 0));


    // vao
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // vbo positions
    glGenBuffers(1, &m_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, m_buffer);
    //std::vector<float> points({-0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.0, 0.5, 0.0});
    //glBufferData(GL_ARRAY_BUFFER, points.size()*sizeof(float), &points, GL_STATIC_DRAW);
    glBufferData(GL_ARRAY_BUFFER, points.size()*3*sizeof(float), points.data(), GL_STATIC_DRAW);

    GLuint loc = 0; // inside vertex shader, position has layout 1
    glVertexAttribPointer(loc, 3, GL_FLOAT, GL_FALSE,  0, 0);
    glEnableVertexAttribArray(loc);


    // Appli init:
    // glViewport(0, 0, 800, 600);
    glClearColor(0.9, 0.5, 0.1, 1.f);
    //glClearColor(0.15f, 0.15f, 0.15f, 1.f);
    glClearDepth(1.f);
    glDepthFunc(GL_LESS);
    glDisable(GL_DEPTH_TEST);
    glFrontFace(GL_CCW);
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);

    // cleaning
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    return;
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // glLoadIdentity();
   // glScalef(_scale, _scale, _scale);
   // std::cout << "SCALE " << _scale << std::endl;
  //  b_Printer printer;
    //printer(_scalem);

    Mat4 m = Mat4();
    Mat4 v = Mat4(); //m_camera.view();
    Mat4 p = Mat4(); // m_camera.projection();
    Mat4 mvp = _scalem * p * v * m;

    // configurer le pipeline, selectionner le vertex array a utiliser
    glBindVertexArray(m_vao);
    assert(m_vao!=0);
    // configurer le pipeline, selectionner le shader program a utiliser
    glUseProgram(m_program);
    assert(m_program != 0);

    GLuint loc = glGetUniformLocation(m_program, "mvp");
    glUniformMatrix4fv(loc, 1, GL_TRUE, mvp.data());
    glDrawArrays(GL_TRIANGLES, 0, 3);


    // cleaning
    glUseProgram(0);
    glBindVertexArray(0);

    return;
}
