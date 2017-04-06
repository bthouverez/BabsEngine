#ifndef GL_VERSION_H
#define GL_VERSION_H


// Settings for OpenGL version
#define OpenGLMajor 3
#define OpenGLMinor 3

#if OpenGLMajor == 3 && OpenGLMinor == 3
#include <QtGui/QOpenGLFunctions_3_3_Core>
typedef QOpenGLFunctions_3_3_Core OpenGLFunctionCore;
#elif OpenGLMajor == 4 && OpenGLMinor == 3
#include <QtGui/QOpenGLFunctions_4_3_Core>
typedef QOpenGLFunctions_4_3_Core OpenGLFunctionCore;
#endif


#endif
