cmake_minimum_required(VERSION 2.8)
project(BabsEngine)

#add_library(BabsEngine ./src/Algebra/vec3.h)

add_executable(

    BabsEngine
	./src/main.cpp
	./src/Algebra/mat3.cpp
	./src/Algebra/mat3.h
	./src/Algebra/mat4.cpp
	./src/Algebra/mat4.h
	./src/Algebra/vec3.cpp
	./src/Algebra/vec3.h
	./src/Algebra/vec4.cpp
	./src/Algebra/vec4.h
	./src/Viewer/camera.cpp
	./src/Viewer/camera.h
	./src/Viewer/gl_version.h
	./src/Viewer/glwidget.cpp
	./src/Viewer/glwidget.h
	./src/Viewer/shader.cpp
	./src/Viewer/shader.h

)


#find_package(OpenGL REQUIRED)

#if(OPENGL_FOUND)
	#message(STATUS "OpenGL found")
#	include_directories( ${OPENGL_INCLUDE_DIRS} )
#	target_link_libraries(BabsEngine ${OPENGL_LIBRARIES} )
#else()
#	message(ERROR "OpenGL not found")
#endif()


#find_package(GLUT REQUIRED)

#if(OPENGL_GLU_FOUND)
#	#message(STATUS "GLU found")
#	include_directories( ${GLUT_INCLUDE_DIRS} )
#	target_link_libraries(BabsEngine ${GLUT_LIBRARY} )
#else()
#	message(ERROR "GLU not found")
#endif()

