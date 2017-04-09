cmake_minimum_required(VERSION 2.8)
project(BabsEngine)

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

	./src/Viewer/window.cpp
	./src/Viewer/window.h	

)

find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
	include_directories( ${OPENGL_INCLUDE_DIRS} )
	target_link_libraries(BabsEngine ${OPENGL_LIBRARIES} )
else()
	message(ERROR "OpenGL not found")
endif()

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
    target_link_libraries(BabsEngine ${GLEW_LIBRARIES} )
else()
    message(ERROR "Glew not found")
endif()



find_package(GLUT REQUIRED)
if(OPENGL_GLU_FOUND)
	include_directories( ${GLUT_INCLUDE_DIRS} )
    link_libraries(${GLUT_LIBRARY})
	target_link_libraries(BabsEngine ${GLUT_LIBRARY} )
else()
	message(ERROR "GLUT not found")
endif()



find_package(SDL2 REQUIRED)
if(SDL2_FOUND)
	message(STATUS "Found SDL2")
	include_directories(${SDL2_INCLUDE_DIR})
    link_libraries(${SDL2_LIBRARY})
	target_link_libraries(BabsEngine ${SDL2_LIBRARY})
else()
	message(ERROR "SDL2 not found")
endif()