cmake_minimum_required(VERSION 2.8)
project(RenderBoy)

######## C++ 14 ##################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
# eventually
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
###################################################################################

add_executable(

    RenderBoy
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

        # TO REMOVE
	./src/Viewer/window.cpp
	./src/Viewer/window.h	

)

find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
        message(STATUS "Found OpenGL")
        include_directories( ${OPENGL_INCLUDE_DIRS} )
        target_link_libraries(RenderBoy ${OPENGL_LIBRARIES} )
else()
        message(ERROR "OpenGL not found")
endif()


find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    message(STATUS "Found GLEW")
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
    target_link_libraries(RenderBoy ${GLEW_LIBRARIES} )
else()
    message(ERROR "GLEW not found")
endif()



#find_package(GLUT REQUIRED)
#if(OPENGL_GLU_FOUND)
#       message(STATUS "Found GLUT")
#	include_directories( ${GLUT_INCLUDE_DIRS} )
#    link_libraries(${GLUT_LIBRARY})
#	target_link_libraries(RenderBoy ${GLUT_LIBRARY} )
#else()
#	message(ERROR "GLUT not found")
#endif()



find_package(SDL2 REQUIRED)
if(SDL2_FOUND)
	message(STATUS "Found SDL2")
        include_directories(${SDL2_INCLUDE_DIRS})
    link_libraries(${SDL2_LIBRARIES})
        target_link_libraries(RenderBoy ${SDL2_LIBRARIES})
else()
	message(ERROR "SDL2 not found")
endif()
