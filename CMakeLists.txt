#Déclaration du projet

project(BabsEngine)



#Déclaration de l'exécutable

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
	./src/Viewer/glwidget_babs.cpp
	./src/Viewer/glwidget_babs.h
	./src/Viewer/glwidget.cpp
	./src/Viewer/glwidget.h
	./src/Viewer/mainwindow.cpp
	./src/Viewer/mainwindow.h
	./src/Viewer/shader.cpp
	./src/Viewer/shader.h

)

