#include "window.h"
#include <iostream>

/**
 * @brief b_Window::b_Window
 */
b_Window::b_Window() : _fullscreen(true), _width(800), _height(600) {
}

/**
 * @brief b_Window::b_Window
 * @param w
 * @param l
 */
b_Window::b_Window(const int &w, const int &h) : _width(w), _height(h) {

}

/**
 * @brief b_Window::init
 * @return
 */
bool b_Window::init() {
    if(SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL: Initialization failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    // cf gl_version.h
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, OpenGLMajor);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, OpenGLMinor);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

    _window = SDL_CreateWindow(
        "Renderboy",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        _width,
        _height,
        SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

    if(_window == nullptr){
        std::cerr <<  "SDL: Window creation failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    if(_fullscreen) SDL_SetWindowFullscreen(_window, SDL_WINDOW_FULLSCREEN);


    _context = SDL_GL_CreateContext(_window);
    if(_context == 0) {
        std::cerr << "ERROR: openGL context: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(_window);
        SDL_Quit();
       return -1;
    }

   // _widget.initializeGL();

}

/**
 * @brief b_Window::events
 */
void b_Window::events() {
    SDL_Event event;

    while(SDL_PollEvent( &event ) != 0) {

        //User requests quit
        if( event.type == SDL_QUIT ){
            std::cout << "Red cross quit" << std::endl;
            _quit = true;

        } else if(event.type == SDL_KEYDOWN) {

            switch(event.key.keysym.sym) {
                case SDLK_UP:
                    std::cout<<"SDLK_UP gedrückt!"<<std::endl;
                break;

                case SDLK_DOWN:
                    std::cout<<"SDLK_DOWN gedrückt!"<<std::endl;
                break;

                case SDLK_LEFT:
                    std::cout<<"SDLK_LEFT gedrückt!"<<std::endl;
                break;

                case SDLK_RIGHT:
                    std::cout<<"SDLK_RIGHT gedrückt!"<<std::endl;
                break;

                case SDLK_q:
                    std::cout<<"SDLK_q gedrückt!"<<std::endl;
                    _quit = true;
                break;

                case SDLK_f:
                    std::cout<<"SDLK_f gedrückt!"<<std::endl;
                    if(_fullscreen) {
                        SDL_SetWindowFullscreen(_window, 0);
                        _fullscreen =  false;
                    } else {
                        SDL_SetWindowFullscreen(_window, SDL_WINDOW_FULLSCREEN);
                        _fullscreen = true;
                    }
                break;

                case SDLK_ESCAPE:
                    std::cout << "Escape quit" << std::endl;
                    _quit = true;
                break;

                default:
                    std::cout<<"Another key pressed"<<std::endl;
                break;
            }
        }
    }

}


float vertices[] = { -0.5, -0.5, 0.5, -0.5, 0.0, 0.5};

/**
 * @brief b_Window::run
 * @return
 */
bool b_Window::run() {
    while(!_quit){
        events();

        glClearColor(0.9, 0.5, 0.1, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);
        //_widget.paintGL();
        /*
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, vertices);
        glEnableVertexAttribArray(0);

        glDrawArrays(GL_TRIANGLES, 0, 3);
        glDisableVertexAttribArray(0);
        */
        SDL_GL_SwapWindow(_window);



    }
}

/**
 * @brief b_Window::release
 * @return
 */
void b_Window::release() {
    SDL_GL_DeleteContext(_context);
    SDL_DestroyWindow(_window);
    SDL_Quit();
}
