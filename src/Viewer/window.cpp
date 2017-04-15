#include "window.h"
#include <iostream>

/**
 * @brief b_Window::b_Window
 */
b_Window::b_Window() {

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
        "Babs Engine",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        _width,
        _height,
        SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);

    if(_window == nullptr){
        std::cerr <<  "SDL: Window creation failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    _context = SDL_GL_CreateContext(_window);

    if(_context == 0) {
        std::cout << SDL_GetError() << std::endl;
        SDL_DestroyWindow(_window);
        SDL_Quit();
       return -1;
    }

    _renderer = SDL_CreateRenderer(_window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
}

/**
 * @brief b_Window::events
 */
void b_Window::events() {
    SDL_Event event;

    while(SDL_PollEvent( &event ) != 0) {

        //User requests quit
        if( event.type == SDL_QUIT ){
            _quit = true;

        } else if(event.type == SDL_KEYDOWN){

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

                case SDLK_ESCAPE:
                    _quit = true;
                break;

                default:
                    std::cout<<"Irgendeine Taste wurde gedrückt"<<std::endl;
                break;
            }
        }
    }

}

/**
 * @brief b_Window::run
 * @return
 */
bool b_Window::run() {
    while(!_quit){
        events();
        SDL_SetRenderDrawColor(_renderer, 255, 255, 255, 255);
        if (SDL_RenderClear(_renderer) != 0) {
            std::cout << "Fehler: " << SDL_GetError() << std::endl;
        }
        SDL_SetRenderDrawColor(_renderer, 255, 0, 0, 255);
        SDL_RenderPresent(_renderer);
    }
}
