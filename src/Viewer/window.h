#ifndef WINDOW_H_
#define WINDOW_H_

#include "gl_version.h"
#include "glwidget.h"
#include <SDL2/SDL.h>

class b_Window {
public:
    b_Window();
    b_Window(const int &w, const int &h);

    bool init();
    void events();
    bool run();
    void release();
private:


    SDL_Window* _window = nullptr;
    SDL_GLContext _context;

    //GLWidget _widget;

    int _width;
    int _height;
    bool _quit = false;
    bool _fullscreen = false;
};

#endif
