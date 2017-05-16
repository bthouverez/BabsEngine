#ifndef SHELL_COLORS_H_
#define SHELL_COLORS_H_

#include <iostream>

class ShellColors {

private:
    std::ostream stream;
    bool b_bold = false; // 1
    bool b_underline = false; // 4
    bool b_blink = false; // 5
    bool b_highlight = false; // 7
    bool b_background = false;
    int color_text;
    int color_background;

public:
    static const int TEXT_COLOR_BLACK = 30;
    static const int TEXT_COLOR_RED = 31;
    static const int TEXT_COLOR_GREEN = 32;
    static const int TEXT_COLOR_BROWN = 33;
    static const int TEXT_COLOR_BLUE = 34;
    static const int TEXT_COLOR_MAGENTA = 35;
    static const int TEXT_COLOR_CYAN = 36;
    static const int TEXT_COLOR_LIGHT_GRAY = 37;

    static const int BACKGROUND_COLOR_BLACK = 40;
    static const int BACKGROUND_COLOR_RED = 41;
    static const int BACKGROUND_COLOR_GREEN = 42;
    static const int BACKGROUND_COLOR_YELLOW = 43;
    static const int BACKGROUND_COLOR_BLUE = 44;
    static const int BACKGROUND_COLOR_PURPLE = 45;
    static const int BACKGROUND_COLOR_CYAN = 46;
    static const int BACKGROUND_COLOR_WHITE = 47;


    /////// Constructors ///////
    ShellColors() : stream(std::cout.rdbuf()), color_text(TEXT_COLOR_BLACK) {

    }

    ShellColors(std::ostream os) : stream(os.rdbuf()), color_text(TEXT_COLOR_BLACK) {

    }

    /////// Text color ///////
    // Unbolded
    void text_black() {
        color_text = TEXT_COLOR_BLACK;
        b_bold = false;
        update();
    }

    void text_red() {
        color_text = TEXT_COLOR_RED;
        b_bold = false;
        update();
    }

    void text_green() {
        color_text = TEXT_COLOR_GREEN;
        b_bold = false;
        update();
    }

    void text_brown() {
        color_text = TEXT_COLOR_BROWN;
        b_bold = true;
        update();
    }

    void text_blue() {
        color_text = TEXT_COLOR_BLUE;
        b_bold = false;
        update();
    }

    void text_magenta() {
        color_text = TEXT_COLOR_MAGENTA;
        b_bold = false;
        update();
    }

    void text_cyan() {
        color_text = TEXT_COLOR_CYAN;
        b_bold = false;
        update();
    }

    void text_light_gray() {
        color_text = TEXT_COLOR_LIGHT_GRAY;
        b_bold = false;
        update();
    }

    // Bolded
    void text_dark_gray() {
        color_text = TEXT_COLOR_BLACK;
        b_bold = true;
        update();
    }

    void text_dark_red() {
        color_text = TEXT_COLOR_RED;
        b_bold = true;
        update();
    }

    void text_dark_green() {
        color_text = TEXT_COLOR_GREEN;
        b_bold = true;
        update();
    }

    void text_yellow() {
        color_text = TEXT_COLOR_BROWN;
        b_bold = true;
        update();
    }

    void text_dark_blue() {
        color_text = TEXT_COLOR_BLUE;
        b_bold = true;
        update();
    }

    void text_dark_magenta() {
        color_text = TEXT_COLOR_MAGENTA;
        b_bold = true;
        update();
    }

    void text_dark_cyan() {
        color_text = TEXT_COLOR_CYAN;
        b_bold = true;
        update();
    }

    void text_white() {
        color_text = TEXT_COLOR_LIGHT_GRAY;
        b_bold = true;
        update();
    }


    /////// background color ///////
    void background_black() {
        color_background = BACKGROUND_COLOR_BLACK;
        b_background = true;
        update();
    }

    void background_red() {
        color_background = BACKGROUND_COLOR_RED;
        b_background = true;
        update();
    }

    void background_green() {
        color_background = BACKGROUND_COLOR_GREEN;
        b_background = true;
        update();
    }

    void background_yellow() {
        color_background = BACKGROUND_COLOR_YELLOW;
        b_background = true;
        update();
    }

    void background_blue() {
        color_background = BACKGROUND_COLOR_BLUE;
        b_background = true;
        update();
    }

    void background_purple() {
        color_background = BACKGROUND_COLOR_PURPLE;
        b_background = true;
        update();
    }

    void background_cyan() {
        color_background = BACKGROUND_COLOR_CYAN;
        b_background = true;
        update();
    }

    void background_white() {
        color_background = BACKGROUND_COLOR_WHITE;
        b_background = true;
        update();
    }

    /////// Text style ///////
    void underline() {
        b_underline = true;
        update();
    }

    void blink() {
        b_blink = true;
        update();
    }

    void highlight() {
        b_highlight = true;
        update();
    }

    void reset() {
        b_bold = false;
        b_background = false;
        b_blink = false;
        b_highlight = false;
        b_underline = false;
        stream << "\e[0m";
    }

    void update() {
       stream << "\e[";
       if(b_bold) stream << "1;";
       stream << color_text;
       if(b_background) stream << ";" << color_background;
       if(b_underline) stream << ";4";
       if(b_blink) stream << ";5";
       if(b_highlight) stream << ";7";
       stream<< "m";
    }
};

#endif
