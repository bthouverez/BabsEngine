/**
  @file vec4.h
  @authors  B.Thouverez <bastien.thouverez -at- etu.univ-lyon1.fr>
 */

/* TODO
 - comparaison de vectors? < > <= >= ?
 */

#ifndef VEC4_H_
#define VEC4_H_

#include <iostream>
#include <cmath>
#include <cassert>
#include "vec3.h"



/**
 * @brief The Vec4 class
 */
class Vec4 {
protected:
    float _x, _y, _z, _w;

public:
    /////// Constructors ///////
    Vec4();
    Vec4(const Vec3 &f);
    Vec4(const Vec3 &f, const float &w);
    Vec4(const float& val);
    Vec4(const float& vx, const float& vy, const float& vz, const float& vw);
    Vec4(const float values[4]);

    /////// Accessors ///////
    float& operator[] (int i);
    float operator[] (int i) const;

    float x() const;
    float y() const;
    float z() const;
    float w() const;

    float* data();

    /////// Unary operators ///////
    Vec4 operator+ () const;
    Vec4 operator- () const;

    /////// Assignment unary operators ///////
    Vec4& operator+= (const Vec4 &v);
    Vec4& operator-= (const Vec4 &v);
    Vec4& operator*= (const float &s);
    Vec4& operator/= (const float &s);


    /////// Binary operators ///////
    friend Vec4 operator+ (const Vec4 &u, const Vec4 &v);
    friend Vec4 operator- (const Vec4 &u, const Vec4 &v);
    friend Vec4 operator* (const Vec4 &v, const float &s);
    friend Vec4 operator* (const float &s, const Vec4 &v);
    friend Vec4 operator/ (const Vec4 &v, const float &s);

    /////// Boolean functions ///////
    friend int operator== (const Vec4 &u, const Vec4 &v);
    friend int operator!= (const Vec4 &u, const Vec4 &v);

    friend std::ostream& operator<< (std::ostream &os, const Vec4 &v);

    ////// Vector operations //////
};


#endif

