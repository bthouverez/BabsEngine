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
    double _x, _y, _z, _w;

public:
    /////// Constructors ///////
    Vec4();
    Vec4(const Vec3 &f);
    Vec4(const Vec3 &f, const double &w);
    Vec4(const double& val);
    Vec4(const double& vx, const double& vy, const double& vz, const double& vw);
    Vec4(const double values[4]);

    /////// Accessors ///////
    double& operator[] (int i);
    double operator[] (int i) const;

    double x() const;
    double y() const;
    double z() const;
    double w() const;

    /////// Unary operators ///////
    Vec4 operator+ () const;
    Vec4 operator- () const;

    /////// Assignment unary operators ///////
    Vec4& operator+= (const Vec4 &v);
    Vec4& operator-= (const Vec4 &v);
    Vec4& operator*= (const double &s);
    Vec4& operator/= (const double &s);


    /////// Binary operators ///////
    friend Vec4 operator+ (const Vec4 &u, const Vec4 &v);
    friend Vec4 operator- (const Vec4 &u, const Vec4 &v);
    friend Vec4 operator* (const Vec4 &v, const double &s);
    friend Vec4 operator* (const double &s, const Vec4 &v);
    friend Vec4 operator/ (const Vec4 &v, const double &s);

    /////// Boolean functions ///////
    friend int operator== (const Vec4 &u, const Vec4 &v);
    friend int operator!= (const Vec4 &u, const Vec4 &v);

    friend std::ostream& operator<< (std::ostream &os, const Vec4 &v);

    ////// Vector operations //////
};


#endif

