/**
  @file mat3.h
  @authors B.Thouverez <bastien.thouverez -at- etu.univ-lyon1.fr>
 */

/* TODO
 - comparaison de vectors? < > <= >= ?
 */

#ifndef VEC3_H_
#define VEC3_H_

#include <iostream>
#include <cmath>
#include <cassert>

class Vec3;
typedef Vec3 Point;

/**
 * @brief The Vec3 class
 */
class Vec3 {
protected:
    float _x, _y, _z;

public:
    /////// Constructors ///////
    Vec3();
    Vec3(const Point &f, const Point &t);
    Vec3(const float& val);
    Vec3(const float& vx, const float& vy, const float& vz);
    Vec3(const float values[3]);

    /////// Accessors ///////
    float& operator[] (int i);
    float operator[] (int i) const;

    float x() const;
    float y() const;
    float z() const;

    float* data();

    /////// Unary operators ///////
    Vec3 operator+ () const;
    Vec3 operator- () const;

    /////// Assignment unary operators ///////
    Vec3& operator+= (const Vec3 &v);
    Vec3& operator-= (const Vec3 &v);
    Vec3& operator*= (const float &s);
    Vec3& operator/= (const float &s);

    // Binary operators
    // Does it make sense to compare vector on coordinates?
    // Maybe on length ?
    friend int operator> (const Vec3&, const Vec3&);
    friend int operator< (const Vec3&, const Vec3&);
    friend int operator>= (const Vec3&, const Vec3&);
    friend int operator<= (const Vec3&, const Vec3&);

    /////// Binary operators ///////
    friend Vec3 operator+ (const Vec3 &u, const Vec3 &v);
    friend Vec3 operator- (const Vec3 &u, const Vec3 &v);
    friend Vec3 operator* (const Vec3 &v, const float &s);
    friend Vec3 operator* (const float &s, const Vec3 &v);
    friend Vec3 operator/ (const Vec3 &v, const float &s);

    /////// Boolean functions ///////
    friend int operator== (const Vec3 &u, const Vec3 &v);
    friend int operator!= (const Vec3 &u, const Vec3 &v);

    friend std::ostream& operator<< (std::ostream &os, const Vec3 &v);

    ////// Vector operations //////
    friend float dotProduct(const Vec3 &u, const Vec3 &v);
    friend Vec3 crossProduct(const Vec3 &u, const Vec3 &v);
    friend float length(const Vec3 &v);
    float length();
    void normalize();
    friend Vec3 normalized(const Vec3 &v);
    friend float cosine(const Vec3 &u,const Vec3 &v);
    friend bool collinear(const Vec3 &u, const Vec3 &v);
    friend bool coplanar(const Vec3 &u, const Vec3 &v, const Vec3 &w);
    friend float distance(const Vec3 &u, const Vec3 &v);
};


#endif

