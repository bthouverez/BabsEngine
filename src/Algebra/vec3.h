/**
  @file mat3.h
  @authors ??, refactored by B.Thouverez <bastien.thouverez -at- etu.univ-lyon1.fr>
 */

/* TODO
 - comparaison de vectors? < > <= >= ?
 */

#ifndef VEC3_H_
#define VEC3_H_

#include <iostream>
#include <cmath>
#include <cassert>



/**
 * @brief The Vec3 class
 */
class Vec3 {
protected:
    double _x, _y, _z;

public:
    /////// Constructors ///////
    Vec3();
    Vec3(const Vec3 &f, const Vec3 &t);
    Vec3(const double& val);
    Vec3(const double& vx, const double& vy, const double& vz);
    Vec3(const double values[3]);

    /////// Accessors ///////
    double& operator[] (int i);
    double operator[] (int i) const;

    double x() const;
    double y() const;
    double z() const;

    /////// Unary operators ///////
    Vec3 operator+ () const;
    Vec3 operator- () const;

    /////// Assignment unary operators ///////
    Vec3& operator+= (const Vec3 &v);
    Vec3& operator-= (const Vec3 &v);
    Vec3& operator*= (const double &s);
    Vec3& operator/= (const double &s);

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
    friend Vec3 operator* (const Vec3 &v, const double &s);
    friend Vec3 operator* (const double &s, const Vec3 &v);
    friend Vec3 operator/ (const Vec3 &v, const double &s);

    /////// Boolean functions ///////
    friend int operator== (const Vec3 &u, const Vec3 &v);
    friend int operator!= (const Vec3 &u, const Vec3 &v);

    friend std::ostream& operator<< (std::ostream &os, const Vec3 &v);

    ////// Vector operations //////
    friend double dotProduct (const Vec3 &u, const Vec3 &v);
    friend Vec3 crossProduct (const Vec3 &u, const Vec3 &v);
    friend double length(const Vec3 &v);
    double length();
    void normalize();
    friend Vec3 normalized(const Vec3 &v);
    friend double cosine(const Vec3 &u,const Vec3 &v);
    friend bool collinear(const Vec3 &u, const Vec3 &v);
    friend bool coplanar(const Vec3 &u, const Vec3 &v, const Vec3 &w);
    friend double distance(const Vec3 &u, const Vec3 &v);
};


#endif

