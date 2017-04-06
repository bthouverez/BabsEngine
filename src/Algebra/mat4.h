/**
  @file mat4.h
  @authors ??, refactored by B.Thouverez <bastien.thouverez -at- etu.univ-lyon1.fr>
 */

// Is this class really necessary? Is not a QMatrix4x4 enough?
#ifndef MAT4_H_
#define MAT4_H_

#include "vec3.h"
#include "mat3.h"

/**
 * @brief The Mat4 class
 */
class Mat4 {
protected:
  double values[4][4];

public:
  Mat4();
  Mat4(const double &s);
  Mat4(const Vec3 &v);
  Mat4(const double &m00, const double &m01, const double &m02, const double &m03, const double &m10, const double &m11, const double &m12, const double &m13, const double &m20, const double &m21, const double &m22, const double &m23, const double &m30, const double &m31, const double &m32, const double &m33);

  static const Mat4 Zero;
  static const Mat4 Identity;

  double& operator() (const int &i,const int &j);
  const double& operator() (const int &i,const int &j) const;

  // Unary operators
  Mat4 operator+ () const;
  Mat4 operator- () const;

  // Assignment operators
  Mat4& operator+= (const Mat4 &m);
  Mat4& operator-= (const Mat4 &m);
  Mat4& operator*= (const Mat4 &m);
  Mat4& operator*= (const double &s);
  Mat4& operator/= (const double &s);

  // Binary operators
  friend Mat4 operator+ (const Mat4 &m1, const Mat4 &m2);
  friend Mat4 operator- (const Mat4 &m1, const Mat4 &m2);
  friend Mat4 operator* (const Mat4 &m1, const Mat4 &m2);

  friend std::ostream& operator<< (std::ostream &os, const Mat4 &m);

  double determinant() const;
  friend double determinant(const Mat4 &m);
  double trace() const;
  friend double trace(const Mat4 &m);
  Mat4 inverse() const;

};

float radians(const float deg);
float degrees(const float rad);
Mat4 Scale(const double &x, const double &y, const double &z);
Mat4 Scale(const Vec3 &v);
Mat4 Scale(const double &s);
Mat4 Translation(const Vec3 &v);
Mat4 Translation(const double &x, const double &y, const double &z);
Mat4 RotationX(const double &a);
Mat4 RotationY(const double &a);
Mat4 RotationZ(const double &a);
Mat4 Rotation(const Vec3 &ax, const double &an);

Mat4 LookAt(const Vec3 &from, const Vec3 &to, const Vec3 &up);
Mat4 Perspective(const float &fov, const float &aspect, const float &znear,	const float &zfar);


#endif
