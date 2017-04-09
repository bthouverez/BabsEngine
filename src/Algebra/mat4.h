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
  float _values[4][4];

public:
  Mat4();
  Mat4(const float &s);
  Mat4(const Vec3 &v);
  Mat4(const float &m00, const float &m01, const float &m02, const float &m03, const float &m10, const float &m11, const float &m12, const float &m13, const float &m20, const float &m21, const float &m22, const float &m23, const float &m30, const float &m31, const float &m32, const float &m33);

  static const Mat4 Zero;
  static const Mat4 Identity;

  float& operator() (const int &i,const int &j);
  const float& operator() (const int &i,const int &j) const;

  float* data();

  // Unary operators
  Mat4 operator+ () const;
  Mat4 operator- () const;

  // Assignment operators
  Mat4& operator+= (const Mat4 &m);
  Mat4& operator-= (const Mat4 &m);
  Mat4& operator*= (const Mat4 &m);
  Mat4& operator*= (const float &s);
  Mat4& operator/= (const float &s);

  // Binary operators
  friend Mat4 operator+ (const Mat4 &m1, const Mat4 &m2);
  friend Mat4 operator- (const Mat4 &m1, const Mat4 &m2);
  friend Mat4 operator* (const Mat4 &m1, const Mat4 &m2);

  friend std::ostream& operator<< (std::ostream &os, const Mat4 &m);

  float determinant() const;
  friend float determinant(const Mat4 &m);
  float trace() const;
  friend float trace(const Mat4 &m);
  Mat4 inverse() const;

  void lookAt(const Vec3 &from, const Vec3 &to, const Vec3 &up);
  void perspective(const float &fov, const float &aspect, const float &znear, const float &zfar);

};

float radians(const float deg);
float degrees(const float rad);
Mat4 Scale(const float &x, const float &y, const float &z);
Mat4 Scale(const Vec3 &v);
Mat4 Scale(const float &s);
Mat4 Translation(const Vec3 &v);
Mat4 Translation(const float &x, const float &y, const float &z);
Mat4 RotationX(const float &a);
Mat4 RotationY(const float &a);
Mat4 RotationZ(const float &a);
Mat4 Rotation(const Vec3 &ax, const float &an);

Mat4 LookAt(const Vec3 &from, const Vec3 &to, const Vec3 &up);
Mat4 Perspective(const float &fov, const float &aspect, const float &znear,	const float &zfar);


#endif
