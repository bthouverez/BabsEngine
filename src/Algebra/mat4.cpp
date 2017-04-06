#include "mat4.h"

/*!
  \class Mat4 matrix.h
  \brief Mat4 implements a 4x4 double matrix. Components are stored in a single
  dimension array, starting from element a<SUB>00</SUB> and sorting components
  by column.
*/

//! Empty matrix, full of 0
const Mat4 Mat4::Zero(0.0);

//! Identity matrix
const Mat4 Mat4::Identity(1.0);


/**
 * @brief Construct empty matrix
 * @details all values are set to 0
 */
Mat4::Mat4() {
    values[0][0] = values[0][1] = values[0][2] = values[0][3] = values[1][0] = values[1][1] = values[1][2] = values[1][3] = values[2][0] = values[2][1] = values[2][2] = values[2][3] = values[3][0] = values[3][1] = values[3][2] = values[3][3] = 0.0;
}

/**
 * @brief Construct matrix with given value
 * @details Given value is used for third first diagonal points, last point of diagonal set to 1
 * @param s value of diagonal
 */
Mat4::Mat4(const double &s) {
    values[0][0] = values[1][1] = values[2][2] = s;
    values[3][3] = 1.0;
    values[0][1] =
            values[0][2] =
            values[0][3] =
            values[1][0] =
            values[1][2] =
            values[1][3] =
            values[2][0] =
            values[2][1] =
            values[2][3] =
            values[3][0] =
            values[3][1] =
            values[3][2] = 0.0;
}

/**
 * @brief Construct matrix with given vector
 * @details Given vector is used for third first diagonal points, last point of diagonal set to 1
 * @param v values of diagonal
 */
Mat4::Mat4(const Vec3 &v) {
    Mat4(1.0);
    values[0][0] = v.x();
    values[1][1] = v.y();
    values[2][2] = v.z();
}

/**
 * @brief Construct matrix with given values
 * @param m00 value at row 0 and column 0
 * @param m01 value at row 0 and column 1
 * @param m02 value at row 0 and column 2
 * @param m03 value at row 0 and column 3
 * @param m10 value at row 1 and column 0
 * @param m11 value at row 1 and column 1
 * @param m12 value at row 1 and column 2
 * @param m13 value at row 1 and column 3
 * @param m20 value at row 2 and column 0
 * @param m21 value at row 2 and column 1
 * @param m22 value at row 2 and column 2
 * @param m23 value at row 2 and column 3
 * @param m30 value at row 3 and column 0
 * @param m31 value at row 3 and column 1
 * @param m32 value at row 3 and column 2
 * @param m33 value at row 3 and column 3
 */
Mat4::Mat4(const double &m00, const double &m01, const double &m02, const double &m03, const double &m10, const double &m11, const double &m12, const double &m13, const double &m20, const double &m21, const double &m22, const double &m23, const double &m30, const double &m31, const double &m32, const double &m33) {
    values[0][0] = m00;
    values[0][1] = m01;
    values[0][2] = m02;
    values[0][3] = m03;
    values[1][0] = m10;
    values[1][1] = m11;
    values[1][2] = m12;
    values[1][3] = m13;
    values[2][0] = m20;
    values[2][1] = m21;
    values[2][2] = m22;
    values[2][3] = m23;
    values[3][0] = m30;
    values[3][1] = m31;
    values[3][2] = m32;
    values[3][3] = m33;
}

/**
 * @brief Write access to matrix values
 * @param i row in matrix
 * @param j column in matrix
 * @return value in matrix at i, j
 */
double& Mat4::operator() (const int &i,const int &j) {
    return values[i][j];
}

/**
 * @brief Read access to matrix values
 * @param i row in matrix
 * @param j column in matrix
 * @return value in matrix at i, j
 */
const double& Mat4::operator() (const int &i,const int &j) const {
    return values[i][j];
}

/**
 * @brief Get positive value of current matrix
 * @return positive value of current matrix
 */
Mat4 Mat4::operator+ () const {
    return *this;
}

/**
 * @brief Get negative value of current matrix
 * @return negative value of current matrix
 */
Mat4 Mat4::operator- () const {
    Mat4 m;
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            m(ii, jj) = -values[ii][jj];
    return m;
}

/**
 * @brief Add a matrix to current one
 * @param m matrix to add
 * @return addition of curernt matrix and m
 */
Mat4& Mat4::operator+= (const Mat4 &m) {
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            (*this)(ii, jj) = values[ii][jj] + m(ii, jj);
    return *this;
}

/**
 * @brief Subtract a matrix to current one
 * @param m matrix to subtract with
 * @return subtraction of current matrix and m
 */
Mat4& Mat4::operator-= (const Mat4 &m) {
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            (*this)(ii, jj) = values[ii][jj] - m(ii, jj);
    return *this;
}

/**
 * @brief Multiply a matrix to current one
 * @param m matrix to multiply with
 * @return multiplication of current matrix with m
 */
Mat4& Mat4::operator*= (const Mat4 &m) {
    Mat4 r(*this);
    *this = r * m;
    return *this;
}

/**
 * @brief Multiply current matrix with scalar
 * @param s scalar to multiply with
 * @return multiplication of current matrix with s
 */
Mat4& Mat4::operator*= (const double &s) {
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            (*this)(ii, jj) *= s;
    return *this;
}

/**
 * @brief Divide current matrix with scalar
 * @param s scalar to divide with
 * @return division of current matrix with s
 */
Mat4& Mat4::operator/= (const double &s) {
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            (*this)(ii, jj) /= s;
    return *this;
}

/**
 * @brief Add two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return addition of m1 + m2
 */
Mat4 operator+ (const Mat4 &m1, const Mat4 &m2) {
    Mat4 r;
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            r = m1(ii, jj) + m2(ii, jj);
    return r;
}

/**
 * @brief Subtract two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return subtraction of m1 - m2
 */
Mat4 operator- (const Mat4 &m1, const Mat4 &m2) {
    Mat4 r;
    for(int ii=0; ii<4; ii++)
        for(int jj=0; ii<4; jj++)
            r = m1(ii, jj) - m2(ii, jj);
    return r;
}

/**
 * @brief Multiply two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return multiplication of m1 * m2
 */
Mat4 operator* (const Mat4 &m1, const Mat4 &m2) {
    Mat4 m;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            m.values[i][j]= m1(i, 0) * m2(0, j) + m1(i, 1) * m2(1, j) + m1(i, 2) * m2(2, j) + m1(i, 3) * m2(3, j);
    return m;
}

/**
 * @brief display matrix
 * @param os stream to print in
 * @param m matrix to display
 * @return stream
 */
std::ostream& operator<< (std::ostream &os, const Mat4 &m) {
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            os << m(i,j) << " ";
        }
        os << std::endl;
    }
    return os;
}


/**
 * @brief Calculate determinant of current matrix
 * @details Decomposition of current Mat4 in four Mat3 and compute determinant for each Mat3
 * @return determinant of current matrix
 */
double Mat4::determinant() const {
    Mat3 m1(values[1][1], values[1][2], values[1][3], values[2][1], values[2][2], values[2][3], values[3][1], values[3][2], values[3][3]);
    Mat3 m2(values[1][0], values[1][2], values[1][3], values[2][0], values[2][2], values[2][3], values[3][0], values[3][2], values[3][3]);
    Mat3 m3(values[1][0], values[1][1], values[1][3], values[2][0], values[2][1], values[2][3], values[3][0], values[3][1], values[3][3]);
    Mat3 m4(values[1][0], values[1][1], values[1][2], values[2][0], values[2][1], values[2][2], values[3][0], values[3][1], values[3][2]);
    return values[0][0]*m1.determinant() - values[0][1]*m2.determinant() + values[0][2]*m3.determinant() - values[0][3]*m4.determinant();
}

/**
 * @brief Calculate determinant of given matrix
 * @param m matrix to calculate determinant on
 * @return determiannt of m
 */
double determinant(const Mat4 &m) {
    return m.determinant();
}

/**
 * @brief Calculate trace of current matrix
 * @details Trace of matrix is the sum of its diagonal elements
 * @return trace of current matrix
 */
double Mat4::trace() const {
    return values[0][0]+values[1][1]+values[2][2]+values[3][3];
}

/**
 * @brief Calculate trace of given matrix
 * @param m matrix to calculate trace on
 * @return trace of m
 */
double trace(const Mat4 &m) {
    return m.trace();
}

/**
 * @brief compute inverse of current matrix
 * @details based on gKit2light framework, thank you J.C.Iehl
 * @return current matrix inverted
 */
Mat4 Mat4::inverse() const {
  Mat4 minv= *this;

  int indxc[4], indxr[4];
  int ipiv[4] = { 0, 0, 0, 0 };

  for (int i = 0; i < 4; i++) {
      int irow = -1, icol = -1;
      float big = 0.f;

      // Choose pivot
      for (int j = 0; j < 4; j++) {
          if (ipiv[j] != 1) {
              for (int k = 0; k < 4; k++) {
                  if (ipiv[k] == 0) {
                      if (fabsf(minv(j, k)) >= big) {
                          big = std::abs(minv(j, k));
                          irow = j;
                          icol = k;
                      }
                  }
                  else if (ipiv[k] > 1)
                      printf("singular matrix in make_inverse()\n");
              }
          }
      }

      assert(irow >= 0 && irow < 4);
      assert(icol >= 0 && icol < 4);

      ++ipiv[icol];
      // Swap rows _irow_ and _icol_ for pivot
      if (irow != icol) {
          for (int k = 0; k < 4; ++k)
              std::swap(minv(irow, k), minv(icol, k));
      }

      indxr[i] = irow;
      indxc[i] = icol;
      if (minv(icol, icol) == 0.)
          printf("singular matrix in make_inverse()\n");

      // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
      float pivinv = 1.f / minv(icol, icol);
      minv(icol, icol) = 1.f;
      for (int j = 0; j < 4; j++)
          minv(icol, j) *= pivinv;

      // Subtract this row from others to zero out their columns
      for (int j = 0; j < 4; j++) {
          if (j != icol) {
              float save = minv(j, icol);
              minv(j, icol) = 0;
              for (int k = 0; k < 4; k++)
                  minv(j, k) -= minv(icol, k)*save;
          }
      }
  }

  // Swap columns to reflect permutation
  for (int j = 3; j >= 0; j--) {
      if (indxr[j] != indxc[j]) {
          for (int k = 0; k < 4; k++)
              std::swap(minv(k, indxr[j]), minv(k, indxc[j]));
      }
  }

  return minv;
}



////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Conversion from degree to radian
 * @param deg angle in degree to convert
 * @return deg in radian
 */
float radians(const float deg) {
    const float pi = 3.1415926535f;
    return (pi  / 180.f) * deg;
}

/**
 * @brief Conversion from radian to degree
 * @param rad angle in radian to convert
 * @return rad in degree
 */
float degrees(const float rad) {
    const float pi = 3.1415926535f;
    return (180.f / pi) * rad;
}


/**
 * @brief Create a scaling matrix
 * @details \f$\displaystyle S = \begin{bmatrix}
 * x & 0 & 0 & 0 \\
 * 0 & y & 0 & 0\\
 * 0 & 0 & z & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param x scale on x axis
 * @param y scale on y axis
 * @param z scale on z axis
 * @return scaling matrix
 */
Mat4 Scale(const double &x, const double &y, const double &z) {
    return Mat4(Vec3(x, y, z));
}

/**
 * @brief Create a scaling matrix
 * @details \f$\displaystyle S = \begin{bmatrix}
 * v_x & 0 & 0 & 0 \\
 * 0 & v_y & 0 & 0\\
 * 0 & 0 & v_z & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param v vector of scaling
 * @return scaling matrix
 */
Mat4 Scale(const Vec3 &v) {
    return Mat4(v);
}

/**
 * @brief Create a scaling matrix
 * @details \f$\displaystyle S = \begin{bmatrix}
 * s & 0 & 0 & 0 \\
 * 0 & s & 0 & 0\\
 * 0 & 0 & s & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param s scaling rate
 * @return scaling matrix
 */
Mat4 Scale(const double &s) {
   return  Mat4(s);
}

/**
 * @brief Create a translation matrix
 * @details \f$\displaystyle T = \begin{bmatrix}
 * 1 & 0 & 0 & v_x \\
 * 0 & 1 & 0 & v_y\\
 * 0 & 0 & 1 & v_z\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param v vector of translation
 * @return translation matrix for v
 */
Mat4 Translation(const Vec3 &v) {
    Mat4 r(Mat4::Identity);
    r(0, 3) = v.x();
    r(1, 3) = v.y();
    r(2, 3) = v.z();
    return r;
}

/**
 * @brief Create a translation matrix
 * @details \f$\displaystyle T = \begin{bmatrix}
 * 1 & 0 & 0 & x \\
 * 0 & 1 & 0 & y\\
 * 0 & 0 & 1 & z\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param x translation along X axis
 * @param y translation along Y axis
 * @param z translation along Z axis
 * @return translation matrix for vector (x, y, z)
 */
Mat4 Translation(const double &x, const double &y, const double &z) {
    return Translation(Vec3(x, y, z));
}

/**
 * @brief Create a rotation matrix around X axis
 * @details \f$\displaystyle R_x = \begin{bmatrix}
 * 1 & 0 & 0 & 0 \\
 * 0 & cos(a) & -sin(a) & 0\\
 * 0 & sin(a) & cos(a) & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param a angle or rotation in degree
 * @return rotation matrix around X axis
 */
Mat4 RotationX(const double &a) {
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));
    return Mat4(
        1,     0,      0, 0,
        0, cos_t, -sin_t, 0,
        0, sin_t,  cos_t, 0,
        0,     0,      0, 1 );

}

/**
 * @brief Create a rotation matrix around y axis
 * @details \f$\displaystyle R_y = \begin{bmatrix}
 * cos(a) & 0 & sin(a) & 0 \\
 * 0 & 1 & 0 & 0\\
 * -sin(a) & 0 & cos(a) & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param a angle or rotation in degree
 * @return rotation matrix around y axis
 */
Mat4 RotationY(const double &a) {
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));
    return Mat4(
         cos_t, 0, sin_t, 0,
             0, 1,     0, 0,
        -sin_t, 0, cos_t, 0,
             0, 0,     0, 1 );
}

/**
 * @brief Create a rotation matrix around Z axis
 * @details \f$\displaystyle R_z = \begin{bmatrix}
 * cos(a) & -sin(a) & 0 & 0 \\
 * sin(a) & cos(a) & 0 & 0\\
 * 0 & 0 & 1 & 0\\
 * 0 & 0 & 0 & 1
 * \end{bmatrix}\f$
 * @param a angle or rotation in degree
 * @return rotation matrix around Z axis
 */
Mat4 RotationZ(const double &a) {
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));
    return Mat4(
        cos_t, -sin_t, 0, 0,
        sin_t,  cos_t, 0, 0,
            0,      0, 1, 0,
            0,      0, 0, 1 );
}

/**
 * @brief Create a rotation matrix around given axis
 * @param ax axis of rotation
 * @param an angle of rotation in degree
 * @return rotation matrix around ax axis of an degree
 */
Mat4 Rotation(const Vec3 &ax, const double &an) {
    Vec3 a= normalized(ax);
    float s= sinf(radians(an));
    float c= cosf(radians(an));

    return Mat4(
        a.x() * a.x() + (1 - a.x() * a.x() ) * c,
        a.x() * a.y() * (1 - c ) - a.z() * s,
        a.x() * a.z() * (1 - c ) + a.y() * s,
        0,

        a.x() * a.y() * (1 - c ) + a.z() * s,
        a.y() * a.y() + (1 - a.y() * a.y() ) * c,
        a.y() * a.z() * (1 - c ) - a.x() * s,
        0,

        a.x() * a.z() * (1 - c ) - a.y() * s,
        a.y() * a.z() * (1 - c ) + a.x() * s,
        a.z() * a.z() + (1 - a.z() * a.z() ) * c,
        0,

        0, 0, 0, 1);
}


/**
 * @brief Create a Lookat matrix
 * @details based on gKit2light framework, thank you J.C.Iehl
 * @param from point of view
 * @param to point to look
 * @param up vector representing up direction
 * @return matrix representing camera placement from from, looking point to with up direction
 */
Mat4 LookAt(const Vec3 &from, const Vec3 &to, const Vec3 &up) {
    Vec3 dir   = normalized( Vec3(from, to) );
    Vec3 right = normalized( crossProduct(dir, normalized(up)) );
    Vec3 newUp = normalized( crossProduct(right, dir) );

    Mat4 m(
        right.x(), newUp.x(), -dir.x(), from.x(),
        right.y(), newUp.y(), -dir.y(), from.y(),
        right.z(), newUp.z(), -dir.z(), from.z(),
        0,       0,        0,     1);

    return m.inverse();
}

/**
 * @brief Create a perspective matrix
 * @details based on gKit2light framework, thank you J.C.Iehl
 * @param fov opening angle
 * @param aspect image size rapport
 * @param znear near depth limit
 * @param zfar far depth limit
 * @return
 */
Mat4 Perspective(const float &fov, const float &aspect, const float &znear,	const float &zfar) {
    // perspective, openGL version
    float itan= 1 / tanf(radians(fov) * 0.5f);
    float id= 1 / (znear - zfar);

    return Mat4(
        itan/aspect,    0,      0,               0,
        0,              itan,   0,               0,
        0,              0,      (zfar+znear)*id, 2.f*zfar*znear*id,
        0,              0,      -1,              0);
}





