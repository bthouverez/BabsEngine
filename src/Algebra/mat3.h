/**
  @file mat3.h
  @authors ??, refactored by B.Thouverez <bastien.thouverez -at- etu.univ-lyon1.fr>
 */

#ifndef MAT3_H_
#define MAT3_H_

#include "vec3.h"

/**
 * @brief 3x3 matrix class
 */
class Mat3 {

protected:
    //! Array containing values
    float _values[9];

public:
    //! Empty matrix, full of zeros
    static const Mat3 Zero;

    //! Identity matrix, ones pn diagonal and rest fill with zeros
    static const Mat3 Identity;

    /////// Constructors ///////
    Mat3();
    Mat3(const float &s);
    Mat3(const Vec3 &v);
    Mat3(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);
    Mat3(const float &v00, const float &v01, const float &v02, const float &v10, const float &v11, const float &v12, const float &v20, const float &v21, const float &v22);


    /////// Accessors ///////
    float& operator[] (const int &i);
    const float& operator[] (const int &i) const;
    float& operator() (const int &i, const int &j);
    const float& operator() (const int &i, const int &j) const;

    float* data();


    //////// Unary operators ///////
    Mat3 operator+ () const;
    Mat3 operator- () const;
    Mat3& operator+= (const Mat3 &m);
    Mat3& operator-= (const Mat3 &m);
    Mat3& operator*= (const Mat3 &m);
    Mat3& operator*= (const float &s);
    Mat3& operator/= (const float &s);

    /////// Binary operators ///////
    friend Mat3 operator+ (const Mat3 &m, const Mat3 &n);
    friend Mat3 operator- (const Mat3 &m, const Mat3 &n);
    friend Mat3 operator* (const Mat3 &m, const Mat3 &n);
    friend Vec3 operator* (const Vec3 &v, const Mat3 &m);
    friend Mat3 operator* (const Mat3 &m, const float &s);
    friend Mat3 operator* (const float &s, const Mat3 &m);
    friend Mat3 operator/ (const Mat3 &m, const float &s);

    friend std::ostream& operator<< (std::ostream &os, const Mat3 &m);

    /////// Matrix operations ///////
    void transpose();
    friend Mat3 transposed(const Mat3 &m);
    float determinant() const;
    friend float determinant(const Mat3 &m);
    float trace() const;
    friend float trace(const Mat3 &m);

    /////// TO CHECK ///////
    static Mat3 Rotation(const Vec3&);
    static Mat3 Rotation(const Vec3&,const float&);
    static Mat3 Rotation(const Vec3&,const Vec3&);

    friend Mat3 Inverse(const Mat3&);
    friend Mat3 Adjoint(const Mat3&);
    friend float Trace(const Mat3&);
    // Algebra
    float SpectralNorm() const;
    // Singular values methods
    void SingularValueDecomposition (Mat3& L, Vec3& S,Mat3& R) const;
    void SingularValueComposition (const Mat3& L, const Vec3& S,const Mat3& R);
    void QDU(Mat3&,Vec3&,Vec3&) const;
    void ExtractAngleAxis(float&, Vec3&) const;
    // Eigensolver, matrix must be symmetric
    void EigenSolveSymmetric (float eigenvalue[3], Vec3 eigenvectovalues[3]) const;
    Vec3 Eigen() const;
    void Tridiagonal (float diag[3], float subd[2]);
    int QLAlgorithm (float diag[3], float subd[3]);

private:
  float MCR2 (float afCoeff[3]) const;
};

#endif
