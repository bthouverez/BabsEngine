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
    double values[9];

public:
    //! Empty matrix, full of zeros
    static const Mat3 Zero;

    //! Identity matrix, ones pn diagonal and rest fill with zeros
    static const Mat3 Identity;

    /////// Constructors ///////
    Mat3();
    Mat3(const double &s);
    Mat3(const Vec3 &v);
    Mat3(const Vec3 &v1, const Vec3 &v2, const Vec3 &v3);
    Mat3(const double &v00, const double &v01, const double &v02, const double &v10, const double &v11, const double &v12, const double &v20, const double &v21, const double &v22);


    /////// Accessors ///////
    double& operator[] (const int &i);
    const double& operator[] (const int &i) const;
    double& operator() (const int &i, const int &j);
    const double& operator() (const int &i, const int &j) const;


    //////// Unary operators ///////
    Mat3 operator+ () const;
    Mat3 operator- () const;
    Mat3& operator+= (const Mat3 &m);
    Mat3& operator-= (const Mat3 &m);
    Mat3& operator*= (const Mat3 &m);
    Mat3& operator*= (const double &s);
    Mat3& operator/= (const double &s);

    /////// Binary operators ///////
    friend Mat3 operator+ (const Mat3 &m, const Mat3 &n);
    friend Mat3 operator- (const Mat3 &m, const Mat3 &n);
    friend Mat3 operator* (const Mat3 &m, const Mat3 &n);
    friend Vec3 operator* (const Vec3 &v, const Mat3 &m);
    friend Mat3 operator* (const Mat3 &m, const double &s);
    friend Mat3 operator* (const double &s, const Mat3 &m);
    friend Mat3 operator/ (const Mat3 &m, const double &s);

    friend std::ostream& operator<< (std::ostream &os, const Mat3 &m);

    /////// Matrix operations ///////
    void transpose();
    friend Mat3 transposed(const Mat3 &m);
    double determinant() const;
    friend double determinant(const Mat3 &m);
    double trace() const;
    friend double trace(const Mat3 &m);

    /////// TO CHECK ///////
    static Mat3 Rotation(const Vec3&);
    static Mat3 Rotation(const Vec3&,const double&);
    static Mat3 Rotation(const Vec3&,const Vec3&);

    friend Mat3 Inverse(const Mat3&);
    friend Mat3 Adjoint(const Mat3&);
    friend double Trace(const Mat3&);
    // Algebra
    double SpectralNorm() const;
    // Singular values methods
    void SingularValueDecomposition (Mat3& L, Vec3& S,Mat3& R) const;
    void SingularValueComposition (const Mat3& L, const Vec3& S,const Mat3& R);
    void QDU(Mat3&,Vec3&,Vec3&) const;
    void ExtractAngleAxis(double&, Vec3&) const;
    // Eigensolver, matrix must be symmetric
    void EigenSolveSymmetric (double eigenvalue[3], Vec3 eigenvectovalues[3]) const;
    Vec3 Eigen() const;
    void Tridiagonal (double diag[3], double subd[2]);
    int QLAlgorithm (double diag[3], double subd[3]);

private:
  double MCR2 (double afCoeff[3]) const;
};

#endif
