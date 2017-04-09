#include "mat3.h"


const Mat3 Mat3::Zero(0.0);

const Mat3 Mat3::Identity(1.0);

/**
 * @brief Empty constructor, create Identity matrix
 */
Mat3::Mat3() {
    _values[1]=_values[2]=_values[3]=_values[5]=_values[6]=_values[7]=0.0;
    _values[0]=_values[4]=_values[8]=1.0;
}

/**
 * @brief Create a diagonal matrix
 * @param s value of diagonal
 */
Mat3::Mat3(const float &s) {
    _values[1]=_values[2]=_values[3]=_values[5]=_values[6]=_values[7]=0.0;
    _values[0]=_values[4]=_values[8]=s;
}

/**
 * @brief Create a diagonal matrix
 * @param v vector containing _values for diagonal
 */
Mat3::Mat3(const Vec3 &v) {
    _values[1]=_values[2]=_values[3]=_values[5]=_values[6]=_values[7]=0.0;
    _values[0]=v[0];
    _values[4]=v[1];
    _values[8]=v[2];
}

/**
 * @brief Create and fill a matrix with 3 given vectors
 * @param v1 first vector -> first row
 * @param v2 second vector -> second row
 * @param v3 third vector -> third row
 */
Mat3::Mat3(const Vec3 &u, const Vec3 &v, const Vec3 &w) {
    _values[0]=u[0];
    _values[1]=u[1];
    _values[2]=u[2];

    _values[3]=v[0];
    _values[4]=v[1];
    _values[5]=v[2];

    _values[6]=w[0];
    _values[7]=w[1];
    _values[8]=w[2];
}

/**
 * @brief Create and fill matrix with 9 given float _values
 * @param v00 row 0 column 0 value
 * @param v01 row 0 column 1 value
 * @param v02 row 0 column 2 value
 * @param v10 row 1 column 0 value
 * @param v11 row 1 column 1 value
 * @param v12 row 1 column 2 value
 * @param v20 row 2 column 0 value
 * @param v21 row 2 column 1 value
 * @param v22 row 2 column 2 value
 */
Mat3::Mat3(const float &v00, const float &v01, const float &v02, const float &v10, const float &v11, const float &v12, const float &v20, const float &v21, const float &v22) {
    _values[0]=v00;
    _values[1]=v01;
    _values[2]=v02;

    _values[3]=v10;
    _values[4]=v11;
    _values[5]=v12;

    _values[6]=v20;
    _values[7]=v21;
    _values[8]=v22;
}

/**
 * @brief Write access to the i<SUP>th</SUP> matrix value
 * @param i rank of matrix value
 * @return i<SUP>th</SUP> matrix value
 */
float& Mat3::operator[] (const int &i) {
    assert(i >= 0 && i <= 9);
    return _values[i];
}

/**
 * @brief Read access to the i<SUP>th</SUP> matrix value
 * @param i rank of matrix value
 * @return i<SUP>th</SUP> matrix value
 */
const float& Mat3::operator[] (const int &i) const {
    assert(i >= 0 && i <= 9);
    return _values[i];
}

/**
 * @brief Write access to the i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 * @param i row of matrix value
 * @param j column of matrix value
 * @return i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 */
float& Mat3::operator() (const int &i, const int &j) {
    assert(i >= 0 && i <= 2);
    assert(j >= 0 && j <= 2);
    return _values[3*j+i];
}

/**
 * @brief Read access to the i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 * @param i row of matrix value
 * @param j column of matrix value
 * @return i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 */
const float& Mat3::operator() (const int &i, const int &j) const {
    assert(i >= 0 && i <= 2);
    assert(j >= 0 && j <= 2);
    return _values[3*j+i];
}


/**
 * @brief Acces to address of first data
 * @return Pointer on x
 */
float* Mat3::data() {
    return &_values[0];
}

/**
 * @brief Get positive value of current matrix
 * @return current matrix
 */
Mat3 Mat3::operator+ () const {
    return *this;
}

/**
 * @brief Get negative value of current matrix
 * @return negative current matrix
 */
Mat3 Mat3::operator-() const {
    return Mat3(-_values[0], -_values[1], -_values[2], -_values[3], -_values[4], -_values[5], -_values[6], -_values[7], -_values[8]);
}

/**
 * @brief Add a matrix to current one
 * @param m matrix to add
 * @return addition of current matrix and m
 */
Mat3& Mat3::operator+=(const Mat3 &u) {
    for (int i=0;i<9;i++) {
        _values[i]+=u._values[i];
    }
    return *this;
}

/**
 * @brief Subtract a matrix from current one
 * @param m matrix to subtract with
 * @return subtraction of current matrix with m
 */
Mat3& Mat3::operator-=(const Mat3 &u) {
    for (int i=0;i<9;i++) {
        _values[i]-=u._values[i];
    }
    return *this;
}

/**
 * @brief Multiply a matrix to current one
 * @param m matrix to multiply with
 * @return multiplication of current matrix and m
 */
Mat3& Mat3::operator*= (const Mat3 &m) {
    *this = Mat3(*this*m);
    return *this;
}

/**
 * @brief Multiply current matrix with scalar
 * @param s scalar to multiply with
 * @return multiplication of current matrix with s
 */
Mat3& Mat3::operator*=(const float &a) {
    for (int i=0;i<9;i++) {
        _values[i]*=a;
    }
    return *this;
}

/**
 * @brief Divide current matrix with scalar
 * @param s scalar to divide with
 * @return division of current matrix with s
 */
Mat3& Mat3::operator/=(const float &a) {
    for (int i=0; i<9; i++) {
        _values[i]/=a;
    }
    return *this;
}

/**
 * @brief Add two matrix
 * @param m first matrix
 * @param n second matrix
 * @return result of addition m + n
 */
Mat3 operator+(const Mat3 &u, const Mat3 &v) {
    return Mat3(u[0]+v[0], u[1]+v[1], u[2]+v[2], u[3]+v[3], u[4]+v[4], u[5]+v[5], u[6]+v[6], u[7]+v[7], u[8]+v[8]);
}

/**
 * @brief Subtract two matrix
 * @param m first matrix
 * @param n second matrix
 * @return result of subtraction m - n
 */
Mat3 operator-(const Mat3 &u, const Mat3 &v) {
    return Mat3(u[0]-v[0], u[1]-v[1], u[2]-v[2], u[3]-v[3], u[4]-v[4], u[5]-v[5], u[6]-v[6], u[7]-v[7], u[8]-v[8]);
}

/**
 * @brief Multiply a vector with a matrix
 * @param v vector to use
 * @param m matrix to use
 * @return vector result of multiplication v * m
 */
Vec3 operator*(const Vec3 &v, const Mat3 &A) {
    return Vec3 (v[0]*A[0]+v[1]*A[1]+v[2]*A[2], v[0]*A[3]+v[1]*A[4]+v[2]*A[5], v[0]*A[6]+v[1]*A[7]+v[2]*A[8]);
}

/**
 * @brief Multiply two matrix
 * @param m first matrix
 * @param n second matrix
 * @return result of multiplication m * n
 */
Mat3 operator* (const Mat3 &u, const Mat3 &v) {
    Mat3 a;
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
        a._values[3*i+j]=u._values[j] * v._values[3*i] + u._values[3+j] * v._values[3*i+1] + u._values[6+j] * v._values[3*i+2];
    }
    return a;
}

/**
 * @brief Multiply a matrix with a scalar
 * @param m matrix to use
 * @param s scalar to use
 * @return result of multiplication m * s
 */
Mat3 operator* (const Mat3 &m, const float &a) {
    return Mat3(a*m[0],a*m[1],a*m[2],a*m[3],a*m[4],a*m[5],a*m[6],a*m[7],a*m[8]);
}

/**
 * @brief Multiply a matrix with a scalar
 * @param s scalar to use
 * @param m matrix to use
 * @return result of multiplication s * m
 */
Mat3 operator* (const float &a, const Mat3 &m) {
    return m * a;
}

/**
 * @brief Divide a matrix by a scalar
 * @param m matrix to use
 * @param s scalar to use
 * @return result of division m / s
 */
Mat3 operator/ (const Mat3 &m, const float &s) {
    return Mat3(m[0]/s, m[1]/s, m[2]/s, m[3]/s, m[4]/s, m[5]/s, m[6]/s, m[7]/s, m[8]/s);
}

/**
 * @brief Display matrix
 * @param os stream to display
 * @param m matrix to display
 * @return
 */
std::ostream& operator<< (std::ostream& os, const Mat3 &m) {
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            os << m(i,j) << " ";
        }
        os << std::endl;
    }
    return os;
}

/////// Matrix utilities ///////
/**
 * @brief Transpose current matrix
 */
void Mat3::transpose() {
  float n__values[9] = { _values[0], _values[3], _values[6], _values[1], _values[4], _values[7], _values[2], _values[5], _values[8] };
  for(int ii = 0; ii < 9; ii++)
      _values[ii] = n__values[ii];
}

/**
 * @brief Return given matrix transposed
 * @param m matrix to transpose
 * @return m transposed
 */
Mat3 transposed(const Mat3 &m) {
  return Mat3(m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]);
}

/**
 * @brief Calculate determinant of current matrix
 * @return determinant of current matrix
 */
float Mat3::determinant() const {
    return  _values[0] * (_values[4]*_values[8] - _values[7]*_values[5]) - _values[1] * (_values[3]*_values[8] - _values[6]*_values[5]) + _values[2] * (_values[3]*_values[7] - _values[6]*_values[4]);
}

/**
 * @brief Calculate determinant of given matrix
 * @param m matrix to calculate determinant
 * @return determinant of m
 */
float determinant(const Mat3 &m) {
    return m.determinant();
}

/**
 * @brief Calculate trace of current matrix
 * @details Trace of diagonal matrix is the sum of its diagonal _values
 * @return trace of current matrix
 */
float Mat3::trace() const {
    return _values[0] + _values[4] + _values[8];
}

/**
 * @brief Calculate trace of given matrix
 * @details Trace of diagonal matrix is the sum of its diagonal _values
 * @param m matrix to calculate trace
 * @return trace of m
 */
float trace(const Mat3 &m) {
    return m.trace();
}



/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
////////////////////// TOCHECK //////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////



/*!
  \brief Returns the opposite of current matrix
*/

/*!
  \brief Computes the inverse of a matrix A<SUP>-1</SUP>. Recall that A<SUP>-1</SUP>
  can be defined as the transposed adjoint matrix divided by the determinant.

  Returns the null matrix as a warn case if A cannot be inversed.

  The threshold value involved in the singular matrix detection
  is set to 10<SUP>-18</SUP>.
*//*
Mat3 Inverse(const Mat3& M)
{
  float e=M[0]*M[4]*M[8]+M[1]*M[5]*M[6]+M[2]*M[3]*M[7]-M[2]*M[4]*M[6]-M[3]*M[1]*M[8]-M[0]*M[5]*M[7];

  if (fabs(e)<1e-18)
  {
   return Mat3::Zero;
  }

  return T(Adjoint(M))/e;
}
*/

/*!
  \brief Multiplication.
*/

/*!
  \brief Create a rotation matrix given a vector of angles
  that specifies the rotation around each world coordinate axis.
  \param vector Vector of angles in radian.
*/
/*
Mat3 Mat3::Rotation(const Vector& vector)
{
  Mat3 matrix;
  Mat3 temp;

  float cosx=cos(vecto_values[0]);
  float sinx=sin(vecto_values[0]);
  float cosy=cos(vecto_values[1]);
  float siny=sin(vecto_values[1]);
  float cosz=cos(vecto_values[2]);
  float sinz=sin(vecto_values[2]);

  matrix.Id();
  matrix[4]=cosx;
  matrix[8]=cosx;
  matrix[5]=sinx;
  matrix[7]=0.0 - sinx;

  temp.Id();
  temp[0]=cosy;
  temp[8]=cosy;
  temp[2]=0.0 - siny;
  temp[6]=siny;

  matrix=matrix*temp;

  temp.Id();
  temp[0]=cosz;
  temp[4]=cosz;
  temp[1]=sinz;
  temp[3]=0.0 - sinz;

  matrix=matrix*temp;

  return matrix;
}

*/
//! Compute the trace (sum of diagonal terms) of a matrix.
/*
float Trace(const Mat3& a)
{
  return a[0]+a[4]+a[8];
}
*/

/*!
  \brief Create a rotation matrix that rotates a normalized
  vector into another one.
  \param a, b Initial and final vector (should be normalized).
*/
/*
Mat3 Mat3::Rotation(const Vector& a,const Vector& b)
{
  #define EPSILON 1e-8
  Mat3 matrix;

  Vector v=a/b;
  float e=a*b;

  // Almost identical vectors
  if(e>1.0-EPSILON)
  {
   return Mat3::Identity;
  }
  // Almost opposite vectors
  else if(e<-1.0+EPSILON)
  {
   Vector up,left;
   float invlen;
   float fxx,fyy,fzz,fxy,fxz,fyz;
   float uxx,uyy,uzz,uxy,uxz,uyz;
   float lxx,lyy,lzz,lxy,lxz,lyz;

   left[0]=0.0; left[1]=a[2]; left[2]=-a[1];
   if(left*left<EPSILON)
   {
    left[0]=-a[2]; left[1]=0.0; left[2]=a[0];
   }
   invlen=1.0/sqrt(left*left);
   left*=invlen;
   up=left/a;
   // now we have a coordinate system, i.e., a basis
   /* M=(a, up, left), and we want to rotate to:    */
   /* N=(-a, up, -left). This is done with the matrix:*/
   /* N*M^T where M^T is the transpose of M          */
/*
   fxx=-a[0]*a[0]; fyy=-a[1]*a[1]; fzz=-a[2]*a[2];
   fxy=-a[0]*a[1]; fxz=-a[0]*a[2]; fyz=-a[1]*a[2];

   uxx=up[0]*up[0]; uyy=up[1]*up[1]; uzz=up[2]*up[2];
   uxy=up[0]*up[1]; uxz=up[0]*up[2]; uyz=up[1]*up[2];

   lxx=-left[0]*left[0]; lyy=-left[1]*left[1]; lzz=-left[2]*left[2];
   lxy=-left[0]*left[1]; lxz=-left[0]*left[2]; lyz=-left[1]*left[2];
   /* symmetric matrix */
/*
   matrix(0,0)=fxx+uxx+lxx; matrix(0,1)=fxy+uxy+lxy; matrix(0,2)=fxz+uxz+lxz;
   matrix(1,0)=matrix(0,1); matrix(1,1)=fyy+uyy+lyy; matrix(1,2)=fyz+uyz+lyz;
   matrix(2,0)=matrix(0,2); matrix(2,1)=matrix(1,2); matrix(2,2)=fzz+uzz+lzz;
  }
  else
  {
   float h,hvx,hvz,hvxy,hvxz,hvyz;
   h=(1.0-e)/(v*v);
   hvx=h*v[0];
   hvz=h*v[2];
   hvxy=hvx*v[1];
   hvxz=hvx*v[2];
   hvyz=hvz*v[1];
   matrix(0,0)=e+hvx*v[0]; matrix(0,1)=hvxy-v[2];    matrix(0,2)=hvxz+v[1];
   matrix(1,0)=hvxy+v[2];  matrix(1,1)=e+h*v[1]*v[1]; matrix(1,2)=hvyz-v[0];
   matrix(2,0)=hvxz-v[1];  matrix(2,1)=hvyz+v[0];    matrix(2,2)=e+hvz*v[2];
  }
  return matrix;
}

*/
/*!
  \brief Create a rotation matrix about an arbitrary axis.
  \param vector Rotation axis, which should be normalized.
  \param angle Rotation angle (should be in radian).
*/
/*
Mat3 Mat3::Rotation(const Vector& vector,const float& angle)
{
  Vector v=Normalized(vector);

/*   float cosx=cos(angle); */
/*   float sinx=sin(angle); */
/*   _values[0]=v[0] * v[0] + cosx * (1.0 - v[0] * v[0]); */
/*   _values[1]=v[0] * v[1] * (1.0 - cosx) + v[2] * sinx; */
/*   _values[2]=v[0] * v[2] * (1.0 - cosx) - v[1] * sinx; */
/*   _values[3]=v[0] * v[1] * (1.0 - cosx) - v[2] * sinx; */
/*   _values[4]=v[1] * v[1] + cosx * (1.0 - v[1] * v[1]); */
/*   _values[5]=v[1] * v[2] * (1.0 - cosx) + v[0] * sinx; */
/*   _values[6]=v[0] * v[2] * (1.0 - cosx) + v[1] * sinx; */
/*   _values[7]=v[1] * v[2] * (1.0 - cosx) - v[0] * sinx; */
/*   _values[8]=v[2] * v[2] + cosx * (1.0 - v[2] * v[2]); */
/*
  v *= sin(0.5*angle);
  float w = cos(0.5*angle);
  float tx  = 2.0*v[0];
  float ty  = 2.0*v[1];
  float tz  = 2.0*v[2];
  float twx = tx*w;
  float twy = ty*w;
  float twz = tz*w;
  float txx = tx*v[0];
  float txy = ty*v[0];
  float txz = tz*v[0];
  float tyy = ty*v[1];
  float tyz = tz*v[1];
  float tzz = tz*v[2];

  Mat3 r;
  _values[0] = 1.0-(tyy+tzz);
  _values[1] = txy-twz;
  _values[2] = txz+twy;
  _values[3] = txy+twz;
  _values[4] = 1.0-(txx+tzz);
  _values[5] = tyz-twx;
  _values[6] = txz-twy;
  _values[7] = tyz+twx;
  _values[8] = 1.0-(txx+tyy);
  return r;
}
*/
/*!
  \brief Compute the adjoint of the argument matrix.
*/
/*
Mat3 Adjoint(const Mat3& M)
{
  Mat3 A;

  A[0]=M[4]*M[8]-M[7]*M[5];
  A[1]=-(M[3]*M[8]-M[6]*M[5]);
  A[2]=M[3]*M[7]-M[6]*M[4];

  A[3]=-(M[1]*M[8]-M[7]*M[2]);
  A[4]=M[0]*M[8]-M[6]*M[2];
  A[5]=-(M[0]*M[7]-M[6]*M[1]);

  A[6]=M[1]*M[5]-M[4]*M[2];
  A[7]=-(M[0]*M[5]-M[3]*M[2]);
  A[8]=M[0]*M[4]-M[3]*M[1];

  return A;
}
*/
