#include "mat3.h"


const Mat3 Mat3::Zero(0.0);

const Mat3 Mat3::Identity(1.0);

/**
 * @brief Empty constructor, create Identity matrix
 */
Mat3::Mat3() {
    values[1]=values[2]=values[3]=values[5]=values[6]=values[7]=0.0;
    values[0]=values[4]=values[8]=1.0;
}

/**
 * @brief Create a diagonal matrix
 * @param s value of diagonal
 */
Mat3::Mat3(const double &s) {
    values[1]=values[2]=values[3]=values[5]=values[6]=values[7]=0.0;
    values[0]=values[4]=values[8]=s;
}

/**
 * @brief Create a diagonal matrix
 * @param v vector containing values for diagonal
 */
Mat3::Mat3(const Vec3 &v) {
    values[1]=values[2]=values[3]=values[5]=values[6]=values[7]=0.0;
    values[0]=v[0];
    values[4]=v[1];
    values[8]=v[2];
}

/**
 * @brief Create and fill a matrix with 3 given vectors
 * @param v1 first vector -> first row
 * @param v2 second vector -> second row
 * @param v3 third vector -> third row
 */
Mat3::Mat3(const Vec3 &u, const Vec3 &v, const Vec3 &w) {
    values[0]=u[0];
    values[1]=u[1];
    values[2]=u[2];

    values[3]=v[0];
    values[4]=v[1];
    values[5]=v[2];

    values[6]=w[0];
    values[7]=w[1];
    values[8]=w[2];
}

/**
 * @brief Create and fill matrix with 9 given double values
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
Mat3::Mat3(const double &v00, const double &v01, const double &v02, const double &v10, const double &v11, const double &v12, const double &v20, const double &v21, const double &v22) {
    values[0]=v00;
    values[1]=v01;
    values[2]=v02;

    values[3]=v10;
    values[4]=v11;
    values[5]=v12;

    values[6]=v20;
    values[7]=v21;
    values[8]=v22;
}

/**
 * @brief Write access to the i<SUP>th</SUP> matrix value
 * @param i rank of matrix value
 * @return i<SUP>th</SUP> matrix value
 */
double& Mat3::operator[] (const int &i) {
    assert(i >= 0 && i <= 9);
    return values[i];
}

/**
 * @brief Read access to the i<SUP>th</SUP> matrix value
 * @param i rank of matrix value
 * @return i<SUP>th</SUP> matrix value
 */
const double& Mat3::operator[] (const int &i) const {
    assert(i >= 0 && i <= 9);
    return values[i];
}

/**
 * @brief Write access to the i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 * @param i row of matrix value
 * @param j column of matrix value
 * @return i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 */
double& Mat3::operator() (const int &i, const int &j) {
    assert(i >= 0 && i <= 2);
    assert(j >= 0 && j <= 2);
    return values[3*j+i];
}

/**
 * @brief Read access to the i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 * @param i row of matrix value
 * @param j column of matrix value
 * @return i<SUP>th</SUP> row and j<SUP>th</SUP> column matrix value
 */
const double& Mat3::operator() (const int &i, const int &j) const {
    assert(i >= 0 && i <= 2);
    assert(j >= 0 && j <= 2);
    return values[3*j+i];
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
    return Mat3(-values[0], -values[1], -values[2], -values[3], -values[4], -values[5], -values[6], -values[7], -values[8]);
}

/**
 * @brief Add a matrix to current one
 * @param m matrix to add
 * @return addition of current matrix and m
 */
Mat3& Mat3::operator+=(const Mat3 &u) {
    for (int i=0;i<9;i++) {
        values[i]+=u.values[i];
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
        values[i]-=u.values[i];
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
Mat3& Mat3::operator*=(const double &a) {
    for (int i=0;i<9;i++) {
        values[i]*=a;
    }
    return *this;
}

/**
 * @brief Divide current matrix with scalar
 * @param s scalar to divide with
 * @return division of current matrix with s
 */
Mat3& Mat3::operator/=(const double &a) {
    for (int i=0; i<9; i++) {
        values[i]/=a;
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
        a.values[3*i+j]=u.values[j] * v.values[3*i] + u.values[3+j] * v.values[3*i+1] + u.values[6+j] * v.values[3*i+2];
    }
    return a;
}

/**
 * @brief Multiply a matrix with a scalar
 * @param m matrix to use
 * @param s scalar to use
 * @return result of multiplication m * s
 */
Mat3 operator* (const Mat3 &m, const double &a) {
    return Mat3(a*m[0],a*m[1],a*m[2],a*m[3],a*m[4],a*m[5],a*m[6],a*m[7],a*m[8]);
}

/**
 * @brief Multiply a matrix with a scalar
 * @param s scalar to use
 * @param m matrix to use
 * @return result of multiplication s * m
 */
Mat3 operator* (const double &a, const Mat3 &m) {
    return m * a;
}

/**
 * @brief Divide a matrix by a scalar
 * @param m matrix to use
 * @param s scalar to use
 * @return result of division m / s
 */
Mat3 operator/ (const Mat3 &m, const double &s) {
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
  double n_values[9] = { values[0], values[3], values[6], values[1], values[4], values[7], values[2], values[5], values[8] };
  for(int ii = 0; ii < 9; ii++)
      values[ii] = n_values[ii];
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
double Mat3::determinant() const {
    return  values[0] * (values[4]*values[8] - values[7]*values[5]) - values[1] * (values[3]*values[8] - values[6]*values[5]) + values[2] * (values[3]*values[7] - values[6]*values[4]);
}

/**
 * @brief Calculate determinant of given matrix
 * @param m matrix to calculate determinant
 * @return determinant of m
 */
double determinant(const Mat3 &m) {
    return m.determinant();
}

/**
 * @brief Calculate trace of current matrix
 * @details Trace of diagonal matrix is the sum of its diagonal values
 * @return trace of current matrix
 */
double Mat3::trace() const {
    return values[0] + values[4] + values[8];
}

/**
 * @brief Calculate trace of given matrix
 * @details Trace of diagonal matrix is the sum of its diagonal values
 * @param m matrix to calculate trace
 * @return trace of m
 */
double trace(const Mat3 &m) {
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
  double e=M[0]*M[4]*M[8]+M[1]*M[5]*M[6]+M[2]*M[3]*M[7]-M[2]*M[4]*M[6]-M[3]*M[1]*M[8]-M[0]*M[5]*M[7];

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

  double cosx=cos(vectovalues[0]);
  double sinx=sin(vectovalues[0]);
  double cosy=cos(vectovalues[1]);
  double siny=sin(vectovalues[1]);
  double cosz=cos(vectovalues[2]);
  double sinz=sin(vectovalues[2]);

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
double Trace(const Mat3& a)
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
  double e=a*b;

  // Almost identical vectors
  if(e>1.0-EPSILON)
  {
   return Mat3::Identity;
  }
  // Almost opposite vectors
  else if(e<-1.0+EPSILON)
  {
   Vector up,left;
   double invlen;
   double fxx,fyy,fzz,fxy,fxz,fyz;
   double uxx,uyy,uzz,uxy,uxz,uyz;
   double lxx,lyy,lzz,lxy,lxz,lyz;

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
   double h,hvx,hvz,hvxy,hvxz,hvyz;
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
Mat3 Mat3::Rotation(const Vector& vector,const double& angle)
{
  Vector v=Normalized(vector);

/*   double cosx=cos(angle); */
/*   double sinx=sin(angle); */
/*   values[0]=v[0] * v[0] + cosx * (1.0 - v[0] * v[0]); */
/*   values[1]=v[0] * v[1] * (1.0 - cosx) + v[2] * sinx; */
/*   values[2]=v[0] * v[2] * (1.0 - cosx) - v[1] * sinx; */
/*   values[3]=v[0] * v[1] * (1.0 - cosx) - v[2] * sinx; */
/*   values[4]=v[1] * v[1] + cosx * (1.0 - v[1] * v[1]); */
/*   values[5]=v[1] * v[2] * (1.0 - cosx) + v[0] * sinx; */
/*   values[6]=v[0] * v[2] * (1.0 - cosx) + v[1] * sinx; */
/*   values[7]=v[1] * v[2] * (1.0 - cosx) - v[0] * sinx; */
/*   values[8]=v[2] * v[2] + cosx * (1.0 - v[2] * v[2]); */
/*
  v *= sin(0.5*angle);
  double w = cos(0.5*angle);
  double tx  = 2.0*v[0];
  double ty  = 2.0*v[1];
  double tz  = 2.0*v[2];
  double twx = tx*w;
  double twy = ty*w;
  double twz = tz*w;
  double txx = tx*v[0];
  double txy = ty*v[0];
  double txz = tz*v[0];
  double tyy = ty*v[1];
  double tyz = tz*v[1];
  double tzz = tz*v[2];

  Mat3 r;
  values[0] = 1.0-(tyy+tzz);
  values[1] = txy-twz;
  values[2] = txz+twy;
  values[3] = txy+twz;
  values[4] = 1.0-(txx+tzz);
  values[5] = tyz-twx;
  values[6] = txz-twy;
  values[7] = tyz+twx;
  values[8] = 1.0-(txx+tyy);
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
