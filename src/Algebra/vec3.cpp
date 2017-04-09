#include "vec3.h"

/**
 * @brief Empty constructor
 */
Vec3::Vec3() {}

/**
 * @brief Create vector from two points
 * @details Create vector from point f to point t
 * @param f start point of vector
 * @param t end of vector
 */
Vec3::Vec3(const Vec3 &f, const Vec3 &t) {
    _x = t.x() - f.x();
    _y = t.y() - f.y();
    _z = t.z() - f.z();
}

/**
 * @brief Create a vector with same coordinates
 * @param val value given to each coordinate
 */
Vec3::Vec3(const float& val) {
    _x=_y=_z=val;
}

/**
 * @brief Create a vector with specified coordinates in 3 floats
 * @param vx x value of vector
 * @param vy y value of vector
 * @param vz z value of vector
 */
Vec3::Vec3(const float& vx, const float& vy, const float& vz) {
    _x=vx; _y=vy; _z=vz;
}

/**
 * @brief Create a vector with specified coordinates in an array
 * @param values array of floats
 */
Vec3::Vec3(const float values[3]) {
    _x=values[0]; _y=values[1]; _z=values[2];
}

/**
* @brief Write access to the i<SUP>th</SUP> coordinate
* @param i rank of coordinate (0->x, 1->y, 2->z)
* @return i<SUP>th</SUP> coordinate
*/
float& Vec3::operator[] (int i) {
    assert(i>=0 && i<=2);
    if(i == 0) return _x;
    else if (i == 1) return _y;
    else return _z;
}

/**
* @brief Read access to the i<SUP>th</SUP> coordinate
* @param i rank of coordinate (0->x, 1->y, 2->z)
* @return i<SUP>th</SUP> coordinate
*/
float Vec3::operator[] (int i) const {
    assert(i>=0 && i<=2);
    if(i == 0) return _x;
    else if (i == 1) return _y;
    else return _z;
}

/**
 * @brief Access to x coordinate
 * @return x coordinate
 */
float Vec3::x() const {
    return _x;
}

/**
 * @brief Access to y coordinate
 * @return y coordinate
 */
float Vec3::y() const {
    return _y;
}

/**
 * @brief Access to z coordinate
 * @return z coordinate
 */
float Vec3::z() const {
    return _z;
}

/**
 * @brief Acces to address of first data
 * @return Pointer on x
 */
float* Vec3::data() {
    return &_x;
}

/**
 * @brief Get positive value of current vector
 * @return current vector
 */
Vec3 Vec3::operator+ () const {
    return *this;
}

/**
 * @brief Get negative value of current vector
 * @return negative current vector
 */
Vec3 Vec3::operator- () const {
    return Vec3(-_x,-_y,-_z);
}


/**
 * @brief Add a vector to current one
 * @param v vector to add
 * @return current vector added with v
 */
Vec3& Vec3::operator+= (const Vec3 &v) {
    _x+=v._x; _y+=v._y; _z+=v._z;
    return *this;
}

/**
 * @brief Subtract a vector to current one
 * @param v vector to subtract with
 * @return current vector subtracted with v
 */
Vec3& Vec3::operator-= (const Vec3 &v) {
    _x-=v._x; _y-=v._y; _z-=v._z;
    return *this;
}

/**
 * @brief Multiply current vector with scalar
 * @param s scalar to multiply with
 * @return current vector multiplied by s
 */
Vec3& Vec3::operator*= (const float &s) {
    _x*=s; _y*=s; _z*=s;
    return *this;
}

/**
 * @brief Divide current vector with scalar
 * @param s scalar to divide with
 * @return current vector divided by s
 */
Vec3& Vec3::operator/= (const float &s) {
    _x/=s; _y/=s; _z/=s;
    return *this;
}

/*
//! Compare two vectors.
int operator> (const Vec3& u, const Vec3& v) {
    return ((u.x>v.x) && (u.y>v.y) && (u.z>v.z));
}

//! Compare two vectors.
int operator< (const Vec3& u, const Vec3& v) {
    return ((u.x<v.x) && (u.y<v.y) && (u.z<v.z));
}
//! Overloaded
int operator>= (const Vec3& u, const Vec3& v) {
    return ((u.x>=v.x) && (u.y>=v.y) && (u.z>=v.z));
}
//! Overloaded
int operator<= (const Vec3& u, const Vec3& v) {
    return ((u.x<=v.x) && (u.y<=v.y) && (u.z<=v.z));
}
*/

/**
 * @brief Add two vectors
 * @param u first vector to add
 * @param v second vector to add
 * @return addition of both vectors
 */
Vec3 operator+ (const Vec3 &u, const Vec3 &v) {
    return Vec3(u._x+v._x, u._y+v._y, u._z+v._z);
}

/**
 * @brief Subtract two vectors
 * @param u first vector
 * @param v vector to subtract
 * @return subtraction u-v
 */
Vec3 operator- (const Vec3 &u, const Vec3 &v) {
    return Vec3(u._x-v._x, u._y-v._y, u._z-v._z);
}

/**
 * @brief Multiply vector with scalar
 * @param v vector to multiply
 * @param s scalar to multiply
 * @return multiplication member to member of v and s
 */
Vec3 operator* (const Vec3 &u, const float &a) {
    return Vec3(u._x*a, u._y*a, u._z*a);
}

/**
 * @brief Multiply vector with scalar
 * @param s scalar to multiply
 * @param v vector to multiply
 * @return multiplication member to member of v and s
 */
Vec3 operator* (const float &a, const Vec3 &v) {
    return v*a;
}

/**
 * @brief divide vector with scalar
 * @param v start vector
 * @param s scalar to divide with
 * @return division member to member of v with s
 */
Vec3 operator/ (const Vec3 &u, const float &a) {
    return Vec3(u._x/a, u._y/a, u._z/a);
}


/**
 * @brief Check if two vectors are equal
 * @param u first vector
 * @param v second vector
 * @return true if u and v are the same, false else
 */
int operator== (const Vec3 &u,const Vec3 &v) {
    return u._x==v._x && u._y==v._y && u._z==v._z;
}

/**
 * @brief Check if two vectors are different
 * @param u first vector
 * @param v second vector
 * @return true if u and v are different, false else
 */
int operator!= (const Vec3 &u,const Vec3 &v) {
    return !(u==v);
}


/**
 * @brief Display vector
 * @param os stream to print in
 * @param v vector to display
 * @return stream with printed vector
 */
std::ostream& operator<<(std::ostream &os, const Vec3 &v) {
    os  << "(" << v._x << ", " << v._y << ", "  << v._z << ")";
    return os;
}


/**
 * @brief Compute scalar product between two vectors
 * @param u first vector
 * @param v second vector
 * @return float resulting the dot product of u and v
 */
float dotProduct(const Vec3 &u, const Vec3 &v) {
    return u._x*v._x + u._y*v._y + u._z*v._z;
}

/**
 * @brief Compute cross product between two vectors
 * @param u first vector
 * @param v second vector
 * @return vector resulting the cross product of u and v
 */
Vec3 crossProduct(const Vec3 &u, const Vec3 &v) {
    return Vec3(u._y*v._z-u._z*v._y, u._z*v._x-u._x*v._z, u._x*v._y-u._y*v._x);
}

/**
 * @brief Calculate the length of vector
 * @param v vector to compute norm
 * @return length of v
 */
float length(const Vec3 &v) {
    return sqrt(v._x*v._x + v._y*v._y + v._z*v._z);
}

/**
 * @brief Calculate length of current vector
 * @return length of current vector
 */
float Vec3::length() {
    return sqrt(_x*_x + _y*_y + _z*_z);
}

/**
 * @brief Normalize current vector
 */
void Vec3::normalize() {
    float n = length();
    _x /= n; _y /= n; _z /= n;
}

/**
 * @brief Get normalized vector
 * @param v vector to normalize
 * @return v normalized
 */
Vec3 normalized(const Vec3 &v) {
    float n = length(v);
    return Vec3(v._x/n, v._y/n, v._z/n);
}

/**
 * @brief Calculate cosine between two vectors
 * @details Based on dot product propriety, cos(u, v) = ||u||.||v||
 * @param u first vector
 * @param v second vector
 * @return cosine of u and v
 */
float cosine(const Vec3 &u,const Vec3 &v) {
    return dotProduct(normalized(u), normalized(v));
}

/**
 * @brief Check if two vectors are collinear
 * @details Collinearity is verified if cross product between both vectors results in zero vector
 * @param u first vector
 * @param v second vector
 * @return true if u and v are collinear, false else
 */
bool collinear(const Vec3 &u, const Vec3 &v) {
    return crossProduct(u, v) == Vec3(0.0, 0.0, 0.0);
}

/**
 * @brief Check if two vectors are coplanar
 * @details Using scalar triple product, if u.(v*w) = 0 then u, v and w are coplanar
 * @param u first vector
 * @param v second vector
 * @param w third vector
 * @return true if u, v and w are coplanar, false else
 */
bool coplanar(const Vec3 &u, const Vec3 &v, const Vec3 &w) {
    return dotProduct(u, crossProduct(v, w)) == 0;
}


/**
 * @brief Calculate L2-distance bewteen two points
 * @param u first point
 * @param v second point
 * @return distance bewteen u and v
 */
float distance(const Vec3 &u, const Vec3 &v) {
    return sqrt( (u.x()-v.x())*(u.x()-v.x()) + (u.y()-v.y())*(u.y()-v.y()) + (u.z()-v.z())*(u.z()-v.z()) );
}

// Are these functions usefull? std can make it
/*
//! Return the minimum between two integers
inline int min(const int a, const int b) {
    return a<b ? a : b;
}

//! Return the maximum between two integers
inline int max(const int a, const int b) {
    return a>b ? a : b;
}

//! Return the minimum between two floats
inline float min(const float& a, const float& b) {
    return a<b ? a : b;
}

//! Return the maximum between two floats
inline float max(const float& a, const float& b) {
    return a>b ? a : b;
}

//! Return the minimum bexteen three floats
inline float min(const float& a, const float& b, const float& c) {
    return min(min(a, c), min(b, c));
}

//! Return the maximum bexteen three floats
inline float max(const float& a, const float& b, const float& c) {
    return max(max(a, c), max(b, c));
}
*/
