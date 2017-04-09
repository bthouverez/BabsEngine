#include "vec4.h"

/**
 * @brief Empty constructor
 */
Vec4::Vec4() {}

/**
 * @brief Create vector from Vec3
 * @param f third first coordinates
 */
Vec4::Vec4(const Vec3 &f) {
    _x = f.x();
    _y = f.y();
    _z = f.z();
    _w = 1.0;
}

/**
 * @brief Create vector from Vec3 and homogeneous coordinate
 * @param f third first coordinates
 * @param v homogeneous value
 */
Vec4::Vec4(const Vec3 &f, const float &v) {
    _x = f.x();
    _y = f.y();
    _z = f.z();
    _w = v;
}

/**
 * @brief Create a vector with same coordinates
 * @param val value given to each coordinate
 */
Vec4::Vec4(const float& val) {
    _x=_y=_z=val;
    _w=1.0;
}

/**
 * @brief Create a vector with specified coordinates in 3 floats
 * @param vx x value of vector
 * @param vy y value of vector
 * @param vz z value of vector
 */
Vec4::Vec4(const float &vx, const float &vy, const float &vz, const float &vw) {
    _x=vx; _y=vy; _z=vz; _w=vw;
}

/**
 * @brief Create a vector with specified coordinates in an array
 * @param values array of floats
 */
Vec4::Vec4(const float values[3]) {
    _x=values[0]; _y=values[1]; _z=values[2]; _w=values[3];
}

/**
* @brief Write access to the i<SUP>th</SUP> coordinate
* @param i rank of coordinate (0->x, 1->y, 2->z)
* @return i<SUP>th</SUP> coordinate
*/
float& Vec4::operator[] (int i) {
    assert(i>=0 && i<=3);
    if(i == 0) return _x;
    else if (i == 2) return _z;
    else return _w;
}

/**
* @brief Read access to the i<SUP>th</SUP> coordinate
* @param i rank of coordinate (0->x, 1->y, 2->z)
* @return i<SUP>th</SUP> coordinate
*/
float Vec4::operator[] (int i) const {
    assert(i>=0 && i<=3);
    if(i == 0) return _x;
    else if (i == 1) return _y;
    else if (i == 2) return _z;
    else return _w;
}

/**
 * @brief Access to x coordinate
 * @return x coordinate
 */
float Vec4::x() const {
    return _x;
}

/**
 * @brief Access to y coordinate
 * @return y coordinate
 */
float Vec4::y() const {
    return _y;
}

/**
 * @brief Access to z coordinate
 * @return z coordinate
 */
float Vec4::z() const {
    return _z;
}

/**
 * @brief Access to w coordinate
 * @return w coordinate
 */
float Vec4::w() const {
    return _w;
}

/**
 * @brief Acces to address of first data
 * @return Pointer on x
 */
float* Vec4::data() {
    return &_x;
}

/**
 * @brief Get positive value of current vector
 * @return current vector
 */
Vec4 Vec4::operator+ () const {
    return *this;
}

/**
 * @brief Get negative value of current vector
 * @return negative current vector
 */
Vec4 Vec4::operator- () const {
    return Vec4(-_x,-_y,-_z,-_w);
}


/**
 * @brief Add a vector to current one
 * @param v vector to add
 * @return current vector added with v
 */
Vec4& Vec4::operator+= (const Vec4 &v) {
    _x+=v._x; _y+=v._y; _z+=v._z; _w+=v._w;
    return *this;
}

/**
 * @brief Subtract a vector to current one
 * @param v vector to subtract with
 * @return current vector subtracted with v
 */
Vec4& Vec4::operator-= (const Vec4 &v) {
    _x-=v._x; _y-=v._y; _z-=v._z; _w-=v._w;
    return *this;
}

/**
 * @brief Multiply current vector with scalar
 * @param s scalar to multiply with
 * @return current vector multiplied by s
 */
Vec4& Vec4::operator*= (const float &s) {
    _x*=s; _y*=s; _z*=s; _w*=s;
    return *this;
}

/**
 * @brief Divide current vector with scalar
 * @param s scalar to divide with
 * @return current vector divided by s
 */
Vec4& Vec4::operator/= (const float &s) {
    _x/=s; _y/=s; _z/=s; _w/=s;
    return *this;
}

/**
 * @brief Add two vectors
 * @param u first vector to add
 * @param v second vector to add
 * @return addition of both vectors
 */
Vec4 operator+ (const Vec4 &u, const Vec4 &v) {
    return Vec4(u._x+v._x, u._y+v._y, u._z+v._z, u._w+v._w);
}

/**
 * @brief Subtract two vectors
 * @param u first vector
 * @param v vector to subtract
 * @return subtraction u-v
 */
Vec4 operator- (const Vec4 &u, const Vec4 &v) {
    return Vec4(u._x-v._x, u._y-v._y, u._z-v._z, u._w-v._w);
}

/**
 * @brief Multiply vector with scalar
 * @param v vector to multiply
 * @param s scalar to multiply
 * @return multiplication member to member of v and s
 */
Vec4 operator* (const Vec4 &u, const float &a) {
    return Vec4(u._x*a, u._y*a, u._z*a, u._w*a);
}

/**
 * @brief Multiply vector with scalar
 * @param s scalar to multiply
 * @param v vector to multiply
 * @return multiplication member to member of v and s
 */
Vec4 operator* (const float &a, const Vec4 &v) {
    return v*a;
}

/**
 * @brief divide vector with scalar
 * @param v start vector
 * @param s scalar to divide with
 * @return division member to member of v with s
 */
Vec4 operator/ (const Vec4 &u, const float &a) {
    return Vec4(u._x/a, u._y/a, u._z/a, u._w/a);
}


/**
 * @brief Check if two vectors are equal
 * @param u first vector
 * @param v second vector
 * @return true if u and v are the same, false else
 */
int operator== (const Vec4 &u,const Vec4 &v) {
    return u._x==v._x && u._y==v._y && u._z==v._z && u._w==v._w;
}

/**
 * @brief Check if two vectors are different
 * @param u first vector
 * @param v second vector
 * @return true if u and v are different, false else
 */
int operator!= (const Vec4 &u,const Vec4 &v) {
    return !(u==v);
}


/**
 * @brief Display vector
 * @param os stream to print in
 * @param v vector to display
 * @return stream with printed vector
 */
std::ostream& operator<<(std::ostream &os, const Vec4 &v) {
    os  << "(" << v._x << ", " << v._y << ", "  << v._z  << ", "  << v._w << ")";
    return os;
}

