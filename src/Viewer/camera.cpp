#include "camera.h"

/**
 * @brief Default constructor
 * @details Camera looking center, in direction X
 */
b_Camera::b_Camera() {
    v = Mat4();
    v.lookAt(Vec3(0.f, 0.f, -5.f), Vec3(0.f, 0.f, 0.f), Vec3(0.f, 1.f, 0.f));
    p = Mat4();
}

/**
 * @brief Camera looking at a point
 * @param o center of the eye
 * @param t point to look at
 * @param u up vector
 */
b_Camera::b_Camera(const Vec3 &o, const Vec3 &t, const Vec3 &u) {
    v = Mat4();
    v.lookAt(o, t, u);
    p = Mat4();
}

/**
 * @brief Camera destructor
 */
b_Camera::~b_Camera() {}

/**
 * @brief Set view matrix with lookAt matrix
 * @param f center of the eye
 * @param t point to look at
 * @param u up vector
 */
void b_Camera::lookAt(const Vec3 &f, const Vec3 &t, const Vec3 &u) {
    v = Mat4();
    v.lookAt(f, t, u);
}

/**
 * @brief Get projection set of image of aspect width/height, with a fov degrees opening
 * @param width width of image
 * @param height height of image
 * @param fov opening angle
 * @return projection matrix
 */
void b_Camera::projection(const float &width, const float &height, const float &fov) {
    p = Mat4();
    p.perspective(fov, width / height, 0.05f, 250.0f);
}


