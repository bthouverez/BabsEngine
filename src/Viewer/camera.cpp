#include "camera.h"

/**
 * @brief Default constructor
 * @details Camera looking center, in direction X
 */
b_Camera::b_Camera() {
    v = QMatrix4x4();
    v.lookAt(QVector3D(0.f, 0.f, -5.f), QVector3D(0.f, 0.f, 0.f), QVector3D(0.f, 1.f, 0.f));
    p = QMatrix4x4();
}

/**
 * @brief Camera looking at a point
 * @param o center of the eye
 * @param t point to look at
 * @param u up vector
 */
b_Camera::b_Camera(const QVector3D &o, const QVector3D &t, const QVector3D &u) {
    v = QMatrix4x4();
    v.lookAt(o, t, u);
    p = QMatrix4x4();
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
void b_Camera::lookAt(const QVector3D &f, const QVector3D &t, const QVector3D &u) {
    v = QMatrix4x4();
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
    p = QMatrix4x4();
    p.perspective(fov, width / height, 0.05f, 250.0f);
}


