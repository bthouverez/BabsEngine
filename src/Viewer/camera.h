#ifndef CAMERA_H_
#define CAMERA_H_

#include "../Algebra/mat4.h"
#include "../Algebra/vec3.h"

class b_Camera {

public:
    b_Camera();
    b_Camera(const Vec3 &o, const Vec3 &d, const Vec3 &u);

    ~b_Camera();

    Mat4 view() { return v; }
    Mat4 projection() { return p; }

    void projection(const float &width, const float &height, const float &fov);
    void lookAt(const Vec3 &f, const Vec3 &t, const Vec3 &u);


private:
    Mat4 v;
    Mat4 p;
};

#endif // CAMERA_H_
