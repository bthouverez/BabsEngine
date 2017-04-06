#ifndef CAMERA_H_
#define CAMERA_H_

#include <QMatrix4x4>
#include <QVector3D>

class b_Camera {

public:
    b_Camera();
    b_Camera(const QVector3D &o, const QVector3D &d, const QVector3D &u);

    ~b_Camera();

    QMatrix4x4 view() { return v; }
    QMatrix4x4 projection() { return p; }

    void projection(const float &width, const float &height, const float &fov);
    void lookAt(const QVector3D &f, const QVector3D &t, const QVector3D &u);


private:
    QMatrix4x4 v;
    QMatrix4x4 p;
};

#endif // CAMERA_H_
