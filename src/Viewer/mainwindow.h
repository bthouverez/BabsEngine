#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_rotY_clicked();

    void on_pushButton_rotZ_clicked();

    void on_pushButton_rotX_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
