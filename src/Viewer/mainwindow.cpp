#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_rotX_clicked()
{
    //ui->widget->rotateX(5.0);
}

void MainWindow::on_pushButton_rotY_clicked()
{
    //ui->widget->rotateY(5.0);
}

void MainWindow::on_pushButton_rotZ_clicked()
{
    //ui->widget->rotateZ(5.0);
}
