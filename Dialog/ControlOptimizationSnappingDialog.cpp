#include "controloptimizationsnappingdialog.h"
#include "ui_controloptimizationsnappingdialog.h"
#include "pclviewer.h"

ControlOptimizationSnappingDialog::ControlOptimizationSnappingDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ControlOptimizationSnappingDialog)
{
    ui->setupUi(this);
    m_count = 0;
    pclviewer = NULL;
}

ControlOptimizationSnappingDialog::~ControlOptimizationSnappingDialog()
{
    delete ui;
}

void ControlOptimizationSnappingDialog::resetCount()
{
    m_count = 0;
}

void ControlOptimizationSnappingDialog::setPolyID(int id)
{
    this->poly_id = id;
}

// 执行优化
void ControlOptimizationSnappingDialog::on_pushButton_clicked()
{
    // perform optimise
    ((PCLViewer*)pclviewer)->performOptimize();

    // increase count
    ++m_count;
    QString str;
    str = "迭代第";
    QString tmp;
    char buf[10];
    itoa(m_count, buf, 10);
    tmp = buf;
    str += tmp;
    str += "次";
    this->ui->label->setText(str);
}

void ControlOptimizationSnappingDialog::on_buttonBox_accepted()
{
    QString str;
    str = "迭代第0次";
    this->ui->label->setText(str);
}
