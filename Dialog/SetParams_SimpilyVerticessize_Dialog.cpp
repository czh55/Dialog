#include "setparams_simpilyverticessize_dialog.h"
#include "ui_setparams_simpilyverticessize_dialog.h"
#include <QMessageBox>
#include <QStringList>

SetParams_SimpilyVerticesSize_Dialog::SetParams_SimpilyVerticesSize_Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetParams_SimpilyVerticesSize_Dialog)
{
    ui->setupUi(this);
    model = new QStandardItemModel();
    initTable();
}

SetParams_SimpilyVerticesSize_Dialog::~SetParams_SimpilyVerticesSize_Dialog()
{
    delete ui;
    delete model;
    delete configIniRead;
}

void SetParams_SimpilyVerticesSize_Dialog::initTable()
{
    // 设置表格每列的宽度
    ui->tableView->setColumnWidth(0, 100);
    ui->tableView->setColumnWidth(1, 600);
    ui->tableView->setColumnWidth(2, 50);

    // 设置表头
    model->setHorizontalHeaderItem(0, new QStandardItem(QObject::tr("variable name")));
    model->setHorizontalHeaderItem(1, new QStandardItem(QObject::tr("role")));
    model->setHorizontalHeaderItem(2, new QStandardItem(QObject::tr("value")));
    ui->tableView->setModel(model);

    configIniRead = new QSettings("config.ini", QSettings::IniFormat);

    // 设置表项
    setTableItem(0, "T_angle_betPointAndPoint", "选取种子点时，种子点法向量与相邻顶点法向量的最大夹角");
    setTableItem(1, "T_angle_betPointAndLine", "点的法向量与直线法向量的最大夹角");
    setTableItem(2, "T_dist_betPointAndLine", "点到直线的最大距离");
    setTableItem(3, "T_lenRatio", "确定可以被剪掉的最小线段占平均线段长度的最大比例");
    setTableItem(4, "T_angle_betLines", "可被合并的两个线段的最大夹角");
}

void SetParams_SimpilyVerticesSize_Dialog::setTableItem(int row, QString str1, QString str2)
{
    model->setItem(row,0, new QStandardItem(str1));
    model->setItem(row,1, new QStandardItem(str2));
    QString str3 = "simplify vertices size/";
    str3.append(str1);
    model->setItem(row,2, new QStandardItem(configIniRead->value(str3).toString()));
}

// 设置参数
void SetParams_SimpilyVerticesSize_Dialog::on_buttonBox_accepted()
{
    QStringList strList;
    strList << "T_angle_betPointAndPoint"
            << "T_angle_betPointAndLine"
            << "T_dist_betPointAndLine"
            << "T_lenRatio"
            << "T_angle_betLines";
    int size = strList.size();
    for(int i = 0; i < size; ++i)
    {
        QString key, value;
        key = model->data(model->index(i, 0), Qt::DisplayRole).toString();
        value = model->data(model->index(i, 2), Qt::DisplayRole).toString();
        key = "simplify vertices size/" + key;
        configIniRead->setValue(key, value);
    }

    QMessageBox::about(NULL, "提示", "参数已设置好");
}
