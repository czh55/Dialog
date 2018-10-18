#include "setparams_mergingvertices_dialog.h"
#include "ui_setparams_mergingvertices_dialog.h"
#include <QMessageBox>
#include <QStringList>

SetParams_MergingVertices_Dialog::SetParams_MergingVertices_Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetParams_MergingVertices_Dialog)
{
    ui->setupUi(this);
    model = new QStandardItemModel();
    initTable();
}

SetParams_MergingVertices_Dialog::~SetParams_MergingVertices_Dialog()
{
    delete ui;
    delete model;
    delete configIniRead;
}

void SetParams_MergingVertices_Dialog::initTable()
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
    setTableItem(0, "ratio_of_scale", "尺度的比例因子，用来设置邻域半径");
    setTableItem(1, "T_ratio_lineSeg_middle_area", "线段中部区域占据线段长度的比例");
    setTableItem(2, "T_dist_proj_to_lineSegEnd", "投影点与线段端点的最大距离");
    setTableItem(3, "T_maxAngle_bet_two_polys", "两个平行多边形的法向量最大夹角");
    setTableItem(4, "osnap_delta", "数值计算偏导数的步长");
    setTableItem(5, "osnap_max_k_count", "共轭梯度的最大迭代次数");
	setTableItem(6, "f1_weight", "更新后顶点的位置与原来位置的距离");
	setTableItem(7, "f2_weight", "更新后顶点位置与目标顶点位置的距离");
	setTableItem(8, "f3_1_weight", "dist[更新后顶点，目标顶点的前驱节点] - dist[目标顶点，目标顶点的前驱节点]");
	setTableItem(9, "f3_2_weight", "dist[更新后顶点，目标顶点的后继节点] - dist[目标顶点，目标顶点的后继节点]");
	setTableItem(10, "f8_weight", "f5的对偶版本(点-点)，更新后顶点位置与两个平面交线的距离");
	setTableItem(11, "f4_weight", "更新后顶点位置到匹配边的距离");
	setTableItem(12, "f5_weight", "更新后顶点位置与两个平面交线的距离");
	setTableItem(13, "f6_weight", "更新后顶点到匹配边两个端点距离之和 - 匹配边长度");
	setTableItem(14, "f7_weight", "更新后顶点到匹配关系所属支撑平面的距离");
}

void SetParams_MergingVertices_Dialog::setTableItem(int row, QString str1, QString str2)
{
    model->setItem(row,0, new QStandardItem(str1));
    model->setItem(row,1, new QStandardItem(str2));
    QString str3 = "merge vertices/";
    str3.append(str1);
    model->setItem(row,2, new QStandardItem(configIniRead->value(str3).toString()));
}

// 设置参数
void SetParams_MergingVertices_Dialog::on_buttonBox_accepted()
{
    QStringList strList;
    strList << "ratio_of_scale"
            << "T_ratio_lineSeg_middle_area"
            << "T_dist_proj_to_lineSegEnd"
            << "T_maxAngle_bet_two_polys"
            << "osnap_delta"
            << "osnap_max_k_count"
			<< "f1_weight"
			<< "f2_weight"
			<< "f3_1_weight"
			<< "f3_2_weight"
			<< "f8_weight"
			<< "f4_weight"
			<< "f5_weight"
			<< "f6_weight"
			<< "f7_weight";
    int size = strList.size();
    for(int i = 0; i < size; ++i)
    {
        QString key, value;
        key = model->data(model->index(i, 0), Qt::DisplayRole).toString();
        value = model->data(model->index(i, 2), Qt::DisplayRole).toString();
        key = "merge vertices/" + key;
        configIniRead->setValue(key, value);
    }

    QMessageBox::about(NULL, "提示", "参数已设置好");
}
