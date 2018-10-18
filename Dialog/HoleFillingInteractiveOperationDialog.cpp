#include "holefillinginteractiveoperationdialog.h"
#include "ui_holefillinginteractiveoperationdialog.h"
#include <iostream>
#include "pclviewer.h"
#include "earclip2.h"
#include <qmessagebox.h>
#include "earclip.h"
#include <pcl/common/centroid.h>
#include "ConstrainedDelaunayTriangulation.h"

using namespace std;

HoleFillingInteractiveOperationDialog::HoleFillingInteractiveOperationDialog(QWidget *parent) :
	QDialog(parent),
	ui(new Ui::HoleFillingInteractiveOperationDialog)
{
	ui->setupUi(this);
	m_vertex1 = NULL;
	m_vertex2 = NULL;
	selected_poly_id = -1;
	m_vi_one_vertex.poly_id = -1;
	m_vi_one_vertex.v = NULL;
	m_camera_pos << 0, 0, 0;
	m_min_x = m_max_x = m_min_y = m_min_z = DBL_MAX;
	m_bottom_left_vertex = m_bottom_right_vertex = NULL;
	m_four_vertices.resize(4);
}

HoleFillingInteractiveOperationDialog::~HoleFillingInteractiveOperationDialog()
{
	delete ui;
	m_CMD = NULL;
	m_pclviewer = NULL;
}

// 选择多边形
void HoleFillingInteractiveOperationDialog::on_pushButton_select_poly_clicked()
{
	memset(m_CMD, 0, 128);
	strcpy(m_CMD, "select poly");
	cout << "进入选择多边形模式" << endl;
}

// 更新操作面板上的信息
void HoleFillingInteractiveOperationDialog::updateInfo()
{
	if (strcmp(m_CMD, "select poly") == 0) {
		QString str;
		str = QString::number(this->selected_poly_id, 10);
		ui->lineEdit_poly_id->setText(str);
		cout << "selected_poly_id = " << selected_poly_id << endl;
	}

	if (strcmp(m_CMD, "select vertex1") == 0) {
		if (!m_vertex1) return;
		QString str;
		int pos = (int)m_vertex1;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex1->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[2]));
		str.append("]");
		ui->lineEdit_p1->setText(str);
	}

	if (strcmp(m_CMD, "select vertex2") == 0) {
		if (!m_vertex2) return;
		QString str;
		int pos = (int)m_vertex2;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex2->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[2]));
		str.append("]");
		ui->lineEdit_p2->setText(str);
	}
}

// 取消对当前多边形的选择
void HoleFillingInteractiveOperationDialog::on_pushButton_simplify_cancel_poly_selection_clicked()
{
	if (selected_poly_id == -1) return;

	// 删除viewer中的多边形
	char buf[64];
	string name;
	//POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	for (int i = 0; i < 1024; ++i) {
		itoa(i, buf, 16);
		name = "display_";
		name.append(buf);
		name.append("_l");
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}

	selected_poly_id = -1;
	memset(m_CMD, 0, 128);
	ui->lineEdit_poly_id->setText("-1");

	// 删除viewer中的顶点1
	((PCLViewer*)m_pclviewer)->viewer->removeShape("vertex1");
	ui->lineEdit_p1->setText("");
	// 删除viewer中的顶点2
	((PCLViewer*)m_pclviewer)->viewer->removeShape("vertex2");
	ui->lineEdit_p2->setText("");

	//  删除viewer中的三角化点集边界
	int max_poly_size = 1024;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}
}

// 选择顶点1
void HoleFillingInteractiveOperationDialog::on_pushButton_select_p1_clicked()
{
	memset(m_CMD, 0, 128);
	strcpy(m_CMD, "select vertex1");
	cout << "准备选择第一个顶点" << endl;
	//is_first_call = true;
}

// 选择顶点2
void HoleFillingInteractiveOperationDialog::on_pushButton_select_p2_clicked()
{
	memset(m_CMD, 0, 128);
	strcpy(m_CMD, "select vertex2");
	cout << "准备选择第二个顶点" << endl;
}

// 生成线段
void HoleFillingInteractiveOperationDialog::on_pushButton_gener_line_seg_clicked()
{
	on_pushButton_del_lineseg_clicked();

	if (selected_poly_id == -1) {
		cerr << "没有选中多边形!" << endl;
		return;
	}

	if (!m_vertex1 || !m_vertex2) {
		cerr << "顶点不能为空!" << endl;
		return;
	}

	POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	if ((POLYGON2 *)m_vertex1->relevant_polys[pPoly->getPolyIndex(m_vertex1)].poly != pPoly ||
		(POLYGON2 *)m_vertex2->relevant_polys[pPoly->getPolyIndex(m_vertex2)].poly != pPoly) {
		cerr << "顶点必须都属于同一个多边形!" << endl;
		return;
	}

	if (m_vertex1 == m_vertex2) {
		cerr << "两个顶点不能相等!" << endl;
		return;
	}

	// 默认从m_vertex1找m_vertex2
	line_segs.push_back(m_vertex1);
	int poly_index = pPoly->getPolyIndex(m_vertex1);
	Vertex2* pv2 = (Vertex2*)m_vertex1->relevant_polys[poly_index].next_point;
	while (pv2 != m_vertex2) {
		line_segs.push_back(pv2);
		poly_index = pPoly->getPolyIndex(pv2);
		pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
	}
	line_segs.push_back(pv2);

	// 显示线段
	char buf[64];
	string name;
	int i = 0;
	Vertex2* pv2_next;
	pv2 = m_vertex1;
	pcl::PointXYZ p1, p2;
	while (pv2 != m_vertex2) {
		poly_index = pPoly->getPolyIndex(pv2);
		pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
		p1.x = pv2->pos[0];
		p1.y = pv2->pos[1];
		p1.z = pv2->pos[2];
		p2.x = pv2_next->pos[0];
		p2.y = pv2_next->pos[1];
		p2.z = pv2_next->pos[2];
		itoa(i++, buf, 16);
		name = "l_";
		name.append(buf);

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 0, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		pv2 = pv2_next;
	}
}

// 清空线段，并移除在viewer中的显示
void HoleFillingInteractiveOperationDialog::on_pushButton_del_lineseg_clicked()
{
	char buf[64];
	string name;
	int line_segs_size = line_segs.size();
	for (int i = 0; i < line_segs_size - 1; ++i) {
		itoa(i, buf, 16);
		name = "l_";
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}

	line_segs.clear();

	// 删除viewer中的顶点1
	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex1");
	ui->lineEdit_p1->setText("");
	// 删除viewer中的顶点2
	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex2");
	ui->lineEdit_p2->setText("");
}

// 切换线段
void HoleFillingInteractiveOperationDialog::on_pushButton_switch_ling_seg_clicked()
{
	if (!m_vertex1 || !m_vertex2) return;

	Vertex2* v_tmp = m_vertex1;
	m_vertex1 = m_vertex2;
	m_vertex2 = v_tmp;

	on_pushButton_gener_line_seg_clicked();

	QString str;
	int pos = (int)m_vertex1;
	str.append(QString::number(pos));
	str.append(" [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_p1->setText(str);

	pos = (int)m_vertex2;
	str = "";
	str.append(QString::number(pos));
	str.append(" [");
	str.append(QString::number(m_vertex2->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex2->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex2->pos[2]));
	str.append("]");
	ui->lineEdit_p2->setText(str);
}

// 线段化简
// 将m_vertex1与m_vertex2之间的顶点全部删除(保留m_vertex1与m_vertex2)
void HoleFillingInteractiveOperationDialog::on_pushButton_simplify_line_segs_clicked()
{
	if (selected_poly_id == -1) {
		cerr << "线段化简 selected_poly_id == -1" << endl;
		return;
	}

	if (line_segs.size() < 2) {
		cerr << "线段化简 line_segs.size() == " << line_segs.size() << endl;
		return;
	}

	POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];

	int size = line_segs.size();
	for (int i = 1; i < size - 1; ++i) {
		pPoly->delVertex2(line_segs[i]);
	}
	pPoly->is_tri = true;	// 当前多边形需要重新执行三角化

							// 更新多边形顶点信息
	((PCLViewer*)m_pclviewer)->m_hf.updateVerCloud();

	// 对更新后的多边形重新进行三角化
	on_pushButton_clicked();
}

// 找当前顶点的前一个顶点
void HoleFillingInteractiveOperationDialog::on_pushButton_find_prev_point_clicked()
{
	if (selected_poly_id == -1) {
		cerr << "selected_poly_id == -1" << endl;
		return;
	}
	POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];

	if (strcmp(m_CMD, "select vertex1") == 0) {
		if (!m_vertex1) {
			m_vertex1 = pPoly->start_point;
		}

		int poly_index = pPoly->getPolyIndex(m_vertex1);

		if (poly_index == -1) {
			m_vertex1 = pPoly->start_point;
		}
		else {
			m_vertex1 = (Vertex2*)m_vertex1->relevant_polys[poly_index].prev_point;
		}

		QString str;
		int pos = (int)m_vertex1;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex1->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[2]));
		str.append("]");
		ui->lineEdit_p1->setText(str);

		pcl::PointXYZ point;
		point.x = m_vertex1->pos[0];
		point.y = m_vertex1->pos[1];
		point.z = m_vertex1->pos[2];
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex1");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "vertex1");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex1");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex1");
	}

	if (strcmp(m_CMD, "select vertex2") == 0) {
		if (!m_vertex2) {
			m_vertex2 = pPoly->start_point;
		}

		int poly_index = pPoly->getPolyIndex(m_vertex2);

		if (poly_index == -1) {
			m_vertex2 = pPoly->start_point;
		}
		else {
			m_vertex2 = (Vertex2*)m_vertex2->relevant_polys[poly_index].prev_point;
		}

		QString str;
		int pos = (int)m_vertex2;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex2->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[2]));
		str.append("]");
		ui->lineEdit_p2->setText(str);

		pcl::PointXYZ point;
		point.x = m_vertex2->pos[0];
		point.y = m_vertex2->pos[1];
		point.z = m_vertex2->pos[2];
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex2");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "vertex2");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex2");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex2");
	}
}

// 找当前顶点的后一个顶点
void HoleFillingInteractiveOperationDialog::on_pushButton_find_next_point_clicked()
{
	if (selected_poly_id == -1) {
		cerr << "selected_poly_id == -1" << endl;
		return;
	}
	POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];

	if (strcmp(m_CMD, "select vertex1") == 0) {
		if (!m_vertex1) {
			m_vertex1 = pPoly->start_point;
		}

		int poly_index = pPoly->getPolyIndex(m_vertex1);

		if (poly_index == -1) {
			m_vertex1 = pPoly->start_point;
		}
		else {
			m_vertex1 = (Vertex2*)m_vertex1->relevant_polys[poly_index].next_point;
		}

		QString str;
		int pos = (int)m_vertex1;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex1->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex1->pos[2]));
		str.append("]");
		ui->lineEdit_p1->setText(str);

		pcl::PointXYZ point;
		point.x = m_vertex1->pos[0];
		point.y = m_vertex1->pos[1];
		point.z = m_vertex1->pos[2];
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex1");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "vertex1");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex1");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex1");
	}

	if (strcmp(m_CMD, "select vertex2") == 0) {
		if (!m_vertex2) {
			m_vertex2 = pPoly->start_point;
		}

		int poly_index = pPoly->getPolyIndex(m_vertex2);

		if (poly_index == -1) {
			m_vertex2 = pPoly->start_point;
		}
		else {
			m_vertex2 = (Vertex2*)m_vertex2->relevant_polys[poly_index].next_point;
		}

		QString str;
		int pos = (int)m_vertex2;
		str.append(QString::number(pos));
		str.append(" [");
		str.append(QString::number(m_vertex2->pos[0]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[1]));
		str.append(" ");
		str.append(QString::number(m_vertex2->pos[2]));
		str.append("]");
		ui->lineEdit_p2->setText(str);

		pcl::PointXYZ point;
		point.x = m_vertex2->pos[0];
		point.y = m_vertex2->pos[1];
		point.z = m_vertex2->pos[2];
		((PCLViewer*)m_pclviewer)->viewer->removePointCloud("vertex2");
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(point);
		((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "vertex2");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex2");
		((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex2");
	}
}

// 寻找具有角对角关系的顶点
void HoleFillingInteractiveOperationDialog::on_pushButton_search_angle2angle_vertices_clicked()
{
	// 测试代码: 针对m_vertex1
	//cout << ((PCLViewer*)m_pclviewer)->m_hf.getSingleRelationCount(m_vertex1) << endl;
	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("angle2angle");

	pcl::PointXYZ p;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size(); ++i) {
		if (((PCLViewer*)m_pclviewer)->m_hf.getSingleRelationCount(((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]) > 2) {
			p = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices_cloud->points[i];
			cloud->push_back(p);
		}
	}

	if (cloud->size() == 0) {
		cout << "不存在具有角对角关系的顶点" << endl;
		return;
	}

	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "angle2angle");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "angle2angle");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 8, "angle2angle");
	cout << "共存在" << cloud->size() << "个具有角对角关系的顶点" << endl;
}

///////////////////////////////////////构造封闭三角化顶点边界////////////////////////////////////////////
// 模式1
// 单击只选择一个顶点模式
void HoleFillingInteractiveOperationDialog::on_pushButton_select_one_vertex_clicked()
{
	ui->lineEdit_one_vertex->setText("");
	ui->radioButton_select_one_vertex->setChecked(true);

	this->m_v_begin_mode2.poly_id = -1;
	this->m_v_begin_mode2.v = NULL;
	this->m_v_end_mode2.poly_id = -1;
	this->m_v_end_mode2.v = NULL;
	this->m_v1_mode3.poly_id = -1;
	this->m_v1_mode3.v = NULL;
	this->m_v2_mode3.poly_id = -1;
	this->m_v2_mode3.v = NULL;
	this->m_v3_mode3.poly_id = -1;
	this->m_v3_mode3.v = NULL;
	this->m_v4_mode3.poly_id = -1;
	this->m_v4_mode3.v = NULL;

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "m_vertex1 == NULL");
		return;
	}

	if (((PCLViewer*)m_pclviewer)->m_hf.getSingleRelationCount(m_vertex1) > 2) {
		QMessageBox::information(this, "提示", "不能选择角对角顶点！");
		return;
	}

	bool is_vertex_belongto_poly = false;
	for (int i = 0; i < m_vertex1->relevant_polys.size(); ++i) {
		if ((POLYGON2*)m_vertex1->relevant_polys[i].poly == ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]) {
			is_vertex_belongto_poly = true;
			break;
		}
	}
	if (!is_vertex_belongto_poly) {
		QMessageBox::information(this, "提示", "顶点不属于该多边形！");
		return;
	}

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_one_vertex->setText(str);
}

// 获取当前camera位置
void HoleFillingInteractiveOperationDialog::on_pushButton_get_camera_pos_clicked()
{
	vector<pcl::visualization::Camera> cameras;
	((PCLViewer*)m_pclviewer)->viewer->getCameras(cameras);
	// cout << "cameras.size = " << cameras.size() << endl;
	m_camera_pos << cameras[0].pos[0], cameras[0].pos[1], cameras[0].pos[2];

	QString str;
	str.append("[");
	str.append(QString::number(m_camera_pos[0]));
	str.append(" ");
	str.append(QString::number(m_camera_pos[1]));
	str.append(" ");
	str.append(QString::number(m_camera_pos[2]));
	str.append("]");

	ui->lineEdit_camera_pos->setText(str);
}

// 构造三角化点集
void HoleFillingInteractiveOperationDialog::on_pushButton_gen_points_set_for_tri_clicked()
{
	m_points_set_for_tri.clear();
	if (ui->radioButton_select_one_vertex->isChecked()) {
		genPointsSetWithOneVertex();
	}

	if (ui->radioButton_select_two_vertices->isChecked()) {
		genPointsSetWithTwoVertices();
	}

	if (ui->radioButton_select_four_vertices->isChecked()) {
		genPointsSetWithFourVertices();
	}
	cout << "三角化点集size = " << m_points_set_for_tri.size() << endl;
}

// 模式3：利用四个顶点生成孔洞边界点集，适用于连接两个分离多边形的情况
bool HoleFillingInteractiveOperationDialog::genPointsSetWithFourVertices()
{
	/*
	* 选择顶点的方式：
	* 2--------------------------------->3
	* ↑                                 ↓
	* ↑                                 ↓
	* 1<---------------------------------4
	*/
	if (m_v1_mode3.poly_id == -1 || !m_v1_mode3.v ||
		m_v2_mode3.poly_id == -1 || !m_v2_mode3.v ||
		m_v3_mode3.poly_id == -1 || !m_v3_mode3.v ||
		m_v4_mode3.poly_id == -1 || !m_v4_mode3.v) {
		QMessageBox::information(this, "提示", "顶点或多边形未设置好，请重新设置!");
		return false;
	}

	if (m_v1_mode3.poly_id == m_v2_mode3.poly_id &&
		m_v1_mode3.v == m_v2_mode3.v) {
		QMessageBox::information(this, "提示", "顶点1与顶点2不能重合!");
		return false;
	}

	if (m_v3_mode3.poly_id == m_v4_mode3.poly_id &&
		m_v3_mode3.v == m_v4_mode3.v) {
		QMessageBox::information(this, "提示", "顶点3与顶点4不能重合!");
		return false;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	m_points_set_for_tri.clear();   // 清理孔洞边界点集数据结构

									// step1: 线段[v2----------->v3]
	Vertex2 *cur_v, *next_v;
	POLYGON2 *cur_ply, *next_ply;
	cur_v = m_v2_mode3.v;
	cur_ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[m_v2_mode3.poly_id];
	m_points_set_for_tri.push_back(cur_v);

	while (true) {
		// 获取下一个顶点
		if (!findNextPoint(cur_v, cur_ply, next_v, next_ply)) {
			m_points_set_for_tri.clear();
			cerr << "获取下一个顶点失败" << endl;
			return false;
		}

		m_points_set_for_tri.push_back(next_v);

		// 找到线段终点v3
		if (next_v == m_v3_mode3.v) {
			break;
		}

		// 更新顶点与多边形
		cur_v = next_v;
		cur_ply = next_ply;
	}

	// step2: 线段[v4----------->v1]
	//is_first_call = true;   //!!!!!!!!!!!!!!!
	cur_v = m_v4_mode3.v;
	cur_ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[m_v4_mode3.poly_id];
	m_points_set_for_tri.push_back(cur_v);

	while (true) {
		// 获取下一个顶点
		if (!findNextPoint(cur_v, cur_ply, next_v, next_ply)) {
			m_points_set_for_tri.clear();
			cerr << "获取下一个顶点失败" << endl;
			return false;
		}

		m_points_set_for_tri.push_back(next_v);

		// 找到线段终点v1
		if (next_v == m_v1_mode3.v) {
			break;
		}

		// 更新顶点与多边形
		cur_v = next_v;
		cur_ply = next_ply;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 显示孔洞边界
	char buf[64];
	string name;

	// 清空之前的多边形
	int max_poly_size = 1024;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}

	pcl::PointXYZ p1, p2;
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		int i_next = i == m_points_set_for_tri.size() - 1 ? 0 : i + 1;

		p1.x = m_points_set_for_tri[i]->pos[0];
		p1.y = m_points_set_for_tri[i]->pos[1];
		p1.z = m_points_set_for_tri[i]->pos[2];

		p2.x = m_points_set_for_tri[i_next]->pos[0];
		p2.y = m_points_set_for_tri[i_next]->pos[1];
		p2.z = m_points_set_for_tri[i_next]->pos[2];

		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1.0, 0, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, name);
	}
}

// 模式2： 利用一个起点与一个终点生成孔洞边界点集
bool HoleFillingInteractiveOperationDialog::genPointsSetWithTwoVertices()
{
	if (m_v_begin_mode2.poly_id == -1 || !m_v_begin_mode2.v ||
		m_v_end_mode2.poly_id == -1 || !m_v_end_mode2.v) {
		QMessageBox::information(this, "提示", "顶点或多边形未设置好，请重新设置!");
		return false;
	}

	if (m_v_begin_mode2.poly_id == m_v_end_mode2.poly_id &&
		m_v_begin_mode2.v == m_v_end_mode2.v) {
		QMessageBox::information(this, "提示", "起点与终点不能重合，请重新设置!");
		return false;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	m_points_set_for_tri.clear();   // 清理孔洞边界点集数据结构
	m_points_set_for_tri.push_back(m_v_begin_mode2.v);
	Vertex2 *cur_v, *next_v;
	POLYGON2 *cur_ply, *next_ply;
	cur_v = m_v_begin_mode2.v;
	cur_ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[m_v_begin_mode2.poly_id];

	while (true) {
		// 获取下一个顶点
		if (!findNextPoint(cur_v, cur_ply, next_v, next_ply)) {
			m_points_set_for_tri.clear();
			cerr << "获取下一个顶点失败" << endl;
			return false;
		}

		m_points_set_for_tri.push_back(next_v);

		if (next_v == m_v_end_mode2.v) {
			break;
		}

		// 更新顶点与多边形
		cur_v = next_v;
		cur_ply = next_ply;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 显示孔洞边界
	char buf[64];
	string name;

	// 清空之前的多边形
	int max_poly_size = 1024;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}

	pcl::PointXYZ p1, p2;
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		int i_next = i == m_points_set_for_tri.size() - 1 ? 0 : i + 1;

		p1.x = m_points_set_for_tri[i]->pos[0];
		p1.y = m_points_set_for_tri[i]->pos[1];
		p1.z = m_points_set_for_tri[i]->pos[2];

		p2.x = m_points_set_for_tri[i_next]->pos[0];
		p2.y = m_points_set_for_tri[i_next]->pos[1];
		p2.z = m_points_set_for_tri[i_next]->pos[2];

		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1.0, 0, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, name);
	}
}

// 模式1：仅利用一个顶点生成孔洞边界点集
bool HoleFillingInteractiveOperationDialog::genPointsSetWithOneVertex()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return false;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return false;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return false;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return false;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	m_points_set_for_tri.clear();
	m_points_set_for_tri.push_back(m_vertex1);

	POLYGON2* ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	POLYGON2* next_ply = NULL;
	Vertex2* v = m_vertex1;
	Vertex2* next_v = NULL;

	while (true) {
		if (!findNextPoint(v, ply, next_v, next_ply))
			break;

		if (next_v != m_vertex1) {
			m_points_set_for_tri.push_back(next_v);
		}
		else {
			break;
		}

		// 更新顶点与多边形
		v = next_v;
		ply = next_ply;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 显示孔洞边界
	char buf[64];
	string name;

	// 清空之前的多边形
	int max_poly_size = 1024;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}

	pcl::PointXYZ p1, p2;
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		int i_next = i == m_points_set_for_tri.size() - 1 ? 0 : i + 1;

		p1.x = m_points_set_for_tri[i]->pos[0];
		p1.y = m_points_set_for_tri[i]->pos[1];
		p1.z = m_points_set_for_tri[i]->pos[2];

		p2.x = m_points_set_for_tri[i_next]->pos[0];
		p2.y = m_points_set_for_tri[i_next]->pos[1];
		p2.z = m_points_set_for_tri[i_next]->pos[2];

		itoa(i, buf, 16);
		name = "hb";    // hole border
		name.append(buf);
		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1.0, 0, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, name);
	}
}

// Input: v: 当前顶点；poly：v的隶属多边形
// Output: next_v: 下一个顶点；next_poly：next_v的隶属多边形
// 应用该方法的前提是需要保证目标多边形的v的前驱节点务必在边界上！
bool HoleFillingInteractiveOperationDialog::findNextPoint(Vertex2* v, POLYGON2* poly, Vertex2* &next_v, POLYGON2* &next_poly)
{
	/////////////////////////////////保证v是poly中的v///////////////////////////////////////////////////
	if (v->relevant_polys.size() == 0) {
		cerr << "in HoleFillingInteractiveOperationDialog::findNextPoint: v->relevant_polys.size() == 0" << endl;
		next_v = NULL;
		next_poly = NULL;
		return false;
	}

	if (!v || !poly) {
		cerr << "in HoleFillingInteractiveOperationDialog::findNextPoint: v==NULL or poly==NULL" << endl;
		next_v = NULL;
		next_poly = NULL;
		return false;
	}

	int poly_index_in_v = poly->getPolyIndex(v);
	if (poly_index_in_v == -1) {
		cerr << "in HoleFillingInteractiveOperationDialog::findNextPoint: v的隶属多边形不是当前多边形poly" << endl;
		return false;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	/***********************************************************************************************/
	// 只关联当前多边形：返回下一个顶点
	if (v->relevant_polys.size() == 1) {
		next_v = (Vertex2*)v->relevant_polys[0].next_point;
		next_poly = poly;
		return true;
	}
	/***********************************************************************************************/
	// 关联两个多边形：首先切换至另一个多边形，然后返回下一个顶点以及另一个多边形的位置
	if (v->relevant_polys.size() == 2) {
		// 切换到另一个多边形
		int the_other_poly_index_in_v = 1 - poly_index_in_v;
		next_v = (Vertex2*)v->relevant_polys[the_other_poly_index_in_v].next_point;
		next_poly = (POLYGON2*)v->relevant_polys[the_other_poly_index_in_v].poly;
		return true;
	}
	/***********************************************************************************************/
	// 关联多个多边形
	vector<int> single_relations;
	for (int i = 0; i < v->relevant_polys.size(); ++i) {
		if (i == poly_index_in_v) continue;

		// 当前多边形:cur_ply
		POLYGON2* cur_ply = (POLYGON2*)v->relevant_polys[i].poly;
		// 当前多边形的当前顶点是v，下一个顶点是next_v
		Vertex2* next_v = (Vertex2*)v->relevant_polys[i].next_point;

		// 检查下一个顶点在别的多边形(非cur_ply)中是否有指向v的关系
		bool is_bi_relation = false;
		for (int j = 0; j < next_v->relevant_polys.size(); ++j) {
			POLYGON2* pPoly = (POLYGON2*)next_v->relevant_polys[j].poly;
			if (pPoly == cur_ply) continue;

			if ((Vertex2*)next_v->relevant_polys[j].next_point == v) {
				is_bi_relation = true;
				break;
			}
		}

		if (!is_bi_relation) {    // 保存单向的关联关系
			single_relations.push_back(i);
		}
	}

	// 如果只有一个，那就是你了
	if (single_relations.size() == 1) {
		next_v = (Vertex2*)v->relevant_polys[single_relations[0]].next_point;
		next_poly = (POLYGON2*)v->relevant_polys[single_relations[0]].poly;
		return true;
	}

	// 如果有多个，需要找选择角度最小的那一个
	// 构造局部法向量
	Eigen::Vector3d local_normal;
	local_normal.setZero();
	for (int i = 0; i < v->relevant_polys.size(); ++i) {
		local_normal[0] += ((POLYGON2*)v->relevant_polys[i].poly)->coeff[0];
		local_normal[1] += ((POLYGON2*)v->relevant_polys[i].poly)->coeff[1];
		local_normal[2] += ((POLYGON2*)v->relevant_polys[i].poly)->coeff[2];
	}
	local_normal /= v->relevant_polys.size();
	local_normal.normalize();

	double d = v->pos.dot(local_normal) * -1;
	Eigen::Vector4d local_coeff;
	local_coeff << local_normal[0],
		local_normal[1],
		local_normal[2],
		d;

	// 构造向量
	vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> e;
	Eigen::Vector3d proj_pos;
	getProjPointOnPlane(((Vertex2*)v->relevant_polys[poly_index_in_v].prev_point)->pos, local_coeff, proj_pos);
	Eigen::Vector3d dir;
	dir = proj_pos - v->pos;
	dir.normalize();
	e.push_back(dir);
	for (int i = 0; i < single_relations.size(); ++i) {
		getProjPointOnPlane(((Vertex2*)v->relevant_polys[single_relations[i]].next_point)->pos, local_coeff, proj_pos);
		dir = proj_pos - v->pos;
		dir.normalize();
		e.push_back(dir);
	}

	int index = -1;
	float min_angle = FLT_MAX;
	for (int i = 1; i < e.size(); ++i) {
		double angle = 0;
		if (e[0].cross(e[i]).dot(local_normal) > 0) {
			angle = acos(e[0].dot(e[i]));
		}
		else {
			angle = 2.0 * 3.141592653589793238 - acos(e[0].dot(e[i]));
		}
		if (angle < min_angle) {
			index = i - 1;
			min_angle = angle;
		}
	}

	next_v = (Vertex2*)v->relevant_polys[single_relations[index]].next_point;
	next_poly = (POLYGON2*)v->relevant_polys[single_relations[index]].poly;
	return true;
}

// 点向平面的投影
void HoleFillingInteractiveOperationDialog::getProjPointOnPlane(Eigen::Vector3d &v_pos, Eigen::Vector4d &coeff, Eigen::Vector3d &proj_pos)
{
	if (coeff.norm() < 0.0000001) {
		cerr << "in HoleFillingInteractiveOperationDialog::getProjPointOnPlane: coeff.norm() < 0.0000001" << endl;
		cerr << "平面方程可能没有设置好" << endl;
		proj_pos.setZero();
		return;
	}
	// 平面法向量
	Eigen::Vector3d normal;
	normal << coeff[0], coeff[1], coeff[2];

	// 首先计算点与平面的符号距离
	double dist_v2plane = v_pos[0] * coeff[0] +
		v_pos[1] * coeff[1] +
		v_pos[2] * coeff[2] +
		coeff[3];
	// 计算投影点的位置
	proj_pos = v_pos - dist_v2plane * normal;
}

// e2顺时针旋转到e1转动的角度
double HoleFillingInteractiveOperationDialog::getRotationAngle(Eigen::Vector3d &e1, Eigen::Vector3d &e2, Eigen::Vector3d &normal)
{
	e1.normalize();
	e2.normalize();

	if (e1.cross(e2).dot(normal) >= 0) {
		return acos(e1.dot(e2));
	}
	else {
		double PI = 3.141592653589793;
		return 2.0*PI - acos(e1.dot(e2));
	}
}

// 测试功能按钮：找到当前顶点的下一个顶点
void HoleFillingInteractiveOperationDialog::on_pushButton_test_findNextPoint_clicked()
{
	/*
	* test code: findNextPoint
	*/
	Vertex2* v = m_vertex1;
	POLYGON2* poly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[this->selected_poly_id];
	Vertex2* next_v = NULL;
	POLYGON2* next_ply = NULL;

	findNextPoint(v, poly, next_v, next_ply);
	m_vertex1 = next_v;	// 更新顶点
						// 更新多边形索引
	selected_poly_id = std::find(((PCLViewer*)m_pclviewer)->m_hf.m_polygons.begin(), ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.end(), next_ply) -
		((PCLViewer*)m_pclviewer)->m_hf.m_polygons.begin();

	char buf[64];
	string name;

	// 显示
	// 清空之前的多边形
	int max_poly_size = 128;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "display_";
		name.append(buf);
		name.append("_l");
		//cout << name << endl;
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}
	// cout << "//////////////////////////////////////////////////////////////////////////////////" << endl;
	Vertex2 *pv2, *pv2_next;
	int poly_index;
	pcl::PointXYZ p1, p2;
	for (int j = 0; j < next_ply->getSize(); ++j) {
		if (j == 0) {
			pv2 = next_ply->start_point;
		}
		else {
			poly_index = next_ply->getPolyIndex(pv2);
			pv2 = (Vertex2 *)pv2->relevant_polys[poly_index].next_point;
		}

		itoa(j, buf, 16);
		name = "display_";
		name.append(buf);
		name.append("_l");
		//cout << name << endl;
		p1.x = pv2->pos[0];
		p1.y = pv2->pos[1];
		p1.z = pv2->pos[2];

		poly_index = next_ply->getPolyIndex(pv2);
		pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
		p2.x = pv2_next->pos[0];
		p2.y = pv2_next->pos[1];
		p2.z = pv2_next->pos[2];

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 0.3, 0.8, 0.7, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
	}

	((PCLViewer*)m_pclviewer)->viewer->removeShape("vertex1");
	m_vertex1 = next_v;
	pcl::PointXYZ point;
	point.x = m_vertex1->pos[0];
	point.y = m_vertex1->pos[1];
	point.z = m_vertex1->pos[2];
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	cloud->push_back(point);
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "vertex1");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertex1");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "vertex1");
}
///////////////////////////////////////////mode2///////////////////////////////////////////////////
// 选择两个顶点构造封闭边界
void HoleFillingInteractiveOperationDialog::on_pushButton_select_two_vertices_clicked()
{
	ui->radioButton_select_two_vertices->setChecked(true);
	m_v_end_mode2.poly_id = -1;
	m_v_end_mode2.v = NULL;
	m_v_begin_mode2.poly_id = -1;
	m_v_begin_mode2.v = NULL;
	ui->lineEdit_two_vertices_1->setText("");
	ui->lineEdit_two_vertices_2->setText("");


	this->m_v_begin_mode2.poly_id = -1;
	this->m_v_begin_mode2.v = NULL;
	this->m_v_end_mode2.poly_id = -1;
	this->m_v_end_mode2.v = NULL;
	this->m_v1_mode3.poly_id = -1;
	this->m_v1_mode3.v = NULL;
	this->m_v2_mode3.poly_id = -1;
	this->m_v2_mode3.v = NULL;
	this->m_v3_mode3.poly_id = -1;
	this->m_v3_mode3.v = NULL;
	this->m_v4_mode3.poly_id = -1;
	this->m_v4_mode3.v = NULL;
}

// 选择第一个顶点，第一个顶点是顺时针旋转的终点
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v1_for_mode2_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v_end_mode2.poly_id = selected_poly_id;
	m_v_end_mode2.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_two_vertices_1->setText(str);
}

// 选择第二个顶点，第二个顶点是顺时针选择的起点
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v2_for_mode2_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v_begin_mode2.poly_id = selected_poly_id;
	m_v_begin_mode2.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_two_vertices_2->setText(str);
}
//////////////////////////////////////////mode3////////////////////////////////////////////////////////
// 选择四个顶点构造封闭边界（适用于连接两个分离的多边形）
void HoleFillingInteractiveOperationDialog::on_pushButton_select_four_vertices_clicked()
{
	ui->radioButton_select_four_vertices->setChecked(true);
	m_v1_mode3.poly_id = m_v2_mode3.poly_id = m_v3_mode3.poly_id = m_v4_mode3.poly_id = -1;
	m_v1_mode3.v = m_v2_mode3.v = m_v3_mode3.v = m_v4_mode3.v = NULL;

	ui->lineEdit_three_vertices_1->setText("");
	ui->lineEdit_three_vertices_2->setText("");
	ui->lineEdit_three_vertices_3->setText("");
	ui->lineEdit_three_vertices_4->setText("");


	this->m_v_begin_mode2.poly_id = -1;
	this->m_v_begin_mode2.v = NULL;
	this->m_v_end_mode2.poly_id = -1;
	this->m_v_end_mode2.v = NULL;
	this->m_v1_mode3.poly_id = -1;
	this->m_v1_mode3.v = NULL;
	this->m_v2_mode3.poly_id = -1;
	this->m_v2_mode3.v = NULL;
	this->m_v3_mode3.poly_id = -1;
	this->m_v3_mode3.v = NULL;
	this->m_v4_mode3.poly_id = -1;
	this->m_v4_mode3.v = NULL;
}
/*
* 选择顶点的方式：
* 2--------------------------------->3
* ↑                                 ↓
* ↑                                 ↓
* 1<---------------------------------4
*/
// 顶点1
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v1_for_mode3_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v1_mode3.poly_id = selected_poly_id;
	m_v1_mode3.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_three_vertices_1->setText(str);
}
// 顶点2
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v2_for_mode3_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v2_mode3.poly_id = selected_poly_id;
	m_v2_mode3.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_three_vertices_2->setText(str);
}
// 顶点3
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v3_for_mode3_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v3_mode3.poly_id = selected_poly_id;
	m_v3_mode3.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_three_vertices_3->setText(str);
}
// 顶点4
void HoleFillingInteractiveOperationDialog::on_pushButton_sel_v4_for_mode3_clicked()
{
	if (m_camera_pos.norm() < 0.0000001) {
		QMessageBox::information(this, "提示", "为局部投影还需设置合适的摄像机位置！");
		return;
	}

	if (selected_poly_id == -1) {
		QMessageBox::information(this, "提示", "还没选择多边形: selected_poly_id == -1");
		return;
	}

	if (!m_vertex1) {
		QMessageBox::information(this, "提示", "还没选择顶点: m_vertex1 == NULL");
		return;
	}

	int poly_index_in_v = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) {
		QMessageBox::information(this, "提示", "当前顶点不隶属指定多边形!");
		return;
	}

	m_v4_mode3.poly_id = selected_poly_id;
	m_v4_mode3.v = m_vertex1;

	QString str;
	str = "poly_id = ";
	str.append(QString::number(selected_poly_id));
	str.append("; point = [");
	str.append(QString::number(m_vertex1->pos[0]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[1]));
	str.append(" ");
	str.append(QString::number(m_vertex1->pos[2]));
	str.append("]");
	ui->lineEdit_three_vertices_4->setText(str);
}

// 利用三个顶点构造新的三角形多边形
void HoleFillingInteractiveOperationDialog::genNewTriPly(Vertex2* v0, Vertex2* v1, Vertex2* v2, Eigen::Vector3d &normal)
{
	Eigen::Vector3d e01, e02;
	e01 = v1->proj_pos - v0->proj_pos;
	e02 = v2->proj_pos - v0->proj_pos;

	bool is012 = e01.cross(e02).dot(normal) > 0;

	if (!is012) {
		Vertex2* tmp = v1;
		v1 = v2;
		v2 = tmp;
	}
	// 用真实位置而不是投影位置计算新三角形的法向量
	e01 = v1->pos - v0->pos;
	e02 = v2->pos - v0->pos;
	//Eigen::Vector3d cross_product = e01.cross(e02);
	Eigen::Vector3d ply_normal = e01.cross(e02);
	ply_normal.normalize();
	double d = v0->pos.dot(ply_normal) * -1;

	POLYGON2* ply = new POLYGON2();
	ply->coeff << ply_normal[0], ply_normal[1], ply_normal[2], d;

	ply->insertFirstVertex2(v0);
	ply->insertMiddleVertex2(v1);
	ply->insertLastVertex2(v2);

	Triangle2 tri;
	tri.v[0] = v0;
	tri.v[1] = v1;
	tri.v[2] = v2;
	ply->triangles.push_back(tri);

	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);
}

// 获取点到线段的距离
double HoleFillingInteractiveOperationDialog::getDistP2LS(Vertex2* l_begin, Vertex2* l_end, Vertex2* v)
{
	Eigen::Vector3d p_base, line_dir;
	p_base = l_begin->proj_pos;
	line_dir = l_end->proj_pos - l_begin->proj_pos;
	line_dir.normalize();

	Eigen::Vector3d delta;
	delta = p_base - v->proj_pos;
	double lambda = -1.0 * line_dir.dot(delta);

	Eigen::Vector3d p_proj;
	p_proj = p_base + lambda * line_dir;

	if ((p_proj - l_begin->proj_pos).dot(line_dir) < 0 ||	// 线段外侧
		(p_proj - l_end->proj_pos).dot(line_dir) > 0)		// 线段外侧
	{
		double dist1 = (v->proj_pos - l_begin->proj_pos).norm();
		double dist2 = (v->proj_pos - l_end->proj_pos).norm();
		if (dist1 < dist2)
			return dist1;
		else
			return dist2;
	}
	else {	// 线段上
		return (v->proj_pos - p_proj).norm();
	}
}

// 获取指定线段到某个点集中每个顶点的距离，并按照从小到大的顺序进行排列
void HoleFillingInteractiveOperationDialog::getDistInfo(Vertex2* l_begin, Vertex2 *l_end, vector<Vertex2*> &points, vector<int> &indices, vector<float> &dists)
{
	indices.clear();
	dists.clear();

	vector<double> d(points.size(), 0);
	for (int i = 0; i < points.size(); ++i) {
		d[i] = getDistP2LS(l_begin, l_end, points[i]);
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ p;
	cloud->reserve(d.size());
	for (int i = 0; i < d.size(); ++i) {
		p.x = d[i];
		p.y = p.z = 0;
		cloud->push_back(p);
	}
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(cloud);
	p.x = p.y = p.z = 0;
	kdtree.nearestKSearch(p, d.size(), indices, dists);
}

// 三角形合法性检查
// 除了边界上，内部不能有其它顶点
// 除了平行，线段不能与其它线段相交
// 必须落在m_points_set_for_tri围成的封闭多边形的内部
bool HoleFillingInteractiveOperationDialog::isValidTriangle(Vertex2* l_begin, Vertex2* l_end, Vertex2* v,
	vector<Vertex2*> &m_points_set_for_tri, int poly_start_index,
	Eigen::Vector3d &normal)
{
	// 除了三角形的三个顶点，其余顶点要么在线段上要么在其外部，不能位于三角形内部
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		if (m_points_set_for_tri[i] == l_begin) continue;
		if (m_points_set_for_tri[i] == l_end) continue;
		if (m_points_set_for_tri[i] == v) continue;

		if (isVertexInTriangle(l_begin, l_end, v, m_points_set_for_tri[i]) &&
			!isVertexOnLingSeg(l_begin, l_end, m_points_set_for_tri[i]) &&
			!isVertexOnLingSeg(l_begin, v, m_points_set_for_tri[i]) &&
			!isVertexOnLingSeg(l_end, v, m_points_set_for_tri[i])) {
			return false;
		}
	}

	// 除了平行，线段不能与其它线段相交
	// 其它线段包括：m_points_set_for_tri中点集构成的封闭边界以及新生成的三角形的边界
	// 首先检查m_points_set_for_tri中点集构成的封闭边界
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		int i_next = i == m_points_set_for_tri.size() - 1 ? 0 : i + 1;

		if (!isTwoLineSegsParal(l_begin, l_end, m_points_set_for_tri[i], m_points_set_for_tri[i_next]) &&
			isTwoLineSegsIntersect(l_begin, l_end, m_points_set_for_tri[i], m_points_set_for_tri[i_next]))
			return false;

		if (!isTwoLineSegsParal(l_begin, v, m_points_set_for_tri[i], m_points_set_for_tri[i_next]) &&
			isTwoLineSegsIntersect(l_begin, v, m_points_set_for_tri[i], m_points_set_for_tri[i_next]))
			return false;

		if (!isTwoLineSegsParal(v, l_end, m_points_set_for_tri[i], m_points_set_for_tri[i_next]) &&
			isTwoLineSegsIntersect(v, l_end, m_points_set_for_tri[i], m_points_set_for_tri[i_next]))
			return false;
	}

	// 然后检查位于m_points_set_for_tri中点集构成的封闭边界内的其它多边形的边界
	for (int i = poly_start_index; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		POLYGON2* ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i];
		int poly_index = -1;
		Vertex2 *v1, *v2, *v3;
		v1 = ply->start_point;
		poly_index = ply->getPolyIndex(v1);
		v2 = (Vertex2*)v1->relevant_polys[poly_index].next_point;
		v3 = (Vertex2*)v1->relevant_polys[poly_index].prev_point;

		// v1,v2
		if (!isTwoLineSegsParal(l_begin, l_end, v1, v2) &&
			isTwoLineSegsIntersect(l_begin, l_end, v1, v2))
			return false;

		if (!isTwoLineSegsParal(l_begin, v, v1, v2) &&
			isTwoLineSegsIntersect(l_begin, v, v1, v2))
			return false;

		if (!isTwoLineSegsParal(v, l_end, v1, v2) &&
			isTwoLineSegsIntersect(v, l_end, v1, v2))
			return false;

		// v1,v3
		if (!isTwoLineSegsParal(l_begin, l_end, v1, v3) &&
			isTwoLineSegsIntersect(l_begin, l_end, v1, v3))
			return false;

		if (!isTwoLineSegsParal(l_begin, v, v1, v3) &&
			isTwoLineSegsIntersect(l_begin, v, v1, v3))
			return false;

		if (!isTwoLineSegsParal(v, l_end, v1, v3) &&
			isTwoLineSegsIntersect(v, l_end, v1, v3))
			return false;

		// v2,v3
		if (!isTwoLineSegsParal(l_begin, l_end, v2, v3) &&
			isTwoLineSegsIntersect(l_begin, l_end, v2, v3))
			return false;

		if (!isTwoLineSegsParal(l_begin, v, v2, v3) &&
			isTwoLineSegsIntersect(l_begin, v, v2, v3))
			return false;

		if (!isTwoLineSegsParal(v, l_end, v2, v3) &&
			isTwoLineSegsIntersect(v, l_end, v2, v3))
			return false;
	}

	// 最后检查v是否位于l_begin, l_end的左侧，如果在左侧则意味着构成的三角形在封闭点集的外部
	Eigen::Vector3d e1, e2;
	e1 = l_end->proj_pos - l_begin->proj_pos;
	e2 = v->proj_pos - l_begin->proj_pos;
	if (e1.cross(e2).dot(normal) > 0) {
		return false;
	}

	return true;
}

// 判断顶点是否在线段上
bool HoleFillingInteractiveOperationDialog::isVertexOnLingSeg(Vertex2* l_begin, Vertex2* l_end, Vertex2* v)
{
	if ((v->proj_pos - l_begin->proj_pos).norm() < 0.0001 ||
		(v->proj_pos - l_end->proj_pos).norm() < 0.0001) {
		return true;
	}

	Eigen::Vector3d p_base, line_dir;
	p_base = l_begin->proj_pos;
	line_dir = l_end->proj_pos - l_begin->proj_pos;
	line_dir.normalize();

	Eigen::Vector3d delta;
	delta = p_base - v->proj_pos;
	double lambda = -1.0 * line_dir.dot(delta);

	Eigen::Vector3d p_proj;
	p_proj = p_base + lambda * line_dir;

	if ((p_proj - l_begin->proj_pos).dot(line_dir) < 0 ||	// 线段外侧
		(p_proj - l_end->proj_pos).dot(line_dir) > 0)		// 线段外侧
	{
		return false;
	}
	else {	// 线段上
		return (v->proj_pos - p_proj).norm() <= 0.0001;
	}
}

// 判断顶点是否在三角形内部
bool HoleFillingInteractiveOperationDialog::isVertexInTriangle(Vertex2* l_begin, Vertex2* l_end, Vertex2* v, Vertex2* test_point)
{
	Eigen::Vector3d v1, v2;
	v1 = l_begin->proj_pos - v->proj_pos;
	v2 = l_end->proj_pos - v->proj_pos;

	// 三角形面积
	double S = abs((v1.cross(v2)).norm())*0.5;

	Eigen::Vector3d t1, t2, t3;
	t1 = l_begin->proj_pos - test_point->proj_pos;
	t2 = l_end->proj_pos - test_point->proj_pos;
	t3 = v->proj_pos - test_point->proj_pos;
	// 三个小三角形的面积
	double S1 = abs((t1.cross(t2)).norm())*0.5;
	double S2 = abs((t1.cross(t3)).norm())*0.5;
	double S3 = abs((t2.cross(t3)).norm())*0.5;

	return abs(S1 + S2 + S3 - S) <= 0.00001;
}

// 两个线段是否平行
bool HoleFillingInteractiveOperationDialog::isTwoLineSegsParal(Vertex2* l1_begin, Vertex2* l1_end, Vertex2* l2_begin, Vertex2* l2_end)
{
	Eigen::Vector3d e1, e2;
	e1 = l1_end->proj_pos - l1_begin->proj_pos;
	e2 = l2_end->proj_pos - l2_begin->proj_pos;

	return (e1.cross(e2)).norm() <= 0.0001;
}

// 两个线段是否相交(首位相连不算相交)
bool HoleFillingInteractiveOperationDialog::isTwoLineSegsIntersect(Vertex2* l1_begin, Vertex2* l1_end, Vertex2* l2_begin, Vertex2* l2_end)
{
	// 注意：若两个线段的端点重合则不认为是相交!
	/*if((l1_begin->pos-l2_begin->pos).norm() < 0.0000000001) return false;
	if((l1_begin->pos-l2_end->pos).norm() < 0.0000000001) return false;
	if((l1_end->pos  -l2_begin->pos).norm() < 0.0000000001) return false;
	if((l1_end->pos  -l2_end->pos).norm() < 0.0000000001) return false;*/

	Eigen::Vector3d p1, p2, p3, p4, p12, p34, p31;
	p1 = l1_begin->proj_pos;
	p2 = l1_end->proj_pos;
	p3 = l2_begin->proj_pos;
	p4 = l2_end->proj_pos;
	p12 = p2 - p1;
	p34 = p4 - p3;
	p31 = p1 - p3;

	Eigen::Matrix2d A;
	Eigen::Vector2d x, b;
	A << p12.dot(p12), -1 * p12.dot(p34),
		-1 * p12.dot(p34), p34.dot(p34);
	b << -1 * p31.dot(p12), p31.dot(p34);
	x = A.colPivHouseholderQr().solve(b);
	double alpha1, alpha2;
	alpha1 = x[0];
	alpha2 = x[1];

	Eigen::Vector3d p_l1, p_l2;
	p_l1 = p1 + alpha1*p12;
	p_l2 = p3 + alpha2*p34;

	double dist_p2p = (p_l1 - p_l2).norm();

	// alpha1在0，1之间意味着离着直线2最近的点在直线1的内部，同理对alpha2
	// 两点的距离充分小意味着两条线段确实是相交的
	if (alpha1 >= 0.00000001 && alpha1 <= 0.999999999 &&
		alpha2 >= 0.00000001 && alpha2 <= 0.999999999 &&
		dist_p2p < 0.0002) {
		return true;
	}
	else {
		return false;
	}
}

// 测试按钮
void HoleFillingInteractiveOperationDialog::on_pushButton_test_clicked()
{
	// 利用点集重心构造局部投影平面
	Eigen::Vector3d centroid;
	centroid << 0, 0, 0;
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		centroid += m_points_set_for_tri[i]->pos;
	}
	centroid = 1.0 / m_points_set_for_tri.size() * centroid;

	Eigen::Vector3d normal;	// 平面法向量
	normal = m_camera_pos - centroid;
	normal.normalize();
	double d;
	d = centroid.dot(normal) * -1;
	Eigen::Vector4d plane_coeff;	// 局部平面方程
	plane_coeff << normal[0], normal[1], normal[2], d;

	// 将m_points_set_for_tri中所有的点都投影到二维平面上
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		this->getProjPointOnPlane(m_points_set_for_tri[i]->pos, plane_coeff, m_points_set_for_tri[i]->proj_pos);
	}

	//((PCLViewer*)m_pclviewer)

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ p;
	p.x = centroid[0];
	p.y = centroid[1];
	p.z = centroid[2];

	cloud->push_back(p);

	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		p.x = m_points_set_for_tri[i]->proj_pos[0];
		p.y = m_points_set_for_tri[i]->proj_pos[1];
		p.z = m_points_set_for_tri[i]->proj_pos[2];
		cloud->push_back(p);
		cout << "p[" << i << "] = \n" << p << endl;
	}

	((PCLViewer*)m_pclviewer)->viewer->removeAllPointClouds();
	((PCLViewer*)m_pclviewer)->viewer->removeAllShapes();
	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("plane");
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "plane");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.6, 0.8, 0.9, "plane");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "plane");

	EIGEN_ALIGN16 Eigen::Matrix3f cov_mat;
	Eigen::Vector4f xyz_centroid;

	pcl::computeMeanAndCovarianceMatrix(*cloud, cov_mat, xyz_centroid);

	EIGEN_ALIGN16 Eigen::Vector3f eigen_values;
	EIGEN_ALIGN16 Eigen::Matrix3f eigen_vectors;
	pcl::eigen33(cov_mat, eigen_vectors, eigen_values);
	cout << "normal = \n" << normal << endl;
	cout << "eigen_vectors = \n" << eigen_vectors << endl;
	cout << "eigen_values = \n" << eigen_values << endl;

	Eigen::Vector3d p1, p2, p3, p4, p12, p34, p31;
	p1 << cloud->points[4].x, cloud->points[4].y, cloud->points[4].z;
	p2 << cloud->points[2].x, cloud->points[2].y, cloud->points[2].z;
	p3 << cloud->points[3].x, cloud->points[3].y, cloud->points[3].z;
	p4 << cloud->points[1].x, cloud->points[1].y, cloud->points[1].z;
	p12 = p2 - p1;
	p34 = p4 - p3;
	p31 = p1 - p3;

	Eigen::Matrix2d A;
	Eigen::Vector2d x, b;
	A << p12.dot(p12), -1 * p12.dot(p34),
		-1 * p12.dot(p34), p34.dot(p34);
	b << -1 * p31.dot(p12), p31.dot(p34);
	x = A.colPivHouseholderQr().solve(b);
	double alpha1, alpha2;
	alpha1 = x[0];
	alpha2 = x[1];
	cout << "alpha1 = " << alpha1 << endl;
	cout << "alpha2 = " << alpha2 << endl;

	Eigen::Vector3d p_l1, p_l2;
	p_l1 = p1 + alpha1*p12;
	p_l2 = p3 + alpha2*p34;
	cout << "交点1 = \n" << p_l1 << endl;
	cout << "交点2 = \n" << p_l2 << endl;

	double dist_p2p = (p_l1 - p_l2).norm();
	cout << "dist_p2p = " << dist_p2p << endl;
	cout << "(p_l1-p_l2) = \n" << (p_l1 - p_l2) << endl;

	// 显示p1,p2
	((PCLViewer*)m_pclviewer)->viewer->addLine(cloud->points[4], cloud->points[2], 1, 1, 0, "1");
	// 显示p3,p4
	((PCLViewer*)m_pclviewer)->viewer->addLine(cloud->points[3], cloud->points[1], 1, 1, 0, "2");

	pcl::PointXYZ pl1, pl2;
	pl1.x = p_l1[0];
	pl1.y = p_l1[1];
	pl1.z = p_l1[2];
	pl2.x = p_l2[0];
	pl2.y = p_l2[1];
	pl2.z = p_l2[2];

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, cloud2;
	cloud1.reset(new pcl::PointCloud<pcl::PointXYZ>);
	cloud2.reset(new pcl::PointCloud<pcl::PointXYZ>);
	cloud1->push_back(pl1);
	cloud2->push_back(pl2);

	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud1, "plane1");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "plane1");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "plane1");

	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud2, "plane2");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 0, "plane2");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "plane2");/**/
}

// 设置摄像机位置
void HoleFillingInteractiveOperationDialog::on_pushButton_set_camera_pos_clicked()
{
	//m_camera_pos << cameras[0].pos[0], cameras[0].pos[1], cameras[0].pos[2];
	cout << "m_camera_pos.x = ";
	cin >> m_camera_pos[0];
	cout << "m_camera_pos.y = ";
	cin >> m_camera_pos[1];
	cout << "m_camera_pos.z = ";
	cin >> m_camera_pos[2];

	QString str;
	str.append("[");
	str.append(QString::number(m_camera_pos[0]));
	str.append(" ");
	str.append(QString::number(m_camera_pos[1]));
	str.append(" ");
	str.append(QString::number(m_camera_pos[2]));
	str.append("]");

	ui->lineEdit_camera_pos->setText(str);
}

// 获取当前多边形数量
void HoleFillingInteractiveOperationDialog::on_pushButton_get_poly_num_clicked()
{
	cout << "当前多边形数量：" << ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size() << "个" << endl;
}

// 获取当前顶点数量
void HoleFillingInteractiveOperationDialog::on_pushButton_get_ver_num_clicked()
{
	cout << "当前顶点数量：" << ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size() << "个" << endl;
}

// 取消显示多边形线框
void HoleFillingInteractiveOperationDialog::on_pushButton_cancel_poly_wireframe_display_clicked()
{
	char buf[64];
	string name;
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		for (int j = 0; j < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->getSize(); ++j) {
			itoa(i, buf, 16);
			name = buf;
			name.append("_");
			itoa(j, buf, 16);
			name.append(buf);
			name.append("_l");
			((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
		}
	}
}

// 显示多边形线框
void HoleFillingInteractiveOperationDialog::on_pushButton_display_poly_wireframe_clicked()
{
	on_pushButton_cancel_poly_wireframe_display_clicked();

	char buf[64];
	string name;

	Vertex2* pv2, *pv2_next;
	int poly_index;
	pcl::PointXYZ p1, p2;
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		for (int j = 0; j < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->getSize(); ++j) {
			if (j == 0) {
				pv2 = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->start_point;
			}
			else {
				pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
			}

			poly_index = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->getPolyIndex(pv2);
			pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;

			p1.x = pv2->pos[0];
			p1.y = pv2->pos[1];
			p1.z = pv2->pos[2];
			p2.x = pv2_next->pos[0];
			p2.y = pv2_next->pos[1];
			p2.z = pv2_next->pos[2];

			itoa(i, buf, 16);
			name = buf;
			name.append("_");
			itoa(j, buf, 16);
			name.append(buf);
			name.append("_l");
			((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 1, 0, name);
			((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 内部三角化完成后需要构造向下，向右和向左的投影边界，然后再构造向后的投影边界，最后将最后的投影边界给封装起来从而构造一个封闭的岩体三维流形
// 进行投影时首先需要指定一段边界
// 使用的数据结构与进行线段化简的数据结构相同：vector<Vertex2*> line_segs
void HoleFillingInteractiveOperationDialog::on_pushButton_gen_proj_line_segs_clicked()
{
	on_pushButton_del_lineseg_clicked();

	if (!m_vertex1 || !m_vertex2) {
		cerr << "顶点不能为空!" << endl;
		return;
	}

	if (!isVertexHasSingleConnRelation(m_vertex1)) {
		cerr << "顶点1不能选在模型内部!" << endl;
		return;
	}

	if (!isVertexHasSingleConnRelation(m_vertex2)) {
		cerr << "顶点2不能选在模型内部!" << endl;
		return;
	}

	if (m_vertex1 == m_vertex2) {
		cerr << "两个顶点不能相等!" << endl;
		return;
	}

	// 默认从m_vertex1找m_vertex2
	Vertex2 *cur_v, *next_v;

	//this->is_first_call = true;	// findNextPoint

	cur_v = m_vertex1;
	//cur_ply = pPoly;

	line_segs.clear();
	line_segs.push_back(cur_v);

	while (true) {
		if (!findNextPointVertexOnModelBorder(cur_v, next_v))
		{
			cerr << "查找下一个顶点失败!" << endl;
			line_segs.clear();
			return;
		}

		line_segs.push_back(next_v);

		if (next_v == m_vertex2)
			break;

		cur_v = next_v;
	}

	cout << "选择的线段序列中，有" << line_segs.size() << "个顶点." << endl;

	// 显示线段
	char buf[64];
	string name;
	pcl::PointXYZ p1, p2;
	for (int i = 0; i < line_segs.size() - 1; ++i) {
		p1.x = line_segs[i]->pos[0];
		p1.y = line_segs[i]->pos[1];
		p1.z = line_segs[i]->pos[2];

		p2.x = line_segs[i + 1]->pos[0];
		p2.y = line_segs[i + 1]->pos[1];
		p2.z = line_segs[i + 1]->pos[2];

		itoa(i, buf, 16);
		name = "l_";
		name.append(buf);

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 0, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
*           ↑y(垂直)
*           |
*           |
*           |------>x(水平)
*          /
*         /
*        z(前后)
*/
// 向下延伸构建底座
void HoleFillingInteractiveOperationDialog::on_pushButton_extend_to_bottom_clicked()
{
	performExtention(0);
}
// 向左延伸
void HoleFillingInteractiveOperationDialog::on_pushButton_extend_to_left_clicked()
{
	performExtention(1);
}
// 向右延伸
void HoleFillingInteractiveOperationDialog::on_pushButton_extend_to_right_clicked()
{
	performExtention(2);
}
// 向后延伸
void HoleFillingInteractiveOperationDialog::on_pushButton_extend_to_backward_clicked()
{
	performExtention(3);
}
//// 封闭模型{必须执行完向下，向左，向右，和向后的延伸操作之后才能执行!}
//void HoleFillingInteractiveOperationDialog::on_pushButton_close_the_model_clicked()
//{
//    EarClip2 ec2;

//    // 构造底面
//    if(plane_type == 0){
//        m_four_vertices[0].push_back(m_bottom_left_vertex);
//        m_four_vertices[0].push_back(m_bottom_right_vertex);
//        POLYGON2* bottom_ply = new POLYGON2();
//        bottom_ply->coeff << 0, -1, 0, -1*m_min_y;
//        for(int i = 0; i < m_four_vertices[0].size(); ++i){
//            if(i == 0){
//                bottom_ply->insertFirstVertex2(m_four_vertices[0][i]);
//            }else if(i == m_four_vertices[0].size()-1){
//                bottom_ply->insertLastVertex2(m_four_vertices[0][i]);
//            }else{
//                bottom_ply->insertMiddleVertex2(m_four_vertices[0][i]);
//            }
//        }
//        ec2.setInputPoly(bottom_ply);
//        ec2.triangulatePoly();
//        ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(bottom_ply);
//    }

//    // 构造左侧立面
//    if(plane_type == 1){
//        m_four_vertices[1].push_back(m_bottom_left_vertex);

//            POLYGON2* left_ply = new POLYGON2();
//            left_ply->coeff << -1, 0, 0, -1*m_min_x;
//            for(int i = 0; i < m_four_vertices[1].size(); ++i){
//                if(i == 0){
//                    left_ply->insertFirstVertex2(m_four_vertices[1][i]);
//                }else if(i == m_four_vertices[1].size()-1){
//                    left_ply->insertLastVertex2(m_four_vertices[1][i]);
//                }else{
//                    left_ply->insertMiddleVertex2(m_four_vertices[1][i]);
//                }
//            }
//            ec2.setInputPoly(left_ply);
//            ec2.triangulatePoly();
//            ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(left_ply);
//    }

//    // 构造右侧立面
//    if(plane_type == 2){
//        m_four_vertices[2].push_back(m_bottom_right_vertex);

//        POLYGON2* right_ply = new POLYGON2();
//        right_ply->coeff << 1, 0, 0, -1*m_max_x;
//        for(int i = 0; i < m_four_vertices[2].size(); ++i){
//            if(i == 0){
//                right_ply->insertFirstVertex2(m_four_vertices[2][i]);
//            }else if(i == m_four_vertices[2].size()-1){
//                right_ply->insertLastVertex2(m_four_vertices[2][i]);
//            }else{
//                right_ply->insertMiddleVertex2(m_four_vertices[2][i]);
//            }
//        }
//        ec2.setInputPoly(right_ply);
//        ec2.triangulatePoly();
//        ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(right_ply);
//    }

//    // 构造后部立面
//    if(plane_type == 3){
//        m_four_vertices[3].push_back(m_bottom_right_vertex);
//        m_four_vertices[3].push_back(m_bottom_left_vertex);

//        POLYGON2* backward_ply = new POLYGON2();
//        backward_ply->coeff << 0, 0, -1, -1*m_min_z;
//        for(int i = 0; i < m_four_vertices[3].size(); ++i){
//            if(i == 0){
//                backward_ply->insertFirstVertex2(m_four_vertices[3][i]);
//            }else if(i == m_four_vertices[3].size()-1){
//                backward_ply->insertLastVertex2(m_four_vertices[3][i]);
//            }else{
//                backward_ply->insertMiddleVertex2(m_four_vertices[3][i]);
//            }
//        }
//        ec2.setInputPoly(backward_ply);
//        ec2.triangulatePoly();
//        ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(backward_ply);
//    }

//    // 更新显示
//    // 设置三角面片重心点云
//    ((PCLViewer*)m_pclviewer)->m_hf.set_cloud_for_indentify_poly_id();

//    // 更新视图
//    ((PCLViewer*)m_pclviewer)->displayWireframe2();	// 线框，顶点点云，三角面片
//    ((PCLViewer*)m_pclviewer)->viewer->addPointCloud(((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id,
//                                                     "triCentroids"); // 三角面片重心
//}
// 执行函数
void HoleFillingInteractiveOperationDialog::performExtention(int extend_type)
{
	if (line_segs.size() <= 0) {
		cerr << "还未选择线段序列!" << endl;
		return;
	}

	// 线段序列沿着投影方向投影后，只准首尾相连，不允许出现线段相交和重合的情况
	// 第一步：输入延伸的空间尺度
	cout << "延伸尺度：";
	double extend_scale;
	cin >> extend_scale;
	if (extend_scale <= 0) {
		cerr << "extend_scale = " << extend_scale << endl;
		return;
	}

	// 第二步：构造投影平面
	Eigen::Vector3d p_base;             // 投影平面穿过的点
	Eigen::Vector3d proj_plane_normal;  // 投影平面法向量
	Eigen::Vector4d proj_plane_coeff;   // 投影平面参数

										// 向下延伸，需找min y(最底部的点)
	if (extend_type == 0) {
		double min_y = DBL_MAX;
		int min_y_index = -1;
		for (int i = 0; i < line_segs.size(); ++i) {
			if (line_segs[i]->pos[1] < min_y) {
				min_y = line_segs[i]->pos[1];
				min_y_index = i;
			}
		}
		//m_min_y = min_y;
		p_base = line_segs[min_y_index]->pos;
		p_base[1] -= extend_scale;
		m_min_y = p_base[1];
		proj_plane_normal << 0, -1, 0;
		double d = p_base.dot(proj_plane_normal) * -1;
		proj_plane_coeff << proj_plane_normal[0], proj_plane_normal[1], proj_plane_normal[2], d;
	}
	// 向左延伸，需要找min x(最左侧的点)
	if (extend_type == 1) {
		double min_x = DBL_MAX;
		int min_x_index = -1;
		for (int i = 0; i < line_segs.size(); ++i) {
			if (line_segs[i]->pos[0] < min_x) {
				min_x = line_segs[i]->pos[0];
				min_x_index = i;
			}
		}
		//m_min_x = min_x;
		p_base = line_segs[min_x_index]->pos;
		p_base[0] -= extend_scale;
		m_min_x = p_base[0];
		proj_plane_normal << -1, 0, 0;
		double d = p_base.dot(proj_plane_normal) * -1;
		proj_plane_coeff << proj_plane_normal[0], proj_plane_normal[1], proj_plane_normal[2], d;
	}
	// 向右延伸，需要找max x(最右侧的点)
	if (extend_type == 2) {
		double max_x = DBL_MIN;
		int max_x_index = -1;
		for (int i = 0; i < line_segs.size(); ++i) {
			if (line_segs[i]->pos[0] > max_x) {
				max_x = line_segs[i]->pos[0];
				max_x_index = i;
			}
		}
		//m_max_x = max_x;
		p_base = line_segs[max_x_index]->pos;
		p_base[0] += extend_scale;
		m_max_x = p_base[0];
		proj_plane_normal << 1, 0, 0;
		double d = p_base.dot(proj_plane_normal) * -1;
		proj_plane_coeff << proj_plane_normal[0], proj_plane_normal[1], proj_plane_normal[2], d;
	}
	// 向后延伸，需要找min z(最靠里的点)
	if (extend_type == 3) {
		double min_z = DBL_MAX;
		int min_z_index = -1;
		for (int i = 0; i < line_segs.size(); ++i) {
			if (line_segs[i]->pos[2] < min_z) {
				min_z = line_segs[i]->pos[2];
				min_z_index = i;
			}
		}
		//m_min_z = min_z;
		p_base = line_segs[min_z_index]->pos;
		p_base[2] -= extend_scale;
		m_min_z = p_base[2];
		proj_plane_normal << 0, 0, -1;
		double d = p_base.dot(proj_plane_normal) * -1;
		proj_plane_coeff << proj_plane_normal[0], proj_plane_normal[1], proj_plane_normal[2], d;
	}

	// 第三步将线段序列中所有的点都投影到投影平面上
	for (int i = 0; i < line_segs.size(); ++i) {
		getProjPointOnPlane(line_segs[i]->pos, proj_plane_coeff, line_segs[i]->proj_pos);
	}

	// 第四步检验投影线段的合法性
	int problem_i, problem_j;
	if (!isProjLineSegsValid(line_segs, problem_i, problem_j)) {
		// 显示有问题的两个线段，然后退出
		cerr << "第" << problem_i << "条线段与第" << problem_j << "条线段的投影出现问题" << endl;

		// 显示投影线段
		((PCLViewer*)m_pclviewer)->viewer->removeShape("prob_1");
		((PCLViewer*)m_pclviewer)->viewer->removeShape("prob_2");
		pcl::PointXYZ p1, p2;
		// 第一条线段
		p1.x = line_segs[problem_i]->proj_pos[0];
		p1.y = line_segs[problem_i]->proj_pos[1];
		p1.z = line_segs[problem_i]->proj_pos[2];

		p2.x = line_segs[problem_i + 1]->proj_pos[0];
		p2.y = line_segs[problem_i + 1]->proj_pos[1];
		p2.z = line_segs[problem_i + 1]->proj_pos[2];

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 0, 0, "prob_1");
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "prob_1");

		// 第二条线段
		p1.x = line_segs[problem_j]->proj_pos[0];
		p1.y = line_segs[problem_j]->proj_pos[1];
		p1.z = line_segs[problem_j]->proj_pos[2];

		p2.x = line_segs[problem_j + 1]->proj_pos[0];
		p2.y = line_segs[problem_j + 1]->proj_pos[1];
		p2.z = line_segs[problem_j + 1]->proj_pos[2];

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 0, 0, "prob_2");
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "prob_2");

		return;
	}
	else {
		cout << "投影线段的位置关系合法" << endl;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 第五步：生成四边形序列
	//EarClip2 ec2;
	int original_ply_size = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size();
	int original_ver_size = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size();
	/*
	* 有多少条线段就会在投影时生成多少个四边形
	* 线段顺序：i, i+1
	* 四边形顶点顺序:i(旧顶点pos),proj_i(新顶点proj_pos）, proj_i+1（新顶点proj_pos）, i+1（旧顶点pos）
	*/
	m_four_vertices[extend_type].clear();   // 构造封闭模型时使用
											// 线段数量
	int line_segs_size = line_segs.size() - 1;
	for (int i = 0; i < line_segs_size; ++i) {
		POLYGON2* ply = new POLYGON2();
		// 多边形法向量（以向下投影为例）
		/*
		*            i-------
		*            ↑       -------i+1
		*            |               |
		*         e1 |               |
		*            |        e2     |
		*          proj_i --------->proj_i+1
		*
		*  ∠proj_i必定为直角（i->proj_i:垂直；proj_i->proj_i+1: 水平)
		*
		*  多边形法向量n = e2×e1
		*/

		Eigen::Vector3d e1, e2, n;
		e1 = line_segs[i]->pos - line_segs[i]->proj_pos;
		e2 = line_segs[i + 1]->proj_pos - line_segs[i]->proj_pos;
		n = e2.cross(e1);
		n.normalize();
		double d = line_segs[i]->pos.dot(n)*-1;
		ply->coeff << n[0], n[1], n[2], d;  // 设置新多边形的法向量

		ply->insertFirstVertex2(line_segs[i]);
		Vertex2* v_proj_i = NULL;
		if (i == 0) {
			v_proj_i = new Vertex2();
			v_proj_i->pos = line_segs[i]->proj_pos;
			ply->insertMiddleVertex2(v_proj_i);
		}
		else {
			v_proj_i = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[
				((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size() - 1];
			ply->insertMiddleVertex2(v_proj_i);
		}

		Vertex2* v_proj_i_plus_one = new Vertex2();
		v_proj_i_plus_one->pos = line_segs[i + 1]->proj_pos;
		ply->insertMiddleVertex2(v_proj_i_plus_one);

		ply->insertLastVertex2(line_segs[i + 1]);

		//ec2.setInputPoly(ply);
		//ec2.triangulatePoly();

		// 顶点
		if (i == 0)
			((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(v_proj_i);
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(v_proj_i_plus_one);
		pcl::PointXYZ p1, p2;
		p1.x = v_proj_i->pos[0];
		p1.y = v_proj_i->pos[1];
		p1.z = v_proj_i->pos[2];
		p2.x = v_proj_i_plus_one->pos[0];
		p2.y = v_proj_i_plus_one->pos[1];
		p2.z = v_proj_i_plus_one->pos[2];
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices_cloud->push_back(p1);
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices_cloud->push_back(p2);

		// 多边形
		((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);
	}

	// 为构造封闭模型，需将新生成的线段序列保存起来original_ver_size
	for (int i = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size() - 1; i >= original_ver_size; --i) {
		m_four_vertices[extend_type].push_back(((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]);
	}
	////////////////////////////////////////////////////////////////////////////////////////
	// 第六步：更新显示
	// 设置三角面片重心点云
	((PCLViewer*)m_pclviewer)->m_hf.set_cloud_for_indentify_poly_id();

	// 更新视图
	((PCLViewer*)m_pclviewer)->displayWireframe2();	// 线框，顶点点云，三角面片
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id, "triCentroids"); // 三角面片重心
}

// 水平方向，垂直方向延伸模型边界
// 投影的线段序列是否合法
// 有问题的线段索引: problem_i, problem_j
bool HoleFillingInteractiveOperationDialog::isProjLineSegsValid(vector<Vertex2*> &line_segs, int &problem_i, int &problem_j)
{
	// 线段数量 = 顶点数量-1
	int size = line_segs.size() - 1;
	problem_i = problem_j = -1;
	for (int i = 0; i < size - 1; ++i) {

		// 检查相邻线段是否具有合法的位置关系
		if (!isAdjacentLineSegsHasValidPosRelation(line_segs, i)) {
			problem_i = i;
			problem_j = i + 1;
			return false;
		}

		for (int j = i + 1; j < size; ++j) {
			// 检查不相邻线段是否具有重合关系
			if (isBothLineSegsOverlap(line_segs, i, j)) {
				problem_i = i;
				problem_j = j;
				return false;
			}

			// 检查不相邻线段是否具有相交关系
			if (isTwoLineSegsIntersect(line_segs[i], line_segs[i + 1], line_segs[j], line_segs[j + 1])) {
				problem_i = i;
				problem_j = j;
				return false;
			}
		}
	}

	return true;
}

// 两个线段的方向向量是否平行(i,j均为线段的起点)
bool HoleFillingInteractiveOperationDialog::isBothLineSegsParal(vector<Vertex2*> &line_segs, int i, int j)
{
	if (abs(i - j) <= 1) return false; // 不考虑首尾相连的线段

	Eigen::Vector3d i_dir, j_dir;
	i_dir = line_segs[i + 1]->proj_pos - line_segs[i]->proj_pos;
	j_dir = line_segs[j + 1]->proj_pos - line_segs[j]->proj_pos;
	i_dir.normalize();
	j_dir.normalize();
	Eigen::Vector3d cross_product = i_dir.cross(j_dir);
	return cross_product.norm() <= 0.000001;
}

// 两条线段是否部分重叠
bool HoleFillingInteractiveOperationDialog::isBothLineSegsOverlap(vector<Vertex2*> &line_segs, int i, int j)
{
	if (!isBothLineSegsParal(line_segs, i, j)) return false;

	Eigen::Vector3d l1_begin, l1_end, l2_begin, l2_end;
	l1_begin = line_segs[i]->proj_pos;
	l1_end = line_segs[i + 1]->proj_pos;
	l2_begin = line_segs[j]->proj_pos;
	l2_end = line_segs[j + 1]->proj_pos;

	// l1_begin是否位于[l2_begin, l2_end]之间
	double l2_len = (l2_end - l2_begin).norm();
	double d1 = (l1_begin - l2_begin).norm();
	double d2 = (l1_begin - l2_end).norm();
	if (d1 + d2 - l2_len < 0.00001) {
		return true;
	}

	// l1_end是否位于[l2_begin, l2_end]之间
	d1 = (l1_end - l2_begin).norm();
	d2 = (l1_end - l2_end).norm();
	if (d1 + d2 - l2_len < 0.00001) {
		return true;
	}
	else {
		return false;
	}
}

// 线段i与线段i+1是否具有合法的位置关系
/*
*           [i+2]
*             **
*               **
*                 **
* [i]******************[i+1]
* ∠[i+1]不能太小!
*/
bool HoleFillingInteractiveOperationDialog::isAdjacentLineSegsHasValidPosRelation(vector<Vertex2*> &line_segs, int i)
{
	double PI = 3.14159265358979323846;
	double min_angle_i_plus_one = 45.0; // degree
	min_angle_i_plus_one = min_angle_i_plus_one / 180 * PI;

	Eigen::Vector3d e1, e2;
	e1 = line_segs[i]->proj_pos - line_segs[i + 1]->proj_pos;
	e2 = line_segs[i + 2]->proj_pos - line_segs[i + 1]->proj_pos;

	e1.normalize();
	e2.normalize();

	double cos_val = e1.dot(e2);

	return cos_val <= cos(min_angle_i_plus_one);
}

//调用本函数时默认模型内部已经不存在孔洞了
// input：v，模型边界上的顶点
// output: next_v, 模型边界上的顶点
bool HoleFillingInteractiveOperationDialog::findNextPointVertexOnModelBorder(Vertex2* v, Vertex2* &next_ver)
{
	next_ver = NULL;

	for (int i = 0; i < v->relevant_polys.size(); ++i) {
		POLYGON2* ply = (POLYGON2*)v->relevant_polys[i].poly;
		Vertex2 *next_v = (Vertex2*)v->relevant_polys[i].next_point;

		// 看next_v在其它多边形的下一个顶点是否含有v
		bool is_v_next_single = true;	// 假设v与next_v的连接是单向的
		for (int j = 0; j < next_v->relevant_polys.size(); ++j) {
			if ((Vertex2*)next_v->relevant_polys[j].next_point == v &&
				next_v->relevant_polys[j].poly != v->relevant_polys[i].poly) {
				is_v_next_single = false;
				break;
			}
		}

		if (is_v_next_single) {
			next_ver = next_v;
			return true;
		}
	}

	return next_ver;
}

// 判断顶点v是否具有单向的连接关系(可用来判断v是否为边界处的顶点)
bool HoleFillingInteractiveOperationDialog::isVertexHasSingleConnRelation(Vertex2*v)
{
	int single_relation_count = 0;
	for (int i = 0; i < v->relevant_polys.size(); ++i) {
		POLYGON2* ply = (POLYGON2*)v->relevant_polys[i].poly;
		Vertex2 *prev_v = (Vertex2*)v->relevant_polys[i].prev_point;
		Vertex2 *next_v = (Vertex2*)v->relevant_polys[i].next_point;
		// 看v在其它多边形的下一个顶点是否含有prev_v
		bool is_v_prev_single = true;	// 假设v与prev_v的连接是单向的
		for (int j = 0; j < v->relevant_polys.size(); ++j) {
			if (j == i) continue;
			if ((Vertex2*)v->relevant_polys[j].next_point == prev_v) {	// 假设不成立
				is_v_prev_single = false;
				break;
			}
		}
		if (is_v_prev_single) ++single_relation_count;

		// 看next_v在其它多边形的下一个顶点是否含有v
		bool is_v_next_single = true;	// 假设v与next_v的连接是单向的
		for (int j = 0; j < next_v->relevant_polys.size(); ++j) {
			if ((Vertex2*)next_v->relevant_polys[j].next_point == v &&
				next_v->relevant_polys[j].poly != v->relevant_polys[i].poly) {
				is_v_next_single = false;
				break;
			}
		}
		if (is_v_next_single) ++single_relation_count;
	}

	if (single_relation_count % 2 != 0) {
		cerr << "single_relation_count = " << single_relation_count << endl;
	}

	return single_relation_count > 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
//// 生成边界顶点
//void HoleFillingInteractiveOperationDialog::on_pushButton_gen_the_two_border_vertices_clicked()
//{
//	m_min_x = DBL_MAX;
//	m_max_x = DBL_MIN;
//	m_min_y = m_min_z = m_min_x;

//	for(int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size(); ++i){
//		Vertex2* v = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i];
//		if(v->pos[0] < m_min_x) m_min_x = v->pos[0];
//		if(v->pos[0] > m_max_x) m_max_x = v->pos[0];
//		if(v->pos[1] < m_min_y) m_min_y = v->pos[1];
//		if(v->pos[2] < m_min_z) m_min_z = v->pos[2];
//	}

//    // 构造左下角顶点
//    m_bottom_left_vertex = new Vertex2();
//    m_bottom_left_vertex->pos << m_min_x, m_min_y, m_min_z;
//    // 构造右下角顶点
//    m_bottom_right_vertex = new Vertex2();
//    m_bottom_right_vertex->pos << m_max_x, m_min_y, m_min_z;

//    ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(m_bottom_left_vertex);
//    ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(m_bottom_right_vertex);
//    cout << "m_bottom_left_vertex->pos = " << m_bottom_left_vertex->pos << endl;
//    cout << "m_bottom_right_vertex->pos = " << m_bottom_right_vertex->pos << endl;
//}

//// 生成底面
//void HoleFillingInteractiveOperationDialog::on_pushButton_gen_bottom_plane_clicked()
//{
//    on_pushButton_gen_proj_line_segs_clicked();
//    plane_type = 0;
//    m_four_vertices[plane_type].clear();
//	int line_segs_size = line_segs.size();
//	for(int i = 0; i < line_segs_size; ++i){
//        m_four_vertices[plane_type].push_back(line_segs[line_segs_size-1-i]);
//	}
//}
//// 生成左侧立面
//void HoleFillingInteractiveOperationDialog::on_pushButton_gen_left_plane_clicked()
//{
//    on_pushButton_gen_proj_line_segs_clicked();
//    plane_type = 1;
//    m_four_vertices[plane_type].clear();
//    int line_segs_size = line_segs.size();
//    for(int i = 0; i < line_segs_size; ++i){
//        m_four_vertices[plane_type].push_back(line_segs[line_segs_size-1-i]);
//    }
//}
//// 生成右侧立面
//void HoleFillingInteractiveOperationDialog::on_pushButton_gen_right_plane_clicked()
//{
//    on_pushButton_gen_proj_line_segs_clicked();
//    plane_type = 2;
//    m_four_vertices[plane_type].clear();
//    int line_segs_size = line_segs.size();
//    for(int i = 0; i < line_segs_size; ++i){
//        m_four_vertices[plane_type].push_back(line_segs[line_segs_size-1-i]);
//    }
//}
//// 生成后部立面
//void HoleFillingInteractiveOperationDialog::on_pushButton_gen_backward_plane_clicked()
//{
//    on_pushButton_gen_proj_line_segs_clicked();
//    plane_type = 3;
//    m_four_vertices[plane_type].clear();
//    int line_segs_size = line_segs.size();
//    for(int i = 0; i < line_segs_size; ++i){
//        m_four_vertices[plane_type].push_back(line_segs[line_segs_size-1-i]);
//    }
//}

// 对当前多边形执行CDT三角化
void HoleFillingInteractiveOperationDialog::on_pushButton_clicked()
{
	ConstrainedDelaunayTriangulation cdt;

	// 对每个多边形执行三角化
	vector<vector<int>> triangles;
	vector<float> min_angles;
	int prob_poly_id = -1;
	float min_angle = FLT_MAX;
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		POLYGON2 *pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i];

		if (pPoly->getSize() < 3) {
			cerr << "多边形顶点数量=" << pPoly->getSize() << " < 3!" << endl;
			pPoly->triangles.clear();
			continue;
		}
		if (!pPoly->is_tri) continue;	// 对已经成功执行过三角化的多边形而言，就没有必要再次执行三角化

										// 多边形顶点序列
										/*
										* vertices: pcl格式的顶点序列
										* vertices_: Vertex2*格式的顶点序列，两者具有一一对应关系
										*/
		pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
		vertices->reserve(pPoly->getSize());
		pcl::PointXYZ pcl_p;
		Vertex2* pv2;
		vector<Vertex2*> vertices_;
		vertices_.reserve(pPoly->getSize());
		int poly_index_in_v;
		for (int j = 0; j < pPoly->getSize(); ++j) {
			if (j == 0) {
				pv2 = pPoly->start_point;
			}
			else {
				pv2 = (Vertex2*)pv2->relevant_polys[poly_index_in_v].next_point;
			}
			poly_index_in_v = pPoly->getPolyIndex(pv2);	// 获取当前多边形在pv2的关联多边形的位置
			vertices_.push_back(pv2);
			pcl_p.x = pv2->pos[0];
			pcl_p.y = pv2->pos[1];
			pcl_p.z = pv2->pos[2];
			vertices->push_back(pcl_p);
		}

		cdt.setInputPointCloud(vertices);
		cdt.triangulate(triangles, min_angles);
		if (triangles.size() != pPoly->getSize() - 2) {
			cerr << "三角化失败，poly id = " << i << endl;
			cerr << "triangles.size() = " << triangles.size() << ", pPoly->getSize()-2 = " << pPoly->getSize() - 2 << endl;
			break;
		}

		pPoly->triangles.clear();	// 清空原有的三角面片集合
		for (int j = 0; j < triangles.size(); ++j) {
			Triangle2 tri;
			tri.v[0] = vertices_[triangles[j][0]];
			tri.v[1] = vertices_[triangles[j][1]];
			tri.v[2] = vertices_[triangles[j][2]];
			if (min_angles[j] >= 5.0)
				pPoly->triangles.push_back(tri);	// 将合法的三角形插入多边形的三角形集合中
			else
				cerr << "当前三角形的最小角度：" << min_angles[j] << endl;
			if (min_angles[j] < min_angle) {
				min_angle = min_angles[j];
			}
		}

		if (min_angle < 5.0) {
			cerr << "poly id = " << i << "; 三角形最小角度为" << min_angle << ", 小于5°" << endl;
			prob_poly_id = i;
			break;
		}
		pPoly->is_tri = false;
	}

	// 更新三角面片重心点云
	((PCLViewer*)m_pclviewer)->m_hf.updateTriCentroids();

	updateViewer();

	if (prob_poly_id == -1) { // 如果三角化没有问题，直接返回
		cout << "cdt 三角化完毕" << endl;
		return;
	}

	// 显示有问题的多边形(如果有的话)
	selected_poly_id = prob_poly_id;
	memset(m_CMD, 0, 128);
	strcpy(m_CMD, "select poly");
	updateInfo();

	/////////////////////////////////////////////////////////////////////////////////
	// 显示选取的多边形
	POLYGON2* pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	char buf[64];
	string name;
	Vertex2* pv2, *pv2_next;
	int poly_index;
	pcl::PointXYZ p1, p2;

	// 清空之前的多边形
	int max_poly_size = 256;
	for (int i = 0; i < max_poly_size; ++i) {
		itoa(i, buf, 16);
		name = "display_";
		name.append(buf);
		name.append("_l");
		//cout << name << endl;
		((PCLViewer*)m_pclviewer)->viewer->removeShape(name);
	}
	// cout << "//////////////////////////////////////////////////////////////////////////////////" << endl;
	for (int j = 0; j < pPoly->getSize(); ++j) {
		if (j == 0) {
			pv2 = pPoly->start_point;
		}
		else {
			poly_index = pPoly->getPolyIndex(pv2);
			pv2 = (Vertex2 *)pv2->relevant_polys[poly_index].next_point;
		}

		itoa(j, buf, 16);
		name = "display_";
		name.append(buf);
		name.append("_l");
		//cout << name << endl;
		p1.x = pv2->pos[0];
		p1.y = pv2->pos[1];
		p1.z = pv2->pos[2];

		poly_index = pPoly->getPolyIndex(pv2);
		pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
		p2.x = pv2_next->pos[0];
		p2.y = pv2_next->pos[1];
		p2.z = pv2_next->pos[2];

		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 0.3, 0.8, 0.7, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
	}
}

// 视图更新函数
// 全局更新
void HoleFillingInteractiveOperationDialog::updateViewer()
{
	((PCLViewer*)m_pclviewer)->viewer->removeAllPointClouds();
	((PCLViewer*)m_pclviewer)->viewer->removeAllShapes();

	// 显示多边形线框
	char buf[64];
	string name;

	Vertex2* pv2, *pv2_next;
	int poly_index;
	pcl::PointXYZ p1, p2;
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		for (int j = 0; j < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->getSize(); ++j) {
			if (j == 0) {
				pv2 = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->start_point;
			}
			else {
				pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
			}

			poly_index = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->getPolyIndex(pv2);
			pv2_next = (Vertex2*)pv2->relevant_polys[poly_index].next_point;

			p1.x = pv2->pos[0];
			p1.y = pv2->pos[1];
			p1.z = pv2->pos[2];
			p2.x = pv2_next->pos[0];
			p2.y = pv2_next->pos[1];
			p2.z = pv2_next->pos[2];

			itoa(i, buf, 16);
			name = buf;
			name.append("_");
			itoa(j, buf, 16);
			name.append(buf);
			name.append("_l");
			((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 1, 0, name);
			((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		}
	}

	// 显示三角面片
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		for (int j = 0; j < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles.size(); ++j) {
			itoa(i, buf, 16);
			name = buf;
			name.append("_");
			itoa(j, buf, 16);
			name.append(buf);
			name.append("_t");

			pcl::PointCloud<pcl::PointXYZ>::Ptr tri(new pcl::PointCloud<pcl::PointXYZ>);
			tri->reserve(3);
			pcl::PointXYZ p0, p1, p2;
			p0.x = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[0]->pos[0];
			p0.y = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[0]->pos[1];
			p0.z = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[0]->pos[2];
			tri->push_back(p0);
			p1.x = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[1]->pos[0];
			p1.y = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[1]->pos[1];
			p1.z = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[1]->pos[2];
			tri->push_back(p1);
			p2.x = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[2]->pos[0];
			p2.y = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[2]->pos[1];
			p2.z = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[j].v[2]->pos[2];
			tri->push_back(p2);

			((PCLViewer*)m_pclviewer)->viewer->addPolygon<pcl::PointXYZ>(tri, name);
		}
	}
	// 插入多边形顶点
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(((PCLViewer*)m_pclviewer)->m_hf.m_vertices_cloud, "vertices");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");

	// 插入三角面片重心
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id,
		"triCentroids");
}

// 删除选定的多边形
void HoleFillingInteractiveOperationDialog::on_pushButton_del_poly_clicked()
{
	//cout << "poly id = " << this->selected_poly_id << endl;
	((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id]->clear();
	auto iter = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.begin();
	for (int i = 0; i <= selected_poly_id; ++i) {
		if (i == 0)
			iter = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.begin();
		else
			++iter;
	}
	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.erase(iter);
	cout << "模型多边形数量：" << ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size() << endl;

	// 更新顶点点云
	((PCLViewer*)m_pclviewer)->m_hf.updateVerCloud();
	// 更新三角面片重心点云
	((PCLViewer*)m_pclviewer)->m_hf.updateTriCentroids();
	// 更新视图
	updateViewer();
}

// 基于CDT的三角化孔洞修复
void HoleFillingInteractiveOperationDialog::on_pushButton_CDT_Hole_Filling_clicked()
{
	if (m_points_set_for_tri.size() == 0) {
		QMessageBox::information(this, "提示", "尚未构造三角化点集！");
		return;
	}
	int poly_start_index = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size();

	// 利用点集重心构造局部投影平面
	Eigen::Vector3d centroid;
	centroid << 0, 0, 0;
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		centroid += m_points_set_for_tri[i]->pos;
	}
	centroid = 1.0 / m_points_set_for_tri.size() * centroid;

	// 估计局部投影平面的法向量
	Eigen::Vector3d normal;	// 平面法向量
	normal = m_camera_pos - centroid;
	normal.normalize();
	double d;
	d = centroid.dot(normal) * -1;
	// 局部投影平面参数
	Eigen::Vector4d plane_coeff;
	plane_coeff << normal[0], normal[1], normal[2], d;

	// 将点投影到局部平面上
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		getProjPointOnPlane(m_points_set_for_tri[i]->pos, plane_coeff, m_points_set_for_tri[i]->proj_pos);
	}

	// 多边形顶点序列
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_vertices(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ pcl_p;
	pcl_vertices->reserve(m_points_set_for_tri.size());
	for (int i = 0; i < m_points_set_for_tri.size(); ++i) {
		pcl_p.x = m_points_set_for_tri[i]->proj_pos[0];
		pcl_p.y = m_points_set_for_tri[i]->proj_pos[1];
		pcl_p.z = m_points_set_for_tri[i]->proj_pos[2];
		pcl_vertices->push_back(pcl_p);
	}
	vector<vector<int>> triangles;
	vector<float> min_angles;
	//float min_angle;
	//int min_index;

	// 三角化点集：m_points_set_for_tri
	ConstrainedDelaunayTriangulation cdt;
	cdt.setInputPointCloud(pcl_vertices);
	cdt.triangulate(triangles, min_angles);

	// 将新生成的三角形插入当前的多边形集合中
	Eigen::Vector3d tri_n, cross_product, e01, e02, e10, e12, e20;
	//min_angle = FLT_MAX;
	for (int i = 0; i < triangles.size(); ++i) {
		if (min_angles[i] < 5.0) {
			cerr << "min_angle = " << min_angles[i] << ", < 5°" << endl;
			continue;
		}

		// 校正三角形的顶点排序2(前驱)->0->1(后继)
		e20 = m_points_set_for_tri[triangles[i][0]]->pos - m_points_set_for_tri[triangles[i][2]]->pos;
		e01 = m_points_set_for_tri[triangles[i][1]]->pos - m_points_set_for_tri[triangles[i][0]]->pos;
		cross_product = e20.cross(e01);
		if (cross_product.dot(normal) < 0) {  // 需调整三角形的顶点顺序
			int tmp = triangles[i][1];
			triangles[i][1] = triangles[i][2];
			triangles[i][2] = tmp;
		}

		// 计算三角形的法向量tri_n
		tri_n = e20.cross(e01);
		if (tri_n.dot(normal) < 0) tri_n *= -1;
		tri_n.normalize();

		// 计算三角形平面参数d
		double d = -1 * (m_points_set_for_tri[triangles[i][0]]->pos.dot(tri_n));

		POLYGON2 * pPoly = new POLYGON2();
		pPoly->coeff << tri_n[0], tri_n[1], tri_n[2], d;

		pPoly->insertFirstVertex2(m_points_set_for_tri[triangles[i][0]]);
		pPoly->insertMiddleVertex2(m_points_set_for_tri[triangles[i][1]]);
		pPoly->insertLastVertex2(m_points_set_for_tri[triangles[i][2]]);

		// 插入三角面片
		Triangle2 tri;
		tri.v[0] = m_points_set_for_tri[triangles[i][0]];
		tri.v[1] = m_points_set_for_tri[triangles[i][1]];
		tri.v[2] = m_points_set_for_tri[triangles[i][2]];
		pPoly->triangles.push_back(tri);

		// 标明不再执行三角化
		pPoly->is_tri = false;

		((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(pPoly);
	}

	// 更新顶点点云，重心点云和视图
	char buf[64];
	string name;
	pcl::PointXYZ p1, p2, p3, p;
	Vertex2* cur_v;
	for (int i = poly_start_index; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		POLYGON2* ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i];

		_itoa(i, buf, 16);

		p1.x = ply->start_point->pos[0];
		p1.y = ply->start_point->pos[1];
		p1.z = ply->start_point->pos[2];

		cur_v = ply->start_point;
		int ply_index = ply->getPolyIndex(cur_v);
		p2.x = ((Vertex2*)cur_v->relevant_polys[ply_index].next_point)->pos[0];
		p2.y = ((Vertex2*)cur_v->relevant_polys[ply_index].next_point)->pos[1];
		p2.z = ((Vertex2*)cur_v->relevant_polys[ply_index].next_point)->pos[2];

		p3.x = ((Vertex2*)cur_v->relevant_polys[ply_index].prev_point)->pos[0];
		p3.y = ((Vertex2*)cur_v->relevant_polys[ply_index].prev_point)->pos[1];
		p3.z = ((Vertex2*)cur_v->relevant_polys[ply_index].prev_point)->pos[2];

		// 添加线框
		name = buf;
		name.append("_0_l");
		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p2, 1, 1, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		name = buf;
		name.append("_1_l");
		((PCLViewer*)m_pclviewer)->viewer->addLine(p1, p3, 1, 1, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		name = buf;
		name.append("_2_l");
		((PCLViewer*)m_pclviewer)->viewer->addLine(p2, p3, 1, 1, 0, name);
		((PCLViewer*)m_pclviewer)->viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);

		// 添加三角形
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		cloud->push_back(p1);
		cloud->push_back(p2);
		cloud->push_back(p3);
		name = buf;
		name.append("_0_t");
		((PCLViewer*)m_pclviewer)->viewer->addPolygon<pcl::PointXYZ>(cloud, name);

		// 更新重心点云
		p.x = (p1.x + p2.x + p3.x) / 3.0;
		p.y = (p1.y + p2.y + p3.y) / 3.0;
		p.z = (p1.z + p2.z + p3.z) / 3.0;
		p.data[3] = i;
		((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id->push_back(p);
	}
	((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id_kdtree.setInputCloud(((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id);
	//((PCLViewer*)m_pclviewer)->viewer->removePointCloud("triCentroids");
	((PCLViewer*)m_pclviewer)->viewer->updatePointCloud(((PCLViewer*)m_pclviewer)->m_hf.m_cloud_for_indentify_poly_id, "triCentroids");

	cout << "三角化孔洞修复完成" << endl;
}

// 计算顶点距离
void HoleFillingInteractiveOperationDialog::on_pushButton_dist_V2V_clicked()
{
	if (!m_vertex1 || !m_vertex2) {
		return;
	}
	cout << (m_vertex1->pos - m_vertex2->pos).norm() << endl;
}

// 测试多边形操作
void HoleFillingInteractiveOperationDialog::on_pushButton_test_ply_op_clicked()
{
	cout << "测试多边形操作:" << endl;
	if (((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size() == 0) {
		cerr << "需要首先加载模型数据" << endl;
	}

	POLYGON2 ply, *pPoly;
	pPoly = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[0];

	int size = pPoly->getSize();
	Vertex2* pv2;
	for (int i = 0; i < size; ++i) {
		pv2 = i == 0 ? pPoly->start_point : pPoly->getNextVertex(pv2);

		if (i == 0)
			ply.insertFirstVertex2(pv2);
		else if (i == size - 1)
			ply.insertLastVertex2(pv2);
		else
			ply.insertMiddleVertex2(pv2);
	}

	cout << "ply size = " << ply.getSize() << endl;

	cout << "验证正向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getNextVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}

	cout << "验证反向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getPrevVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}

	///////////////////////////////////////////////////////////////////
	// 删除第二个顶点与第三个顶点
	Vertex2* v1, *v2;
	cout << "删除第1个顶点:" << endl;
	pv2 = ply.start_point;
	v1 = pv2;
	ply.delVertex2(pv2);
	cout << "删除第2个顶点：" << endl;
	pv2 = ply.start_point;
	pv2 = ply.getNextVertex(pv2);
	v2 = pv2;
	ply.delVertex2(pv2);
	size = ply.getSize();
	cout << "当前多边形size = " << ply.getSize() << endl;
	cout << "验证正向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getNextVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}
	cout << "验证反向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getPrevVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}

	///////////////////////////////////////////////////////////////////
	cout << "插入v1" << endl;
	pv2 = ply.start_point;
	ply.insertVertex2(pv2, v1);
	cout << "插入v2" << endl;
	ply.insertVertex2(v1, v2);
	size = ply.getSize();
	cout << "当前多边形size = " << ply.getSize() << endl;
	cout << "验证正向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getNextVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}
	cout << "验证反向循环" << endl;
	for (int i = 0; i <= size; ++i) {
		pv2 = i == 0 ? ply.start_point : ply.getPrevVertex(pv2);
		cout << i << ": " << pv2 << endl;
	}
}

// 在当前多边形的边上插入新的顶点
// selected_poly_id
// m_vertex1
void HoleFillingInteractiveOperationDialog::on_pushButton_insert_vertices_clicked()
{
	if (selected_poly_id == -1 || !m_vertex1) return;

	POLYGON2* ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	int poly_index_in_v = ply->getPolyIndex(m_vertex1);
	if (poly_index_in_v == -1) return;

	int v_num = -1;
	cout << "input number of vertices: ";
	cin >> v_num;

	Vertex2* v, *next_v;
	v = m_vertex1;
	next_v = (Vertex2*)v->relevant_polys[poly_index_in_v].next_point;

	Eigen::Vector3d dir_seg = next_v->pos - v->pos;

	Vertex2* cur_v = v;
	double ratio = 1.0 / (v_num + 1);
	for (int i = 0; i < v_num; ++i) {
		Vertex2* new_v = new Vertex2();

		new_v->pos = v->pos + (i + 1) * ratio * dir_seg;

		((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(new_v);
		ply->insertVertex2(cur_v, new_v);
		cur_v = new_v;
	}

	ply->is_tri = true;

	on_pushButton_clicked();    // 执行CDT三角化
}

// 在多边形内部添加顶点
// m_vertex1与m_vertex2必须都位于同一个选中的多边形中
void HoleFillingInteractiveOperationDialog::on_pushButton__insert_interior_vertices_clicked()
{
	if (selected_poly_id == -1 || !m_vertex1 || !m_vertex2) return;

	POLYGON2* ply = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[selected_poly_id];
	if (ply->getPolyIndex(m_vertex1) == -1 || ply->getPolyIndex(m_vertex2) == -1) {
		cerr << "m_vertex1与m_vertex2必须都位于同一个选中的多边形中" << endl;
	}

	int v_num;  // 顶点数量
	cout << "输入插入的顶点数量：";
	cin >> v_num;

	double ratio = 1.0 / (v_num + 1);

	int poly_size = ply->getSize();
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	cloud->reserve(poly_size + v_num);
	Vertex2* pv2;
	pcl::PointXYZ pcl_p;
	for (int i = 0; i < ply->getSize(); ++i) {
		pv2 = i == 0 ? ply->start_point : ply->getNextVertex(pv2);
		pcl_p.x = pv2->pos[0];
		pcl_p.y = pv2->pos[1];
		pcl_p.z = pv2->pos[2];
		cloud->push_back(pcl_p);
	}

	Eigen::Vector3d dir_seg = m_vertex2->pos - m_vertex1->pos;
	for (int i = 0; i < v_num; ++i) {
		Vertex2* new_v = new Vertex2();
		new_v->pos = m_vertex1->pos + (i + 1)*ratio*dir_seg;
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(new_v);
		pcl_p.x = new_v->pos[0];
		pcl_p.y = new_v->pos[1];
		pcl_p.z = new_v->pos[2];
		cloud->push_back(pcl_p);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	ConstrainedDelaunayTriangulation cdt;
	vector<vector<int>> triangles;
	vector<float> min_angles;
	int prob_poly_id = -1;
	cdt.setInputPointCloud_and_triangulate(cloud, poly_size, triangles, min_angles);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
// 生成两个关键顶点
void HoleFillingInteractiveOperationDialog::on_pushButton_gen_two_key_vertices_clicked()
{
	m_min_y = m_min_x = m_min_z = DBL_MAX;
	m_max_x = DBL_MIN;

	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size(); ++i) {
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0] < m_min_x) {
			m_min_x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
		}
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0] > m_max_x) {
			m_max_x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
		}

		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1] < m_min_y) {
			m_min_y = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1];
		}
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2] < m_min_z) {
			m_min_z = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2];
		}
	}

	if (!m_bottom_left_vertex) {
		m_bottom_left_vertex = new Vertex2();
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(m_bottom_left_vertex);
	}

	if (!m_bottom_right_vertex) {
		m_bottom_right_vertex = new Vertex2();
		((PCLViewer*)m_pclviewer)->m_hf.m_vertices.push_back(m_bottom_right_vertex);
	}

	m_bottom_left_vertex->pos << m_min_x, m_min_y, m_min_z;
	m_bottom_right_vertex->pos << m_max_x, m_min_y, m_min_z;

	cout << "m_min_x = " << m_min_x << "; m_max_x = " << m_max_x << "; m_min_y = " << m_min_y << "; m_min_z = " << m_min_z << endl;
}

// 模型空间范围
void HoleFillingInteractiveOperationDialog::on_pushButton_2_clicked()
{
	m_min_y = m_min_x = m_min_z = DBL_MAX;
	m_max_x = DBL_MIN;

	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size(); ++i) {
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0] < m_min_x) {
			m_min_x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
		}
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0] > m_max_x) {
			m_max_x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
		}

		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1] < m_min_y) {
			m_min_y = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1];
		}
		if (((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2] < m_min_z) {
			m_min_z = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2];
		}
	}

	cout << "m_min_x = " << m_min_x << "; m_max_x = " << m_max_x << "; m_min_y = " << m_min_y << "; m_min_z = " << m_min_z << endl;
}

// 显示多边形法向量
void HoleFillingInteractiveOperationDialog::on_pushButton_display_ply_normal_clicked()
{
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::Normal pn;
	pcl::PointXYZ pcl_p;
	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_polygons.size(); ++i) {
		pn.normal_x = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->coeff[0];
		pn.normal_y = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->coeff[1];
		pn.normal_z = ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->coeff[2];
		normals->push_back(pn);
		Eigen::Vector3d vec;
		vec.setZero();
		vec += ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[0].v[0]->pos;
		vec += ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[0].v[1]->pos;
		vec += ((PCLViewer*)m_pclviewer)->m_hf.m_polygons[i]->triangles[0].v[2]->pos;
		vec /= 3;
		pcl_p.x = vec[0];
		pcl_p.y = vec[1];
		pcl_p.z = vec[2];
		cloud->push_back(pcl_p);
	}

	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("ply normals");
	((PCLViewer*)m_pclviewer)->viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(cloud, normals, 1, 1, "ply normals");
}

// 生成底面
void HoleFillingInteractiveOperationDialog::on_pushButton_gen_bottom_plane_clicked()
{
	// 数据结构：line_segs(需首先生成投影线段)
	POLYGON2* ply = new POLYGON2();
	ply->coeff << 0, -1, 0, m_min_y;
	for (int i = 0; i < line_segs.size(); ++i) {
		if (i == 0) {
			ply->insertFirstVertex2(line_segs[line_segs.size() - 1 - i]);
		}
		else {
			ply->insertMiddleVertex2(line_segs[line_segs.size() - 1 - i]);
		}
	}
	ply->insertMiddleVertex2(m_bottom_left_vertex);
	ply->insertLastVertex2(m_bottom_right_vertex);

	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);

	on_pushButton_clicked();
}

// 生成左侧立面(最好先点击“模型空间范围”按钮）
void HoleFillingInteractiveOperationDialog::on_pushButton_4_clicked()
{
	// 数据结构：line_segs(需首先生成投影线段)
	POLYGON2* ply = new POLYGON2();
	ply->coeff << -1, 0, 0, m_min_x;

	int size = line_segs.size();
	for (int i = 0; i < size; ++i) {
		if (i == 0) {
			ply->insertFirstVertex2(line_segs[size - 1 - i]);
		}
		else if (i == size - 1) {
			ply->insertLastVertex2(line_segs[size - 1 - i]);
		}
		else {
			ply->insertMiddleVertex2(line_segs[size - 1 - i]);
		}
	}

	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);

	on_pushButton_clicked();
}

// 生成右侧立面
void HoleFillingInteractiveOperationDialog::on_pushButton_5_clicked()
{
	// 数据结构：line_segs(需首先生成投影线段)
	POLYGON2* ply = new POLYGON2();
	ply->coeff << 1, 0, 0, -1 * m_max_x;

	int size = line_segs.size();
	for (int i = 0; i < size; ++i) {
		if (i == 0) {
			ply->insertFirstVertex2(line_segs[size - 1 - i]);
		}
		else if (i == size - 1) {
			ply->insertLastVertex2(line_segs[size - 1 - i]);
		}
		else {
			ply->insertMiddleVertex2(line_segs[size - 1 - i]);
		}
	}

	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);

	on_pushButton_clicked();
}
// 在生成前三个立面时均手动选取顶点1和顶点2构成投影线段
// 在生成后侧立面时，手动选取顶点1和对应的多边形（保证在该多边形的前提下，顶点1的前驱仍位于孔洞边界上）
// 生成后侧立面
void HoleFillingInteractiveOperationDialog::on_pushButton_6_clicked()
{
	// 数据结构：m_points_set_for_tri(模式一)
	POLYGON2* ply = new POLYGON2();
	ply->coeff << 0, 0, -1, m_min_y;

	int size = m_points_set_for_tri.size();
	for (int i = 0; i < size; ++i) {
		if (i == 0) {
			ply->insertFirstVertex2(m_points_set_for_tri[size - 1 - i]);
		}
		else if (i == size - 1) {
			ply->insertLastVertex2(m_points_set_for_tri[size - 1 - i]);
		}
		else {
			ply->insertMiddleVertex2(m_points_set_for_tri[size - 1 - i]);
		}
	}

	((PCLViewer*)m_pclviewer)->m_hf.m_polygons.push_back(ply);

	on_pushButton_clicked();
}

// 模型封闭性检测
void HoleFillingInteractiveOperationDialog::on_pushButton_3_clicked()
{
	// 首先检测模型是否存在孔洞
	cout << "孔洞检测..." << endl;
	pcl::PointXYZ p, pcl_p;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);

	for (int i = 0; i < ((PCLViewer*)m_pclviewer)->m_hf.m_vertices.size(); ++i) {
		if (((PCLViewer*)m_pclviewer)->m_hf.getSingleRelationCount(((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]) > 0) {
			cerr << "模型不封闭，问题顶点id：" << i << endl;

			p.x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
			p.y = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1];
			p.z = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2];
			cloud->push_back(p);
		}

		pcl_p.x = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[0];
		pcl_p.y = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[1];
		pcl_p.z = ((PCLViewer*)m_pclviewer)->m_hf.m_vertices[i]->pos[2];
		vertices->push_back(pcl_p);
	}

	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("borer_v");
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "borer_v");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "borer_v");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "borer_v");

	if (cloud->size() == 0)
		cout << "模型封闭" << endl;
	else
		cout << "边界顶点数量：" << cloud->size() << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////
	// 然后检测冗余顶点
	cout << "冗余顶点检测..." << endl;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(vertices);
	vector<int> indices;
	vector<float> dists;
	cloud->clear();
	for (int i = 0; i < vertices->size(); ++i) {
		kdtree.radiusSearch(vertices->points[i], 0.0001, indices, dists);
		if (indices.size() > 1) {
			cloud->push_back(vertices->points[i]);
		}
	}
	if (cloud->size() == 0) {
		cout << "没有冗余顶点" << endl;
	}
	else {
		cout << "冗余顶点数量：" << cloud->size() << endl;
	}

	((PCLViewer*)m_pclviewer)->viewer->removePointCloud("rongyu_v");
	((PCLViewer*)m_pclviewer)->viewer->addPointCloud(cloud, "rongyu_v");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1, 0, "rongyu_v");
	((PCLViewer*)m_pclviewer)->viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "rongyu_v");
}
