#include "pclviewer.h"
#include "ui_pclviewer.h"

#include <QFileDialog>
#include <QString>
#include <QSettings>
#include <QStringList>
#include <QFile>
#include <QTextStream>

#include <queue>
#include <algorithm>

#include "setbgcolordialog.h"
#include "removepointclouddialog.h"
#include "setparams_simpilyverticessize_dialog.h"
#include "setparams_mergingvertices_dialog.h"

#include "Event.h"

/****************************************前半段：开始**************************************************/
#include <QLayout>
#include <QFileDialog>
#include <QSettings>
#include <QMessageBox>
#include <omp.h>

#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/convolution_3d.h>
#include <pcl/filters/convolution.h>
#include <pcl/filters/fast_bilateral.h>
#include <pcl/filters/random_sample.h>
#include <pcl/surface/mls.h> 
#include <pcl/common/transforms.h>
#include <pcl/filters/filter.h>

#include "PlaneDetect.h"
#include "ui_MainWindow.h"
#include "Registration.h"
#include "HoleFilling_PLY.h"
#include "TriangularMeshing.h"
#include "PlaneDetectSetParamDialog.h"
#include "RemovePointCloudDialog.h"
#include "SetBGColorDialog.h"
#include "SetPointCloudProDialog.h"
/****************************************前半段：结束**************************************************/

PCLViewer::PCLViewer (QWidget *parent) :
  QMainWindow (parent),
  ui (new Ui::PCLViewer)
{
  ui->setupUi (this);
  this->setWindowTitle ("Rockmass 3D Reconstruction");

  // Set up the QVTK window
  viewer.reset (new pcl::visualization::PCLVisualizer ("viewer", false));
  ui->qvtkWidget->SetRenderWindow (viewer->getRenderWindow ());
  viewer->setupInteractor (ui->qvtkWidget->GetInteractor (), ui->qvtkWidget->GetRenderWindow ());
  ui->qvtkWidget->update();

  // 加载参数
  loadParams();

  // 设置顶点化简器的各项参数
  m_svs.setParams(T_angle_betPointAndPoint, T_angle_betPointAndLine, T_dist_betPointAndLine,
                  T_lenRatio, T_angle_betLines);  
  // 设置顶点合并的各项参数
  m_sv.ratio_of_scale = this->ratio_of_scale;

  m_vertices.reset(new pcl::PointCloud<pcl::PointXYZ>);
}

PCLViewer::~PCLViewer ()
{
  delete ui;
}

/****************************************前半段：开始**************************************************/
// 打开文件
void PCLViewer::on_openPointCloudAction_triggered()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("open file"), " ", tr("pcdFiles(*.pcd)"));
	if (fileName == "") return;
	pcl::io::loadPCDFile(fileName.toStdString(), *source_cloud);
	cout << "loaded " << source_cloud->size() << " points." << endl;
	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();
}
void PCLViewer::on_openTxtAction_triggered()
{
	int number_Txt;
	FILE *fp_txt;
	tagPOINT_3D TxtPoint;
	vector<tagPOINT_3D> m_vTxtPoints;
	QString fileName = QFileDialog::getOpenFileName(this, tr("open file"), " ", tr("txtFiles(*.txt)"));
	if (fileName == "") return;
	//QApplication app(argc, argv);

	QByteArray ba = fileName.toLatin1();
	char *c_str2 = ba.data();
	//printf("str2: %s", c_str2);
	//return app.exec();
	fp_txt = fopen(c_str2, "r");

	if (fp_txt)
	{
		while (fscanf(fp_txt, "%lf %lf %lf", &TxtPoint.x, &TxtPoint.y, &TxtPoint.z) != EOF)
		{
			m_vTxtPoints.push_back(TxtPoint);
		}
	}
	else
		cout << "txt数据加载失败！" << endl;
	number_Txt = m_vTxtPoints.size();
	pcl::PointCloud<pcl::PointXYZ> cloud_txt;


	// Fill in the cloud data  
	cloud_txt.width = number_Txt;
	cloud_txt.height = 1;
	cloud_txt.is_dense = false;
	cloud_txt.points.resize(cloud_txt.width * cloud_txt.height);


	for (size_t i = 0; i < cloud_txt.points.size(); ++i)
	{
		cloud_txt.points[i].x = m_vTxtPoints[i].x;
		cloud_txt.points[i].y = m_vTxtPoints[i].y;
		cloud_txt.points[i].z = m_vTxtPoints[i].z;
	}
	*source_cloud = cloud_txt;
	cout << "loaded " << source_cloud->size() << " points." << endl;
	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();
}
//define by czh
// 打开配准点云file_1  target_cloud_registration
void PCLViewer::on_openFile1RegistrationAction_triggered()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("open file"), " ", tr("pcdFiles(*.pcd)"));
	if (fileName == "") return;
	pcl::io::loadPCDFile(fileName.toStdString(), *target_cloud_registration);
	cout << "loaded " << target_cloud_registration->size() << " points." << endl;
	viewer->removePointCloud("source");
	viewer->addPointCloud(target_cloud_registration, "source");
	ui->qvtkWidget->update();
}

//define by czh
// 打开配准点云file_2  source_cloud_registration
void PCLViewer::on_openFile2RegistrationAction_triggered()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("open file"), " ", tr("pcdFiles(*.pcd)"));
	if (fileName == "") return;
	pcl::io::loadPCDFile(fileName.toStdString(), *source_cloud_registration);
	cout << "loaded " << source_cloud_registration->size() << " points." << endl;
	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud_registration, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(source_cloud_registration, cloud_result);
}

//*******************************************开始：移除NAN点和体素滤波直接覆盖原点云对象*********************************************************************

//define by czh
// 移除NaN点
//输入变量：source_cloud_registration  target_cloud_registration
//输出变量：source_cloud_registration  target_cloud_registration
void PCLViewer::on_removeNanAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();
	//变量转换器
	verb_transform(cloud_result, source_cloud_registration);

	//去除NAN点 source_cloud_registration
	std::vector<int> indices_src; //保存去除的点的索引
	pcl::removeNaNFromPointCloud(*source_cloud_registration, *source_cloud_registration, indices_src);
	std::cout << "remove *source_cloud_registration nan" << endl;
	//去除NAN点 target_cloud_registration
	std::vector<int> indices_tgt;
	pcl::removeNaNFromPointCloud(*target_cloud_registration, *target_cloud_registration, indices_tgt);
	std::cout << "remove *target_cloud_registration nan" << endl;

	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud_registration, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(source_cloud_registration, cloud_result);
}

//define by czh
//移除source_cloud的NAN点
void PCLViewer::on_removeNan1Action_triggered() {
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	std::vector<int> indices;
	PointCloudT::Ptr new_cloud(new PointCloudT);
	pcl::removeNaNFromPointCloud<pcl::PointXYZ>(*source_cloud, *new_cloud, indices);
	source_cloud = new_cloud;

	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();

	cout << "After nan points been removed, points size =" << source_cloud->size() << endl;
}


//define by czh
/*
将两个点云进行合并后的滤波
*/
void PCLViewer::on_voxelGridFiltMergeCloudAction_triggered() {
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr mergeCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	PointCloudT::Ptr cloud_merge(new PointCloudT);

	//变量转换器
	//verb_transform(cloud_result, cloud_result);

	//清除界面上所有点云
	viewer->removeAllPointClouds();

	savePointCloudFile(target_cloud_registration, cloud_result, "double_shadow");

	//合并两个点云
	mergePointCloud(target_cloud_registration, cloud_result, mergeCloud);

	//将 PointXYZRGB to PointXYZ
	xyzrgbTransformxyz(mergeCloud, cloud_merge);

	//下采样滤波 source_cloud_registration
	pcl::VoxelGrid<pcl::PointXYZ> voxel_grid;
	float delta{};
	cout << "set  leafsize :" << endl;
	cin >> delta;
	voxel_grid.setLeafSize(delta, delta, delta);
	voxel_grid.setInputCloud(cloud_merge);
	int pointCloud_num1 = mergeCloud->size();
	PointCloudT::Ptr cloud_tr(new PointCloudT);
	voxel_grid.filter(*cloud_merge);
	std::cout << "down size *cloud_tr_o from " << pointCloud_num1 << "to" << cloud_merge->size() << endl;

	//显示点云
	viewer->removePointCloud("source");
	viewer->addPointCloud(cloud_merge, "source");

	//将局部变量cloud_merge赋值给全局source_cloud，进行平面提取
	*source_cloud = *cloud_merge;
}



//define by czh
// 两个配准点云，进行体素滤波
//输入变量：source_cloud_registration  target_cloud_registration
//输出变量：source_cloud_registration  target_cloud_registration
void PCLViewer::on_voxelGridFiltAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	//变量转换器
	verb_transform(cloud_result, source_cloud_registration);

	//清除界面上所有点云
	viewer->removeAllPointClouds();

	//下采样滤波 source_cloud_registration
	pcl::VoxelGrid<pcl::PointXYZ> voxel_grid;
	float delta{};
	cout << "set both leafsize :" << endl;
	cin >> delta;
	voxel_grid.setLeafSize(delta, delta, delta);
	voxel_grid.setInputCloud(source_cloud_registration);
	int pointCloud_num1 = source_cloud_registration->size();
	PointCloudT::Ptr cloud_tr(new PointCloudT);
	voxel_grid.filter(*source_cloud_registration);
	std::cout << "down size *cloud_tr_o from " << pointCloud_num1 << "to" << source_cloud_registration->size() << endl;
	//下采样滤波 cloud_tgt_o
	pcl::VoxelGrid<pcl::PointXYZ> voxel_grid_2;
	voxel_grid_2.setLeafSize(delta, delta, delta);
	voxel_grid_2.setInputCloud(target_cloud_registration);
	int pointCloud_num2 = target_cloud_registration->size();
	PointCloudT::Ptr cloud_tgt(new PointCloudT);
	voxel_grid_2.filter(*target_cloud_registration);
	std::cout << "down size *cloud_tgt_o.pcd from " << pointCloud_num2 << "to" << target_cloud_registration->size() << endl;

	//直接显示在主界面
	viewer->removePointCloud("source");
	viewer->addPointCloud(source_cloud_registration, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(source_cloud_registration, cloud_result);
}

//define by czh
//对source_cloud进行体素滤波
void PCLViewer::on_voxelGridFilt1Action_triggered() {
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	PointCloudT::Ptr cloud_filtered(new PointCloudT);
	pcl::VoxelGrid<pcl::PointXYZ> sor;
	sor.setInputCloud(source_cloud);
	float delta{};
	cout << "set leafsize :" << endl;
	cin >> delta;
	sor.setLeafSize(delta, delta, delta);
	sor.filter(*cloud_filtered);
	cout << "after voxel filterd, cloud size = " << cloud_filtered->size() << endl;
	source_cloud = cloud_filtered;

	viewer->removeAllPointClouds();
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();
}
//统计滤波
void PCLViewer::on_statisticalOutlierRemovalFiltAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
	PointCloudT::Ptr cloud_filtered(new PointCloudT);
	sor.setInputCloud(source_cloud);
	float delta{};
	cout << "设置在进行统计时考虑查询点临近点数，推荐50-150:" << endl;
	cin >> delta;
	sor.setMeanK(delta);
	float felta{};
	cout << "设置判断是否为离群点的阀值，推荐0.5-1.5:" << endl;
	cin >> felta;
	sor.setStddevMulThresh(felta);
	sor.filter(*cloud_filtered);
	cout << "after statisticalOutlierRemoval filt, cloud size = " << cloud_filtered->size() << endl;
	source_cloud = cloud_filtered;

	viewer->removeAllPointClouds();
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();
}
void PCLViewer::on_setNumberforpointAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	//PointCloudT::Ptr cloud_in(new PointCloudT), cloud_out(new PointCloudT);
	//boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Viewer"));
	//pcl::io::loadPCDFile("unknown_object.pcd", *cloud_in);
	//std::cerr << *cloud_in << std::endl;
	cout << "after statisticalOutlierRemoval filt, cloud size = " << endl;
	pcl::RandomSample<pcl::PointXYZ> rs;
	PointCloudT::Ptr cloud_filtered(new PointCloudT);
	rs.setInputCloud(source_cloud);
	float felta{};
	cout << "设置点云点数" << endl;
	cin >> felta;

	rs.setSample(felta);
	rs.filter(*cloud_filtered);
	cout << "after statisticalOutlierRemoval filt, cloud size = " << cloud_filtered->size() << endl;
	source_cloud = cloud_filtered;

	viewer->removeAllPointClouds();
	viewer->addPointCloud(source_cloud, "source");
	ui->qvtkWidget->update();

}
void PCLViewer::on_rebuildPlaneAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	pcl::search::KdTree<pcl::PointXYZ>::Ptr treerebuild(new pcl::search::KdTree<pcl::PointXYZ>);
	PointCloudT::Ptr mls_cloud(new PointCloudT);
	pcl::PointCloud<pcl::PointNormal> mls_points;
	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
	mls.setComputeNormals(true);
	mls.setInputCloud(source_cloud);
	mls.setPolynomialFit(true);
	mls.setSearchMethod(treerebuild);
	float felta{};
	cout << "设置搜索范围 提示0.3" << endl;
	cin >> felta;
	mls.setSearchRadius(felta);
	cout << "start" << endl;
	mls.process(mls_points);
	mls_cloud->resize(mls_points.size());

	for (size_t i = 0; i < mls_points.points.size(); ++i)
	{
		mls_cloud->points[i].x = mls_points.points[i].x; //error 
		mls_cloud->points[i].y = mls_points.points[i].y; //error 
		mls_cloud->points[i].z = mls_points.points[i].z; //error 
	}
	cout << "after rebuile, cloud size = " << mls_cloud->size() << endl;
	source_cloud = mls_cloud;
	viewer->removeAllPointClouds();
	viewer->addPointCloud(source_cloud, "source");
	cout << "done" << endl;
	ui->qvtkWidget->update();
	
	//pcl::copyPointCloud(*mls_points,*source_cloud)
}


//*******************************************结束：移除NAN点和体素滤波直接覆盖原点云对象*********************************************************************

//define by czh
//保存当前点云
void PCLViewer::on_savePointCloudAction_triggered() {
	pcl::io::savePCDFile("now_source_cloud.pcd", *source_cloud);
	cout << "save source_cloud success" << endl;
}

//define by czh
//默认对【source_cloud】进行网格化
void  PCLViewer::on_TriangularMeshingAction_triggered() {
	pcl::PolygonMesh mesh_result;
	//执行网格化
	triangularMeshing(source_cloud, mesh_result);

	//保存网格化的结果
	//保存一个具体的.ply文件，用于下一步的孔洞修复，孔洞修复是基于cgal中的openMesh结构，由于PCL的PointCloud类型不能直接转换成openMesh的mesh结构
	pcl::io::savePLYFile("result_mesh/result_mesh.ply", mesh_result);
}

//define by czh
//孔洞修复，保存结果为.off或者ply格式都行文件
//注：只能读取"result_mesh/result_mesh.ply"，进行孔洞修复
void PCLViewer::on_repairHolesOFFAction_triggered() {
	const char* fileNameIn = "result_mesh/result_mesh.ply";
	std::string fileNameOut = "filled_OM.ply";
	//执行孔洞修复
	holeFilling_PLY(fileNameIn, fileNameOut);

	//将结果读取出来，并且从ply转换为pcd格式
	pcl::PLYReader reader;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	reader.read<pcl::PointXYZ>("result_HoleFilling/filled_OM.ply", *cloud);
	pcl::io::savePCDFile("result_pcd/filed_OM.pcd", *cloud);
	cout << "结果从filled_OM.ply转换为了filled_OM.pcd,并保存在result_pcd文件夹下" << endl;
}


//define by czh
//旋转原点云
//输入数据：source_cloud_registration
//输出数据：cloud_tr
void PCLViewer::on_rotatePointCloudAction_triggered() {
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	//变量转换器
	verb_transform(cloud_result, source_cloud_registration);

	//清除界面上所有点云
	viewer->removeAllPointClouds();

	//tools->transformation(*source_cloud_registration, *cloud_icp, transformation_matrix);
	transformation(*source_cloud_registration, *cloud_tr);

	//直接显示在主界面
	viewer->removePointCloud("source");
	viewer->addPointCloud(cloud_tr, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(cloud_tr, cloud_result);
}
//define by czh 
//SAC配准
//输入数据：cloud_tr target_cloud_registration 
//输出数据：cloud_icp
void PCLViewer::on_registrationSACAction_triggered() {
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	//变量转换器
	verb_transform(cloud_result, cloud_sac);

	//清除界面上所有点云
	viewer->removeAllPointClouds();

	//计算表面法线1
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_src;
	ne_src.setInputCloud(cloud_sac);
	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_src(new pcl::search::KdTree< pcl::PointXYZ>());
	ne_src.setSearchMethod(tree_src);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_src_normals(new pcl::PointCloud< pcl::Normal>);
	ne_src.setRadiusSearch(0.02);
	ne_src.compute(*cloud_src_normals);
	//计算表面法线2
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne_tgt;
	ne_tgt.setInputCloud(target_cloud_registration);
	pcl::search::KdTree< pcl::PointXYZ>::Ptr tree_tgt(new pcl::search::KdTree< pcl::PointXYZ>());
	ne_tgt.setSearchMethod(tree_tgt);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_tgt_normals(new pcl::PointCloud< pcl::Normal>);
	//ne_tgt.setKSearch(20);
	ne_tgt.setRadiusSearch(0.02);
	ne_tgt.compute(*cloud_tgt_normals);

	//计算FPFH 1
	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_src;
	fpfh_src.setInputCloud(cloud_sac);
	fpfh_src.setInputNormals(cloud_src_normals);
	pcl::search::KdTree<PointT>::Ptr tree_src_fpfh(new pcl::search::KdTree<PointT>);
	fpfh_src.setSearchMethod(tree_src_fpfh);
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_src(new pcl::PointCloud<pcl::FPFHSignature33>());
	fpfh_src.setRadiusSearch(0.05);
	fpfh_src.compute(*fpfhs_src);
	std::cout << "compute *cloud_tr fpfh" << endl;
	//计算FPFH 2
	pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_tgt;
	fpfh_tgt.setInputCloud(target_cloud_registration);
	fpfh_tgt.setInputNormals(cloud_tgt_normals);
	pcl::search::KdTree<PointT>::Ptr tree_tgt_fpfh(new pcl::search::KdTree<PointT>);
	fpfh_tgt.setSearchMethod(tree_tgt_fpfh);
	pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_tgt(new pcl::PointCloud<pcl::FPFHSignature33>());
	fpfh_tgt.setRadiusSearch(0.05);
	fpfh_tgt.compute(*fpfhs_tgt);
	std::cout << "compute *cloud_tgt fpfh" << endl;

	//SAC配准
	pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> scia;
	scia.setInputSource(cloud_sac);
	scia.setInputTarget(target_cloud_registration);
	scia.setSourceFeatures(fpfhs_src);
	scia.setTargetFeatures(fpfhs_tgt);
	//scia.setMinSampleDistance(1);
	//scia.setNumberOfSamples(2);
	//scia.setCorrespondenceRandomness(20);
	scia.align(*cloud_sac);
	//输出SAC旋转矩阵
	std::cout << "sac has converged:" << scia.hasConverged() << "  score: " << scia.getFitnessScore() << endl;
	Eigen::Matrix4f sac_trans;
	sac_trans = scia.getFinalTransformation();
	std::cout << sac_trans << endl;
	clock_t sac_time = clock();
	std::cout << "SAC done" << endl;
	//新建页面显示SAC效果
	visualization(target_cloud_registration, cloud_tr, cloud_sac);

	//直接显示SAC效果-单点云
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr mergeCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	mergePointCloud(target_cloud_registration, cloud_sac, mergeCloud);
	viewer->removePointCloud("source");
	viewer->addPointCloud(mergeCloud, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(cloud_sac, cloud_result);

}

//define by czh
//配准
//输入数据：target_cloud_registration cloud_icp
//输出数据：cloud_icp
void PCLViewer::on_registrationICPAction_triggered()
{
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	//变量转换器
	verb_transform(cloud_result, cloud_icp);

	//清除界面上所有点云
	viewer->removeAllPointClouds();

	pcl::console::TicToc time;
	time.tic();

	pcl::IterativeClosestPoint<PointT, PointT> icp;
	icp.setMaximumIterations(iterations);
	icp.setInputSource(cloud_icp);
	icp.setInputTarget(target_cloud_registration);
	icp.align(*cloud_icp);
	icp.setMaximumIterations(1);  // We set this variable to 1 for the next time we will call .align () function
	std::cout << "Applied " << iterations << " ICP iteration(s) in " << time.toc() << " ms" << std::endl;

	if (icp.hasConverged())
	{
		std::cout << "\nICP has converged, score is " << icp.getFitnessScore() << std::endl;
		std::cout << "\nICP transformation " << iterations << " : cloud_icp -> cloud_in" << std::endl;
		//transformation_matrix = icp.getFinalTransformation().cast<double>();
		//print4x4Matrix(transformation_matrix);

		//update file
		//tools->savePointCloudFile(cloud_in_1, cloud_icp, iterations);
	}
	else
	{
		//检测ICP算法是否收敛，不然就退出程序。
		PCL_ERROR("\nICP has not converged.\n");
		system("pause");
		return;
	}

	//双通道显示点云
	visualization(target_cloud_registration, cloud_tr, cloud_icp);
	//保存点云
	savePointCloudFile(target_cloud_registration, cloud_icp, "icp-run" + iterations);

	//直接显示ICP效果-单点云
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr mergeCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	mergePointCloud(target_cloud_registration, cloud_icp, mergeCloud);
	viewer->removePointCloud("source");
	viewer->addPointCloud(mergeCloud, "source");
	ui->qvtkWidget->update();

	//变量转换器
	verb_transform(cloud_icp, cloud_result);
}

//define by czh
/*基于平面的配准
平面配准同样使用了变量转换器，只不过是写了第一个和最后一个函数的内部！
*/
void PCLViewer::on_registrationPlaneAction_triggered() {
	
	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	/*在执行平面配准时，要事先求得当前输入的两个点云的平面数据*.pcd *.txt
	*其实下面三个函数都是对PlaneDetect.h中函数的组合使用而已。
	*/

	/*原程序对全局变量source_cloud的平面提取*/
	//临时保存source_cloud到temp_cloud，使用完之后，再恢复
	PointCloudT::Ptr temp_cloud(new PointCloudT);
	*temp_cloud = *source_cloud;
	//变量转换器：将source_cloud赋值为source_cloud_registration
	verb_transform(source_cloud_registration, source_cloud);
	//为平面提取加载参数
	on_load_param_action_triggered();
	//估计法向量
	on_normalEstimateAction_triggered();
	//自动执行
	on_autoPerformAction_triggered();
	//保存多边形数据
	//on_savePolyDataAction_triggered();
	source_polygon = polygon_transform();

	//变量转换器：将source_cloud赋值为target_cloud_registration
	verb_transform(target_cloud_registration, source_cloud);
	//为平面提取加载参数
	on_load_param_action_triggered();
	//估计法向量
	on_normalEstimateAction_triggered();
	//自动执行
	on_autoPerformAction_triggered();
	//保存多边形数据
	//on_savePolyDataAction_triggered();
	target_polygon = polygon_transform();

	//恢复source_cloud
	*source_cloud = *temp_cloud;

	start = std::clock();
	cout << "auto run......." << endl;
	//loadpolygon();
	while (source_index < source_polygon.size())
	{
		FindSimilarPoly();
		cout << "match continue..." << endl;
		source_index++;
		IsMatch = false;
	}
	registration();
	registration_cloud();
	//icp();
	finish = std::clock();
	cout << "finally cost " << finish - start << " ms." << endl;
	memset(state, 0, CMD_SIZE);
}

// 改变背景颜色
void PCLViewer::on_bgColorMenu_triggered()
{

	//清除界面
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	SetBGColorDialog dlg;

	if (dlg.exec() != QDialog::Accepted) return;

	int r, g, b;
	dlg.getRGB(r, g, b);

	viewer->setBackgroundColor(r, g, b);
	viewer->spinOnce();
	ui->qvtkWidget->update();
	cout << "修改完背景颜色" << endl;
}
//change by czh
// 改变点云颜色
void PCLViewer::on_pointCloudColorMenu_triggered()
{
	SetPointCloudProDialog *dlg = new SetPointCloudProDialog();
	if (dlg->exec() != QDialog::Accepted) return;

	int r, g, b, size;
	dlg->getData(r, g, b, size);
	delete dlg;

	//清除界面上所有点云
	viewer->removeAllPointClouds();
	pcl::visualization::PointCloudColorHandlerCustom<PointT> cloud_color_h(source_cloud, r, g, b);
	viewer->addPointCloud(source_cloud, cloud_color_h);
	ui->qvtkWidget->update();
	cout << "修改完点云颜色" << endl;
}

//define by czh
//清除屏幕中的点云
void PCLViewer::on_cleanPointCloudAction_triggered() {
	//清除界面上所有点云
	viewer->removeAllPointClouds();
	ui->qvtkWidget->update();
	cout << "清除屏幕中的点云" << endl;
}
// 把点云移动至重心
void PCLViewer::on_translateToCentroidAction_triggered()
{
	if (source_cloud->size() == 0) return;
	pcl::PointXYZ p;
	p.x = p.y = p.z = 0;
	for (int i = 0; i < source_cloud->size(); ++i)
	{
		p.x += source_cloud->points[i].x;
		p.y += source_cloud->points[i].y;
		p.z += source_cloud->points[i].z;
	}
	p.x /= float(source_cloud->size());
	p.y /= float(source_cloud->size());
	p.z /= float(source_cloud->size());

	Eigen::Affine3f transform = Eigen::Affine3f::Identity();
	transform.translation() << -1.0*p.x,
		-1.0*p.y,
		-1.0*p.z;
	pcl::transformPointCloud(*source_cloud, *source_cloud, transform);
	viewer->updatePointCloud(source_cloud, "source");
	viewer->resetCamera();
	ui->qvtkWidget->update();
	cout << "把点云移动至重心" << endl;
}

// 设置平面提取参数
void PCLViewer::on_plane_detect_set_param_Action_triggered()
{
	PlaneDetectSetParamDialog dlg;
	dlg.exec();
}
// 移除冗余点
void PCLViewer::on_removeRedundantPointsAction_triggered()
{
	std::vector<bool> is_redundant(source_cloud->size(), false);
	kdtree_source.setInputCloud(source_cloud);
	std::vector<int> indices;
	std::vector<float> dists;
	for (int i = 0; i < source_cloud->size(); ++i)
	{
		if (is_redundant[i]) continue;
		kdtree_source.radiusSearch(source_cloud->points[i], min_dist_between_points, indices, dists);

		for (int j = 1; j < indices.size(); ++j)
			is_redundant[indices[j]] = true;
	}
	PointCloudT::Ptr new_cloud(new PointCloudT);
	new_cloud->reserve(source_cloud->size());
	for (int i = 0; i < source_cloud->size(); ++i)
	{
		if (!is_redundant[i])
			new_cloud->push_back(source_cloud->points[i]);
	}
	source_cloud = new_cloud;
	cout << "after remove redundant points, cloud size = " << source_cloud->size() << endl;

}


// 修正坐标
// 默认y轴朝上；x轴水平；z轴指向岩体方向
// 简单处理：swap y与z的值，z值再乘以-1
void PCLViewer::on_regulateCoorAction_triggered()
{
	float tmp;
	for (int i = 0; i < source_cloud->size(); ++i)
	{
		tmp = source_cloud->points[i].y;
		source_cloud->points[i].y = source_cloud->points[i].z;
		source_cloud->points[i].z = tmp;
		source_cloud->points[i].z *= -1.0f;
	}
	viewer->updatePointCloud(source_cloud, "source");
	ui->qvtkWidget->update();
}

// 为平面提取加载参数
void PCLViewer::on_load_param_action_triggered()
{
	QSettings configIni("config.txt", QSettings::IniFormat);

	isConductTranslate = configIni.value("PlaneDetect/isConductTranslate").toBool();
	min_dist_between_points = configIni.value("PlaneDetect/min_dist_between_points").toFloat();
	r_for_estimate_normal = configIni.value("PlaneDetect/r_for_estimate_normal").toFloat();
	interval_level = configIni.value("PlaneDetect/interval_level").toInt();
	normal_length_for_display = configIni.value("PlaneDetect/normal_length_for_display").toFloat();
	normal_length_for_selected = configIni.value("PlaneDetect/normal_length_for_selected").toFloat();
	is_norm_direction_valid = configIni.value("PlaneDetect/is_norm_direction_valid").toBool();
	r_for_regulate_normal = configIni.value("PlaneDetect/r_for_regulate_normal").toFloat();
	num_res = configIni.value("PlaneDetect/num_res").toFloat();
	color_gain = configIni.value("PlaneDetect/color_gain").toFloat();
	std_dev_gaussian = configIni.value("PlaneDetect/std_dev_gaussian").toFloat();
	search_radius_create_gradient = configIni.value("PlaneDetect/search_radius_create_gradient").toFloat();
	T_vote_num = configIni.value("PlaneDetect/T_vote_num").toInt();
	radius_base = configIni.value("PlaneDetect/radius_base").toFloat();
	delta_radius = configIni.value("PlaneDetect/delta_radius").toFloat();
	S_threshold = configIni.value("PlaneDetect/S_threshold").toFloat();
	dev_threshold = configIni.value("PlaneDetect/dev_threshold").toFloat();
	search_radius_for_plane_seg = configIni.value("PlaneDetect/search_radius_for_plane_seg").toFloat();
	T_num_of_single_plane = configIni.value("PlaneDetect/T_num_of_single_plane").toInt();
	radius_border = configIni.value("PlaneDetect/radius_border").toFloat();
	radius_local = configIni.value("PlaneDetect/radius_local").toFloat();
	T_angle_bias_integrate = configIni.value("PlaneDetect/T_angle_bias_integrate").toFloat();
	r_local = configIni.value("PlaneDetect/r_local").toFloat();
	T_curvature_bias = configIni.value("PlaneDetect/T_curvature_bias").toFloat();
	T_point_normal_bias = configIni.value("PlaneDetect/T_point_normal_bias").toFloat();
	T_ratio_merge_planes = configIni.value("PlaneDetect/T_ratio_merge_planes").toFloat();
	alpha_poly = configIni.value("PlaneDetect/alpha_poly").toFloat();
	T_dist_point_plane = configIni.value("PlaneDetect/T_dist_point_plane").toFloat();
	T_cluster_num = configIni.value("PlaneDetect/T_cluster_num").toInt();
	cout << "参数加载完成" << endl;

	kdtree_source.setInputCloud(source_cloud);
}
// 估计法向量
void PCLViewer::on_normalEstimateAction_triggered()
{
	estimateNormal();
	memset(state, 0, CMD_SIZE);
	strcpy(state, "estimate normal");
}
// 校正点云法向量
void PCLViewer::on_regulateNormalAction_triggered()
{
	memset(state, 0, CMD_SIZE);
	strcpy(state, "regulate normal");
	RegulateNormalDialog.show();
}
// 执行法向量校正操作
void PCLViewer::on_performRegulateAction_triggered()
{
	is_norm_direction_valid = RegulateNormalDialog.is_norm_direction_valid;
	regulateNormal();
}
// 创建参数空间
void PCLViewer::on_createPSAction_triggered()
{
	start = std::clock();
	createPS();
	voteForPS();
	finish = std::clock();
	cout << "parameter space has been constructed, ps has " << ps_cloud->size() << " elements." << endl;
	cout << "time used: " << finish - start << "ms." << endl;

	memset(state, 0, CMD_SIZE);
	strcpy(state, "create ps");
	processStateMsg();
}
// 保存法向量点云文件
void PCLViewer::on_savePointNormalFileAction_triggered()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Open PCD Files"),
		"",
		tr("PCD Files (*.pcd)"));

	if (!fileName.isNull())
	{
		pcl::PointCloud<pcl::PointNormal>::Ptr cloud(new pcl::PointCloud<pcl::PointNormal>);
		pcl::PointNormal p;
		for (int i = 0; i < source_cloud->size(); ++i)
		{
			p.x = source_cloud->points[i].x;
			p.y = source_cloud->points[i].y;
			p.z = source_cloud->points[i].z;
			p.normal_x = source_normal->points[i].normal_x;
			p.normal_y = source_normal->points[i].normal_y;
			p.normal_z = source_normal->points[i].normal_z;
			cloud->push_back(p);
		}
		pcl::io::savePCDFile(fileName.toStdString(), *cloud);
		cout << "file:" << fileName.toStdString() << " has been saved." << endl;
	}
}
// 打开带有法向量的点云文件
void PCLViewer::on_openPointCloudNormalFileAction_triggered()
{
	on_load_param_action_triggered();
	QString fileName = QFileDialog::getOpenFileName(this, tr("open file"), " ", tr("pcdFiles(*.pcd)"));
	if (fileName == "") return;
	strcpy(pcdFileName, fileName.toStdString().data());
	loadPointNormal();
	cout << "loaded " << source_cloud->size() << " points." << endl;
}
// 平滑参数空间
void PCLViewer::on_filtPSAction_triggered()
{
	start = std::clock();
	filtPS();
	finish = std::clock();
	cout << "filt used: " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "filt ps");
	processStateMsg();
}
// 分割参数空间
void PCLViewer::on_segPSAction_triggered()
{
	start = std::clock();
	createGradientField();
	rectifySinkFields();
	finish = std::clock();
	cout << "segment ps used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "segment ps");
	processStateMsg();
}
// 分割平面
void PCLViewer::on_segPlaneAction_triggered()
{
	start = std::clock();
	segmentPlanes();
	finish = std::clock();
	cout << "plane segmentation used " << finish - start << "ms." << endl;
	cout << "plane number = " << plane_clouds.size() << endl;

	// 设置增长单元集合
	for (int i = 0; i < plane_clouds.size(); ++i)
	{
		Plane plane;
		PointCloudT::Ptr cloud(new PointCloudT);
		*cloud = *plane_clouds[i].points_set;
		plane.points_set = cloud;
		growth_unit_set.push_back(plane);
	}

	memset(state, 0, CMD_SIZE);
	strcpy(state, "segment planes");
	processStateMsg();
}
// 区域生长
void PCLViewer::on_regionGrowingAction_triggered()
{
	start = std::clock();
	growPlaneArea();
	finish = std::clock();
	cout << "plane growing used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "grow planes");
	processStateMsg();
}
// 平面合并
void PCLViewer::on_mergePlanesAction_triggered()
{
	start = std::clock();
	mergePlanes();
	finish = std::clock();
	cout << "plane merging used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "merge planes");
	processStateMsg();
}

//define by czh
//读入平面，方便进行平面多边形化
void PCLViewer::on_OpenPlanesAction_triggered()
{
	QString QfileName = QFileDialog::getOpenFileName(this, tr("Open PCD Files"), "/", tr("PCD Files (*.pcd)"));
	// 存有所有的平面
	std::string vertices_fileName = QfileName.toStdString();
	pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::io::loadPCDFile(vertices_fileName, *vertices);
	cout << "total vertices size = " << vertices->size() << endl;

	string fileName = vertices_fileName;
	for (int i = 0; i < 4; ++i) {
		fileName.pop_back();
	}
	// 平面法向量
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	string normals_fileName = fileName;
	normals_fileName.append("_polyNormal.pcd");
	pcl::io::loadPCDFile(normals_fileName, *normals);
	int poly_size = normals->size();
	cout << "poly size = " << poly_size << endl;

	// 多边形顶点数量
	char buf[64];
	memset(buf, 0, 64);
	std::ifstream file;
	string v_size_fileName = fileName;
	v_size_fileName.append("_polySize.txt");
	memset(buf, 0, 64);
	vector<int> poly_v_size;
	poly_v_size.reserve(poly_size);
	memset(buf, 0, 64);
	file.open(v_size_fileName, std::ios::in);
	while (file.getline(buf, 64)) {
		int size = atoi(buf);
		poly_v_size.push_back(size);
		memset(buf, 0, 64);
	}
	file.clear();
	file.close();

	// 初始化多边形数据
	plane_clouds_final.resize(poly_size);
	int base_count = 0;
	for (int i = 0; i < poly_size; ++i) {
		plane_clouds_final[i].points_set->reserve(poly_v_size[i]);
		for (int j = 0; j < poly_v_size[i]; ++j) {
			plane_clouds_final[i].points_set->push_back(vertices->points[base_count + j]);
		}

		plane_clouds_final[i].coeff.values.push_back(normals->points[i].normal_x);
		plane_clouds_final[i].coeff.values.push_back(normals->points[i].normal_y);
		plane_clouds_final[i].coeff.values.push_back(normals->points[i].normal_z);
		base_count += poly_v_size[i];
	}
}
// 平面多边形化
void PCLViewer::on_polyPlanesAction_triggered()
{
	start = std::clock();
	polyPlanes();
	finish = std::clock();
	cout << "plane polygonization used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "poly planes");
	processStateMsg();
}
// 点云后处理
void PCLViewer::on_postProcessAction_triggered()
{
	start = std::clock();
	postProcessPlanes();
	finish = std::clock();
	cout << "plane post process used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "postProcess planes");
	processStateMsg();
}
// 继续运行
void PCLViewer::on_runAgainAction_triggered()
{
	// 读参数
	on_load_param_action_triggered();

	// 重新估计法向量
	estimateNormal();
	// 校正法向量
	regulateNormal();

	// 霍夫变换
	// 创建参数空间
	start = std::clock();
	createPS();
	voteForPS();
	finish = std::clock();
	cout << "parameter space has been constructed, ps has " << ps_cloud->size() << " elements." << endl;
	cout << "time used: " << finish - start << "ms." << endl;
	// 过滤参数空间
	start = std::clock();
	filtPS();
	finish = std::clock();
	cout << "filt used: " << finish - start << "ms." << endl;
	// 分割参数空间
	start = std::clock();
	createGradientField();
	rectifySinkFields();
	finish = std::clock();
	cout << "segment ps used " << finish - start << "ms." << endl;
	// 分割增长单元
	start = std::clock();
	segmentPlanes();
	finish = std::clock();
	cout << "plane segmentation used " << finish - start << "ms." << endl;
	cout << "plane number = " << plane_clouds.size() << endl;
	// 区域生长
	start = std::clock();
	growPlaneArea();
	finish = std::clock();
	cout << "plane growing used " << finish - start << "ms." << endl;
	// 平面合并
	start = std::clock();
	mergePlanes();
	finish = std::clock();
	cout << "plane merging used " << finish - start << "ms." << endl;
	// 多边形化
	start = std::clock();
	polyPlanes();
	finish = std::clock();
	cout << "plane polygonization used " << finish - start << "ms." << endl;
	// 后处理
	start = std::clock();
	postProcessPlanes();
	finish = std::clock();
	cout << "plane post process used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "postProcess planes");
	processStateMsg();
}
// 自动执行
void PCLViewer::on_autoPerformAction_triggered()
{
	// 重新估计法向量
	estimateNormal();

	// 霍夫变换
	// 创建参数空间
	start = std::clock();
	// 创建参数空间
	createPS();
	// 向参数空间进行投票
	voteForPS();
	finish = std::clock();
	cout << "parameter space has been constructed, ps has " << ps_cloud->size() << " elements." << endl;
	cout << "time used: " << finish - start << "ms." << endl;
	// 过滤参数空间
	start = std::clock();
	filtPS();
	finish = std::clock();
	cout << "filt used: " << finish - start << "ms." << endl;
	// 分割参数空间
	start = std::clock();
	createGradientField();
	rectifySinkFields();
	finish = std::clock();
	cout << "segment ps used " << finish - start << "ms." << endl;
	// 分割增长单元
	start = std::clock();
	segmentPlanes();
	finish = std::clock();
	cout << "plane segmentation used " << finish - start << "ms." << endl;
	cout << "plane number = " << plane_clouds.size() << endl;
	// 区域生长
	start = std::clock();
	growPlaneArea();
	finish = std::clock();
	cout << "plane growing used " << finish - start << "ms." << endl;
	// 平面合并
	start = std::clock();
	mergePlanes();
	finish = std::clock();
	cout << "plane merging used " << finish - start << "ms." << endl;
	// 多边形化
	start = std::clock();
	polyPlanes();
	finish = std::clock();
	cout << "plane polygonization used " << finish - start << "ms." << endl;
	// 后处理
	start = std::clock();
	postProcessPlanes();
	finish = std::clock();
	cout << "plane post process used " << finish - start << "ms." << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "postProcess planes");
	processStateMsg();
}
// 删除指定的多边形
// 处理对象：plane_clouds_final, g_selected_point, poly_centroids_cloud g_is_poly_del kdtree_poly_centroids_cloud
void PCLViewer::on_editPolyAction_triggered()
{
	// step0: 提示用户
	cout << "进入多边形删除模式" << endl;
	// step1: 提示点选取回调函数
	memset(state, 0, CMD_SIZE);
	strcpy(state, "del poly");
	// step2: 初始化数据
	poly_centroids_cloud->clear();
	g_is_poly_del.clear();
	poly_centroids_cloud->reserve(plane_clouds_final.size());
	g_is_poly_del.reserve(plane_clouds_final.size());
	// step3: 设置全局数据
	pcl::PointXYZ p;
	int size;
	for (int i = 0; i < plane_clouds_final.size(); ++i)
	{
		g_is_poly_del.push_back(false);
		p.x = p.y = p.z = 0;
		size = plane_clouds_final[i].border->size();
		for (int j = 0; j < size; ++j)
		{
			p.x += plane_clouds_final[i].border->points[j].x;
			p.y += plane_clouds_final[i].border->points[j].y;
			p.z += plane_clouds_final[i].border->points[j].z;
		}
		p.x /= float(size);
		p.y /= float(size);
		p.z /= float(size);
		poly_centroids_cloud->push_back(p);
	}
	kdtree_poly_centroids_cloud.setInputCloud(poly_centroids_cloud);
	// step4: 显示数据
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	// 重心显示为绿色
	viewer->addPointCloud(poly_centroids_cloud, "centroids");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, "centroids");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "centroids");
	for (int i = 0; i < plane_clouds_final.size(); ++i)
	{
		QString cloud_id;
		cloud_id = QString::number(i, 10);
		viewer->addPointCloud(plane_clouds_final[i].border, cloud_id.toStdString());
	}

	g_selected_poly_id = -1;
}
// 确定删除多边形
// g_selected_poly_id
void PCLViewer::on_delPolyAction_triggered()
{
	g_is_poly_del[g_selected_poly_id] = true;
	cout << "被删除的多边形id = " << g_selected_poly_id << endl;
	for (int i = 0; i < g_is_poly_del.size(); ++i)
	{
		if (g_is_poly_del[i])
		{
			QString cloud_id = QString::number(i, 10);
			viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, cloud_id.toStdString());
			viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, cloud_id.toStdString());
			ui->qvtkWidget->update();
			// cout << cloud_id.toStdString() << endl;
		}
	}
}
// 执行删除操作
void PCLViewer::on_performDelAction_triggered()
{
	std::vector<int> offset(plane_clouds_final.size(), 0);
	int count = 0;
	for (int i = 0; i < plane_clouds_final.size(); ++i)
	{
		if (g_is_poly_del[i])
		{
			++count;
			offset[i] = -1;
			continue;
		}
		offset[i] = count;
	}

	for (int i = 0; i < offset.size(); ++i)
	{
		if (offset[i] == -1) continue;
		if (offset[i] > 0)
		{
			plane_clouds_final[i - offset[i]] = plane_clouds_final[i];
		}
	}

	for (int i = 0; i < count; ++i)
	{
		plane_clouds_final.pop_back();
	}

	// 更新视图
	on_editPolyAction_triggered();
}

// 保存多边形数据
void PCLViewer::on_savePolyDataAction_triggered()
{
	QString fileName = QFileDialog::getSaveFileName(this,tr("Open PCD Files"),"",tr("PCD Files (*.pcd)"));

	if (!fileName.isNull())
	{
		// 保存顶点序列
		PointCloudT::Ptr vertices(new PointCloudT);
		for (int ii = 0; ii < plane_clouds_final.size(); ++ii)
		{
			for (int i = 0; i < plane_clouds_final[ii].border->size(); ++i)
			{
				vertices->push_back(plane_clouds_final[ii].border->points[i]);
			}
		}
		pcl::io::savePCDFileASCII(fileName.toStdString(), *vertices);
		// 保存每个多边形的size
		int length = fileName.length();
		fileName = fileName.mid(0, length - 3);
		fileName += "_polySize.txt";
		std::ofstream file(fileName.toStdString(), std::ios::out);

		for (int ii = 0; ii < plane_clouds_final.size(); ++ii)
		{
			file << plane_clouds_final[ii].border->size() << endl;
		}
		// 保存每个多边形的法向量
		pcl::PointCloud<pcl::Normal>::Ptr plane_normals(new pcl::PointCloud<pcl::Normal>);
		pcl::Normal pn;
		for (int i = 0; i < plane_clouds_final.size(); ++i)
		{
			Plane plane = plane_clouds_final[i];
			pn.normal_x = plane.coeff.values[0];
			pn.normal_y = plane.coeff.values[1];
			pn.normal_z = plane.coeff.values[2];
			plane_normals->push_back(pn);
		}
		int length1 = fileName.length();
		fileName = fileName.mid(0, length - 3);
		fileName += "_polyNormal.txt";
		pcl::io::savePCDFileASCII(fileName.toStdString(), *plane_normals);
		cout << "all " << plane_clouds_final.size() << " polys' data have been saved." << endl;
		//保存每个多边形的Scale
		int length2 = fileName.length();
		fileName = fileName.mid(0, length - 3);
		fileName += "_polyScale.txt";
		std::ofstream file1(fileName.toStdString(), std::ios::out);
		for (int scale = 0; scale < plane_clouds_final.size(); ++scale)
		{
			file1 << r_for_estimate_normal << endl;
		}
	}
}
// 进入多边形修剪模式
void PCLViewer::on_enterPruneModeAction_triggered()
{
	// step0: 提示用户
	cout << "进入多边形修剪模式" << endl;
	// step1: 提示点选取回调函数
	memset(state, 0, CMD_SIZE);
	strcpy(state, "prune poly");
	// step2: 初始化数据
	poly_centroids_cloud->clear();
	g_is_poly_del.clear();
	poly_centroids_cloud->reserve(plane_clouds_final.size());
	g_is_poly_del.reserve(plane_clouds_final.size());
	// step3: 设置全局数据
	pcl::PointXYZ p;
	int size;
	for (int i = 0; i < plane_clouds_final.size(); ++i)
	{
		g_is_poly_del.push_back(false);
		p.x = p.y = p.z = 0;
		size = plane_clouds_final[i].border->size();
		for (int j = 0; j < size; ++j)
		{
			p.x += plane_clouds_final[i].border->points[j].x;
			p.y += plane_clouds_final[i].border->points[j].y;
			p.z += plane_clouds_final[i].border->points[j].z;
		}
		p.x /= float(size);
		p.y /= float(size);
		p.z /= float(size);
		poly_centroids_cloud->push_back(p);
	}
	kdtree_poly_centroids_cloud.setInputCloud(poly_centroids_cloud);
	// step4: 显示数据
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	// 重心显示为绿色
	viewer->addPointCloud(poly_centroids_cloud, "centroids");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, "centroids");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "centroids");
	for (int i = 0; i < plane_clouds_final.size(); ++i)
	{
		QString cloud_id;
		cloud_id = QString::number(i, 10);
		viewer->addPointCloud(plane_clouds_final[i].border, cloud_id.toStdString());
	}

	g_selected_poly_id = -1;
	// 提示用户选择某个多边形
	cout << "通过选中重心选择特定多边形" << endl;
}
// 选中当前多边形
// 仅显示当前多边形
void PCLViewer::on_selCurPolyAction_triggered()
{
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	viewer->addPointCloud(plane_clouds_final[g_selected_poly_id].border, "border");
	// 提示用户选择顶点
	cout << "可以开始选择第一个顶点" << endl;
	memset(state, 0, CMD_SIZE);
	strcpy(state, "prune poly select points");
}
// 设置选中点为第一个顶点
// 选中的点为绿色
void PCLViewer::on_setFirstPointAction_triggered()
{
	first_point = g_selected_point;
	PointCloudT::Ptr cloud(new PointCloudT);
	cloud->push_back(first_point);
	viewer->removePointCloud("point");
	viewer->removePointCloud("first point");
	viewer->addPointCloud(cloud, "first point");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, "first point");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "first point");
	cout << "第一个点已选中" << endl;
}
// 设置选中点为第二个顶点
void PCLViewer::on_setSecondPointAction_triggered()
{
	second_point = g_selected_point;
	PointCloudT::Ptr cloud(new PointCloudT);
	cloud->push_back(second_point);
	viewer->removePointCloud("point");
	viewer->removePointCloud("second point");
	viewer->addPointCloud(cloud, "second point");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 255, 0, "second point");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "second point");
	cout << "第二个点已选中" << endl;
}
// 执行多边形切分操作
// 即，将一个多边形通过选择的两点构成的直线段把原多边形切分成两个多边形
void PCLViewer::on_performPolyCutAction_triggered()
{
	// step0: 分别找到第一个点和第二个点在原来多边形中的索引
	int first_index = -1;
	int second_index = -1;
	for (int i = 0; i < plane_clouds_final[g_selected_poly_id].border->size(); ++i)
	{
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == first_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == first_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == first_point.z) {
			first_index = i;
		}
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == second_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == second_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == second_point.z) {
			second_index = i;
		}
		if (first_index != -1 && second_index != -1)
			break;
	}
	cout << "first_index = " << first_index << "; second_index = " << second_index << endl;
	// step1: 做好备份
	PointCloudT::Ptr cloud_backup(new PointCloudT);
	*cloud_backup = *plane_clouds_final[g_selected_poly_id].border;
	plane_clouds_final[g_selected_poly_id].border->clear();

	// step2: 原有多边形：first_index->second_index
	int cur_index = first_index;
	while (true)
	{
		plane_clouds_final[g_selected_poly_id].border->push_back(cloud_backup->points[cur_index]);

		if (cur_index == second_index) break;

		cur_index = cur_index == cloud_backup->size() - 1 ? 0 : cur_index + 1;
	}

	cout << "ok1" << endl;

	// step3: 增加一个多边形:second_index -> first_index
	Plane plane;
	plane.border.reset(new PointCloudT);
	plane.border->reserve(cloud_backup->size() - plane_clouds_final[g_selected_poly_id].border->size());
	cur_index = second_index;
	while (true)
	{
		plane.border->push_back(cloud_backup->points[cur_index]);
		if (cur_index == first_index) break;
		cur_index = cur_index == cloud_backup->size() - 1 ? 0 : cur_index + 1;
	}
	plane.coeff = plane_clouds_final[g_selected_poly_id].coeff;
	plane_clouds_final.push_back(plane);

	cout << "ok2" << endl;

	// step4: 更新视图
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();
	on_enterPruneModeAction_triggered();
}
// 显示要删除的线段
void PCLViewer::on_displayLineSegAction_triggered()
{
	// step0: 分别找到第一个点和第二个点在原来多边形中的索引
	int first_index = -1;
	int second_index = -1;
	for (int i = 0; i < plane_clouds_final[g_selected_poly_id].border->size(); ++i)
	{
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == first_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == first_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == first_point.z) {
			first_index = i;
		}
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == second_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == second_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == second_point.z) {
			second_index = i;
		}
		if (first_index != -1 && second_index != -1)
			break;
	}
	cout << "first_index = " << first_index << "; second_index = " << second_index << endl;

	// step1: 找到从first_point 到 second_point之间的点
	PointCloudT::Ptr line_points(new PointCloudT);
	int cur_index = first_index == plane_clouds_final[g_selected_poly_id].border->size() - 1 ?
		0 : first_index + 1;
	while (cur_index != second_index)
	{
		line_points->push_back(plane_clouds_final[g_selected_poly_id].border->points[cur_index]);
		cur_index = cur_index == plane_clouds_final[g_selected_poly_id].border->size() - 1 ?
			0 : cur_index + 1;
	}
	viewer->removePointCloud("line seg");
	viewer->addPointCloud(line_points, "line seg");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 255, 0, 0, "line seg");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "line seg");
}
// 切换线段
void PCLViewer::on_switchLineSegAction_triggered()
{
	pcl::PointXYZ tmp_point = first_point;
	first_point = second_point;
	second_point = tmp_point;
	on_displayLineSegAction_triggered();
}
// 执行删除
void PCLViewer::on_performLineSegDelAction_triggered()
{
	// step0: 分别找到第一个点和第二个点在原来多边形中的索引
	int first_index = -1;
	int second_index = -1;
	for (int i = 0; i < plane_clouds_final[g_selected_poly_id].border->size(); ++i)
	{
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == first_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == first_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == first_point.z) {
			first_index = i;
		}
		if (plane_clouds_final[g_selected_poly_id].border->points[i].x == second_point.x &&
			plane_clouds_final[g_selected_poly_id].border->points[i].y == second_point.y &&
			plane_clouds_final[g_selected_poly_id].border->points[i].z == second_point.z) {
			second_index = i;
		}
		if (first_index != -1 && second_index != -1)
			break;
	}

	// step1: 将[second_index first_index]闭区间的点保存起来
	std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>> line_points;
	int cur_index = second_index;
	while (true)
	{
		line_points.push_back(plane_clouds_final[g_selected_poly_id].border->points[cur_index]);
		if (cur_index == first_index) break;
		cur_index = cur_index == plane_clouds_final[g_selected_poly_id].border->size() - 1 ?
			0 : cur_index + 1;
	}
	plane_clouds_final[g_selected_poly_id].border->clear();
	for (int i = 0; i < line_points.size(); ++i)
	{
		plane_clouds_final[g_selected_poly_id].border->push_back(line_points[i]);
	}

	on_selCurPolyAction_triggered();
}
// 加载多边形数据
void PCLViewer::on_loadPolyDataAction_triggered()
{

}
/****************************************前半段：结束**************************************************/

// 加载参数
void PCLViewer::loadParams(){
    QSettings configIniRead("config.ini", QSettings::IniFormat);
    T_angle_betPointAndPoint = configIniRead.value("simplify vertices size/T_angle_betPointAndPoint").toFloat();
    T_angle_betPointAndLine = configIniRead.value("simplify vertices size/T_angle_betPointAndLine").toFloat();
    T_dist_betPointAndLine = configIniRead.value("simplify vertices size/T_dist_betPointAndLine").toFloat();
    T_lenRatio = configIniRead.value("simplify vertices size/T_lenRatio").toFloat();
    T_angle_betLines = configIniRead.value("simplify vertices size/T_angle_betLines").toFloat();
    m_svs.setParams(T_angle_betPointAndPoint, T_angle_betPointAndLine, T_dist_betPointAndLine,
                    T_lenRatio, T_angle_betLines);

    ratio_of_scale=configIniRead.value("merge vertices/ratio_of_scale").toFloat();
    T_ratio_lineSeg_middle_area=configIniRead.value("merge vertices/T_ratio_lineSeg_middle_area").toFloat();
    T_dist_proj_to_lineSegEnd=configIniRead.value("merge vertices/T_dist_proj_to_lineSegEnd").toFloat();
    T_maxAngle_bet_two_polys=configIniRead.value("merge vertices/T_maxAngle_bet_two_polys").toFloat();
    m_sv.setParams(ratio_of_scale, T_ratio_lineSeg_middle_area, T_dist_proj_to_lineSegEnd,
                   T_maxAngle_bet_two_polys);

    osnap_delta=configIniRead.value("merge vertices/osnap_delta").toDouble();
    osnap_max_k_count=configIniRead.value("merge vertices/osnap_max_k_count").toInt();
    m_osnap.delta = osnap_delta;
    m_osnap.max_k_count = osnap_max_k_count;

	m_osnap.f1_weight=configIniRead.value("merge vertices/f1_weight").toDouble();
	m_osnap.f2_weight=configIniRead.value("merge vertices/f2_weight").toDouble();
	m_osnap.f3_1_weight=configIniRead.value("merge vertices/f3_1_weight").toDouble();
	m_osnap.f3_2_weight=configIniRead.value("merge vertices/f3_2_weight").toDouble();
	m_osnap.f8_weight=configIniRead.value("merge vertices/f8_weight").toDouble();
	m_osnap.f4_weight=configIniRead.value("merge vertices/f4_weight").toDouble();
	m_osnap.f5_weight=configIniRead.value("merge vertices/f5_weight").toDouble();
	m_osnap.f6_weight=configIniRead.value("merge vertices/f6_weight").toDouble();
	m_osnap.f7_weight=configIniRead.value("merge vertices/f7_weight").toDouble();
}

// 打开文件
void PCLViewer::on_OpenAction_triggered()
{
    initial_polys.clear();
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    QString QfileName = QFileDialog::getOpenFileName(this, tr("Open PCD Files"), "/", tr("PCD Files (*.pcd)"));
    // std_fileName: 存有所有的多边形顶点
    std::string vertices_fileName = QfileName.toStdString();
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::io::loadPCDFile(vertices_fileName, *vertices);
    cout << "total vertices size = " << vertices->size() << endl;

    string fileName = vertices_fileName;
    for(int i = 0; i < 4; ++i){
        fileName.pop_back();
    }
    // 多边形法向量
    pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
    string normals_fileName = fileName;
    normals_fileName.append("_polyNormal.pcd");
    pcl::io::loadPCDFile(normals_fileName, *normals);
    int poly_size = normals->size();
    cout << "poly size = " << poly_size << endl;

    // 多边形尺度
    string scale_fileName = fileName;
    scale_fileName.append("_polyScale.txt");
    vector<float> poly_scales;
    poly_scales.reserve(poly_size);
    char buf[64];
    memset(buf, 0, 64);
    std::ifstream file;
    file.open(scale_fileName, std::ios::in);
    while(file.getline(buf, 64)){
        float scale = atof(buf);
        poly_scales.push_back(scale);
        memset(buf, 0, 64);
    }
    file.clear();
    file.close();
    // 整理多边形尺度
    this->ply_scales.clear();
    for(int i = 0; i < poly_scales.size(); ++i){
        if(i == 0){
            ply_scales.push_back(poly_scales[i]);
        }else{
            float diff = poly_scales[i] - ply_scales[ply_scales.size()-1];
            if(abs(diff) > 0.0001f){
                ply_scales.push_back(poly_scales[i]);
            }
        }
    }
    // 打印多边形的所有尺度：
    cout << "number of polygon scales: " << ply_scales.size() << endl;
    for(int i = 0; i < ply_scales.size(); ++i){
        cout << i << ": " << ply_scales[i] << endl;
    }

    // 多边形顶点数量
    string v_size_fileName = fileName;
    v_size_fileName.append("_polySize.txt");
    memset(buf, 0, 64);
    vector<int> poly_v_size;
    poly_v_size.reserve(poly_size);
    memset(buf, 0, 64);
    file.open(v_size_fileName, std::ios::in);
    while(file.getline(buf, 64)){
        int size = atoi(buf);
        poly_v_size.push_back(size);
        memset(buf, 0, 64);
    }
    file.clear();
    file.close();

    /////////////////////////////////////////////////////////////
    // 初始化多边形数据
    initial_polys.resize(poly_size);
    int base_count = 0;
    for(int i = 0; i < poly_size; ++i){
        initial_polys[i].vertices->reserve(poly_v_size[i]);
        for(int j = 0; j < poly_v_size[i]; ++j){
            initial_polys[i].vertices->push_back(vertices->points[base_count + j]);
        }
        initial_polys[i].normal << normals->points[i].normal_x,
                                   normals->points[i].normal_y,
                                   normals->points[i].normal_z;
        initial_polys[i].scale = poly_scales[i];

        base_count += poly_v_size[i];
    }

    // test code:
    srand(time(NULL));
    for(int i = 0; i < poly_size; ++i){
        memset(buf, 0, 64);
        itoa(i, buf, 10);
        string id = buf;
        // cout << "id = " << id << endl;
        viewer->addPointCloud(initial_polys[i].vertices, id);
        double r, g, b;
        r = rand() % 255;
        g = rand() % 255;
        b = rand() % 255;
        r /= 255;
        g /= 255;
        b /= 255;
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, id);
    }    
}

void PCLViewer::init()
{
	viewer->registerKeyboardCallback(keyboardEventOccurred, (void *)this);
	viewer->registerPointPickingCallback(pointPickingEventOccurred, (void *)this);
	viewer->registerAreaPickingCallback(areaPickingEventOccurred, (void *)this);
}

// 显示背景色
void PCLViewer::on_setBgColorAction_triggered()
{
    SetBGColorDialog dlg;
	if (dlg.exec() != QDialog::Accepted) return;

    int r, g, b;
	dlg.getRGB(r, g, b);
    viewer->setBackgroundColor(r, g, b);
}

// 显示或隐藏坐标轴
bool is_show_coor = true;
void PCLViewer::on_setCoordinateAxesAction_triggered()
{
    if(is_show_coor){
        viewer->addCoordinateSystem(1.0f);
    }else{
        viewer->removeCoordinateSystem();
    }
    is_show_coor = !is_show_coor;
}

// 移除指定点云
void PCLViewer::on_removePointCloudAction_triggered()
{
}

// 清屏
void PCLViewer::on_clearScreenAction_triggered()
{
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
}

// 估计顶点法向量
void PCLViewer::on_verticesNormalEstimationAction_triggered()
{
    // computeVerticesNormals(poly.vertices_cloud, poly.vertices_normals, poly.poly_normal);
    for(int i = 0; i < initial_polys.size(); ++i){
        this->m_svs.computeVerticesNormals(initial_polys[i].vertices, initial_polys[i].vertices_normals, initial_polys[i].normal);
    }
    pcl::PointCloud<pcl::Normal>::Ptr normal_cloud (new pcl::PointCloud<pcl::Normal>);
    for(int i = 0; i < initial_polys.size(); ++i){
        for(int j = 0; j < initial_polys[i].vertices_normals.size(); ++j){
            pcl::Normal pn;
            pn.normal_x = initial_polys[i].vertices_normals[j][0];
            pn.normal_y = initial_polys[i].vertices_normals[j][1];
            pn.normal_z = initial_polys[i].vertices_normals[j][2];
            normal_cloud->push_back(pn);
        }
    }
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i = 0; i < initial_polys.size(); ++i){
        for(int j = 0; j < initial_polys[i].vertices->size(); ++j){
            cloud->push_back(initial_polys[i].vertices->points[j]);
        }
    }
    viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(cloud, normal_cloud, 1, 0.2);
}

// 查看特定多边形的顶点法向量
void PCLViewer::on_ObservePolyVerNormalAction_triggered()
{
    RemovePointCloudDialog dlg;
    dlg.exec();
    // 多边形编号: poly_id
    int poly_id = dlg.id.toInt();

    // 平移到坐标原点
    float sum_x = 0;
    float sum_y = 0;
    float sum_z = 0;
    InitialPoly *pPoly = &initial_polys[poly_id];
    for(int i = 0; i < pPoly->vertices->size(); ++i){
        sum_x += pPoly->vertices->points[i].x;
        sum_y += pPoly->vertices->points[i].y;
        sum_z += pPoly->vertices->points[i].z;
    }
    float translate_x = sum_x / pPoly->vertices->size();
    float translate_y = sum_y / pPoly->vertices->size();
    float translate_z = sum_z / pPoly->vertices->size();

    // 平移后的顶点点云
    pcl::PointXYZ p;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i = 0; i < pPoly->vertices->size(); ++i){
        p.x = pPoly->vertices->points[i].x - translate_x;
        p.y = pPoly->vertices->points[i].y - translate_y;
        p.z = pPoly->vertices->points[i].z - translate_z;
        cloud->push_back(p);
    }
    // 法向量点云
    pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
    pcl::Normal pn;
    for(int i = 0; i < pPoly->vertices_normals.size(); ++i){
        pn.normal_x = pPoly->vertices_normals[i][0];
        pn.normal_y = pPoly->vertices_normals[i][1];
        pn.normal_z = pPoly->vertices_normals[i][2];
        normals->push_back(pn);
    }

    //pcl::visualization::PCLVisualizer viewer;
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
    viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(cloud, normals, 1, 0.5f, "cloud normal");
    viewer->addPointCloud(cloud, "cloud");
    // viewer->addCoordinateSystem(1.0f);
}

// 化简顶点数量参数配置
void PCLViewer::on_simplifyVerticesSizeAction_triggered()
{
    SetParams_SimpilyVerticesSize_Dialog dlg;
    dlg.exec();
    this->loadParams(); // 更新参数
}

// 区域生长提取连续线段
void PCLViewer::on_RGtoDetectLineSegsAction_triggered()
{
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly *pPoly = & initial_polys[i];
        // 首先基于区域生长计算初始的线段序列
        m_svs.constructInitLineSegs(pPoly->vertices, pPoly->vertices_normals,
                                    pPoly->normal, pPoly->line_segs);
        // 优化多边形的顶点序列
        m_svs.refineLineSegs(pPoly->vertices, pPoly->simplified_vertices,
                             pPoly->line_segs, pPoly->final_line_segs);

        //cout << pPoly->vertices->size() - pPoly->simplified_vertices->size() << endl;
    }

    // 显示提取效果
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly *pPoly = & initial_polys[i];
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            *vertices += *pPoly->simplified_vertices;
        }
    }
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");

    // 线段：
    string line_name;
    int count = 0;
    char buf[64];
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly *pPoly = & initial_polys[i];
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            // 线段名：str(count)
            memset(buf, 0, 64);
            _itoa(count++, buf, 10);
            line_name.clear();
            line_name.append(buf);

            int next = j == pPoly->simplified_vertices->size()-1 ? 0 : j+1;
            viewer->addLine(pPoly->simplified_vertices->points[j], pPoly->simplified_vertices->points[next],
                            1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }
}

// 显示特定多边形的顶点化简结果
void PCLViewer::on_RGlineSegDisplayAction_triggered()
{
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    RemovePointCloudDialog dlg;
    dlg.exec();
    // 多边形编号: poly_id
    int poly_id = dlg.id.toInt();

    // 平移到坐标原点
    float sum_x = 0;
    float sum_y = 0;
    float sum_z = 0;
    InitialPoly *pPoly = &initial_polys[poly_id];
    for(int i = 0; i < pPoly->vertices->size(); ++i){
        sum_x += pPoly->vertices->points[i].x;
        sum_y += pPoly->vertices->points[i].y;
        sum_z += pPoly->vertices->points[i].z;
    }
    float translate_x = sum_x / pPoly->vertices->size();
    float translate_y = sum_y / pPoly->vertices->size();
    float translate_z = sum_z / pPoly->vertices->size();

    // 平移后的顶点点云
    pcl::PointXYZ p;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i = 0; i < pPoly->vertices->size(); ++i){
        p.x = pPoly->vertices->points[i].x - translate_x;
        p.y = pPoly->vertices->points[i].y - translate_y;
        p.z = pPoly->vertices->points[i].z - translate_z;
        cloud->push_back(p);
    }

    // 用不同颜色显示各条线段
    char buf[64];
    string line_name;
    int count = 0;
    srand(time(NULL));
    for(int i = 0; i < pPoly->final_line_segs.size(); ++i){
        int next = i == pPoly->final_line_segs.size()-1 ? 0 : i+1;

        int cur_index = pPoly->final_line_segs[i];
        int next_index = pPoly->final_line_segs[next];

        pcl::PointCloud<pcl::PointXYZ>::Ptr line_cloud (new pcl::PointCloud<pcl::PointXYZ>);
        int j = cur_index;
        while(j != next_index){
            line_cloud->push_back(cloud->points[j]);
            j = j == pPoly->vertices->size()-1 ? 0 : j+1;
        }
        memset(buf, 0, 64);
        _itoa(count++, buf, 10);
        line_name = buf;

        viewer->addPointCloud(line_cloud, line_name);
        double r = rand()%255/255.0;
        double g = rand()%255/255.0;
        double b = rand()%255/255.0;
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g, b, line_name);
    }

    cout << "poly id = " << poly_id << endl;
    cout << "original vertices size = " << pPoly->vertices->size() << endl;
    cout << "current vertices size = " << pPoly->simplified_vertices->size() << endl;
}

// 执行顶点化简
void PCLViewer::on_doSimplifyVerticesSizeAction_triggered()
{

    for(int i = 0; i < initial_polys.size(); ++i){
        // 第一步：估计顶点法向量
        this->m_svs.computeVerticesNormals(initial_polys[i].vertices, initial_polys[i].vertices_normals, initial_polys[i].normal);

        // 第二步：区域生长提取线段
        InitialPoly *pPoly = & initial_polys[i];
        m_svs.constructInitLineSegs(pPoly->vertices, pPoly->vertices_normals,
                                            pPoly->normal, pPoly->line_segs);
        // 第三步：优化线段序列
        m_svs.refineLineSegs(pPoly->vertices, pPoly->simplified_vertices,
                                     pPoly->line_segs, pPoly->final_line_segs);
    }

    // 显示化简之后的多边形顶点数量
    int v_count = 0;
    for(int i = 0; i < initial_polys.size(); ++i){
        v_count += initial_polys[i].simplified_vertices->size();
    }
    cout << "化简之后的模型顶点数量：" << v_count;

    // 显示线段
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly *pPoly = & initial_polys[i];
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            *vertices += *pPoly->simplified_vertices;
        }
    }
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, "vertices");

    // 线段：
    string line_name;
    int count = 0;
    char buf[64];
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly *pPoly = & initial_polys[i];
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            // 线段名：str(count)
            memset(buf, 0, 64);
            _itoa(count++, buf, 10);
            line_name.clear();
            line_name.append(buf);

            int next = j == pPoly->simplified_vertices->size()-1 ? 0 : j+1;
            viewer->addLine(pPoly->simplified_vertices->points[j], pPoly->simplified_vertices->points[next],
                            1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 移除所有点云
void PCLViewer::on_removeAllPointCloudsAction_triggered()
{
    viewer->removeAllPointClouds();
}

// 移除所有形状
void PCLViewer::on_removeAllShapesAction_triggered()
{
    viewer->removeAllShapes();
}

// 移除初始顶点集合
void PCLViewer::on_removeAllInitialVerticesAction_triggered()
{
    char buf[64];
    string id;
    for(int i = 0; i < initial_polys.size(); ++i){
        itoa(i, buf, 10);
        id = buf;
        viewer->removePointCloud(id);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 顶点合并参数设置
void PCLViewer::on_mergeVerticesParamsAction_triggered()
{
    SetParams_MergingVertices_Dialog dlg;
    dlg.exec();
    loadParams();
}

// 查看顶点的邻域半径
void PCLViewer::on_showNeighborRegionAction_triggered()
{
    int count = 0;
    char buf[64];
    string sphere_name;
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly* pPoly = &initial_polys[i];
        float r = this->ratio_of_scale * pPoly->scale;
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            sphere_name = "sphere";
            itoa(count++, buf, 10);
            sphere_name += buf;

            pcl::ModelCoefficients coeff;
            coeff.values.push_back(pPoly->simplified_vertices->points[j].x);
            coeff.values.push_back(pPoly->simplified_vertices->points[j].y);
            coeff.values.push_back(pPoly->simplified_vertices->points[j].z);
            coeff.values.push_back(r);
            viewer->addSphere(coeff, sphere_name);
        }
    }
}

// 获取点向直线的投影点
void getProjPoint(pcl::PointXYZ &p, pcl::PointXYZ &e1, pcl::PointXYZ &e2, pcl::PointXYZ &p_proj){
    Eigen::Vector3f p_base, line_dir;
    p_base << e1.x, e1.y, e1.z;
    line_dir << e2.x - e1.x,
                e2.y - e1.y,
                e2.z - e1.z;
    line_dir.normalize();

    Eigen::Vector3f v_p, delta;
    v_p << p.x, p.y, p.z;
    delta = p_base - v_p;
    float lambda = -1.0 * line_dir.dot(delta);

    Eigen::Vector3f p_proj_;
    p_proj_ = p_base + lambda * line_dir;
    p_proj.x = p_proj_[0];
    p_proj.y = p_proj_[1];
    p_proj.z = p_proj_[2];
}

// 查看匹配关系
bool is_first_run = true;
void PCLViewer::on_showMatchingAction_triggered()
{
    // 删除所有的圆球
    int count = 0;
    char buf[64];
    string sphere_name;
    for(int i = 0; i < initial_polys.size(); ++i){
        InitialPoly* pPoly = &initial_polys[i];
        for(int j = 0; j < pPoly->simplified_vertices->size(); ++j){
            sphere_name = "sphere";
            itoa(count++, buf, 10);
            sphere_name += buf;

            viewer->removeShape(sphere_name);
        }
    }

    // 用新的数据结构代表多边形
    if(is_first_run){
        this->m_sv.setPolygons(this->initial_polys);
        is_first_run = false;
    }else{
        // 记得调整每个顶点的搜索半径
        POLYGON *pPoly;
        Vertex *pv;
        for(int i = 0; i < m_sv.polygons.size(); ++i){
            pPoly = &m_sv.polygons[i];
            for(int j = 0; j < pPoly->getSize(); ++j){
                pv = j == 0 ? pPoly->start_point : pv->next_point;
                pv->radius = this->ratio_of_scale * initial_polys[i].scale;
            }
        }
    }
    // 为多边形的每个顶点设置邻域关系
    this->m_sv.setNeighborsForEachPoly();

    // 显示邻域关系
    // 点-点：红色线段
    // 点-边：白色线段
    count = 0;
    string line_name;
    pcl::PointXYZ p_proj;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            if(pv->neighbor_poly_id == -1) continue;    // 没有邻域关系

            line_name = "neigh";
            itoa(count++, buf, 10);
            line_name += buf;
            // 邻域关系：点-点
            if(pv->is_vertex_neighbor){
                viewer->addLine(pv->point, pv->neighbor_vertex->point, line_name);
                viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, line_name);
                viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
                //cout << "count = " << count << endl;
            }else{// 邻域关系：点-边
                getProjPoint(pv->point, pv->neighbor_edge->point, pv->neighbor_edge->next_point->point, p_proj);
				pcl::PointXYZ point = p_proj;	// 加一个局部变量，可避免出现abort错误
                viewer->addLine(pv->point, p_proj, line_name);
                viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, line_name);
                viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
                //cout << "count = " << count << endl;
            }
        }
    }
}

// 查看成簇的匹配关系（仅限点-点匹配）
void PCLViewer::on_showMatchingClustersAction_triggered()
{
    // 将单向的匹配关系转换为双向的匹配关系
    m_sv.setMatchingRelations();
    //cout << "ok1" << endl;
    std::vector<std::vector<Vertex *>> *clusters = &(m_sv.clusters);
    //cout << "ok2" << endl;
    m_sv.getMatchingVerticesClusters(*clusters);
    //cout << "ok3" << endl;
    // 显示簇
    cout << "cluster.size() = " << clusters->size() << endl;
    /*for(int i = 0; i < clusters.size(); ++i){
        cout << "cluster " << i << " has " << clusters[i].size() << " vertices." << endl;
    }*/
    char buf[64];
    string cloud_name;
    int count = 0;
    for(int i = 0; i < clusters->size(); ++i){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        for(int j = 0; j < (*clusters)[i].size(); ++j){
            Vertex *pv = (*clusters)[i][j];
            cloud->push_back(pv->point);
        }
        cloud_name = "clu";
        itoa(count++, buf, 10);
        cloud_name += buf;
        viewer->addPointCloud(cloud, cloud_name);
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 7, cloud_name);
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, cloud_name);
    }
}

// 合并点-点匹配关系
void PCLViewer::on_mergeMatchingRelationsAction_triggered()
{
    m_sv.refineMachingRelations();
    m_sv.mergeMatchingVerticesClusters(m_sv.clusters);

    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    // 所有的顶点
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 显示线段
    int count = 0;
    string line_name;
    char buf[64];
    Vertex *pv;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pcl::PointXYZ p1, p2;
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "line";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);

            vertices->push_back(p1);
        }
    }

    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, "vertices");
    // 显示匹配关系
    // “点-点”红色线段
    // “点-边”白色线段
    count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            if(pv->neighbor_poly_id != -1 && !pv->is_updated){
                if(pv->is_vertex_neighbor){
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    p2 = pv->neighbor_vertex->point;
                    //viewer->addLine(p1, p2, 0, 1.0, 1.0, line_name);
                    //viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
                }else{
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    getProjPoint(p1, pv->neighbor_edge->point, pv->neighbor_edge->next_point->point, p2);
                    viewer->addLine(p1, p2, 1.0, 1.0, 1.0, line_name);
                    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
                }
            }
        }
    }
}

// 合并点-边匹配关系
void PCLViewer::on_mergeVertexEdgeMatchingAction_triggered()
{
    // 处理点-边匹配
    m_sv.dealWithVertexEdgeMatchingRelation();

    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    // 所有的顶点
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 显示线段
    int count = 0;
    string line_name;
    char buf[64];
    Vertex *pv;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pcl::PointXYZ p1, p2;
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "line";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);

            vertices->push_back(p1);
        }
    }

    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0,1.0, 1.0, "vertices");
    return;

    // 显示匹配关系
    // “点-点”红色线段
    // “点-边”白色线段
    count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            if(pv->neighbor_poly_id != -1 && !pv->is_updated){
                if(pv->is_vertex_neighbor){
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    p2 = pv->neighbor_vertex->point;
                    viewer->addLine(p1, p2, 0, 1.0, 1.0, line_name);
                    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
                }else{
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    getProjPoint(p1, pv->neighbor_edge->point, pv->neighbor_edge->next_point->point, p2);
                    if(m_sv.distP2P(pv->neighbor_edge->point, pv->neighbor_edge->next_point->point) <= 0.00001f){
                        p2 = pv->neighbor_edge->point;
                    }
                    /*if(i == 21 && pv->neighbor_poly_id == 13){
                        cout << "i == 21 && pv->neighbor_poly_id == 13" << endl;
                        cout << pv->point << endl;
                        cout << pv->neighbor_edge->point << endl;
                        cout << pv->neighbor_edge->next_point->point << endl;
                        viewer->addLine(pv->neighbor_edge->point, pv->neighbor_edge->next_point->point, 0.5, 1.0, 0.5, "line_name");
                        viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, "line_name");
                    }*/
                    viewer->addLine(p1, p2, 1.0, 1.0, 1.0, line_name);
                    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
                }
            }
        }
    }
}

// 查看特定多边形
void PCLViewer::on_viewSpecificPolyAction_triggered()
{
    RemovePointCloudDialog dlg;
    dlg.exec();
    // 多边形编号: poly_id
    int poly_id = dlg.id.toInt();

    // 删除可能的多边形线段
    int size = 1000;
    string poly_line_name;
    char buf[64];
    for(int i = 0; i < size; ++i){
        poly_line_name = "polyLine";
        itoa(i, buf, 16);
        poly_line_name += buf;
        viewer->removeShape(poly_line_name);
    }

    // 添加线段
    POLYGON *pPoly = &m_sv.polygons[poly_id];
    Vertex *pv;
    for(int i = 0; i < pPoly->getSize(); ++i){
        pv = i == 0 ? pPoly->start_point : pv->next_point;
        poly_line_name = "polyLine";
        itoa(i, buf, 16);
        poly_line_name += buf;
        pcl::PointXYZ p1, p2;
        p1 = pv->point;
        p2 = pv->next_point->point;
        viewer->addLine(p1, p2, 1.0, 0.5, 0.5, poly_line_name);
        viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, poly_line_name);
    }
}

// 轮番查看特定多边形
void PCLViewer::on_action_4_triggered()
{
    memset(CMD, 0, CMD_Size);
    strcpy(CMD, "see poly");
    string str = CMD;
    cout << str << endl;
}

// 处理多边形的退化情况
void PCLViewer::on_dealWithPolyDegenerationAction_triggered()
{
    // 处理多边形退化情况
    m_sv.dealWithPolysDegeneration(m_sv.polygons);

    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 加线段
    char buf[64];
    string line_name;
    pcl::PointXYZ p1, p2;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            vertices->push_back(pv->point);
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "l";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }

    // 加顶点点云
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, "vertices");
}

// 执行多边形删除
void PCLViewer::on_action_10_triggered()
{
    RemovePointCloudDialog dlg;
    dlg.exec();
    // 多边形编号: poly_id
    int poly_id = dlg.id.toInt();

    // 首先清理顶点
    m_sv.polygons[poly_id].clear();

    // 然后从多边形集合中将其删除
    cout << "before del, m_sv.polygons.size = " << m_sv.polygons.size() << endl;
    auto iter = std::find(m_sv.polygons.begin(), m_sv.polygons.end(), m_sv.polygons[poly_id]);
    m_sv.polygons.erase(iter);
    cout << "after del, m_sv.polygons.size = " << m_sv.polygons.size() << endl;

    // 最后更新多边形的顶点信息
    if(poly_id == m_sv.polygons.size()-1) return;
    for(int i = poly_id; i < m_sv.polygons.size()-1; ++i){
        Vertex* pv;
        POLYGON* pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pv->point.data[3] = i;
        }
    }

    ////////////////////////////////////////////////////////////////////
    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 加线段
    char buf[64];
    string line_name;
    pcl::PointXYZ p1, p2;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            vertices->push_back(pv->point);
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "l";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }

    // 加顶点点云
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");
    cout << "多边形删除完毕" << endl;
}

// 一步执行顶点合并
void PCLViewer::on_performSnappingVerticesAction_triggered()
{
    m_sv.setParams(ratio_of_scale, T_ratio_lineSeg_middle_area, T_dist_proj_to_lineSegEnd,
                   T_maxAngle_bet_two_polys);
    this->m_sv.setPolygons(this->initial_polys);
    m_sv.perform(initial_polys);

    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 加线段
    char buf[64];
    string line_name;
    pcl::PointXYZ p1, p2;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            vertices->push_back(pv->point);
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "l";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }

    // 加顶点点云
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");
}

// 继续执行顶点合并
void PCLViewer::on_goOnPerformSnappingVerticesAction_triggered()
{
	// 更新参数
	this->loadParams();
    m_sv.perform(initial_polys);

    // 更新视图
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 加线段
    char buf[64];
    string line_name;
    pcl::PointXYZ p1, p2;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            vertices->push_back(pv->point);
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "l";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }

    // 加顶点点云
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");
}

// 在有角对边关系的线段上增加顶点
void PCLViewer::on_action_6_triggered()
{
	m_sv.insertVertexOnEdgeWithAngleToEdgeRelation();
	this->displayWireframe();
    cout << "在有角对边关系的线段上增加顶点完成" << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 查看多边形id
void PCLViewer::on_observePolyIDAction_triggered()
{
    memset(CMD, 0, CMD_Size);
    strcpy(CMD, "observe poly id");
    cout << "enter observe poly id mode" << endl;

    if(m_sv.polygons.size() == 0){
        m_sv.setPolygons(initial_polys);
        if(m_sv.polygons.size() == 0){
            cerr << "系统内没有多边形数据!" << endl;
            return;
        }
    }

    // 重置顶点点云
    this->m_vertices->clear();
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            m_vertices->push_back(pv->point);
        }
    }
    this->m_kdtree.setInputCloud(this->m_vertices);

    viewer->removeAllPointClouds();
    viewer->addPointCloud(m_vertices, "m_vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "m_vertices");
}

void PCLViewer::displayWireframe()
{
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 加线段
    char buf[64];
    string line_name;
    pcl::PointXYZ p1, p2;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            vertices->push_back(pv->point);
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "l";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);
        }
    }

    // 加顶点点云
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");
}

void PCLViewer::performOptimize(){
    cout << "delta = " << m_osnap.delta << endl;
    cout << "max_k_count = " << m_osnap.max_k_count << endl;
	m_osnap.cosntructJacobMat();
	//cout << "J = " << endl;
	//cout << m_osnap.J << endl;	
	m_osnap.performGN();
	//cout << "beta = " << endl;
	//cout << m_osnap.beta << endl;
	// m_osnap.beta[3] = m_osnap.beta[4] = m_osnap.beta[5] = 0;
	/* 匀速螺旋运动测试代码
	m_osnap.beta.setZero();
	m_osnap.beta[0] = 0.1;
	//m_osnap.beta[6] = 0.3;*/
	m_osnap.updateMovingCoor();
	m_osnap.updateVertices();

	// 更新多边形顶点位置
	Vertex *pv = m_osnap.m_poly->start_point;
	int poly_id = floor(pv->point.data[3]+0.5f);
	for(int i = 0; i < m_osnap.m_osnap_poly.vertices.size(); ++i){
		initial_polys[poly_id].simplified_vertices->points[i].x = m_osnap.m_osnap_poly.vertices[i][0];
		initial_polys[poly_id].simplified_vertices->points[i].y = m_osnap.m_osnap_poly.vertices[i][1];
		initial_polys[poly_id].simplified_vertices->points[i].z = m_osnap.m_osnap_poly.vertices[i][2];
	}

	displayWireframe();

	// 显示更新后的多边形
	viewer->removePointCloud("update poly");
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ p;
	for(int i = 0; i < m_osnap.m_osnap_poly.vertices.size(); ++i){
		p.x = m_osnap.m_osnap_poly.vertices[i][0];
		p.y = m_osnap.m_osnap_poly.vertices[i][1];
		p.z = m_osnap.m_osnap_poly.vertices[i][2];
		cloud->push_back(p);
	} 
	viewer->addPointCloud(cloud, "update poly");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "update poly");

	char buf[64];
	string line_name;
	pcl::PointXYZ p1, p2;
	int i_next;
	for(int i = 0; i < cloud->size(); ++i){
		itoa(i, buf, 16);
		line_name = "update";
		line_name += buf;
		viewer->removeShape(line_name);

		i_next = i == cloud->size()-1 ? 0 : i+1;
		p1 = cloud->points[i];
		p2 = cloud->points[i_next];
		viewer->addLine(p1, p2, 0, 1.0, 0, line_name);
		viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
	}
}

// 执行特定多边形的位置优化
void PCLViewer::on_performOptimisePositionAction_triggered()
{
    RemovePointCloudDialog dlg;
    dlg.exec();
    int poly_id = dlg.id.toInt();
    cout << "selected poly id = " << poly_id << endl;

    // 设置多边形poly_id的邻域关系
    if(m_sv.polygons.size() == 0){
        m_sv.setPolygons(this->initial_polys);
    }else{
        // 更新顶点的搜索半径
        m_sv.setSearchRadiusForEachVertex(this->initial_polys);
    }

    if(poly_id >= m_sv.polygons.size()){ cerr << "多边形编号超出范围!" << endl; return;}
    // 设置多边形poly_id中每个顶点的邻域关系
    Vertex *pv;
    for(int i = 0; i < m_sv.polygons[poly_id].getSize(); ++i){
        pv = i == 0 ? m_sv.polygons[poly_id].start_point : pv->next_point;
        m_sv.getNeighborPolygon(poly_id, pv);
    }

    // 优化
	// 优化步骤初始化
	m_osnap.reset();
	// 设置数据结构
	m_osnap.m_polygons = &m_sv.polygons;
	m_osnap.m_poly = &m_sv.polygons[poly_id];
	m_osnap.setOSnapPoly();
	m_osnap.categoriseVertices();
	m_osnap.initBeta();
    //m_osnap.delta = 0.0000001;

    m_cos_dlg.resetCount();
    m_cos_dlg.setPolyID(poly_id);
    m_cos_dlg.pclviewer = (void *)this;
    m_cos_dlg.show();
}

// 生成测试多边形
// 执行优化
bool is_f = true;
void PCLViewer::on_action_triggered()
{
    m_sv.polygons.clear();
    is_f = true;

    // 首先生成四边形
    POLYGON poly1;
    Vertex *v1 = new Vertex();
    v1->point.x = 2;
    v1->point.y = 1;
    v1->point.z = 0;
    poly1.insertFirstVertex(v1);
    Vertex *v2 = new Vertex();
    v2->point.x = 2;
    v2->point.y = 0;
    v2->point.z = 0;
    poly1.insertMiddleVertex(v2);
    Vertex *v3 = new Vertex();
    v3->point.x = 3;
    v3->point.y = 0;
    v3->point.z = 0;
    poly1.insertMiddleVertex(v3);
    Vertex *v4 = new Vertex();
    v4->point.x = 3;
    v4->point.y = 1;
    v4->point.z = 0;
    poly1.insertLastVertex(v4);

    poly1.coeff << 0, 0, 1, 0;
    m_sv.polygons.push_back(poly1);

    // 然后生成三角形
    POLYGON poly2;
    Vertex *v5 = new Vertex();
    v5->point.x = 2.1;
    v5->point.y = 1.15;
    v5->point.z = -0.13;
    poly2.insertFirstVertex(v5);
    Vertex *v6 = new Vertex();
    v6->point.x = 3.15;
    v6->point.y = 1.15;
    v6->point.z = -0.1;
    poly2.insertMiddleVertex(v6);
    Vertex *v7 = new Vertex();
    v7->point.x = 2.65;
    v7->point.y = 1.15;
    v7->point.z = -1.8;
    poly2.insertLastVertex(v7);

    poly2.coeff << 0, 1, 0, 1.15;
    m_sv.polygons.push_back(poly2);

    // 设置第一个多边形的邻域关系
    Vertex* pv;
    for(int i = 0; i < poly1.getSize(); ++i){
        pv = i == 0 ? poly1.start_point : pv->next_point;
        if(i == 0){
            pv->neighbor_poly_id = 1;
            pv->is_vertex_neighbor = true;
            pv->neighbor_vertex = v5;
            pv->neighbor_edge = NULL;
        }else if(i == 3){
            pv->neighbor_poly_id = 1;
            pv->is_vertex_neighbor = true;
            pv->neighbor_vertex = v6;
            pv->neighbor_edge = NULL;
        }else{
            pv->neighbor_poly_id = -1;
            pv->neighbor_edge = pv->neighbor_vertex = NULL;
        }
    }

    // 显示两个多边形
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
    // 所有的顶点
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 显示线段
    int count = 0;
    string line_name;
    char buf[64];
    //Vertex *pv;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pcl::PointXYZ p1, p2;
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "line";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);

            vertices->push_back(p1);
        }
    }

    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0,0.0, 0.0, "vertices");
    return;
}

// 对测试多边形执行位置优化
void PCLViewer::on_action_2_triggered()
{
	if(is_f){
		// 优化
		// 优化步骤初始化
		m_osnap.reset();
		// 设置数据结构
		m_osnap.m_polygons = &m_sv.polygons;
		m_osnap.m_poly = &m_sv.polygons[0];
		m_osnap.setOSnapPoly();
		m_osnap.categoriseVertices();
		m_osnap.initBeta();

        cout << "delta = " << m_osnap.delta << endl;
		cout << "max_k_count = " << m_osnap.max_k_count << endl;

		cout << "f1_weight = " << m_osnap.f1_weight << endl;
		cout << "f2_weight = " << m_osnap.f2_weight << endl;
		cout << "f3_1_weight = " << m_osnap.f3_1_weight << endl;
		cout << "f3_2_weight = " << m_osnap.f3_2_weight << endl;
		cout << "f8_weight = " << m_osnap.f8_weight << endl;
		cout << "f4_weight = " << m_osnap.f4_weight << endl;
		cout << "f5_weight = " << m_osnap.f5_weight << endl;
		cout << "f6_weight = " << m_osnap.f6_weight << endl;
		cout << "f7_weight = " << m_osnap.f7_weight << endl;

        is_f = false;
	}

	// 构造Jacob矩阵
    //m_osnap.delta = 0.000000001;
	m_osnap.cosntructJacobMat();
	//cout << "J = " << endl;
	//cout << m_osnap.J << endl;
	m_osnap.performGN();
	cout << "beta = " << endl;
	cout << m_osnap.beta << endl;

	m_osnap.updateMovingCoor();
	m_osnap.updateVertices();

	// 显示更新后的多边形
	viewer->removePointCloud("update poly");
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointXYZ p;
	for(int i = 0; i < m_osnap.m_osnap_poly.vertices.size(); ++i){
		p.x = m_osnap.m_osnap_poly.vertices[i][0];
		p.y = m_osnap.m_osnap_poly.vertices[i][1];
		p.z = m_osnap.m_osnap_poly.vertices[i][2];
		cloud->push_back(p);
	} 
	viewer->addPointCloud(cloud, "update poly");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "update poly");

	char buf[64];
	string line_name;
	pcl::PointXYZ p1, p2;
	int i_next;
	for(int i = 0; i < cloud->size(); ++i){
		itoa(i, buf, 16);
		line_name = "update";
		line_name += buf;
		viewer->removeShape(line_name);

		i_next = i == cloud->size()-1 ? 0 : i+1;
		p1 = cloud->points[i];
		p2 = cloud->points[i_next];
		viewer->addLine(p1, p2, 0, 1.0, 0, line_name);
		viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
	}
}

// 多边形三角分割显示平面效果
void PCLViewer::on_action_5_triggered()
{
    EarClip earClip;

    for(int i = 0; i < m_sv.polygons.size(); ++i){
        earClip.setInputPoly(&(m_sv.polygons[i]));
        if(!earClip.triangulatePoly()){
            cerr << "poly_id = " << i << ", triangulation failed!" << endl;
        }
    }

    // 显示三角面片
    char buf[64];
    string name;
    int count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        for(int j = 0; j < m_sv.polygons[i].triangles.size(); ++j){
            name = "tri";
            itoa(count++, buf, 16);
            name += buf;

            pcl::PointCloud<pcl::PointXYZ>::Ptr tri (new pcl::PointCloud<pcl::PointXYZ>);
            tri->reserve(3);
            tri->push_back(m_sv.polygons[i].triangles[j].v[0]->point);
            tri->push_back(m_sv.polygons[i].triangles[j].v[1]->point);
            tri->push_back(m_sv.polygons[i].triangles[j].v[2]->point);

            viewer->addPolygon<pcl::PointXYZ>(tri, name);
            // viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.5, 0.5, 0.5, name);
        }
    }
}

// 顶点合并后，检查模型合法性
void PCLViewer::on_action_11_triggered()
{
    Vertex *prob_v = NULL;
    if(m_sv.is_model_valid(prob_v)){
        cout << "当前模型合法" << endl;
    }else{
        cout << "检测到不合法的情况" << endl;
		cout << "poly id = " << floor(prob_v->point.data[3]+0.5) << endl;

        viewer->removePointCloud("prob_v");
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        cloud->push_back(prob_v->point);
		viewer->addPointCloud(cloud, "prob_v");
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1, 0, 0, "prob_v");
        viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "prob_v");
    }

    ///////////////////////////////////////////
    // 多边形顶点数量不得少于三个
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        if(m_sv.polygons[i].getSize() < 3){
            cerr << "poly id = " << i << ", with vertices number = " << m_sv.polygons[i].getSize() << " invalid(shold be removed)!" << endl;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 空洞修复

void PCLViewer::displayWireframe2()
{
	viewer->removeAllPointClouds();
	viewer->removeAllShapes();

	pcl::PointXYZ p1, p2;
	POLYGON2 *poly;

	char buf[64];
	string name;
	Vertex2 *pv;
	int count = 0;
	for(int i = 0; i < m_hf.m_polygons.size(); ++i){
		poly = m_hf.m_polygons[i];
		int poly_index_in_pv;
		for(int j = 0; j < poly->getSize(); ++j){
			/*pv = j == 0 ? poly->start_point : pv->next_point;
			p1 = m_hf.m_vertices->points[pv->index];
			p2 = m_hf.m_vertices->points[pv->next_point->index];*/
			if(j == 0){
				pv = poly->start_point;
				poly_index_in_pv = poly->getPolyIndex(pv);	// 在更新后的pv中找到多边形的索引
			}else{
				pv = (Vertex2*)pv->relevant_polys[poly_index_in_pv].next_point;
				poly_index_in_pv = poly->getPolyIndex(pv);	// 在更新后的pv中找到多边形的索引
			}

			p1.x = pv->pos[0];
			p1.y = pv->pos[1];
			p1.z = pv->pos[2];
            p2.x = ((Vertex2*)pv->relevant_polys[poly_index_in_pv].next_point)->pos[0];
            p2.y = ((Vertex2*)pv->relevant_polys[poly_index_in_pv].next_point)->pos[1];
            p2.z = ((Vertex2*)pv->relevant_polys[poly_index_in_pv].next_point)->pos[2];

            itoa(i, buf, 16);
            name = buf;
            name.append("_");
            itoa(j, buf, 16);
            name.append(buf);
            name.append("_l");           

			viewer->addLine(p1, p2, 1.0, 1.0, 0, name);
			viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, name);
		}
	}

	viewer->addPointCloud(m_hf.m_vertices_cloud, "vertices");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0, 0, "vertices");

	// 显示三角网格
	//char buf[64];
    //string name;
    count = 0;
	for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        for(int j = 0; j < m_hf.m_polygons[i]->triangles.size(); ++j){
            itoa(i, buf, 16);
            name = buf;
            name.append("_");
            itoa(j, buf, 16);
            name.append(buf);
            name.append("_t");

            pcl::PointCloud<pcl::PointXYZ>::Ptr tri (new pcl::PointCloud<pcl::PointXYZ>);
            tri->reserve(3);
			pcl::PointXYZ p0, p1, p2;
			p0.x = m_hf.m_polygons[i]->triangles[j].v[0]->pos[0];
			p0.y = m_hf.m_polygons[i]->triangles[j].v[0]->pos[1];
			p0.z = m_hf.m_polygons[i]->triangles[j].v[0]->pos[2];
            tri->push_back(p0);
			p1.x = m_hf.m_polygons[i]->triangles[j].v[1]->pos[0];
			p1.y = m_hf.m_polygons[i]->triangles[j].v[1]->pos[1];
			p1.z = m_hf.m_polygons[i]->triangles[j].v[1]->pos[2];
            tri->push_back(p1);
			p2.x = m_hf.m_polygons[i]->triangles[j].v[2]->pos[0];
			p2.y = m_hf.m_polygons[i]->triangles[j].v[2]->pos[1];
			p2.z = m_hf.m_polygons[i]->triangles[j].v[2]->pos[2];
            tri->push_back(p2);

            viewer->addPolygon<pcl::PointXYZ>(tri, name);
        }
    }
}

// 构造全局点集
void PCLViewer::on_action_3_triggered()
{
    // 进行数据结构的转换
	m_hf.setInputPolygons(&m_sv.polygons);

	/*/ 三角化多边形
	EarClip2 ec;
	for(int i = 0; i < m_hf.m_polygons.size(); ++i){
		ec.setInputPoly(m_hf.m_polygons[i]);
		if(!ec.triangulatePoly()){
			cerr << "多边形" << i << "未能成功三角化" << endl;
		}
	}

    // 设置三角面片重心点云
    m_hf.set_cloud_for_indentify_poly_id();*/

	// 更新视图
	displayWireframe2();
}

// 进入交互模式
void PCLViewer::on_action_7_triggered()
{
    // 将重心点云加入viewer中
    viewer->removePointCloud("triCentroids");
    viewer->addPointCloud(m_hf.m_cloud_for_indentify_poly_id, "triCentroids");
    //cout << "m_hf.m_cloud_for_indentify_poly_id size = " << m_hf.m_cloud_for_indentify_poly_id->size() << endl;

    m_hfio_dlg.m_pclviewer = this;
    m_hfio_dlg.m_CMD = CMD;
    m_hfio_dlg.show();
}

// 保存模型数据
void PCLViewer::on_action_8_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("Save Model Files"),
                                                    "",
                                                    tr("Model Files (*.model)"));

    if (fileName.isNull()) return;

    QFile file(fileName);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    /////////////////////////////////////////////////////////////////////////
    for(int i = 0; i < this->m_hf.m_vertices.size(); ++i){
        m_hf.m_vertices[i]->index_in_m_vertices = i;
    }
    /////////////////////////////////////////////////////////////////////////
    QTextStream out(&file);
    out.setRealNumberPrecision(16);
    out << "[points info]\n";
    out << m_hf.m_vertices.size() << endl;
    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        out << m_hf.m_vertices[i]->pos[0] << "\t"
            << m_hf.m_vertices[i]->pos[1] << "\t"
            << m_hf.m_vertices[i]->pos[2] << "\n";
    }
    out << "[polygons info]" << endl;
    out << m_hf.m_polygons.size() << endl;
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        out << "[ply" << i << "]" << endl;
        out << m_hf.m_polygons[i]->coeff[0] << "\t"
            << m_hf.m_polygons[i]->coeff[1] << "\t"
            << m_hf.m_polygons[i]->coeff[2] << "\t"
            << m_hf.m_polygons[i]->coeff[3] << "\n";
        out << m_hf.m_polygons[i]->getSize() << endl;
        Vertex2* pv2;
        int poly_index;
        for(int j = 0; j < m_hf.m_polygons[i]->getSize(); ++j){
            if(j == 0){
                pv2 = m_hf.m_polygons[i]->start_point;
            }else{
                pv2 = (Vertex2*)pv2->relevant_polys[poly_index].next_point;
            }
            poly_index = m_hf.m_polygons[i]->getPolyIndex(pv2);
            out << pv2->index_in_m_vertices << endl;
        }
    }
}

// 打开岩体模型数据
void PCLViewer::on_action_9_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open Model Files"),
                                                    "",
                                                    tr("Model Files (*.model)"));

    if (fileName.isNull()) return;

    QFile file(fileName);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    //////////////////////////////////////////////////////////////////////////////////////////
    // 清理当前的数据结构
    m_hf.m_cloud_for_indentify_poly_id->clear();
    m_hf.m_vertices_cloud->clear();
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        m_hf.m_polygons[i]->clear();
        delete m_hf.m_polygons[i];
    }
    m_hf.m_polygons.clear();

    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        m_hf.m_vertices[i]->relevant_polys.clear();
        delete m_hf.m_vertices[i];
    }
    m_hf.m_vertices.clear();
    /////////////////////////////////////////////////////////////////////////////////////////////
    QTextStream in(&file);
    in.setRealNumberPrecision(16);

    QString line = in.readLine();
    // 顶点数量
    line = in.readLine();
    int points_number = line.toInt();
    // 读取顶点数据
    for(int i = 0; i < points_number; ++i){
		line = in.readLine();
        QStringList xyz = line.split("\t");
        Vertex2* v = new Vertex2();
        v->pos[0] = xyz[0].toDouble();
        v->pos[1] = xyz[1].toDouble();
        v->pos[2] = xyz[2].toDouble();
        m_hf.m_vertices.push_back(v);
    }

    // 读取多边形数据
    line = in.readLine();
    line = in.readLine();
    int ply_number = line.toInt();
    //EarClip2 ec2;
    for(int i = 0; i < ply_number; ++i){
        POLYGON2* ply = new POLYGON2();
        line = in.readLine();   // [plyi]
        line = in.readLine();   // coeff
        QStringList abcd = line.split("\t");
        ply->coeff << abcd[0].toDouble(),
                      abcd[1].toDouble(),
                      abcd[2].toDouble(),
                      abcd[3].toDouble();
        line = in.readLine();
        int ply_size = line.toInt();
        for(int j = 0; j < ply_size; ++j){
            line = in.readLine();
            int index = line.toInt();
            if(j == 0){
                ply->insertFirstVertex2(m_hf.m_vertices[index]);
            }else if(j == ply_size -1){
                ply->insertLastVertex2(m_hf.m_vertices[index]);
            }else{
                ply->insertMiddleVertex2(m_hf.m_vertices[index]);
            }
        }
//        ec2.setInputPoly(ply);
//        ec2.triangulatePoly();

        m_hf.m_polygons.push_back(ply);
    }
    ///////////////////////////////////////////////////////////////////////////
    // 显示
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
    // 设置三角面片重心点云
    //m_hf.set_cloud_for_indentify_poly_id();

	// 添加顶点点云
	m_hf.m_vertices_cloud->reserve(m_hf.m_vertices.size());
	pcl::PointXYZ p;
	for(int i = 0; i < m_hf.m_vertices.size(); ++i){
		p.x = m_hf.m_vertices[i]->pos[0];
		p.y = m_hf.m_vertices[i]->pos[1];
		p.z = m_hf.m_vertices[i]->pos[2];
		m_hf.m_vertices_cloud->push_back(p);
	}
	m_hf.m_vertices_kdtree.setInputCloud(m_hf.m_vertices_cloud);

    // 更新视图
    displayWireframe2();

	cout << "loaded " << m_hf.m_vertices.size() << " vertices." << endl;
	cout << "loaded " << m_hf.m_polygons.size() << " polygons." << endl;
}

// 保存成dda格式的数据
void PCLViewer::on_action_dda_triggered()
{
    // 统计三角面片的数量
    int tri_number = 0;
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        POLYGON2* ply = m_hf.m_polygons[i];
        tri_number += ply->triangles.size();
    }

    // 统计固定点的数量
    // 固定点是最底部的点，即具有最小y值的点
    double m_min_y = DBL_MAX;

    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        if(m_hf.m_vertices[i]->pos[1] < m_min_y){
            m_min_y = m_hf.m_vertices[i]->pos[1];
        }
    }

    vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> immovable_points;
    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        if(m_hf.m_vertices[i]->pos[1] == m_min_y){
            immovable_points.push_back(m_hf.m_vertices[i]->pos);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("Save Model Files in 3D-DDA style"),
                                                    "",
                                                    tr("Model Files (*.dda)"));

    if (fileName.isNull()) return;

    QFile file(fileName);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);
    int precision;
    cout << "precision = ";
    cin >> precision;
    out.setRealNumberPrecision(precision);

    // 不连续面的数量
    out << tri_number << endl;

    // 固定点 荷载点 测量点 孔洞点 的数量
    out << immovable_points.size() << "\t"
        << 0 << "\t" << 0 << "\t" << 0 << endl;

    // 锚杆数量
    out << 0 << endl;

    //节理面参数
    // 顶点材料    面材料    使用次数
    for(int i = 0; i < tri_number; ++i){
        out << 1 << "\t" << 1 << "\t" << 1 << endl;
    }

    // 空间位置+平移
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        POLYGON2* ply = m_hf.m_polygons[i];
        for(int j = 0; j < ply->triangles.size(); ++j){
            for(int k = 0; k < 3; ++k){
                out << ply->triangles[j].v[k]->pos[0] << " "
                    << ply->triangles[j].v[k]->pos[1] << " "
                    << ply->triangles[j].v[k]->pos[2] << "\t";
            }
            out << 0 << " " << 0 << " " << 0 << endl;
        }
    }

    // 固定点位置
    for(int i = 0; i < immovable_points.size(); ++i){
        out << immovable_points[i][0] << " "
            << immovable_points[i][1] << " "
            << immovable_points[i][2] << endl;
    }

    file.close();
    cout << "文件已保存" << endl;
}

/*****************************************不连续面相关：*******************************************/
// 角度转弧度
double PCLViewer::degree2rad(double degree){
    return degree/180.0*3.141592653589793;
}

// 将产状转换成法向量
void PCLViewer::convertOrientaion2Normal(Orientation &ori, Eigen::Vector3d &normal){
    double alpha = degree2rad(ori.trend);
    double beta = degree2rad(ori.plunge);

    normal << cos(alpha)*cos(beta),
              sin(alpha)*cos(beta),
              sin(beta);
}

// 计算两个向量的夹角
double PCLViewer::getAngleBtwNormals(Eigen::Vector3d &n1, Eigen::Vector3d &n2){
    n1.normalize();
    n2.normalize();
    double cos_val = n1.dot(n2);
    return acos(cos_val)/3.141592653589793*180.0;
}

// 加载优势结构面数据
void PCLViewer::on_action_12_triggered()
{
    QString QfileName = QFileDialog::getOpenFileName(this, tr("Open Files"),
                                                 "/", tr("Files (*.*)"));
    if(QfileName.isNull()) return;

    vector<Orientation> orientations;

    QFile file(QfileName);

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;

    QTextStream in(&file);

    QString line = in.readLine();

    while(!line.isNull()){
        QStringList str_list = line.split("\t");
        Orientation ori;
        ori.trend = str_list[0].toDouble();
        ori.plunge = str_list[1].toDouble();
        orientations.push_back(ori);

        line = in.readLine();
    }

    cout << "loaded " << orientations.size() << " dominant orientations." << endl;
    cout << "Trend\tPlunge" << endl;
    for(int i = 0; i < orientations.size(); ++i){
        cout << orientations[i].trend << "\t" << orientations[i].plunge << endl;
    }

    // 将优势产状转换为优势法向量
    cout << "将优势产状转换为优势法向量" << endl;

    m_dominant_normals.clear();
    for(int i = 0; i < orientations.size(); ++i){
        Eigen::Vector3d normal;
        convertOrientaion2Normal(orientations[i], normal);
        m_dominant_normals.push_back(normal);
    }

    for(int i = 0; i < m_dominant_normals.size(); ++i){
        cout << "normal " << i << ":\n";
        cout << m_dominant_normals[i] << endl;
    }

    // 打印法向量彼此之间的夹角
    cout << "打印法向量彼此之间的夹角:" << endl;
    vector<vector<double>> angles;
    angles.resize(m_dominant_normals.size());
    for(int i = 0; i < angles.size(); ++i){
        angles[i].resize(m_dominant_normals.size());
    }

    for(int i = 0; i < angles.size(); ++i){
        for(int j = 0; j < angles.size(); ++j){
            angles[i][j] = i != j ? getAngleBtwNormals(m_dominant_normals[i], m_dominant_normals[j]) : 0;
        }
    }

    for(int i = 0; i < angles.size(); ++i){
        cout << "\t" << i;
    }
    cout << endl;
    for(int i = 0; i < angles.size(); ++i){
        cout << i << "\t";
        for(int j = 0; j < angles.size(); ++j){
            cout << angles[i][j] << "\t";
        }
        cout << endl;
    }

    ///////////////////////////////////////////////////////
    // 将法向量的相反方向也添加进去
    int original_size = m_dominant_normals.size();
    for(int i = 0; i < original_size; ++i){
        Eigen::Vector3d n = -1 * m_dominant_normals[i];
        m_dominant_normals.push_back(n);
    }

    m_inserted_discontinuities.clear(); // 清理之前插入的不连续面数据
}

// 加载平面点云
void PCLViewer::on_action_14_triggered()
{
    QString QfileName = QFileDialog::getOpenFileName(this, tr("Open PCD Files"),
                                                 "/", tr("PCD Files (*.pcd)"));
    if(QfileName.isNull()) return;

    string filename = QfileName.toStdString();

    m_plane_point_clouds.reset(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::io::loadPCDFile(filename, *m_plane_point_clouds);

    cout << "loaded " << m_plane_point_clouds->size() << " points." << endl;
}

// 识别不连续面
void PCLViewer::on_action_13_triggered()
{
    if(this->m_hf.m_polygons.size() == 0){
        cout << "还未加载模型数据！" << endl;
    }

    cout << "original poly number = ";
    cin >> m_original_poly_number;

    if(m_original_poly_number > m_hf.m_polygons.size()){
        cout << "m_original_poly_number = " << m_original_poly_number << " > m_hf.m_polygons.size() = " << m_hf.m_polygons.size() << endl;
        return;
    }

    cout << "与优势结构面的最大夹角:";
    double max_angle;
    cin >> max_angle;

    double scale;
    cout << "不连续面的空间尺度：";
    cin >> scale;

    vector<double> angles(m_dominant_normals.size());   // 保存平面法向量与结构面法向量的夹角
    Eigen::Vector3d normal;
    int d_count = 0;
    for(int i = 0; i < m_original_poly_number; ++i){
        POLYGON2* ply = m_hf.m_polygons[i];
        normal << ply->coeff[0], ply->coeff[1], ply->coeff[2];
        double min_angle = DBL_MAX;
        int index = -1;
        for(int j = 0; j < m_dominant_normals.size(); ++j){
            angles[j] = getAngleBtwNormals(m_dominant_normals[j], normal);
            if(min_angle > angles[j]){
                min_angle = angles[j];
                index = j;
            }
        }
        if(min_angle < max_angle){  // 最小角度在阈值以内，则认为该多边形为不连续面
            ply->is_discontinuity = true;
            // 用一个大的正方形表示该多边形
            genBigPlane(scale, ply);
            d_count += 2;
        }
    }
    cout << "不连续面自动识别完毕" << endl;
    cout << "当前的不连续面数量（以三角形计）：" << d_count << endl;
}

/*
 * 根据尺度和多边形id，构造大平面
 * input:
 *      scale: 尺度
 *      ply: 多边形
 * output:
 *      vertices: size为6的一个向量，存有6个顶点，每三个构成一个三角形（将一个大的四方形平面分割成两个三角形）
*/
void PCLViewer::genBigPlane(double scale, POLYGON2* ply)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    Vertex2* pv2;
    pcl::PointXYZ pcl_p;
    for(int i = 0; i < ply->getSize(); ++i){
        pv2 = i == 0 ? ply->start_point : ply->getNextVertex(pv2);
        pcl_p.x = pv2->pos[0];
        pcl_p.y = pv2->pos[1];
        pcl_p.z = pv2->pos[2];
        cloud->push_back(pcl_p);
    }

    // 首先构造局部坐标系（原点是多边形的重心，方向由PCA确定）
    EIGEN_ALIGN16 Eigen::Matrix3f covariance_matrix;
    Eigen::Vector4f xyz_centroid;
    pcl::computeMeanAndCovarianceMatrix(*cloud, covariance_matrix, xyz_centroid);

    // 协方差矩阵的特征值与特征向量
    EIGEN_ALIGN16 Eigen::Vector3f eigen_values;
    EIGEN_ALIGN16 Eigen::Matrix3f eigen_vectors;
    pcl::eigen33(covariance_matrix, eigen_vectors, eigen_values);

    // 重心/原点
    Eigen::Vector3d centroid;
    centroid << xyz_centroid[0],
                xyz_centroid[1],
                xyz_centroid[2];

    // 方向向量
    Eigen::Vector3d e1, e2;
    e1 << eigen_vectors(0,1),
          eigen_vectors(1,1),
          eigen_vectors(2,1);
    e2 << eigen_vectors(0,2),
          eigen_vectors(1,2),
          eigen_vectors(2,2);
    e1.normalize();
    e2.normalize();

    // 正方形边长：2*scale
    // 四个顶点
    // 假设e1水平，e2垂直
    Eigen::Vector3d v0, v1, v2, v3;
    v0 = centroid + -1*scale*e1 + scale*e2;     // 左上角
    v1 = centroid + scale*e1 + scale*e2;        // 右上角
    v2 = centroid + scale*e1 + -1*scale*e2;     // 右下角
    v3 = centroid + -1*scale*e1 + -1*scale*e2;  // 左下角

    /*
     * v0--------------v1
     * | \             |
     * |   \           |
     * |     \         |
     * |       \       |
     * |         \     |
     * |           \   |
     * |             \ |
     * v3--------------v2
    */
    ply->discontinuity_vertices.resize(6);
    ply->discontinuity_vertices[0] = v0;
    ply->discontinuity_vertices[1] = v1;
    ply->discontinuity_vertices[2] = v2;
    ply->discontinuity_vertices[3] = v0;
    ply->discontinuity_vertices[4] = v2;
    ply->discontinuity_vertices[5] = v3;
}

// 准备交互式插入不连续面
void PCLViewer::on_action_16_triggered()
{
    // 清屏
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    // 显示点云数据
    viewer->addPointCloud(m_plane_point_clouds, "cloud");

    strcpy(CMD, "insert dis");
}

// 将当前选中的顶点作为顶点1
void PCLViewer::on_action_1_triggered()
{
    m_p1 = m_selected_p;
    cout << "顶点1：" << endl;
    cout << m_p1 << endl;
}

// 将当前选中的顶点作为顶点2
void PCLViewer::on_action_17_triggered()
{
    m_p2 = m_selected_p;
    cout << "顶点2：" << endl;
    cout << m_p2 << endl;
}

// 生成不连续面
void PCLViewer::on_action_18_triggered()
{
    Eigen::Vector3d v0, v1, v2, v3, v, n;
    v1 << m_p1.x, m_p1.y, m_p1.z;
    v2 << m_p2.x, m_p2.y, m_p2.z;
    v = v1-v2;

    /*
     *           ↑y(垂直)
     *           |
     *           |
     *           |------>x(水平)
     *          /
     *         /
     *        z(前后)
    */
    v3 << 0, 0, -1; // 指向后方

    n = v.cross(v3);
    n.normalize();

    // 用与n夹角最小的优势法向量去代替它
    double min_angle = DBL_MAX;
    int index = -1;
    for(int i = 0; i < m_dominant_normals.size(); ++i){
        double angle = getAngleBtwNormals(n, m_dominant_normals[i]);
        if(min_angle > angle){
            min_angle = angle;
            index = i;
        }
    }

    n = m_dominant_normals[index];

    ////////////////////////////////////////////////////////////////////////
    double scale;
    cout << "input scale = ";
    cin >> scale;

    // 正方形中心
    Eigen::Vector3d centroid = (v1+v2)*0.5;
    // 方向向量
    Eigen::Vector3d e1, e2;
    e1 = v2 - centroid;
    e1.normalize();
    e2 = e1.cross(n);
    e2.normalize();

    v0 = centroid + -1*scale*e1 + scale*e2;     // 左上角
    v1 = centroid + scale*e1 + scale*e2;        // 右上角
    v2 = centroid + scale*e1 + -1*scale*e2;     // 右下角
    v3 = centroid + -1*scale*e1 + -1*scale*e2;  // 左下角

    /*
     * v0--------------v1
     * | \             |
     * |   \           |
     * |     \         |
     * |       \       |
     * |         \     |
     * |           \   |
     * |             \ |
     * v3--------------v2
    */
    DiscontinuityTriangle tr1, tr2;
    tr1.v[0] = v0;
    tr1.v[1] = v1;
    tr1.v[2] = v2;

    tr2.v[0] = v0;
    tr2.v[1] = v2;
    tr2.v[2] = v3;

    m_inserted_discontinuities.push_back(tr1);
    m_inserted_discontinuities.push_back(tr2);
    cout << "不连续面插入完成" << endl;
}

// 显示不连续面
void PCLViewer::on_action_19_triggered()
{
    // 清屏
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();

    char buf[12];
    string name;
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        _itoa(i, buf, 16);

        POLYGON2* ply = m_hf.m_polygons[i];
        if(!ply->is_discontinuity){
            for(int j = 0; j < ply->triangles.size(); ++j){
                pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
                pcl::PointXYZ p;
                char buf2[12];
                _itoa(j, buf2, 16);
                name = buf;
                name.append("_");
                name.append(buf2);
                name.append("_t");// i_j_t
                for(int k = 0; k < 3; ++k){
                    p.x = ply->triangles[j].v[k]->pos[0];
                    p.y = ply->triangles[j].v[k]->pos[1];
                    p.z = ply->triangles[j].v[k]->pos[2];
                    cloud->push_back(p);
                }
                viewer->addPolygon<pcl::PointXYZ>(cloud, name);
            }
        }else{
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointXYZ p;
            for(int j = 0; j < 6; ++j){
                if(j == 3 || j == 4) continue;

                p.x = ply->discontinuity_vertices[j][0];
                p.y = ply->discontinuity_vertices[j][1];
                p.z = ply->discontinuity_vertices[j][2];
                cloud->push_back(p);
            }
            name = buf;
            name.append("_0_t");
            viewer->addPolygon<pcl::PointXYZ>(cloud, name);
        }
    }

    // 显示后插入的不连续面
    _itoa(m_hf.m_polygons.size(), buf, 16);
    for(int i = 0; i < m_inserted_discontinuities.size(); ++i){
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointXYZ p;
        for(int j = 0; j < 3; ++j){
            p.x = m_inserted_discontinuities[i].v[j][0];
            p.y = m_inserted_discontinuities[i].v[j][1];
            p.z = m_inserted_discontinuities[i].v[j][2];
            cloud->push_back(p);
        }
        char buf2[12];
        _itoa(i, buf2, 16);
        name = buf;
        name.append("_");
        name.append(buf2);
        name.append("_t");// i_j_t
        viewer->addPolygon<pcl::PointXYZ>(cloud, name);
    }
}

// 保存岩体模型DDA数据
void PCLViewer::on_action_DDA_triggered()
{
    // 首先统计不连续面的数量
    int tri_number = 0;
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        if(m_hf.m_polygons[i]->is_discontinuity){
            tri_number += 2;
        }else{
            tri_number += m_hf.m_polygons[i]->triangles.size();
        }
    }
    tri_number += m_inserted_discontinuities.size();    // 插入的不连续面

    cout << "不连续面的数量：" << tri_number << endl;

    // 统计固定点的数量
    // 固定点是最底部的点，即具有最小y值的点
    double m_min_y = DBL_MAX;

    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        if(m_hf.m_vertices[i]->pos[1] < m_min_y){
            m_min_y = m_hf.m_vertices[i]->pos[1];
        }
    }

    vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> immovable_points;
    for(int i = 0; i < m_hf.m_vertices.size(); ++i){
        if(m_hf.m_vertices[i]->pos[1] == m_min_y){
            immovable_points.push_back(m_hf.m_vertices[i]->pos);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    tr("Save Model Files in 3D-DDA style"),
                                                    "",
                                                    tr("Model Files (*.dda)"));

    if (fileName.isNull()) return;

    QFile file(fileName);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
        return;

    QTextStream out(&file);
    int precision;
    cout << "precision = ";
    cin >> precision;
    out.setRealNumberPrecision(precision);

    // 不连续面的数量
    out << tri_number << endl;

    // 固定点 荷载点 测量点 孔洞点 的数量
    out << immovable_points.size() << "\t"
        << 0 << "\t" << 0 << "\t" << 0 << endl;

    // 锚杆数量
    out << 0 << endl;

    //节理面参数
    // 顶点材料    面材料    使用次数
    for(int i = 0; i < tri_number; ++i){
        out << 1 << "\t" << 1 << "\t" << 1 << endl;
    }

    // 空间位置+平移
    // 原有多边形
    for(int i = 0; i < m_hf.m_polygons.size(); ++i){
        POLYGON2* ply = m_hf.m_polygons[i];
        if(ply->is_discontinuity){
            for(int j = 0; j < 6; ++j){
                out << ply->discontinuity_vertices[j][0] << " "
                    << ply->discontinuity_vertices[j][1] << " "
                    << ply->discontinuity_vertices[j][2] << "\t";
                if(j == 2 || j == 5){
                    out << 0 << " " << 0 << " " << 0 << endl;
                }
            }
        }else{
            for(int j = 0; j < ply->triangles.size(); ++j){
                for(int k = 0; k < 3; ++k){
                    out << ply->triangles[j].v[k]->pos[0] << " "
                        << ply->triangles[j].v[k]->pos[1] << " "
                        << ply->triangles[j].v[k]->pos[2] << "\t";
                }
                out << 0 << " " << 0 << " " << 0 << endl;
            }
        }
    }
    // 插入的不连续面
    for(int i = 0; i < m_inserted_discontinuities.size(); ++i){
        for(int k = 0; k < 3; ++k){
            out << m_inserted_discontinuities[i].v[k][0] << " "
                << m_inserted_discontinuities[i].v[k][1] << " "
                << m_inserted_discontinuities[i].v[k][2] << "\t";
        }
        out << 0 << " " << 0 << " " << 0 << endl;
    }

    // 固定点位置
    for(int i = 0; i < immovable_points.size(); ++i){
        out << immovable_points[i][0] << " "
            << immovable_points[i][1] << " "
            << immovable_points[i][2] << endl;
    }

    file.close();
    cout << "文件已保存" << endl;
}

// 顶点合并-展示精华匹配关系
void PCLViewer::on_action_20_triggered()
{
    if(is_first_run){
        this->m_sv.setPolygons(this->initial_polys);
        is_first_run = false;
    }else{
        // 记得调整每个顶点的搜索半径
        POLYGON *pPoly;
        Vertex *pv;
        for(int i = 0; i < m_sv.polygons.size(); ++i){
            pPoly = &m_sv.polygons[i];
            for(int j = 0; j < pPoly->getSize(); ++j){
                pv = j == 0 ? pPoly->start_point : pv->next_point;
                pv->radius = this->ratio_of_scale * initial_polys[i].scale;
            }
        }
    }
    // 为多边形的每个顶点设置邻域关系
    this->m_sv.setNeighborsForEachPoly();
    // 将单向的匹配关系转换为双向的匹配关系
    this->m_sv.setMatchingRelations();
    // 精化匹配关系(不添加点-点匹配，但可能会减少点-边匹配)
    this->m_sv.refineMachingRelations();

    // 显示
    viewer->removeAllPointClouds();
    viewer->removeAllShapes();
    // 所有的顶点
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);

    // 显示线段
    int count = 0;
    string line_name;
    char buf[64];
    Vertex *pv;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pcl::PointXYZ p1, p2;
            p1 = pv->point;
            p2 = pv->next_point->point;
            itoa(count++, buf, 16);
            line_name = "line";
            line_name += buf;
            viewer->addLine(p1, p2, 1.0, 1.0, 0.0, line_name);
            viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 2, line_name);

            vertices->push_back(p1);
        }
    }
    viewer->addPointCloud(vertices, "vertices");
    viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0, 1.0, 1.0, "vertices");
    // 显示匹配关系
    // “点-点”红色线段
    // “点-边”白色线段
    count = 0;
    for(int i = 0; i < m_sv.polygons.size(); ++i){
        POLYGON *pPoly = &m_sv.polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            if(pv->neighbor_poly_id != -1 && !pv->is_updated){
                if(pv->is_vertex_neighbor){
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    p2 = pv->neighbor_vertex->point;
                    viewer->addLine(p1, p2, 0, 1.0, 1.0, line_name);
                    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
                }else{
                    itoa(count++, buf, 16);
                    line_name = "ver";
                    line_name += buf;
                    pcl::PointXYZ p1, p2;
                    p1 = pv->point;
                    getProjPoint(p1, pv->neighbor_edge->point, pv->neighbor_edge->next_point->point, p2);
                    viewer->addLine(p1, p2, 1.0, 1.0, 1.0, line_name);
                    viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3, line_name);
                }
            }
        }
    }
}

// 手动指定不连续面
vector<int> plane_id;
void PCLViewer::on_action_21_triggered()
{
    cout << "is clear plane_id?" << endl;
    char is_clear;
    cin >> is_clear;
    if(is_clear == 'y'){
        plane_id.clear();
    }

    cout << "input plane id:";
    int int_plane_id = 0;

    do{
        cin >> int_plane_id;
        if(int_plane_id < 0){
            return;
        }else{
            if(std::find(plane_id.begin(), plane_id.end(), int_plane_id) == plane_id.end()){
                plane_id.push_back(int_plane_id);

                cout << "plane_id.size = " << plane_id.size() << ": ";
                for(int i = 0; i < plane_id.size(); ++i) cout << plane_id[i] << "   ";
                cout << endl;
            }
        }
	} while (true);
}
