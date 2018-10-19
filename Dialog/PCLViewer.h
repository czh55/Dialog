#ifndef PCLVIEWER_H
#define PCLVIEWER_H

// STL
#include <iostream>
#include <fstream>

// Qt
#include <QMainWindow>

// Point Cloud Library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/centroid.h>
/****************************************前半段：开始**************************************************/
#include <pcl/filters/voxel_grid.h>
/****************************************前半段：结束**************************************************/

// Visualization Toolkit (VTK)
#include <vtkRenderWindow.h>

// local
#include "SimplifyVerticesSize.h"
#include "osnap.h"
#include "earclip2.h"

#include "holefillinginteractiveoperationdialog.h"
#include "controloptimizationsnappingdialog.h"
/****************************************前半段：开始**************************************************/
#include <iostream>
#include "RegulateNormalDialog.h"
/****************************************前半段：结束**************************************************/

using namespace std;

struct Orientation{
    double trend;
    double plunge;
};

struct RGB{
    double r;
    double g;
    double b;
};

struct DiscontinuityTriangle{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d v[3];
};

//c++全局变量 不重复包含
//https://blog.csdn.net/a_ran/article/details/39229941?utm_source=blogxgwz1
extern char CMD[128];

namespace Ui
{
  class PCLViewer;
}

class PCLViewer : public QMainWindow
{
  Q_OBJECT

public:
  explicit PCLViewer (QWidget *parent = 0);
  ~PCLViewer ();

  // 最初的多边形集合
  vector<InitialPoly> initial_polys;

  // 多边形的尺度
  vector<float> ply_scales;

  // 多边形顶点点云
  pcl::PointCloud<pcl::PointXYZ>::Ptr m_vertices;

  // 与多边形顶点点云匹配的kdtree
  pcl::KdTreeFLANN<pcl::PointXYZ> m_kdtree;

  /////////////////////////////////////////////////////////////////////////////////////////
  // 不连续面相关：
  pcl::PointCloud<pcl::PointXYZ>::Ptr m_plane_point_clouds;
  vector<pcl::PointCloud<pcl::PointXYZ>::Ptr, Eigen::aligned_allocator<pcl::PointCloud<pcl::PointXYZ>::Ptr>> m_planes;
  // 优势法向量
  vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> m_dominant_normals;

  // 角度转弧度
  double degree2rad(double degree);
  // 将产状转换成法向量
  void convertOrientaion2Normal(Orientation &ori, Eigen::Vector3d &normal);
  // 计算两个向量的夹角
  double getAngleBtwNormals(Eigen::Vector3d &n1, Eigen::Vector3d &n2);
  // 最初的多边形数量
  int m_original_poly_number;
  /*
   * 根据尺度和多边形id，构造大平面
   * input:
   *      scale: 尺度
   *      ply: 多边形
   * output:
   *      vertices: size为6的一个向量，存有6个顶点，每三个构成一个三角形（将一个大的四方形平面分割成两个三角形）
  */
  void genBigPlane(double scale, POLYGON2* ply);
  // 插入的不连续面
  vector<DiscontinuityTriangle> m_inserted_discontinuities;
  // 当前选中的顶点
  pcl::PointXYZ m_selected_p;
  pcl::PointXYZ m_p1;
  pcl::PointXYZ m_p2;
/****************************************前半段：开始**************************************************/
	RegulateNormalDialog RegulateNormalDialog;
/****************************************前半段：结束**************************************************/
private slots:
	void on_openPointCloudAction_triggered();

	void on_openTxtAction_triggered();

	void on_openFile1RegistrationAction_triggered();

	void on_openFile2RegistrationAction_triggered();

	void on_voxelGridFiltMergeCloudAction_triggered();

	void on_rotatePointCloudAction_triggered();

	void on_registrationSACAction_triggered();

	void on_registrationICPAction_triggered();

	void on_registrationPlaneAction_triggered();

	void on_savePointCloudAction_triggered();

	void on_repairHolesOFFAction_triggered();

	void on_TriangularMeshingAction_triggered();

	void on_bgColorMenu_triggered();

	void on_cleanPointCloudAction_triggered();

	void on_removeNan1Action_triggered();

	void on_voxelGridFilt1Action_triggered();

	void on_statisticalOutlierRemovalFiltAction_triggered();

	void on_setNumberforpointAction_triggered();

	void on_rebuildPlaneAction_triggered();

	void on_pointCloudColorMenu_triggered();

	void on_translateToCentroidAction_triggered();

	void on_removeNanAction_triggered();

	void on_plane_detect_set_param_Action_triggered();

	void on_removeRedundantPointsAction_triggered();

	void on_regulateCoorAction_triggered();

	void on_load_param_action_triggered();

	void on_normalEstimateAction_triggered();

	void on_regulateNormalAction_triggered();

	void on_performRegulateAction_triggered();

	void on_createPSAction_triggered();

	void on_savePointNormalFileAction_triggered();

	void on_openPointCloudNormalFileAction_triggered();

	void on_filtPSAction_triggered();

	void on_segPSAction_triggered();

	void on_segPlaneAction_triggered();

	void on_regionGrowingAction_triggered();

	void on_mergePlanesAction_triggered();

	void on_polyPlanesAction_triggered();

	void on_postProcessAction_triggered();

	void on_runAgainAction_triggered();

	void on_autoPerformAction_triggered();

	void on_editPolyAction_triggered();

	void on_delPolyAction_triggered();

	void on_performDelAction_triggered();

	void on_savePolyDataAction_triggered();

	void on_enterPruneModeAction_triggered();

	void on_selCurPolyAction_triggered();

	void on_setFirstPointAction_triggered();

	void on_setSecondPointAction_triggered();

	void on_performPolyCutAction_triggered();

	void on_displayLineSegAction_triggered();

	void on_switchLineSegAction_triggered();

	void on_performLineSegDelAction_triggered();

	void on_loadPolyDataAction_triggered();

	void on_voxelGridFiltAction_triggered();


public:
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;

private slots:
  void on_OpenAction_triggered();

  void on_setBgColorAction_triggered();

  void on_setCoordinateAxesAction_triggered();

  void on_removePointCloudAction_triggered();

  void on_verticesNormalEstimationAction_triggered();

  void on_clearScreenAction_triggered();

  void on_ObservePolyVerNormalAction_triggered();

  void on_simplifyVerticesSizeAction_triggered();

  void on_RGtoDetectLineSegsAction_triggered();

  void on_RGlineSegDisplayAction_triggered();

  void on_doSimplifyVerticesSizeAction_triggered();

  void on_removeAllPointCloudsAction_triggered();

  void on_removeAllShapesAction_triggered();

  void on_removeAllInitialVerticesAction_triggered();

  void on_showNeighborRegionAction_triggered();

  void on_mergeVerticesParamsAction_triggered();

  void on_showMatchingAction_triggered();

  void on_mergeMatchingRelationsAction_triggered();

  void on_showMatchingClustersAction_triggered();

  void on_mergeVertexEdgeMatchingAction_triggered();

  void on_viewSpecificPolyAction_triggered();

  void on_action_4_triggered();

  void on_dealWithPolyDegenerationAction_triggered();

  void on_performSnappingVerticesAction_triggered();

  void on_goOnPerformSnappingVerticesAction_triggered();

  void on_observePolyIDAction_triggered();

  void on_performOptimisePositionAction_triggered();

  void on_action_triggered();

  void on_action_2_triggered();

  void on_action_5_triggered();

  void on_action_6_triggered();

  void on_action_3_triggered();

  void on_action_7_triggered();

  void on_action_8_triggered();

  void on_action_9_triggered();

  void on_action_dda_triggered();

  void on_action_11_triggered();

  void on_action_10_triggered();

  void on_action_12_triggered();

  void on_action_13_triggered();

  void on_action_14_triggered();

  void on_action_16_triggered();

  void on_action_1_triggered();

  void on_action_17_triggered();

  void on_action_18_triggered();

  void on_action_19_triggered();

  void on_action_DDA_triggered();

  void on_action_20_triggered();

  void on_action_21_triggered();

private:
  Ui::PCLViewer *ui;

public:
  SimplifyVerticesSize m_svs;   // 顶点化简
  SnappingVertices m_sv;        // 顶点合并
  OSnap m_osnap;	// 位置优化
  HoleFilling m_hf;				// 孔洞填充

  void init();
  void performOptimize();
  ControlOptimizationSnappingDialog m_cos_dlg;
  HoleFillingInteractiveOperationDialog m_hfio_dlg;

  void displayWireframe();
  void displayWireframe2();

private:

  // 参数
  void loadParams();
  // 顶点化简
  float T_angle_betPointAndPoint;	// 选取种子点时，种子点法向量与相邻顶点法向量的最大夹角（单位：弧度）
  float T_angle_betPointAndLine;	// 点的法向量与直线法向量的最大夹角（单位：弧度）
  float T_dist_betPointAndLine;	// 点到直线的最大距离
  float T_lenRatio;               // 确定可以被剪掉的最小线段占平均线段长度的最大比例
  float T_angle_betLines;         // 可被合并的两个线段的最大夹角（单位：弧度）
  // 顶点合并
  float ratio_of_scale;             // 尺度的比例因子，用来设置邻域半径
  float T_ratio_lineSeg_middle_area; // 线段中部区域占据线段长度的比例
  float T_dist_proj_to_lineSegEnd;   // 投影点与线段端点的最大距离
  float T_maxAngle_bet_two_polys;    // 两个平行多边形的法向量最大夹角
  // osnap
  double osnap_delta;               // 数值计算偏导数的步长
  int osnap_max_k_count;          // 共轭梯度的最大迭代次数
};

#endif // PCLVIEWER_H
