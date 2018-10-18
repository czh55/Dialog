#pragma once

//PCL
#include <pcl\point_cloud.h>
#include <pcl/common/centroid.h>
#include <pcl/visualization/pcl_visualizer.h>

//STL
#include <iostream>
#include <string>
#include <vector>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_conformer_2.h>

struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){ 
    return nesting_level%2 == 1;
  }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;

class ConstrainedDelaunayTriangulation
{
public:
	ConstrainedDelaunayTriangulation(void);
	~ConstrainedDelaunayTriangulation(void);

	Polygon_2 m_ply;

    /***************************************只有多边形边界上的点*****************************************/
	// 输入多边形顶点序列
	void setInputPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud_ply);
	// 执行CDT三角化(不插入新的顶点)
	// triangles: 多边形集合(N-2个，N为多边形顶点个数)
	// tri_index_with_min_angle: 具有最小内角的三角形索引
	// min_angle: 三角形的最小内角
    void triangulate(std::vector<std::vector<int>> &triangles, std::vector<float> &min_angles);
    /**********************************多边形边界上的点+多边形内部的孤立点*********************************/
    // 输入多边形顶点序列以及孤立顶点
    // 所有孤立顶点默认均位于多边形的内部
    void setInputPointCloud_and_triangulate(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, int ply_size,
                                             std::vector<std::vector<int>> &triangles, std::vector<float> &min_angles);

private:
	void mark_domains(CDT& ct, CDT::Face_handle start, int index, std::list<CDT::Edge>& border);
	void mark_domains(CDT& cdt);
	// 执行CDT三角化(不插入新的顶点)
	void triangulate(std::vector<std::vector<int>> &triangles);
};

