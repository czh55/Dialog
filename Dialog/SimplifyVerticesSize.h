#pragma once
#include "initialPoly.h"
#include <pcl/segmentation/sac_segmentation.h>

class SimplifyVerticesSize
{
public:
	SimplifyVerticesSize(void);
    ~SimplifyVerticesSize(void);
    inline
    void setParams(float T1, float T2, float T3, float T4, float T5){
        T_angle_betPointAndPoint = T1;
        T_angle_betPointAndLine = T2;
        T_dist_betPointAndLine = T3;
        T_lenRatio = T4;
        T_angle_betLines = T5;

        // 角度转弧度
        double PI = 3.14159265358979;
        T_angle_betPointAndPoint = T_angle_betPointAndPoint * PI / 180.0f;
        T_angle_betPointAndLine = T_angle_betPointAndLine * PI / 180.0f;
        T_angle_betLines = T_angle_betLines * PI / 180.0f;
    }

    // 构造多边形顶点法向量(默认顶点排序已经满足右手法则)
	// 仅利用前后顶点以及当前顶点确定当前顶点的法向量
	void computeVerticesNormals(pcl::PointCloud<pcl::PointXYZ>::Ptr &contour, 
	                        std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals,
							Eigen::Vector3f &plane_normal);
	// 区域生长构造初始的线段集合
	void constructInitLineSegs( pcl::PointCloud<pcl::PointXYZ>::Ptr &contour,
								std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals,
								Eigen::Vector3f &poly_normal,
								std::vector<std::vector<int>> &line_segs );
	std::vector<int> final_line_segs;
	// 继续优化多边形的顶点序列
	void refineLineSegs( pcl::PointCloud<pcl::PointXYZ>::Ptr &contour,			// 原始顶点
						 pcl::PointCloud<pcl::PointXYZ>::Ptr &final_vertices,	// 最终顶点
						 std::vector<std::vector<int>> &retrieved_line_segs, // 区域生长法提取的线段序列
						 std::vector<int> &final_line_segs );	 // 最终的线段序列

    //void simplifyPolyVertices(MyPoly &poly);

	void vectorSort(std::vector<float> &v, std::vector<int> &indices);

protected:
	void getSeedIndex( pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud,
					   int &seed_index,
					   std::vector<bool> &is_processed,
					   std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> &vertices_normals );
	float distP2Line( pcl::PointXYZ &point,
					  Eigen::Vector3f &p_base,	// 直线模型的基点
					  Eigen::Vector3f &line_dir );
	float distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2);

private:
    float T_angle_betPointAndPoint;	// 选取种子点时，种子点法向量与相邻顶点法向量的最大夹角（单位：弧度）
    float T_angle_betPointAndLine;	// 点的法向量与直线法向量的最大夹角（单位：弧度）
    float T_dist_betPointAndLine;	// 点到直线的最大距离
    float T_lenRatio;               // 确定可以被剪掉的最小线段占平均线段长度的最大比例
    float T_angle_betLines;         // 可被合并的两个线段的最大夹角（单位：弧度）
};

