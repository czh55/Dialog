#pragma once

//#ifndef Q_MOC_RUN
//// Point Cloud Library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
//#endif

class InitialPoly
{
public:
    InitialPoly(){
        vertices.reset(new pcl::PointCloud<pcl::PointXYZ>);
        simplified_vertices.reset(new pcl::PointCloud<pcl::PointXYZ>);
        normal << 0, 0, 0;
        scale = 0;
    }

    // 最初的顶点序列
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices;
    // 多边形法向量
    Eigen::Vector3f normal;
    // 多边形尺度
    float scale;

    // 顶点化简：
    // 化简之后的顶点序列
    pcl::PointCloud<pcl::PointXYZ>::Ptr simplified_vertices;
    // 顶点法向量
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f>> vertices_normals;
    // 区域生长时检测到的线段序列
    std::vector<std::vector<int>> line_segs;
    // 最终的线段序列
    std::vector<int> final_line_segs;
};
