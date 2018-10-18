#ifndef EARCLIP_H
#define EARCLIP_H
#include <pcl/point_cloud.h>
#include <vector>
#include "snappingvertices.h"

class EarClip
{
public:
    EarClip();
    // 设置输入的多边形
    void setInputPoly(POLYGON *pPoly);
    // 对多边形进行三角化
    bool triangulatePoly();

private:
    POLYGON * m_poly;
    vector<Vertex *> vertices;  // 多边形顶点位置顺序存储

    double PI;

    // 获取索引i的前驱结点
    int getPrev(std::vector<int> &indices, int i);
    // 获取索引i的后继结点
    int getNext(std::vector<int> &indices, int i);
    // 删除vector中的第i个元素
    void del_vertex(std::vector<int> &indices, int i);
    // 某点是否在三角形内部
    // 目标顶点的索引: p_index
    // 三角形顶点索引: tri_v1, tri_v2, tri_v3
    bool isPointInTriangle(int p_index, int tri_v1, int tri_v2, int tri_v3);
    // 获取多边形顶点的内角
    // vertices: 顶点序列，逆时针排列
    // normal: 多边形法向量，指向物体外部
    // index: 目标顶点索引
    // prev: index的前驱
    // next: index的后继
    double getInnerAngle(vector<Vertex *> &vertices, Eigen::Vector3f &normal, int index, int prev, int next);
    // 点-点距离函数
    float distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2);
};

#endif // EARCLIP_H
