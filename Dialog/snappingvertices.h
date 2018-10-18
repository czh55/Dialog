#ifndef SNAPPINGVERTICES_H
#define SNAPPINGVERTICES_H

#include "initialPoly.h"
//#include "earclip.h"

using namespace std;

struct MatchingPoly
{
    int poly_id;
    vector<void *> vertices;
};

// 顶点数据结构
struct Vertex
{
    pcl::PointXYZ point;	// 顶点位置

    float radius;			// 搜索半径

    Vertex* next_point;     // 前驱节点
    Vertex* prev_point;     // 后继节点

    int neighbor_poly_id;        // 邻居多边形id
    bool is_vertex_neighbor;    // 距离最近的是顶点还是边
    Vertex* neighbor_vertex;    // 邻居顶点
    Vertex* neighbor_edge;      // 邻居边

    vector<MatchingPoly> maching_polys; // 与该点具有“点-点”匹配关系的多边形顶点(即哪一个多边形的哪几个顶点对其有匹配关系)

    bool is_updated;

	bool is_intersectPoint;	// 线段相交，标记该顶点是否为交点
};

// 三角形
struct Triangle
{
    /*Vertex *v0;
    Vertex *v1;
    Vertex *v2;*/
    Vertex* v[3];
};

class POLYGON
{
private:
    int size;

public:
    Vertex *start_point;        // 多边形起点
    Vertex *cur_point;
    Eigen::Vector4f coeff;      // 多边形平面方程(a,b,c,d: ax+by+cz+d=0)

    vector<Triangle> triangles; // 对多边形进行三角化之后的小三角形集合

    POLYGON(){
        start_point = cur_point = NULL;
        size = 0;
        coeff << 0, 0, 0, 0;
    }

	~POLYGON(){
        //clear();
	}

	// 清空多边形
	void clear(){
		Vertex *pv, *pv_next;
		for(int i = 0; i < size; ++i){
			pv = i == 0 ? this->start_point : pv_next;
			pv_next = pv->next_point;
			delete pv;
		}
		size = 0;
		coeff << 0, 0, 0, 0;
		this->start_point = this->cur_point = NULL;
	}

    int getSize(){return size;}

    // 插入第一个顶点
    void insertFirstVertex(Vertex *v){
        start_point = v;
        cur_point = v;
        ++size;
    }
    // 插入中间顶点
    void insertMiddleVertex(Vertex *v){
        cur_point->next_point = v;
        v->prev_point = cur_point;
        // 把v作为当前顶点
        cur_point = v;
        ++size;
    }
    // 插入最后一个顶点
    void insertLastVertex(Vertex *v){
        cur_point->next_point = v;
        v->prev_point = cur_point;
        v->next_point = start_point;
        start_point->prev_point = v;
        ++size;
    }

    // 删除顶点
    void delVertex(Vertex *v){
        if(isVertexInPoly(v)){
            --size;

            // 如果v是起始顶点，则更新起始顶点
            if(v == this->start_point){
                this->start_point = v->next_point;
            }

            Vertex* prev_v, * next_v;
            prev_v = v->prev_point;
            next_v = v->next_point;
            prev_v->next_point = next_v;
            next_v->prev_point = prev_v;
            delete v;
            v = NULL;
        }else{
            cout << "vertex is not a member of this Polygons" << endl;
        }
    }

    // 在顶点v后面插入一个新的顶点new_v
    bool insertVertex(Vertex *v, Vertex *new_v){
        if(!isVertexInPoly(v)) return false;

		// 设置new_v的搜索半径
		new_v->radius = v->radius;

        Vertex *v_next = v->next_point;
        v->next_point = new_v;
        new_v->prev_point = v;
        new_v->next_point = v_next;
        v_next->prev_point = new_v;
        ++size;
        return true;
    }

    // 判断点是否在多边形内
    bool isVertexInPoly(Vertex *v){
        if(size == 0) return false;

        for(int i = 0; i < size; ++i){
            cur_point = i == 0 ? start_point : cur_point->next_point;
            if(cur_point == v){
                return true;
            }
        }
        return false;
    }

    bool operator==(const POLYGON& poly) const
    {
       return poly.start_point == this->start_point;
    }
};

// 计算不相邻线段交点的数据结构
// 需保存多边形中哪两条线段相交，以及交点信息
struct IntersectingPoint{
    int edge_id1;   // 边id
    int edge_id2;   // 边id
    pcl::PointXYZ point;
	Vertex *edge_id1_begin;
	Vertex *edge_id2_begin;
};

class SnappingVertices
{
public:
    SnappingVertices();

    float ratio_of_scale;              // 尺度的比例因子，用来设置邻域半径
    float T_ratio_lineSeg_middle_area; // 线段中部区域占据线段长度的比例
    float T_dist_proj_to_lineSegEnd;   // 投影点与线段端点的最大距离
    float T_maxAngle_bet_two_polys;    // 两个平行多边形的法向量最大夹角

    void setParams(float &T1, float &T2, float &T3, float &T4);

    // 多边形结构
    vector<POLYGON, Eigen::aligned_allocator<POLYGON>> polygons;
    // 具有匹配关系的顶点簇
    vector<vector<Vertex*>> matchingVerticesClusters;
    // 具有匹配关系的顶点簇
    vector<vector<Vertex*>> clusters;

    // 构造多边形集合
    void setPolygons(vector<InitialPoly> &initial_polys);
    // 寻找多边形P的顶点v距离最近的其它多边形
    void getNeighborPolygon(int poly_id/*P*/, Vertex* v);
    // 为多边形的每个顶点设置邻域关系
    void setNeighborsForEachPoly();
    // 在获取邻居多边形的时候为每个顶点构建了单向的匹配关系
    // 该函数将构建双向的匹配关系
    void setMatchingRelations();
    // 在顶点v的匹配多边形中检查是否存在编号为id的多边形，如果有则设置其位置
    bool findPoly(Vertex *v, int poly_id, int &poly_id_in_v);
    // 匹配点是否存在
    bool findVertex(Vertex *v, int poly_id_in_v, Vertex *neighbor_v);

    // 精化匹配关系，主要处理以下情况：
    // 顶点v具有“点-边”匹配关系，同时别的顶点v_与v又具有"点-点”匹配
    // 处理策略：本着尽可能实现点-点匹配的原则，去除v的“点-边”匹配
    void refineMachingRelations();

    // 提取匹配的顶点簇
    void getMatchingVerticesClusters(vector<vector<Vertex*>> &clusters);

    // 处理点-点匹配
    // 合并具有匹配关系的点集
    void mergeMatchingVerticesClusters(vector<vector<Vertex*>> &clusters);
    void dealWithEachCluster(vector<Vertex*> &cluster);
    // 仅处理分别属于两个多边形的顶点簇
    void dealWithEachCluster2(vector<Vertex*> &cluster);
    // 计算两个多边形的交线
    bool getIntersectLine(Eigen::Vector4f &coeff1, Eigen::Vector4f &coeff2,
                          Eigen::Vector3f &p_base, Eigen::Vector3f &line_dir);
    // 获取点在直线上的投影点
    void getProjPointOnLine(pcl::PointXYZ &p, pcl::PointXYZ &p_proj,
                            Eigen::Vector3f &p_base, Eigen::Vector3f &line_dir);

    // 计算三个多边形的交点
    void getIntersectPointOfTriplePolys(Eigen::Vector4f &coeff1,
                                        Eigen::Vector4f &coeff2,
                                        Eigen::Vector4f &coeff3,
                                        pcl::PointXYZ &p);

    // 计算平面法向量的夹角
    float getIncludedAngleOfTwoPolys(int poly_a, int poly_b);

    // 判断两个多边形是否平行
    bool areTwoPolysParallel(int poly_a, int poly_b);

    // 获取点在线段上的投影点
    // e1, e2: 线段端点
    void getProjPointOnLineSeg(pcl::PointXYZ &p, pcl::PointXYZ &p_proj, pcl::PointXYZ &e1, pcl::PointXYZ &e2);
    // 判断两个顶点是否严格相等
    bool isEqual(pcl::PointXYZ &p1, pcl::PointXYZ &p2){
        return (p1.x == p2.x && p1.y == p2.y && p1.z == p2.z);
    }

    // 距离函数
    // 点-点
    float distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2);
    // 点-边
    // 如果点p的投影点落在线段内部，则返回p到投影点的距离；否则返回FLT_MAX
    float distP2E(pcl::PointXYZ &p, pcl::PointXYZ &e1, pcl::PointXYZ &e2);

    // 处理“点-边”匹配
    void dealWithVertexEdgeMatchingRelation();

    // 处理多边形退化情况
    void dealWithPolysDegeneration(vector<POLYGON, Eigen::aligned_allocator<POLYGON>> &polygons);   // 主体函数
    void dealWithSinglePolyDegeneration(POLYGON &poly);  // 单个多边形的处理函数
    // 处理顶点冗余(相邻顶点位置重合，不相邻顶点位置重合）
    void dealWithVerticesRedundancy(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices);
    void dealWithLineSegCollineation(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices);         // 处理相邻线段共线
    void dealWithLineSegIntersect(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices);        // 处理非近邻线段的相交
    void cleanMatchingRelaion();    // 清理每个顶点的匹配关系

    // 子函数
    // 计算顶点index1到index2经过路径上线段的长度之和
    float lineSegDist(int index1, int index2, vector<float> &line_lens);

    // 计算两条线段的交点
    // 只有交点均落在两条线段内部，才会返回合法的值
    bool getIntersectPointOfTwoLineSeg(pcl::PointXYZ &line1_begin, pcl::PointXYZ &line1_end,
                                       pcl::PointXYZ &line2_begin, pcl::PointXYZ &line2_end,
                                       pcl::PointXYZ &intersect_point);
    // 计算一个多边形中所有不相邻线段的交点
    void getIntersectPointsOfPoly(pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices,
                                  vector<IntersectingPoint, Eigen::aligned_allocator<IntersectingPoint>> &intersect_points);
    // 根据编号判断两条线段是否相连
    bool isConnectTwoLineSegs(int size, int edge_id1, int edge_id2);

	// 多边形周长函数
	float getCircumferenceOfPoly(pcl::PointCloud<pcl::PointXYZ>::Ptr &border);

	// 将交点插入到多边形中
	void insertIPVertexIntoPoly(Vertex *v, POLYGON &poly,const int edge_id, vector<Vertex*> &poly_vertices);

    // 设置顶点的搜索半径
    void setSearchRadiusForEachVertex(vector<InitialPoly> &initial_polys);

    // 一步执行函数
    // 参数设置以及数据准备需提前备好
    void perform(vector<InitialPoly> &initial_polys);

	// 在具有角对边关系的线段上插入顶点
	void insertVertexOnEdgeWithAngleToEdgeRelation();
	void checkLineSegThatPassThroughVertex(int poly_id, Vertex *pv);
	bool isVertexOnEdge(Vertex *pv, Vertex *edge_begin);

    // 检查当前模型是否合法
    // 如果检测到非法情况，则返回问题顶点
    bool is_model_valid(Vertex *&prob_v);
    // 编号为poly_id的多边形的顶点v是否落在其它多边形的边上
    bool is_vertex_on_edge(int poly_id, Vertex* v);
    // 编号为poly_id的多边形的顶点v是否落在其它多边形的内部
    bool is_vertex_in_poly(int poly_id, Vertex* v);
    // 某点是否在三角形内部
    bool is_vertex_in_triangle(pcl::PointXYZ &point, Triangle &tri);
    // 线段[e1, e2]是否穿过其它多边形的内部
    bool is_edge_pass_through_poly(int poly_id, Vertex* e1, Vertex* e2);
};

#endif // SNAPPINGVERTICES_H
