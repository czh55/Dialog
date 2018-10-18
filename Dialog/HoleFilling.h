#ifndef HOLEFILLING_H
#define HOLEFILLING_H

#include "snappingvertices.h"
#include <pcl\point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

struct POLYGON_RELEVANT
{
    void *poly;
	void *prev_point;
	void *next_point;
};

struct Vertex2
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Eigen::Vector3d pos;		// 顶点位置
	Eigen::Vector3d proj_pos;	// 局部投影的位置
	Eigen::Vector4d coeff;		// 局部投影平面方程
	double state_angle;			// 状态角度
	vector<POLYGON_RELEVANT> relevant_polys;	// 关联的多边形信息
    int index_in_m_vertices;
};

struct Triangle2
{
	Vertex2* v[3];
};

class POLYGON2
{
public:
    int size;

public:
	// http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Vertex2 *start_point;        // 多边形起点
    Vertex2 *cur_point;
    Eigen::Vector4d coeff;      // 多边形平面方程(a,b,c,d: ax+by+cz+d=0)

    vector<Triangle2> triangles; // 对多边形进行三角化之后的小三角形集合
    vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> discontinuity_vertices;

	bool is_tri;	// 是否需要执行三角化
    bool is_discontinuity;  // 多边形是否是不连续面

    POLYGON2(){
        start_point = cur_point = NULL;
        size = 0;
        coeff << 0, 0, 0, 0;
		is_tri = true;
        is_discontinuity = false;
    }

	~POLYGON2(){
        //clear();
	}

	// 清空多边形
	void clear(){
        // 首先处理顶点的关联多边形信息
        Vertex2* pv2, *next_pv2;
        int poly_index_in_v;
        for(int i = 0; i < size; ++i){
            if(i == 0){
                pv2 = this->start_point;
            }else{
                pv2 = next_pv2;
            }
            poly_index_in_v = this->getPolyIndex(pv2);
            next_pv2 = (Vertex2*)pv2->relevant_polys[poly_index_in_v].next_point;

            auto iter = pv2->relevant_polys.begin();
            for(int j = 0; j <= poly_index_in_v; ++j){
                if(j == 0)
                    iter = pv2->relevant_polys.begin();
                else
                    ++iter;
            }
            pv2->relevant_polys.erase(iter);
        }

		size = 0;
		coeff << 0, 0, 0, 0;
		this->start_point = this->cur_point = NULL;
        triangles.clear();
        is_tri = true;
        discontinuity_vertices.clear();
	}

    int getSize(){return size;}

	int getPolyIndex(Vertex2 *v){
		int poly_index = -1;
		for(int i = 0; i < v->relevant_polys.size(); ++i){
			if(v->relevant_polys[i].poly == this){
				poly_index = i;
				break;
			}
		}
		return poly_index;
	}

    ///////////////////////////多边形建立阶段-begin////////////////////////////////
    // 插入第一个顶点
    void insertFirstVertex2(Vertex2 *v){
        start_point = v;
        cur_point = v;
        ++size;

		POLYGON_RELEVANT pr;
		pr.poly = this;
		pr.next_point = pr.prev_point = NULL;
		v->relevant_polys.push_back(pr);
    }
    // 插入中间顶点
    void insertMiddleVertex2(Vertex2 *v){
		int poly_index = getPolyIndex(cur_point);
		cur_point->relevant_polys[poly_index].next_point = v;
		POLYGON_RELEVANT pr;
		pr.poly = this;
		pr.prev_point = cur_point;
		pr.next_point = NULL;
		v->relevant_polys.push_back(pr);

        // 把v作为当前顶点
        cur_point = v;
        ++size;
    }
    // 插入最后一个顶点
    void insertLastVertex2(Vertex2 *v){
		int poly_index = getPolyIndex(cur_point);
		cur_point->relevant_polys[poly_index].next_point = v;
		POLYGON_RELEVANT pr;
		pr.poly = this;
		pr.prev_point = cur_point;
		pr.next_point = start_point;
		v->relevant_polys.push_back(pr);
		poly_index = getPolyIndex(start_point);
		start_point->relevant_polys[poly_index].prev_point = v;

        ++size;
    }
    ///////////////////////////多边形建立阶段-end////////////////////////////////

    // 删除顶点
    void delVertex2(Vertex2 *v){
        is_tri = true;
        if(isVertex2InPoly(v)){
            if(size == 3){
                cerr<<"当前多边形size == 3, 无法执行顶点删除操作." << endl;
                return;
            }

            --size;
			int poly_index = getPolyIndex(v);

            // 如果v是起始顶点，则更新起始顶点(用起始顶点的下一个顶点作为起始顶点)
            if(v == this->start_point){
				this->start_point = (Vertex2 *)v->relevant_polys[poly_index].next_point;
            }

            Vertex2* prev_v, * next_v;
			prev_v = (Vertex2 *)v->relevant_polys[poly_index].prev_point;
			next_v = (Vertex2 *)v->relevant_polys[poly_index].next_point;
			int prev_v_poly_index = getPolyIndex(prev_v);
			prev_v->relevant_polys[prev_v_poly_index].next_point = next_v;
			int next_v_poly_index = getPolyIndex(next_v);
			next_v->relevant_polys[next_v_poly_index].prev_point = prev_v;
            /*delete v;
            v = NULL;*/
			for(auto iter = v->relevant_polys.begin(); iter != v->relevant_polys.end(); ++iter){
				if(iter->poly == this){
					v->relevant_polys.erase(iter);
					break;
				}
            }
        }else{
            cout << "vertex is not a member of this Polygon" << endl;
        }
    }

    // 在顶点v后面插入一个新的顶点new_v
	// 前提：已建立好多边形
    bool insertVertex2(Vertex2 *v, Vertex2 *new_v){
        is_tri = true;
        // 首先验证new_v是否位于多边形的支撑平面上
        double d =  new_v->pos[0]*coeff[0] +
                    new_v->pos[1]*coeff[1] +
                    new_v->pos[2]*coeff[2] +
                                  coeff[3];
        if(abs(d) > 0.00001){
            cerr << "new_v 与平面的距离d = " << d << endl;
            return false;
        }

        // 当前顶点必须在已有多边形中
        if(!isVertex2InPoly(v)){
            cerr << "vertex v is not in ply" << endl;
            return false;
        }

        // 新顶点绝不能已经存在于当前多边形中
        if(isVertex2InPoly(new_v)){
            cerr << "vertex new_v is in ply!" << endl;
            return false;
        }

        // 找到当前顶点的下一个顶点
        int ply_index_in_v = getPolyIndex(v);
        Vertex2* next_v = (Vertex2*)v->relevant_polys[ply_index_in_v].next_point;

        // 将当前顶点的后继指向new_v
        v->relevant_polys[ply_index_in_v].next_point = new_v;
        // 将当前下一个顶点的前驱指向new_v
        int ply_index_in_next_v = getPolyIndex(next_v);
        next_v->relevant_polys[ply_index_in_next_v].prev_point = new_v;

        // 为新顶点生成新的关联关系
        POLYGON_RELEVANT ply_rel;
		ply_rel.poly = this;
        ply_rel.prev_point = v;
        ply_rel.next_point = next_v;
        new_v->relevant_polys.push_back(ply_rel);

        ++size;
    }

    // 获取当前顶点的前驱
    Vertex2* getPrevVertex(Vertex2* v){
        int ply_index = this->getPolyIndex(v);
        if(ply_index != -1)
            return (Vertex2*)v->relevant_polys[ply_index].prev_point;
        else
            return NULL;
    }

    // 获取当前顶点的后继
    Vertex2* getNextVertex(Vertex2* v){
        int ply_index = this->getPolyIndex(v);
        if(ply_index != -1)
            return (Vertex2*)v->relevant_polys[ply_index].next_point;
        else
            return NULL;
    }

    // 判断顶点是否在多边形内
    bool isVertex2InPoly(Vertex2 *v){
        if(size == 0) return false;

        for(int i = 0; i < size; ++i){
            //cur_point = i == 0 ? start_point : cur_point->next_point;
			if(i == 0)
				cur_point = start_point;
			else{
				int poly_index = getPolyIndex(cur_point);
				cur_point = (Vertex2 *)cur_point->relevant_polys[poly_index].next_point;
			}

            if(cur_point == v){
                return true;
            }
        }
        return false;
    }
};


class HoleFilling
{
public:
    HoleFilling();

	// 外部的多边形信息
	vector<POLYGON, Eigen::aligned_allocator<POLYGON>> *sv_polygons;

	// 全局顶点序列
	pcl::PointCloud<pcl::PointXYZ>::Ptr m_vertices_cloud;
	pcl::KdTreeFLANN<pcl::PointXYZ> m_vertices_kdtree;
	vector<Vertex2*> m_vertices;
	// 多边形集合
	vector<POLYGON2*> m_polygons;

    // 用来确定多边形id
    pcl::PointCloud<pcl::PointXYZ>::Ptr m_cloud_for_indentify_poly_id;
    pcl::KdTreeFLANN<pcl::PointXYZ> m_cloud_for_indentify_poly_id_kdtree;

	// 清理数据结构
	void clear();

	// step1: 对多边形数据进行格式转换
	void setInputPolygons(vector<POLYGON, Eigen::aligned_allocator<POLYGON>> *p_polygons);

    // step2: 设置m_cloud_for_indentify_poly_id
    void set_cloud_for_indentify_poly_id();

	// 更新顶点数据(顶点与对应的顶点点云同步更新）
	void updateVertices();

    // 获取一个顶点单向匹配关系的数量
    int getSingleRelationCount(Vertex2* v);

	// 基于多边形的三角化剖分结果更新三角面片重心点云
	void updateTriCentroids();

	// 更新顶点点云(包括：删除没有关联关系的顶点，更新顶点点云和对应的kdtree）
	void updateVerCloud();
};

#endif // HOLEFILLING_H
