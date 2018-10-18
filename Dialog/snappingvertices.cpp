#include "snappingvertices.h"
#include "earclip.h"
#include <queue>
#include <set>
#include <pcl/kdtree/kdtree_flann.h>
/*
 SnappingVertices使用方法：
    第一步：
    setParams
    第二步：
    // 用新的数据结构代表多边形
    this->m_sv.setPolygons(this->initial_polys);
    第三步：
    // 为多边形的每个顶点设置邻域关系
    this->m_sv.setNeighborsForEachPoly();
    第四步：
    // 将单向的匹配关系转换为双向的匹配关系
    m_sv.setMatchingRelations();
    第五步：匹配作为一种连通关系，将互相有匹配关系的顶点聚类成一簇
    // 提取具有匹配关系的顶点簇
    m_sv.getMatchingVerticesClusters(clusters);
    第六步：将具有点-边匹配的顶点匹配关系按是否有别的点与其匹配转换为点-点匹配
    // 精化匹配关系
    refineMachingRelations
    第七步：(顶点位置发生改变)
    // 处理每个具有“点-点”匹配关系的顶点簇
    mergeMatchingVerticesClusters
    第八步：处理“点-边”匹配
    dealWithVertexEdgeMatchingRelation()
    第九步：处理多边形的退化情况
    dealWithPolysDegeneration()
*/

SnappingVertices::SnappingVertices()
{
}

void SnappingVertices::setParams(float &T1, float &T2, float &T3, float &T4){
    ratio_of_scale = T1;
    T_ratio_lineSeg_middle_area = T2;
    T_dist_proj_to_lineSegEnd = T3;
    T_maxAngle_bet_two_polys = T4;
    float PI = 3.14159265358979;
    // 角度转弧度
    T_maxAngle_bet_two_polys = T_maxAngle_bet_two_polys * PI / 180.0;
}

// 构造多边形集合
void SnappingVertices::setPolygons(vector<InitialPoly> &initial_polys){
    this->polygons.clear();

    for(int i = 0; i < initial_polys.size(); ++i){
        POLYGON poly;
        for(int j = 0; j < initial_polys[i].simplified_vertices->size(); ++j){
            Vertex *p_vertex = new Vertex();
            p_vertex->is_updated = false;
            p_vertex->radius = ratio_of_scale*initial_polys[i].scale;   // 设置顶点的搜索半径
            // 第i个多边形的第j个顶点
            p_vertex->point = initial_polys[i].simplified_vertices->points[j];
            p_vertex->point.data[3] = i;    // data[3]存储该点对应的多边形的索引
            p_vertex->next_point = p_vertex->prev_point = NULL;

            p_vertex->neighbor_poly_id = -1;    // 邻居多边形编号为-1表明没有匹配关系
            p_vertex->neighbor_vertex = p_vertex->neighbor_edge = NULL;

            // 起点
            if(j == 0){
                poly.insertFirstVertex(p_vertex);
            }else if(j == initial_polys[i].simplified_vertices->size() - 1){
            // 终点
                poly.insertLastVertex(p_vertex);
            }else{
            // 中间的点
                poly.insertMiddleVertex(p_vertex);
            }
        }

        // 多边形平面参数
        pcl::PointXYZ *p;
        InitialPoly *pPoly = &initial_polys[i];
        p = &initial_polys[i].simplified_vertices->points[0];
        float d = -1 * (p->x*pPoly->normal[0] + p->y*pPoly->normal[1] +
                        p->z*pPoly->normal[2]);
        poly.coeff << pPoly->normal[0],
                      pPoly->normal[1],
                      pPoly->normal[2],
                      d;

        polygons.push_back(poly);
    }
}

// 寻找多边形P的一个顶点v距离最近的其它多边形
void SnappingVertices::getNeighborPolygon(int poly_id/*P*/, Vertex* v){
    float dist = FLT_MAX;

    for(int i = 0; i < polygons.size(); ++i){
        if(i == poly_id) continue;

        // 首先看顶点，保存与v距离最近的顶点
        float vertex_dist = FLT_MAX;
        Vertex *neighbor_v = NULL; // 保存与v距离最近的顶点
        Vertex *cur_v;
        for(int j = 0; j < polygons[i].getSize(); ++j){ // 遍历第i个多边形的所有顶点
            cur_v = j == 0 ? polygons[i].start_point : cur_v->next_point;
            float d_ = distP2P(v->point, cur_v->point);
            if(d_ < vertex_dist){
                vertex_dist = d_;
                neighbor_v = cur_v;
            }
        }

        // 然后看顶点与线段
        float edge_dist = FLT_MAX;
        Vertex *neighbor_e = NULL; // 保存与v距离最近的线段
        Vertex *cur_e;
        for(int j = 0; j < polygons[i].getSize(); ++j){ // 遍历第i个多边形的所有边
            cur_e = j == 0 ? polygons[i].start_point : cur_e->next_point;
            float d_ = distP2E(v->point, cur_e->point, cur_e->next_point->point);
            if(d_ < edge_dist){
                edge_dist = d_;
                neighbor_e = cur_e;
            }
        }

        float dist_ = vertex_dist <= edge_dist ? vertex_dist : edge_dist;

        if(dist <= dist_) continue; // 如果保存的最近距离比当前距离更小，则检查下一个多边形

        dist = dist_;               // 否则，更新当前距离
        v->neighbor_poly_id = i;    // 并更新距离最近的多边形编号

        // 与v离得最近的是顶点
        if(neighbor_e == NULL || vertex_dist <= edge_dist){
            v->is_vertex_neighbor = true;
            v->neighbor_vertex = neighbor_v;
			v->neighbor_edge = NULL;
        }else{ // 与v离得最近的是线段
            v->is_vertex_neighbor = false;
            v->neighbor_edge = neighbor_e;
			v->neighbor_vertex = NULL;
		}
    }

    // 如果当前距离超过了顶点v的搜索半径，则顶点v没有任何匹配
    if(dist > v->radius){
        v->neighbor_edge = v->neighbor_vertex = NULL;
        v->neighbor_poly_id = -1;
        return;
    }

    // 当前距离没有超过定点v的搜索半径
    // 匹配关系为点-点
    if(v->is_vertex_neighbor){
        if(dist > v->neighbor_vertex->radius){  // 距离同样不能大于匹配顶点的搜索半径
            v->neighbor_edge = v->neighbor_vertex = NULL;
            v->neighbor_poly_id = -1;
            return;
        }
        return;
    }

    // 匹配关系为点-边，这个要稍微复杂一些，遵循的基本原则：能尽量匹配点就不匹配边
    pcl::PointXYZ e_begin/*线段起点*/, e_end/*线段终点*/, e_middle/*线段中间的点*/;
    e_begin = v->neighbor_edge->point;
    e_end = v->neighbor_edge->next_point->point;
    e_middle.x = (e_begin.x + e_end.x)/2.0;
    e_middle.y = (e_begin.y + e_end.y)/2.0;
    e_middle.z = (e_begin.z + e_end.z)/2.0;

    pcl::PointXYZ p_proj;   // v在e上的投影点
    getProjPointOnLineSeg(v->point, p_proj, e_begin, e_end);

    float lineSeg_len = distP2P(e_begin, e_end);
    float half_lineSeg_len = lineSeg_len * 0.5;

    pcl::PointXYZ e_nearer = distP2P(p_proj, e_begin) <= distP2P(p_proj, e_end) ? e_begin : e_end;

    if(distP2P(p_proj, e_middle) <= T_ratio_lineSeg_middle_area * half_lineSeg_len){    // 投影点落在线段中部区域
        // 点-边匹配
        if(dist > v->neighbor_edge->radius){    // 点到边的距离同样不能超过边的搜索半径
            v->neighbor_poly_id = -1;
            v->neighbor_edge = v->neighbor_vertex = NULL;
        }
        return;
    }else if(distP2P(p_proj, e_nearer) <= T_dist_proj_to_lineSegEnd &&
             distP2P(e_nearer, v->point) <= v->radius &&
             distP2P(e_nearer, v->point) <= v->neighbor_edge->radius){
        // 点-点匹配
        v->is_vertex_neighbor = true;
        v->neighbor_vertex = isEqual(e_nearer, e_begin) ? v->neighbor_edge : v->neighbor_edge->next_point;
        v->neighbor_edge = NULL;
    }else{
        // 点-边匹配
        if(dist > v->neighbor_edge->radius){    // 点到边的距离同样不能超过边的搜索半径
            v->neighbor_poly_id = -1;
            v->neighbor_edge = v->neighbor_vertex = NULL;
        }
        return;
    }
}

// 为多边形的每个顶点设置邻域关系
void SnappingVertices::setNeighborsForEachPoly(){
    for(int i = 0; i < polygons.size(); ++i){
        POLYGON *pPoly = &polygons[i];
        Vertex* pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            getNeighborPolygon(i, pv);
        }
    }
}

// 距离函数
// 点-点
float SnappingVertices::distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2){
    return pow((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z), 0.5f);
}

// 点-边距离
// 如果点p的投影点落在线段内部，则返回p到投影点的距离；否则返回FLT_MAX
float SnappingVertices::distP2E(pcl::PointXYZ &p, pcl::PointXYZ &e1, pcl::PointXYZ &e2){
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

    Eigen::Vector3f p_proj;
    p_proj = p_base + lambda * line_dir;

    Eigen::Vector3f lambda1/*向量：e1e2*/, lambda2/*向量：e1p_proj*/;
    lambda1 << e2.x - e1.x,
               e2.y - e1.y,
               e2.z - e1.z;
    lambda2 << p_proj - p_base;

    bool is_in_lineSeg;
    if(lambda1.dot(lambda2) < 0){
        is_in_lineSeg = false;
    }else{
        if(lambda2.norm()/lambda1.norm() >= 1.0){
            is_in_lineSeg = false;
        }else{
            is_in_lineSeg = true;
        }
    }

    if(is_in_lineSeg){
        return (v_p-p_proj).norm();
    }

    return FLT_MAX;
}

// 获取点在线段上的投影点
void SnappingVertices::getProjPointOnLineSeg(pcl::PointXYZ &p, pcl::PointXYZ &p_proj, pcl::PointXYZ &e1, pcl::PointXYZ &e2){
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

    Eigen::Vector3f p_proj_v;
    p_proj_v = p_base + lambda * line_dir;
    p_proj.x = p_proj_v[0];
    p_proj.y = p_proj_v[1];
    p_proj.z = p_proj_v[2];
}

// 在获取邻居多边形的时候为每个顶点构建了单向的匹配关系
// 该函数将构建双向的匹配关系
void SnappingVertices::setMatchingRelations(){
    for(int i = 0; i < polygons.size(); ++i){
        POLYGON *pPoly = &polygons[i];
        Vertex *pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            // pv: 当前顶点
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            // 如果是”点-边“匹配或没有匹配则检查下一个顶点
            if(!pv->is_vertex_neighbor || pv->neighbor_poly_id == -1) continue;

            // “点-点”匹配
            if(pv->neighbor_poly_id == -1){// 异常处理
                std::cerr << "pv->neighbor_poly_id == -1 in SnappingVertices::setMatchingRelations()" << std::endl;
                continue;
            }

            // 匹配点对应的多边形是否已存在于v的匹配多边形中
            int poly_id_in_v;
            if(findPoly(pv, pv->neighbor_poly_id, poly_id_in_v)){ // 如果多边形存在
                // 匹配点是否已存在
                if(findVertex(pv, poly_id_in_v, pv->neighbor_vertex)){ // 如果顶点存在
                    // do nothing
                }else{ // 如果顶点不存在
                    pv->maching_polys[poly_id_in_v].vertices.push_back((void *)pv->neighbor_vertex);
                }
            }else{ // 多边形不存在
                MatchingPoly mp;
                mp.poly_id = pv->neighbor_poly_id;
                mp.vertices.push_back((void *)pv->neighbor_vertex);
                pv->maching_polys.push_back(mp);
            }

            // 在匹配点的匹配关系中设置v
            int poly_id_in_neighbor;
            if(findPoly(pv->neighbor_vertex, i, poly_id_in_neighbor)){
                if(findVertex(pv->neighbor_vertex, poly_id_in_neighbor, pv)){
                }else{
                    pv->neighbor_vertex->maching_polys[poly_id_in_neighbor].vertices.push_back((void*)pv);
                }
            }else{
                MatchingPoly mp;
                mp.poly_id = i;
                mp.vertices.push_back((void*) pv);
                pv->neighbor_vertex->maching_polys.push_back(mp);
            }
        }
    }
}

// 在顶点v的匹配多边形中检查是否存在编号为id的多边形，如果有则设置其位置
bool SnappingVertices::findPoly(Vertex *v, int poly_id, int &poly_id_in_v){
    for(int i = 0; i < v->maching_polys.size(); ++i){
        if(v->maching_polys[i].poly_id == poly_id){
            poly_id_in_v = i;
            return true;
        }
    }
    poly_id_in_v = -1;
    return false;
}

// 匹配点是否存在
bool SnappingVertices::findVertex(Vertex *v, int poly_id_in_v, Vertex *neighbor_v){
    for(int i = 0; i < v->maching_polys[poly_id_in_v].vertices.size(); ++i){
        if((Vertex *)(v->maching_polys[poly_id_in_v].vertices[i]) == neighbor_v)
            return true;
    }
    return false;
}

// 精化匹配关系，主要处理以下情况：
// 顶点v具有“点-边”匹配关系，同时别的顶点v_与v又具有"点-点”匹配
// 处理策略：本着尽可能实现点-点匹配的原则，去除v的“点-边”匹配
void SnappingVertices::refineMachingRelations(){
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < polygons.size(); ++i){
        pPoly = &polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            // 如果是“点-边”匹配
            if(!pv->is_vertex_neighbor && pv->neighbor_poly_id != -1){
                // 看pv是否有“点-点”匹配关系
                if(pv->maching_polys.size() > 0){
                    pv->neighbor_edge = NULL;
                    pv->neighbor_poly_id = pv->maching_polys[0].poly_id;
                    pv->neighbor_vertex = (Vertex *)pv->maching_polys[0].vertices[0];
                    pv->is_vertex_neighbor = true;
                }
            }
        }
    }
}

// 提取匹配的顶点簇
void SnappingVertices::getMatchingVerticesClusters(vector<vector<Vertex*>> &clusters)
{
    // 清理数据结构
    if(clusters.size() > 0){
        for(int i = 0; i < clusters.size(); ++i){
            clusters[i].clear();
        }
    }
    clusters.clear();

    // 设置每个顶点的被处理状态(is_updated)
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < this->polygons.size(); ++i){
        pPoly = &this->polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pv->is_updated = false;
        }
    }

    // 寻找簇
    for(int i = 0; i < this->polygons.size(); ++i){
        pPoly = &this->polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;

            if(pv->is_updated) continue;

            if(pv->maching_polys.size() == 0){
                pv->is_updated = true;
                continue;
            }

            // pv->maching_polys.size() > 0
            vector<Vertex*> cluster;
            std::queue<Vertex*> Q;
            Q.push(pv);
            while(!Q.empty()){
                Vertex *cur_v = Q.front();
                Q.pop();
                cluster.push_back(cur_v);
                cur_v->is_updated = true;

                Vertex *pv2;
                for(int k = 0; k < cur_v->maching_polys.size(); ++k){
                    for(int l = 0; l < cur_v->maching_polys[k].vertices.size(); ++l){
                        pv2 = (Vertex *)cur_v->maching_polys[k].vertices[l];
                        if(pv2->is_updated) continue;
                        Q.push(pv2);
                    }
                }
            }
            clusters.push_back(cluster);
        }
    }
}

// 合并具有匹配关系的点集
void SnappingVertices::mergeMatchingVerticesClusters(vector<vector<Vertex *> > &clusters){
    // 重置顶点的更新状态
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < polygons.size(); ++i){
        pPoly = &polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pv->is_updated = false;
        }
    }

    for(int i = 0; i < clusters.size(); ++i){
        dealWithEachCluster(clusters[i]);
    }
}
// mergeMatchingVerticesClusters子函数
void SnappingVertices::dealWithEachCluster(vector<Vertex*> &cluster){
    // 计算参与的多边形数量
    set<int> polys_set;
    for(int i = 0; i < cluster.size(); ++i){
        int poly_id = floor(cluster[i]->point.data[3]+0.5f);
        polys_set.insert(poly_id);
    }

    // 关联的多边形数量不得超过3个也不能低于2个，出现此种情况要报错
    if(polys_set.size() > 3 || polys_set.size() < 2){
        cerr << "associative polys' number exceeds 3 or less than 2!" << endl;
        cerr << "polys_set.size() = " << polys_set.size() << endl;
        for(auto iter = polys_set.begin(); iter != polys_set.end(); ++iter){
            cout << *iter << "; ";
        }
        cout << endl;
        return;
    }

    // 关联的多边形有2个
    if(polys_set.size() == 2){
        // 提取多边形id
        int poly_a, poly_b;
        auto iter = polys_set.begin();
        poly_a = *iter;
        ++iter;
        poly_b = *iter;

        // 计算两个多边形的交线
        Eigen::Vector3f p_base, line_dir;
        getIntersectLine(polygons[poly_a].coeff, polygons[poly_b].coeff, p_base, line_dir);
        // 计算这些顶点在交线上的投影
        Vertex*pv;
        pcl::PointXYZ p_proj, p_proj_mean;
        p_proj_mean.x = p_proj_mean.y = p_proj_mean.z = 0;
        for(int i = 0; i < cluster.size(); ++i){
            pv = cluster[i];
            getProjPointOnLine(pv->point, p_proj, p_base, line_dir);
            p_proj_mean.x += p_proj.x;
            p_proj_mean.y += p_proj.y;
            p_proj_mean.z += p_proj.z;
        }
        // 计算投影点的重心
        p_proj_mean.x /= cluster.size();
        p_proj_mean.y /= cluster.size();
        p_proj_mean.z /= cluster.size();

        // 更新顶点位置
        for(int i = 0; i < cluster.size(); ++i){
            pv = cluster[i];
            if(distP2P(pv->point, p_proj_mean) <= pv->radius){
                float pv_poly_index = pv->point.data[3];    // pv所属多边形id
                pv->point = p_proj_mean;
                pv->point.data[3] = pv_poly_index;          // 重置pv的data3属性
                pv->is_updated = true;
            }
        }
        return;
    }

    // 关联的多边形有3个
    // 提取多边形id
    vector<int> poly_indices;
    for(auto iter = polys_set.begin(); iter != polys_set.end(); ++iter) poly_indices.push_back(*iter);
    pcl::PointXYZ p_intersect;
    // 计算交点
    getIntersectPointOfTriplePolys(polygons[poly_indices[0]].coeff,
            polygons[poly_indices[1]].coeff,
            polygons[poly_indices[2]].coeff,
            p_intersect);
    // 尝试更新位置
    bool are_all_dists_within_radius = true;
    Vertex *pv;
    for(int i = 0; i < cluster.size(); ++i){
        pv = cluster[i];
        if(distP2P(pv->point, p_intersect) > pv->radius) {
            are_all_dists_within_radius = false;
            break;
        }
    }
    if(are_all_dists_within_radius){
        for(int i = 0; i < cluster.size(); ++i) {
            pv = cluster[i];
            float pv_poly_index = pv->point.data[3];    // pv所属多边形id
            pv->point = p_intersect;
            pv->point.data[3] = pv_poly_index;          // 重置pv的data[3]属性
            pv->is_updated = true;
        }
        return;
    }

    // 运行到这里说明三个平面的交点不能替换全部的顶点
    // 计算三个多边形法向量的两两夹角
    float n01, n02, n12;
    n01 = getIncludedAngleOfTwoPolys(poly_indices[0], poly_indices[1]);
    n02 = getIncludedAngleOfTwoPolys(poly_indices[0], poly_indices[2]);
    n12 = getIncludedAngleOfTwoPolys(poly_indices[1], poly_indices[2]);
    int i0, i1, i2;// i0与i1保存最小夹角的多边形索引，i2保存与它们均成较大夹角的多边形索引
    float min_angle = FLT_MAX, angle1, angle2;// angle1与angle2保存非最小的两个夹角
    if(min_angle > n01){    // n01是最小的夹角
        min_angle = n01;
        angle1 = n02;
        angle2 = n12;
        i0 = poly_indices[0];
        i1 = poly_indices[1];
        i2 = poly_indices[2];
    }
    if(min_angle > n02){    // n02是最小的夹角
        min_angle = n02;
        angle1 = n01;
        angle2 = n12;
        i0 = poly_indices[0];
        i1 = poly_indices[2];
        i2 = poly_indices[1];
    }
    if(min_angle > n12){    // n12是最小的夹角
        min_angle = n12;
        angle1 = n01;
        angle2 = n02;
        i0 = poly_indices[1];
        i1 = poly_indices[2];
        i2 = poly_indices[0];
    }

    // 分别构造cluster1与cluster2
    vector<Vertex*> cluster1;   // i0的全部顶点+与i0关联更近的i2中的顶点
    vector<Vertex*> cluster2;   // i1的全部顶点+与i1关联更近的i2中的顶点
    for(int i = 0; i < cluster.size(); ++i){
        int poly_id = floor(cluster[i]->point.data[3]+0.5f);
        if(poly_id == i0)   // i0中的顶点
            cluster1.push_back(cluster[i]);
        else if(poly_id == i1)  // i1中的顶点
            cluster2.push_back(cluster[i]);
        else{                   // i2中的顶点
            // 寻找与当前顶点最近的多边形
            float min_dist = FLT_MAX;
            int index = -1; // 保存距离最近的多边形索引
            for(int j = 0; j < cluster[i]->maching_polys.size(); ++j){
                for(int k = 0; k < cluster[i]->maching_polys[j].vertices.size(); ++k){
                    pv = (Vertex*)(cluster[i]->maching_polys[j].vertices[k]);
                    if(min_dist > distP2P(cluster[i]->point, pv->point)){
                        min_dist = distP2P(cluster[i]->point, pv->point);
                        index = j;
                    }
                }
            }
            if(index = i0){
                cluster1.push_back(cluster[i]);
            }else{
                cluster2.push_back(cluster[i]);
            }
        }
    }// 构造cluster1与cluster2结束

    // 处理cluster1:
    dealWithEachCluster2(cluster1);
    // 处理cluster2:
    dealWithEachCluster2(cluster2);
}

// 仅处理分别属于两个多边形的顶点簇
void SnappingVertices::dealWithEachCluster2(vector<Vertex*> &cluster)
{
    // 计算参与的多边形数量
    set<int> polys_set;
    for(int i = 0; i < cluster.size(); ++i){
        int poly_id = floor(cluster[i]->point.data[3]+0.5f);
        polys_set.insert(poly_id);
    }
    if(polys_set.size() == 1){
        return;
    }
    if(polys_set.size() == 2){
        // 提取多边形id
        int poly_a, poly_b;
        auto iter = polys_set.begin();
        poly_a = *iter;
        ++iter;
        poly_b = *iter;

        // 计算两个多边形的交线
        Eigen::Vector3f p_base, line_dir;
        getIntersectLine(polygons[poly_a].coeff, polygons[poly_b].coeff, p_base, line_dir);
        // 计算这些顶点在交线上的投影
        Vertex*pv;
        pcl::PointXYZ p_proj, p_proj_mean;
        p_proj_mean.x = p_proj_mean.y = p_proj_mean.z = 0;
        for(int i = 0; i < cluster.size(); ++i){
            pv = cluster[i];
            getProjPointOnLine(pv->point, p_proj, p_base, line_dir);
            p_proj_mean.x += p_proj.x;
            p_proj_mean.y += p_proj.y;
            p_proj_mean.z += p_proj.z;
        }
        // 计算投影点的重心
        p_proj_mean.x /= cluster.size();
        p_proj_mean.y /= cluster.size();
        p_proj_mean.z /= cluster.size();

        // 更新顶点位置
        for(int i = 0; i < cluster.size(); ++i){
            pv = cluster[i];
            if(distP2P(pv->point, p_proj_mean) <= pv->radius){
                float pv_poly_index = pv->point.data[3];    // pv所属多边形id
                pv->point = p_proj_mean;
                pv->point.data[3] = pv_poly_index;          // 重置pv的data[3]属性
                pv->is_updated = true;
            }
        }
    }else{
        cerr << "SnappingVertices::dealWithEachCluster2(vector<Vertex*> &cluster):" << "error use!" << endl;
    }
}

// 计算三个多边形的交点
void SnappingVertices::getIntersectPointOfTriplePolys(Eigen::Vector4f &coeff1,
                                                      Eigen::Vector4f &coeff2,
                                                      Eigen::Vector4f &coeff3,
                                                      pcl::PointXYZ &p)
{
    Eigen::Matrix3f A;
    Eigen::Vector3f x, b;
    A << coeff1[0], coeff1[1], coeff1[2],
         coeff2[0], coeff2[1], coeff2[2],
         coeff3[0], coeff3[1], coeff3[2];
    b << -1*coeff1[3], -1 *coeff2[3], -1*coeff3[3];
    x = A.colPivHouseholderQr().solve(b);
    p.x = x[0];
    p.y = x[1];
    p.z = x[2];
}

// 计算两个多边形的交线
bool SnappingVertices::getIntersectLine(Eigen::Vector4f &coeff1, Eigen::Vector4f &coeff2,
                                        Eigen::Vector3f &p_base, Eigen::Vector3f &line_dir)
{
    Eigen::Vector3f n1, n2;
    n1 << coeff1[0], coeff1[1], coeff1[2];
    n2 << coeff2[0], coeff2[1], coeff2[2];

    // 两个多边形平行
    if(abs(n1.dot(n2)) >= 0.99999){
        p_base << 0, 0, 0;
        line_dir << 0, 0, 0;
        return false;
    }

    // 用两个法向量的差乘获取交线的方向向量
    line_dir = n1.cross(n2);
    line_dir.normalize();

    // 让系数绝对值最小的数对应的变量值置为0
    int index = 0;
    float min = coeff1[0];
    if(abs(min) > abs(coeff1[1]))
    {
        min = coeff1[1];
        index = 1;
    }
    if(abs(min) > abs(coeff1[2]))
    {
        min = coeff1[2];
        index = 2;
    }
    Eigen::Matrix2d A;
    Eigen::Vector2d x, b;
    p_base[index] = 0;
    switch(index)
    {
    case 0:
        A << coeff1[1], coeff1[2],
                coeff2[1], coeff2[2];
        b << coeff1[3], coeff2[3]; b *= -1;
        x = A.inverse() * b;
        p_base[1] = x[0];
        p_base[2] = x[1];
        break;
    case 1:
        A << coeff1[0], coeff1[2],
                coeff2[0], coeff2[2];
        b << coeff1[3], coeff2[3]; b *= -1;
        x = A.inverse() * b;
        p_base[0] = x[0];
        p_base[2] = x[1];
        break;
    case 2:
        A << coeff1[0], coeff1[1],
                coeff2[0], coeff2[1];
        b << coeff1[3], coeff2[3]; b *= -1;
        x = A.inverse() * b;
        p_base[0] = x[0];
        p_base[1] = x[1];
        break;
    }
    return true;
}

// 获取点在直线上的投影点
void SnappingVertices::getProjPointOnLine(pcl::PointXYZ &p, pcl::PointXYZ &p_proj,
                        Eigen::Vector3f &p_base, Eigen::Vector3f &line_dir){
    Eigen::Vector3f v_p, delta;
    v_p << p.x, p.y, p.z;
    delta = p_base - v_p;
    float lambda = -1.0 * line_dir.dot(delta);

    Eigen::Vector3f p_proj_v;
    p_proj_v = p_base + lambda * line_dir;
    p_proj.x = p_proj_v[0];
    p_proj.y = p_proj_v[1];
    p_proj.z = p_proj_v[2];
}

// 计算平面法向量的夹角
float SnappingVertices::getIncludedAngleOfTwoPolys(int poly_a, int poly_b){
    Eigen::Vector3f n1, n2;
    n1 << polygons[poly_a].coeff[0],
          polygons[poly_a].coeff[1],
          polygons[poly_a].coeff[2];
    n2 << polygons[poly_b].coeff[0],
          polygons[poly_b].coeff[1],
          polygons[poly_b].coeff[2];
    float cos_val = n1.dot(n2);
    float rad_angle = acos(cos_val);
    return rad_angle * 180.0 / 3.141592653589793;
}

// 判断两个多边形是否平行
bool SnappingVertices::areTwoPolysParallel(int poly_a, int poly_b){
    // 计算两个多边形法向量的夹角
    Eigen::Vector3f normal_a, normal_b;
    POLYGON *pPoly = &this->polygons[poly_a];
    normal_a << pPoly->coeff[0], pPoly->coeff[1], pPoly->coeff[2];
    pPoly = &this->polygons[poly_b];
    normal_b << pPoly->coeff[0], pPoly->coeff[1], pPoly->coeff[2];

    float cos_val = normal_a.dot(normal_b);

    cout << "angle:" << acos(cos_val)/3.1415927*180.0 << endl;

    cos_val = abs(cos_val);

    float T_cos_val = cos(T_maxAngle_bet_two_polys);

    return cos_val >= T_cos_val;
}

// 处理“点-边”匹配
void SnappingVertices::dealWithVertexEdgeMatchingRelation()
{
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < polygons.size(); ++i){
        pPoly = &polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;

            if(pv->neighbor_poly_id != -1 && !pv->is_vertex_neighbor){// “点-边匹配”
                // 交线L：p_base, line_dir
                Eigen::Vector3f p_base, line_dir;
                int poly_a, poly_b;
                poly_a = i;
                poly_b = pv->neighbor_poly_id;
                this->getIntersectLine(polygons[poly_a].coeff, polygons[poly_b].coeff, p_base, line_dir);
                // pv在L上的投影点
                pcl::PointXYZ p_proj;
                this->getProjPointOnLine(pv->point, p_proj, p_base, line_dir);
                // 判断pv到投影点的距离
                if(distP2P(pv->point, p_proj) <= pv->radius &&
                        distP2P(pv->point, p_proj) <= pv->neighbor_edge->radius){
                    float pv_poly_index = pv->point.data[3];    // pv所属多边形id
                    pv->point = p_proj; // 更新pv的位置
                    pv->point.data[3] = pv_poly_index;          // 重置pv的data[3]属性
                    pv->is_updated = true;
                    // 向对应边添加一个顶点
                    Vertex* new_v = new Vertex();
                    new_v->point = p_proj;
                    new_v->point.data[3] = poly_b;
                    new_v->radius = pv->neighbor_edge->radius;
                    new_v->neighbor_poly_id = -1;
                    new_v->neighbor_vertex = new_v->neighbor_edge = NULL;
                    if(!polygons[poly_b].insertVertex(pv->neighbor_edge, new_v)){
                        cerr << "顶点插入失败！" << endl;
                    }
                }
            }
        }
    }
}

// 处理多边形退化情况
// 主体函数
void SnappingVertices::dealWithPolysDegeneration(vector<POLYGON, Eigen::aligned_allocator<POLYGON>> &polygons)
{
    for(int i = 0; i < polygons.size(); ++i){
        dealWithSinglePolyDegeneration(polygons[i]);
    }
}
// 单个多边形的处理函数
void SnappingVertices::dealWithSinglePolyDegeneration(POLYGON &poly)
{
    // 构造顶点点云
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices (new pcl::PointCloud<pcl::PointXYZ>);
    Vertex *pv;
    vertices->reserve(poly.getSize());
    for(int i = 0; i < poly.getSize(); ++i){
        pv = i == 0 ? poly.start_point : pv->next_point;
        vertices->push_back(pv->point);
    }

	// 注意vertices中的顶点与poly中的顶点排序具有一对一的关系，在历次更新中这种关系也是得以保持的

    // 首先处理顶点冗余
    dealWithVerticesRedundancy(poly, vertices);
    // 然后处理线段折返
    dealWithLineSegCollineation(poly, vertices);
    // 最后处理不相邻线段的相交问题
    dealWithLineSegIntersect(poly, vertices);
    // 清理每个顶点的匹配关系
    cleanMatchingRelaion();
}

// 清理每个顶点的匹配关系
void SnappingVertices::cleanMatchingRelaion()
{
    for(int i = 0; i < polygons.size(); ++i){
        POLYGON *pPoly = &polygons[i];
        Vertex*pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pv->neighbor_poly_id = -1;
            pv->neighbor_edge = pv->neighbor_vertex = NULL;
            pv->maching_polys.clear();
            pv->is_intersectPoint = pv->is_updated = false;
        }
    }
}

// 子函数
// 计算顶点index1到index2经过路径上线段的长度之和
float SnappingVertices::lineSegDist(int index1, int index2, vector<float> &line_lens){
    float d = line_lens[index1];
    int next = index1 == line_lens.size()-1 ? 0 : index1+1;
    while(next != index2){
        next = next == line_lens.size()-1 ? 0 : next+1;
        d += line_lens[next];
    }
    return d;
}

// 处理顶点冗余(相邻顶点位置重合，不相邻顶点位置重合）
void SnappingVertices::dealWithVerticesRedundancy(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices)
{
    // 建立顶点集合的kdtree
    pcl::KdTreeFLANN<pcl::PointXYZ> tree;
    tree.setInputCloud(vertices);
    vector<int> indices;
    vector<float> dists;

    vector<vector<int>> clusters;   // 保存冗余顶点索引

    vector<bool> is_processed;
	is_processed.resize(vertices->size());
    for(auto iter = is_processed.begin(); iter != is_processed.end(); ++iter) *iter = false;

    // 利用kdtree找到冗余点簇
    for(int i = 0; i < vertices->size(); ++i){
        if(is_processed[i]) continue;

        tree.radiusSearch(vertices->points[i], 0.001f, indices, dists);
        if(indices.size() > 1){
            vector<int> cluster;
            cluster.reserve(indices.size());
            for(int j = 0; j < indices.size(); ++j){
                cluster.push_back(indices[j]);
                is_processed[indices[j]] = true;
            }
            clusters.push_back(cluster);
        }
    }

    // 每条线段长度与多边形周长
    vector<float> line_lens;
    float poly_circumference = 0;
    line_lens.reserve(vertices->size());
    for(int i = 0; i < vertices->size(); ++i){
        int next = i == vertices->size()-1 ? 0 : i+1;
        float d = distP2P(vertices->points[i], vertices->points[next]);
        line_lens.push_back(d);
        poly_circumference += d;
    }

    // 标记顶点是否需要删除
    vector<bool> is_del;
    is_del.resize(vertices->size());
    for(auto iter = is_del.begin(); iter != is_del.end(); ++iter) *iter = false;

    // 处理每个冗余点簇
    for(int i = 0; i < clusters.size(); ++i){
        vector<int>*pCluster = &clusters[i];
        // 升序排序
        std::sort(pCluster->begin(), pCluster->end());
        int j,k;// 注意j,k定义了点簇的边界，看笔记

        j = (*pCluster)[0]; // 第一个
        k = (*pCluster)[pCluster->size()-1];  // 最后一个

        for(int l = 0; l < pCluster->size(); ++l){
            int l_next = l == pCluster->size()-1 ? 0 : l+1;
            // 如果从l走向l_next较长，那么从l_next走向l就更短
            // 因而l_next作为起点；而l作为终点
            if(lineSegDist((*pCluster)[l], (*pCluster)[l_next], line_lens) > 0.5*poly_circumference){
                j = (*pCluster)[l_next];
                k = (*pCluster)[l];
                break;
            }
        }

        // 将顶点序列中j_next到k之间的顶点全部删掉
        int j_next = j == vertices->size()-1 ? 0 : j+1;
        for(int m = j_next; m != k; ){
            is_del[m] = true;
            m = m == vertices->size()-1 ? 0 : m+1;
        }
        is_del[k] = true;
    }

    // 基于is_del确定顶点的删除方案
    vector<int> off_set;
    off_set.resize(is_del.size());
    int del_count = 0;
    for(int i = 0; i < is_del.size(); ++i){
        if(is_del[i]){
            off_set[i] = -1;
            ++del_count;
        }else{
            off_set[i] = del_count;
        }
    }

    // 更新vertices
    for(int i = 0; i < is_del.size(); ++i){
        if(off_set[i] != -1){
            vertices->points[i-off_set[i]] = vertices->points[i];
        }
    }

	for(int i = 0; i < del_count; ++i){ vertices->points.pop_back(); }

    // 更新poly
    Vertex *pv;
    for(int i = 0; i < is_del.size(); ++i){
        pv = i == 0 ? poly.start_point : pv->next_point;

        // 若要删除当前顶点，需要保存当前顶点的前驱节点
        if(is_del[i]){
            Vertex *pv_prev = pv->prev_point;
            poly.delVertex(pv);
            pv = pv_prev;
        }
    }
}
// 处理相邻线段共线
void SnappingVertices::dealWithLineSegCollineation(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices)
{
    vector<bool> is_del; is_del.resize(vertices->size(), false);

    // 当前顶点"前驱线段"与"后继线段"的方向向量
    Eigen::Vector3f prev_line_dir, next_line_dir;

    // 遍历顶点
    for(int i = 0; i < vertices->size(); ++i){
        int prev = i == 0 ? vertices->size()-1 : i-1;
        int next = i == vertices->size()-1 ? 0 : i+1;

        prev_line_dir << vertices->points[i].x - vertices->points[prev].x,
                         vertices->points[i].y - vertices->points[prev].y,
                         vertices->points[i].z - vertices->points[prev].z;

        next_line_dir << vertices->points[next].x - vertices->points[i].x,
                         vertices->points[next].y - vertices->points[i].y,
                         vertices->points[next].z - vertices->points[i].z;

        prev_line_dir.normalize();
        next_line_dir.normalize();

        float cos_val = prev_line_dir.dot(next_line_dir);
        float abs_cos_val = abs(cos_val);

        if(abs(abs_cos_val - 1.0) < 0.000001f){
            is_del[i] = true;
        }
    }

    // 基于is_del确定顶点的删除方案
    vector<int> off_set;
    off_set.resize(is_del.size());
    int del_count = 0;
    for(int i = 0; i < is_del.size(); ++i){
        if(is_del[i]){
            off_set[i] = -1;
            ++del_count;
        }else{
            off_set[i] = del_count;
        }
    }

    // 更新vertices
    for(int i = 0; i < is_del.size(); ++i){
        if(off_set[i] != -1){
            vertices->points[i-off_set[i]] = vertices->points[i];
        }
    }

    for(int i = 0; i < del_count; ++i){ vertices->points.pop_back(); }

    // 更新poly
    Vertex *pv;
    for(int i = 0; i < is_del.size(); ++i){
        pv = i == 0 ? poly.start_point : pv->next_point;

        // 若要删除当前顶点，需要保存当前顶点的前驱节点
        if(is_del[i]){
            Vertex *pv_prev = pv->prev_point;
            poly.delVertex(pv);
            pv = pv_prev;
        }
    }
}

// 将交点v插入到多边形第edge_id个顶点的后面中
void SnappingVertices::insertIPVertexIntoPoly(Vertex *v, POLYGON &poly,const int edge_id, vector<Vertex*> &poly_vertices){
	// 目标线段的起始顶点
	Vertex *begin_v = poly_vertices[edge_id];
	// 定位到目标顶点以及目标顶点的下一个顶点
	Vertex *cur_v = poly_vertices[edge_id];
	Vertex *next_v = cur_v->next_point;
	float d1, d2;
	// 计算被插入顶点与目标线段起点的距离
	d1 = this->distP2P(v->point, begin_v->point);
	do{
		// 计算下一个顶点距离目标线段起点的距离
		d2 = this->distP2P(next_v->point, begin_v->point);
		// v若要能插入到cur_v与next_v之间，必须有d1 < d2
		if(d1 < d2){
			poly.insertVertex(cur_v, v);
			return;
		}else{
			// 如果d1>d2，意味着被插入顶点需要向更远的位置处移动才行
			cur_v = cur_v->next_point;
			next_v = cur_v->next_point;
		}
	}while(true);
}

// 处理非近邻线段的相交
void SnappingVertices::dealWithLineSegIntersect(POLYGON &poly, pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices)
{
    // 第一步：计算所有不相邻线段的交点
	vector<IntersectingPoint, Eigen::aligned_allocator<IntersectingPoint>> intersect_points;	// 用来保存交点信息
	getIntersectPointsOfPoly(vertices, intersect_points);
    if(intersect_points.size() == 0) return;	// 如果没有交点则直接返回

	// 第二步：将交点插入到多边形中
	// 设置当前顶点的"是否为交点"属性为false
	Vertex *pv;
	for(int i = 0; i < poly.getSize(); ++i){
		pv = i == 0 ? poly.start_point : pv->next_point;
		pv->is_intersectPoint = false;
	}

	// 将所有的顶点顺序存储
	vector<Vertex*> poly_vertices;
	for(int i = 0; i < poly.getSize(); ++i){
		pv = i == 0 ? poly.start_point : pv->next_point;
		poly_vertices.push_back(pv);
	}

	// 将交点插入多边形中
	for(int i = 0; i < intersect_points.size(); ++i){
		int edge_id1 = intersect_points[i].edge_id1;
		int edge_id2 = intersect_points[i].edge_id2;
		// 先处理第一条线段
		Vertex *v1 = new Vertex();
		v1->is_intersectPoint = true;
		v1->point = intersect_points[i].point;
        v1->point.data[3] = poly.start_point->point.data[3];
		insertIPVertexIntoPoly(v1, poly, edge_id1, poly_vertices);
		// 再处理第二条线段
		Vertex *v2 = new Vertex();
		v2->is_intersectPoint = true;
		v2->point = intersect_points[i].point;
        v2->point.data[3] = poly.start_point->point.data[3];
		insertIPVertexIntoPoly(v2, poly, edge_id2, poly_vertices);
	}

	// 用多边形中顶点更新vertices
	vertices->clear();
	for(int i = 0; i < poly.getSize(); ++i){
		pv = i == 0 ? poly.start_point : pv->next_point;
		vertices->push_back(pv->point);
		if(pv->is_intersectPoint){
			// 用point中的第四个元素表明该点是否为交点，最大值是，最小值否
			vertices->points[vertices->size()-1].data[3] = FLT_MAX;
		}else{
			vertices->points[vertices->size()-1].data[3] = FLT_MIN;
		}
	}

	// 第三步：将所有的交点集中到一个集合中
	// ip_indices: 保存所有交点的索引
	vector<int> ip_indices;	 // ip: abbv for intersect_point 
	for(int i = 0; i < vertices->size(); ++i){
		if(vertices->points[i].data[3] == FLT_MAX){
			ip_indices.push_back(i);
		}
	}
	// 用来标记该交点是否已被处理过
	vector<bool> is_processed;
	is_processed.resize(ip_indices.size(), false);

	// 第四步：创建新的子多边形
	vector<vector<int>> sub_polys;	// 用来保存新创建的子多边形，索引对应vertices

	// 从每个交点开始遍历
	for(int i = 0; i < ip_indices.size(); ++i){
		if(is_processed[i]) continue;	// 如果当前交点已被处理过，则跳过

		is_processed[i] = true;			// 设置当前交点已被处理过

		int ip_index = ip_indices[i];	// 获取当前交点的索引	

		vector<int> sub_poly;			// 准备创建多边形
		sub_poly.push_back(ip_index);	// 将当前交点作为新多边形的起点

		int next = ip_index == vertices->size()-1 ? 0 : ip_index+1;	// 当前交点在原多边形中的下一个顶点的索引
		do{
			// 普通顶点
			if(vertices->points[next].data[3] == FLT_MIN){
				// 直接加入多边形
				sub_poly.push_back(next);
			}else if(vertices->points[next].data[3] == FLT_MAX &&
				    !isEqual(vertices->points[ip_index], vertices->points[next])){
				// 非ip_index的对偶交点
				// 首先找到next的对偶交点dual_next，就像是火车变轨一样
				int dual_next;
				for(int j = 0; j < ip_indices.size(); ++j){
					if(is_processed[j]) continue;
					if(ip_indices[j] == next) continue;
					if(isEqual(vertices->points[next], vertices->points[ip_indices[j]])){// 找到了对偶点
						dual_next = ip_indices[j];
						break;
					}
				} // forj

				// 然后将dual_next加入多边形
				sub_poly.push_back(dual_next);

				// 接下来标记dual_next已被处理过
				auto iter1 = std::find(ip_indices.begin(), ip_indices.end(), dual_next);
				int dual_next_index = iter1 - ip_indices.begin();
				is_processed[dual_next_index] = true;

				// 最后更新next
				next = dual_next;
			}else{
				// ip_index的对偶顶点
				sub_polys.push_back(sub_poly);
				break;	// 构造完毕一个多边形，退出while循环
			}
			// 指向下一个顶点
			next = next == vertices->size()-1 ? 0 : next+1;
		}while(true);
	}// for i

	// 第五步 计算每个多边形的周长
	vector<float> poly_c;
	float circumference;
	for(int i = 0; i < sub_polys.size(); ++i){
		pcl::PointCloud<pcl::PointXYZ>::Ptr poly_border (new pcl::PointCloud<pcl::PointXYZ>);
		for(int j = 0; j < sub_polys[i].size(); ++j){
			poly_border->push_back(vertices->points[sub_polys[i][j]/*sub_polys中第i个多边形中第j个顶点在vertices中的索引*/]);
		}

		// 计算poly_border的周长
		circumference = this->getCircumferenceOfPoly(poly_border);
		poly_c.push_back(circumference);
	}

	// 第六步：在poly_c/sub_polys中寻找周长最长的子多边形的索引
	float max_c = FLT_MIN;
	int poly_index = -1;
	for(int i = 0; i < poly_c.size(); ++i){
		if(poly_c[i] > max_c){
			max_c = poly_c[i];
			poly_index = i;
		}
	}

	// 第七步：设置周长最长的子多边形中每个顶点的删除状态为保留，而其它顶点则予以删除
	vector<bool> is_del;
	is_del.resize(vertices->size(), true);
	for(int i = 0; i < sub_polys[poly_index].size(); ++i){
		is_del[sub_polys[poly_index][i]] = false;
	}

	// 更新poly_vertices;
	poly_vertices.clear();
	for(int i = 0; i < poly.getSize(); ++i){
		pv = i == 0 ? poly.start_point : pv->next_point;
		poly_vertices.push_back(pv);
	}

    // 根据顶点索引利用poly_vertices找到顶点位置
	for(int i = 0; i < is_del.size(); ++i){
		if(!is_del[i]) continue;
		poly.delVertex(poly_vertices[i]);
	}

	// 再次用多边形中顶点更新vertices
	vertices->clear();
	for(int i = 0; i < poly.getSize(); ++i){
		pv = i == 0 ? poly.start_point : pv->next_point;
		vertices->push_back(pv->point);
	}
}

// 根据编号判断两条线段是否相连
// 即其中一条线段是否是另一条线段的前驱线段或后继线段
bool SnappingVertices::isConnectTwoLineSegs(int size, int edge_id1, int edge_id2)
{
    int prev = edge_id1 == 0 ? size-1 : edge_id1-1;
    int next = edge_id1 == size-1 ? 0 : edge_id1+1;

    return (prev == edge_id2) || (next == edge_id2);
}

// 计算一个多边形中所有不相邻线段的交点
void SnappingVertices::getIntersectPointsOfPoly(pcl::PointCloud<pcl::PointXYZ>::Ptr &vertices,
                              vector<IntersectingPoint, Eigen::aligned_allocator<IntersectingPoint> > &intersect_points)
{
    if(intersect_points.size() > 0) intersect_points.clear();

    int size = vertices->size();

    pcl::PointXYZ intersect_point;

    for(int i = 0; i < size-1; ++i){
        for(int j = i+1; j < size; ++j){
            if(isConnectTwoLineSegs(size, i, j)) continue;

            int i_next = i == size-1 ? 0 : i+1;
            int j_next = j == size-1 ? 0 : j+1;

            if(getIntersectPointOfTwoLineSeg(vertices->points[i], vertices->points[i_next],
                                             vertices->points[j], vertices->points[j_next],
                                             intersect_point)){
                IntersectingPoint ip;
                ip.edge_id1 = i;
                ip.edge_id2 = j;
                ip.point = intersect_point;
                intersect_points.push_back(ip);
            }
        }
    }
}

// 计算两条线段的交点
// 只有交点均落在两条线段内部，才会返回合法的值
bool SnappingVertices::getIntersectPointOfTwoLineSeg(pcl::PointXYZ &line1_begin, pcl::PointXYZ &line1_end,
                                   pcl::PointXYZ &line2_begin, pcl::PointXYZ &line2_end,
                                   pcl::PointXYZ &intersect_point)
{
    Eigen::Vector3f p1, p2, p3, p4, p12, p34, p31;
    p1 << line1_begin.x, line1_begin.y, line1_begin.z;
    p2 << line1_end.x, line1_end.y, line1_end.z;
    p3 << line2_begin.x, line2_begin.y, line2_begin.z;
    p4 << line2_end.x, line2_end.y, line2_end.z;
    p12 = p2-p1;
    p34 = p4-p3;
    p31 = p1-p3;

    Eigen::Matrix2f A;
    Eigen::Vector2f x, b;
    A << p12.dot(p12), -1*p12.dot(p34),
         -1*p12.dot(p34), p34.dot(p34);
    b << -1*p31.dot(p12), p31.dot(p34);
    x = A.colPivHouseholderQr().solve(b);
    float alpha1, alpha2;
    alpha1 = x[0];
    alpha2 = x[1];

    Eigen::Vector3f p_l1, p_l2;
    p_l1 = p1 + alpha1*p12;
    p_l2 = p3 + alpha2*p34;

    float dist_p2p = (p_l1-p_l2).norm();

    // alpha1在0，1之间意味着离着直线2最近的点在直线1的内部，同理对alpha2
    // 两点的距离充分小意味着两条线段确实是相交的
    if(alpha1 > 0 && alpha1 < 1 &&
       alpha2 > 0 && alpha2 < 1 &&
       dist_p2p < 0.0005){
        intersect_point.x = p_l1[0];
        intersect_point.y = p_l1[1];
        intersect_point.z = p_l1[2];
        return true;
    }else{
        intersect_point.x = intersect_point.y = intersect_point.z = FLT_MAX;
        return false;
    }
}

// 多边形周长函数
float SnappingVertices::getCircumferenceOfPoly(pcl::PointCloud<pcl::PointXYZ>::Ptr &border)
{
	float circumference = 0;
	for(int i = 0; i < border->size(); ++i){
		int next = i == border->size()-1 ? 0 : i+1;
		circumference += this->distP2P(border->points[i], border->points[next]);
	}
	return circumference;
}

// 设置顶点的搜索半径
void SnappingVertices::setSearchRadiusForEachVertex(vector<InitialPoly> &initial_polys)
{
    POLYGON *pPoly;
    Vertex *pv;
    for(int i = 0; i < polygons.size(); ++i){
        pPoly = &polygons[i];
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;
            pv->radius = this->ratio_of_scale * initial_polys[i].scale;
        }
    }
}

// 一步执行函数
// 参数设置以及多边形数据需提前备好
//第一步：
//setParams
//第二步：
// 用新的数据结构代表多边形
//this->m_sv.setPolygons(this->initial_polys);
void SnappingVertices::perform(vector<InitialPoly> &initial_polys)
{
    // 更新顶点的搜索半径
    setSearchRadiusForEachVertex(initial_polys);
    // 为多边形的每个顶点设置邻域关系
    setNeighborsForEachPoly();
    // 将单向的匹配关系转换为双向的匹配关系
    setMatchingRelations();
    // 获取具有"点-点"匹配联通关系的顶点簇
    getMatchingVerticesClusters(clusters);
    // 精化匹配关系(不添加点-点匹配，但可能会减少点-边匹配)
    refineMachingRelations();
    // 合并具有"点-点"匹配关系的簇
    mergeMatchingVerticesClusters(clusters);
    // 处理点-边匹配
    dealWithVertexEdgeMatchingRelation();
    // 处理多边形退化情况
    dealWithPolysDegeneration(polygons);
}

// 在具有角对边关系的线段上插入顶点
void SnappingVertices::insertVertexOnEdgeWithAngleToEdgeRelation()
{
	Vertex *pv;
	for(int i = 0; i < polygons.size(); ++i){
		for(int j = 0; j < polygons[i].getSize(); ++j){
			pv = j == 0 ? polygons[i].start_point : pv->next_point;
			// 检测在其它多边形中是否存在穿过点pv的线段，pv需在线段内部
			// 如果有，则在对应线段内部插入与pv在同一个位置的顶点
			checkLineSegThatPassThroughVertex(i, pv);
		}
	}
}

void SnappingVertices::checkLineSegThatPassThroughVertex(int poly_id, Vertex *pv)
{
	Vertex *edge_begin;
	for(int i = 0; i < polygons.size(); ++i){
		if(i == poly_id) continue;
		for(int j = 0; j < polygons[i].getSize(); ++j){
			edge_begin = j == 0 ? polygons[i].start_point : edge_begin->next_point;

			// 检测pv是否在线段(edge_begin, edge_end)的内部
			if(isVertexOnEdge(pv, edge_begin)){
				Vertex * v_new = new Vertex();
				v_new->point = pv->point;
				v_new->point.data[3] = edge_begin->point.data[3];
				v_new->radius = edge_begin->radius;
				polygons[i].insertVertex(edge_begin, v_new);
			}
		}
	}
}

bool SnappingVertices::isVertexOnEdge(Vertex *pv, Vertex *edge_begin)
{
	Eigen::Vector3f p_base, line_dir;
	p_base << edge_begin->point.x, edge_begin->point.y, edge_begin->point.z;
    line_dir << edge_begin->next_point->point.x - edge_begin->point.x,
				edge_begin->next_point->point.y - edge_begin->point.y,
                edge_begin->next_point->point.z - edge_begin->point.z;
    line_dir.normalize();

    Eigen::Vector3f v_p, delta;
	v_p << pv->point.x, pv->point.y, pv->point.z;
    delta = p_base - v_p;
    float lambda = -1.0 * line_dir.dot(delta);

    Eigen::Vector3f p_proj;
    p_proj = p_base + lambda * line_dir;

	float dist = (v_p - p_proj).norm();

	float epsilon = 0.0001;

	if(dist > epsilon) return false;

	float dist1 = this->distP2P(pv->point, edge_begin->point);
	float dist2 = this->distP2P(pv->point, edge_begin->next_point->point);

	if(dist1 <= epsilon*10) return false;				// 落在端点上
	if(dist2 <= epsilon*10) return false;	// 落在端点上

	float dist3 = this->distP2P(edge_begin->point, edge_begin->next_point->point);

	if(dist1+dist2 - dist3 > epsilon*9) return false;	// 落在线段外部
	else return true;						// 落在线段内部
}

// 检查当前模型是否合法
// 如果检测到非法情况，则返回问题顶点
bool SnappingVertices::is_model_valid(Vertex *&prob_v)
{
    // 首先对多边形执行三角化
    EarClip earClip;

    for(int i = 0; i < polygons.size(); ++i){
        earClip.setInputPoly(&(polygons[i]));
        if(!earClip.triangulatePoly()){
            cerr << "poly_id = " << i << ", triangulation failed!" << endl;
        }
    }

    for(int i = 0; i < polygons.size(); ++i){
        POLYGON *pPoly = &polygons[i];

        Vertex* pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;

            // pv 是否落在其它多边形的边上（包括边的起点和终点）
            if(is_vertex_on_edge(i, pv)) continue;
            // pv 是否落在其它多边形内部（不包括落在边界上）
            if(is_vertex_in_poly(i, pv)){
                prob_v = pv;
                return false;
            }
            // [pv pv->next_point]是否穿过其它多边形的内部
            if(is_edge_pass_through_poly(i, pv, pv->next_point)){
                prob_v = pv;
                return false;
            }
        }
    }
    return true;
}

// 编号为poly_id的多边形的顶点v是否落在其它多边形的边上
bool SnappingVertices::is_vertex_on_edge(int poly_id, Vertex* v)
{
    float epsilon = 0.0001;
    for(int i = 0; i < polygons.size(); ++i){
        if(i == poly_id) continue;
        POLYGON *pPoly = &polygons[i];

        Vertex* pv;
        for(int j = 0; j < pPoly->getSize(); ++j){
            pv = j == 0 ? pPoly->start_point : pv->next_point;

            // 点-点
            if(distP2P(v->point, pv->point) <= epsilon)
                return true;

            // 点-边
            if(distP2E(v->point, pv->point, pv->next_point->point) <= epsilon)
                return true;
        }
    }
    return false;
}

// 编号为poly_id的多边形的顶点v是否落在其它多边形的内部
bool SnappingVertices::is_vertex_in_poly(int poly_id, Vertex* v)
{
    float epsilon = 0.1;
    for(int i = 0; i < polygons.size(); ++i){
        if(i == poly_id) continue;
        POLYGON *pPoly = &polygons[i];

        // 首先计算点与平面的符号距离
        float dist = v->point.x*pPoly->coeff[0] +
                     v->point.y*pPoly->coeff[1] +
                     v->point.z*pPoly->coeff[2] +
                                pPoly->coeff[3];
        if(abs(dist) > epsilon) continue;

        // 计算点到平面的投影点
        pcl::PointXYZ proj_p;
        proj_p.x = v->point.x - dist * pPoly->coeff[0];
        proj_p.y = v->point.y - dist * pPoly->coeff[1];
        proj_p.z = v->point.z - dist * pPoly->coeff[2];

        // 判断投影点是否在多边形内部
        for(int j = 0; j < pPoly->triangles.size(); ++j){
            if(is_vertex_in_triangle(v->point, pPoly->triangles[j])){
                return true;
            }
        }
    }
    return false;
}

// 某点是否在三角形内部
bool SnappingVertices::is_vertex_in_triangle(pcl::PointXYZ &point, Triangle &tri)
{
    Eigen::Vector3f e01, e02, e0, e1, e2;
    e01 << tri.v[1]->point.x - tri.v[0]->point.x,
           tri.v[1]->point.y - tri.v[0]->point.y,
           tri.v[1]->point.z - tri.v[0]->point.z;
    e02 << tri.v[2]->point.x - tri.v[0]->point.x,
           tri.v[2]->point.y - tri.v[0]->point.y,
           tri.v[2]->point.z - tri.v[0]->point.z;
    float S = (e01.cross(e02)).norm() * 0.5;
    S = abs(S);

    e0 << tri.v[0]->point.x - point.x,
          tri.v[0]->point.y - point.y,
          tri.v[0]->point.z - point.z;
    e1 << tri.v[1]->point.x - point.x,
          tri.v[1]->point.y - point.y,
          tri.v[1]->point.z - point.z;
    e2 << tri.v[2]->point.x - point.x,
          tri.v[2]->point.y - point.y,
          tri.v[2]->point.z - point.z;
    float S1 = e0.cross(e1).norm() * 0.5;
    float S2 = e0.cross(e2).norm() * 0.5;
    float S3 = e1.cross(e2).norm() * 0.5;
    S1 = abs(S1);
    S2 = abs(S2);
    S3 = abs(S3);

    return abs(S1+S2+S3-S) <= 0.00001;
}

// 线段[e1, e2]是否穿过其它多边形的内部
bool SnappingVertices::is_edge_pass_through_poly(int poly_id, Vertex* e1, Vertex* e2)
{
    //float epsilon = 0.0001;
    for(int i = 0; i < polygons.size(); ++i){
        if(i == poly_id) continue;
        POLYGON *pPoly = &polygons[i];

        // e1到pPoly支撑平面的符号距离
        float d1 = e1->point.x * pPoly->coeff[0] +
                   e1->point.y * pPoly->coeff[1] +
                   e1->point.z * pPoly->coeff[2] +
                                 pPoly->coeff[3];
        float d2 = e2->point.x * pPoly->coeff[0] +
                   e2->point.y * pPoly->coeff[1] +
                   e2->point.z * pPoly->coeff[2] +
                                 pPoly->coeff[3];
		if(d1*d2 >= 0) return false;	// 位于平面的同一方向，不可能相交
        if(d1*d2 >= -0.00001) return false;	// 有一个顶点位于支撑平面上(或其附近)，亦不可能相交

        // 计算线段与多边形支撑平面的交点
        double A, B, lambda;
        A = pPoly->coeff[0] * (e2->point.x - e1->point.x) +
            pPoly->coeff[1] * (e2->point.y - e1->point.y) +
            pPoly->coeff[2] * (e2->point.z - e1->point.z);
        B = pPoly->coeff[0] * e1->point.x +
            pPoly->coeff[1] * e1->point.y +
            pPoly->coeff[2] * e1->point.z +
            pPoly->coeff[3];
        lambda = B/A*-1;
        pcl::PointXYZ p_intersect;
        p_intersect.x = e1->point.x + lambda * (e2->point.x - e1->point.x);
        p_intersect.y = e1->point.y + lambda * (e2->point.y - e1->point.y);
        p_intersect.z = e1->point.z + lambda * (e2->point.z - e1->point.z);

        // 判断交点是否在多边形内部
        for(int j = 0; j < pPoly->triangles.size(); ++j){
            if(is_vertex_in_triangle(p_intersect, pPoly->triangles[j])){
                return true;
            }
        }
    }

    return false;
}
