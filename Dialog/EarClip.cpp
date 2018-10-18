#include "earclip.h"

EarClip::EarClip()
{
    this->m_poly = NULL;
    PI = 3.14159265358979;
}

// 设置输入的多边形
void EarClip::setInputPoly(POLYGON *pPoly)
{
    this->vertices.clear();
    this->m_poly = pPoly;
    Vertex *pv;
    for(int i = 0; i < pPoly->getSize(); ++i){
        pv = i == 0 ? pPoly->start_point : pv->next_point;
        this->vertices.push_back(pv);
    }
    if(vertices.size() != pPoly->getSize()){
        cerr << "in EarClip::setInputPoly(POLYGON *pPoly) vertices.size() != pPoly->getSize()" << endl;
        cerr << "vertices.size() = " << vertices.size() << endl;
        cerr << "pPoly->getSize() = " << pPoly->getSize() << endl;
    }
}

// 对多边形进行三角化
// 前提：平面的法向量指向物体外部
//		平面边界点需要符合右手法则：大拇指指向平面的外侧（这样line(p0,p1)左边的点在多边形内部，右边的点在多边形外部）
bool EarClip::triangulatePoly()
{
    if(!m_poly){
        cerr << "in EarClip::triangulatePoly(), m_poly == NULL" << endl;
        return false;
    }

    // 多边形顶点至少三个
    if(m_poly->getSize() < 3){
        cerr << "in EarClip::triangulatePoly(), m_poly->getSize() = " << m_poly->getSize() << endl;
        return false;
    }

    // 清理多边形的三角形数据
    m_poly->triangles.clear();
    int poly_size = m_poly->getSize();
    // 三角形数量 = 多边形顶点数量 - 2
	m_poly->triangles.reserve(poly_size-2);

    // 处理特殊情况，如果只有三个点，则直接构造三角形，然后退出
    if(poly_size == 3){
        Triangle triangle;
        triangle.v[0] = vertices[0];
        triangle.v[1] = vertices[1];
        triangle.v[2] = vertices[2];

        m_poly->triangles.push_back(triangle);
        return true;
    }

    // 对于顶点数量超过3的一般情况
    // 首先保存多边形的顶点索引
    // indices: 当前的顶点列表
    std::vector<int> indices;
    indices.resize(poly_size);
    int i__ = 0;
    for(auto iter = indices.begin(); iter != indices.end(); ++iter){*iter = i__++;}

    // 然后标记每个顶点在合法性检查中是否需要进行检查（只有内角为PI的顶点才不进行此类检查）
    // 即，判断该点是否在被检测的三角形内部
    std::vector<bool> is_checked;
    is_checked.resize(poly_size);
    for(auto iter = is_checked.begin(); iter != is_checked.end(); ++iter){*iter = true;}    // 默认每个点都需要检查
    // 存储每个顶点的内角
    std::vector<double> inner_angles;
    inner_angles.resize(poly_size);
    for(auto iter = inner_angles.begin(); iter != inner_angles.end(); ++iter){*iter = -1;}
    // 标记该顶点是否已经被砍掉，只有没被砍掉的顶点才可以计算其内角
    std::vector<bool> is_cut;
    is_cut.resize(poly_size);
    for(auto iter = is_cut.begin(); iter != is_cut.end(); ++iter){*iter = false;}   // 一开始都没有被砍掉

    // 在每次迭代中，总是首先处理角度最大的内角
    // 如果内角超过了PI，则直接忽略(因为它不能被砍掉从而构建一个新的三角形)
    // 如果内角为PI，则构建退化三角形，并标记平角顶点为不被检查的顶点
    // 如果内角小于PI，则需进一步检查其内部是否有其它顶点，如果没有其它顶点，则可构造三角形
    // 每构建一个新的三角形，则重新计算受到影响的顶点的内角，再将这些内角从大到小排列
    int vertices_size = poly_size;   // vertices_size: 原始顶点的数量
    Eigen::Vector3f normal;   // normal: 多边形的法向量
    normal << m_poly->coeff[0],
              m_poly->coeff[1],
              m_poly->coeff[2];
    double eps = 0.000001;

    int while_count = 0;  // while_count: 统计while循环的次数，从而规避无限循环的情况
    while(true){
        if(while_count ++ >= 100000){
            cerr << "可能出现了死循环" << endl;
            return false;
        }
        // 如果上一次仅剩下三个顶点，那么切掉一个顶点之后只剩下两个顶点，此时所有的三角形均被获取，可退出该函数
        if(indices.size() < 3){
            return true;
        }

        // 计算顶点内角
        // i: 当前顶点的索引
        // prev: 前驱顶点的索引
        // next: 后继顶点的索引
        for(int i = 0; i < vertices_size; ++i){
            if(is_cut[i]){  // 如果顶点已被切掉，则将其对应的内角赋予特殊值
                inner_angles[i] = -1;
            }else{
                // 找到i在indices中的位置index
                int index = -1;
                for(index = 0; index < indices.size(); ++index){
                    if(indices[index] == i){
                        break;
                    }
                }
                if(index == -1){std::cerr << "index = -1" << endl; return false;}
                // 找到i的前驱与后继
                int prev, next;
                prev = indices[getPrev(indices, index)];
                next = indices[getNext(indices, index)];

                inner_angles[i] = getInnerAngle(vertices, normal, i, prev, next);
            }
        }

        // 寻找内角为PI的顶点
        bool is_find_PI = false;
        for(int i = 0; i < vertices_size; ++i){
            if(abs(inner_angles[i]-PI) < eps){  // 找到平角顶点
                // 找到i在indices中的位置
                int index = -1;
                for(index = 0; index < indices.size(); ++index){
                    if(indices[index] == i){
                        break;
                    }
                }
                if(indices[index]-i != 0) cerr << "indices[index]-i != 0\n" << "i = " << i << "\nindices[index]=" << indices[index] << endl;
                // 找到i的前驱与后继
                int prev, next;
                prev = indices[getPrev(indices, index)];
                next = indices[getNext(indices, index)];

                // 构建三角形
                Triangle triangle;
//                triangle.a = i;
//                triangle.b = prev;
//                triangle.c = next;
                triangle.v[0] = vertices[i];
                triangle.v[1] = vertices[prev];
                triangle.v[2] = vertices[next];
                m_poly->triangles.push_back(triangle);
                del_vertex(indices, index); // 将该顶点从顶点序列中删除
                is_cut[i] = true;           // 切掉该顶点
                is_checked[i] = false;      // 因为是平角，所以消除该点被检测的权利
                is_find_PI = true;          // 标识
                break;
            }
        }

        /*/ test code2:
            cout << "is_find_PI: " << is_find_PI << endl;*/

        if(is_find_PI){
            continue;
        }

        // 在正常角度中寻找最大的合法内角（当然不能超过PI）
        double max_angle = -1;// 最大合法角度
        int max_i = -1;       // 上述角度对应的顶点索引(针对原始的顶点序列)
        int max_prev, max_next;
        for(int i = 0; i < vertices_size; ++i){
            // 首先排除特殊情况，不然会检测已经被删除的顶点，而这些顶点早已不存在于indices中
            // ＞PI的顶点不考虑
            if(inner_angles[i] > PI) continue;
            // 已经被删除的顶点不考虑
            if(is_cut[i]) continue;

            // 寻找该顶点的前驱后继结点
            int index = -1;
            for(index = 0; index < indices.size(); ++index){
                if(indices[index] == i){
                    break;
                }
            }
            if(indices[index]-i != 0) cerr << "indices[index]-i != 0\n" << "i = " << i << "\nindices[index]=" << indices[index] << endl;

            int prev, next;
            prev = indices[getPrev(indices, index)];
            next = indices[getNext(indices, index)];

            // 可能的三角形：<i,prev,next>
            // 该三角形内部是否存在其它顶点
            bool is_other_point_in_triangle = false;
            for(int j = 0; j < vertices_size; ++j){
                if(j == i || j == prev || j == next || !is_checked[j])
                    continue;

                if(isPointInTriangle(j, i, prev, next)){
                    is_other_point_in_triangle = true;
                    break;
                }
            }
            // 该三角形内部没有其它顶点且最大角度小于当前内角，则更新最大角度
            if(!is_other_point_in_triangle && max_angle < inner_angles[i]){
                max_angle = inner_angles[i];
                max_i = i;
                max_prev = prev;
                max_next = next;
            }

        }
        if(max_angle == -1){
            cerr << "max_angle == -1 in trianglatePlane" << endl;
            return false;
        }

        /*/ test code3:
            cout << "max_angle = " << max_angle << endl;
            cout << "max_i = " << max_i << "; max_prev = " << max_prev << "; max_next = " << max_next << endl;
            */

        // 构建三角形
        Triangle triangle;
//        triangle.a = max_i;
//        triangle.b = max_prev;
//        triangle.c = max_next;
        triangle.v[0] = vertices[max_i];
        triangle.v[1] = vertices[max_prev];
        triangle.v[2] = vertices[max_next];
        m_poly->triangles.push_back(triangle);

        is_cut[max_i] = true;   // 标记该顶点已被切掉

        // 从indices中删除该顶点
        // 找到max_i对应的indices中的索引.
        int i_ = -1;
        for(i_ = 0; i_ < indices.size(); ++i_){
            if(indices[i_] == max_i){
                break;
            }
        }
        if(i_ == indices.size()){
            cerr << "i_ == -1" << endl;
            return false;
        }
        /*/ test code4:
            cout << "indices.size = " << indices.size() << endl;
            cout << "i_ = " << i_ << endl;
            cout << "indices[i_] = " << indices[i_] << endl;
            cout << "Sleep(10000)" << endl;
            Sleep(10000);*/
        del_vertex(indices, i_);
        /*cout << "\nindices.size = " << indices.size() << endl;
            cout << "Sleep(10000)" << endl;
            Sleep(10000);
            cout << "nextnextnext" << endl;*/
    }
}

// 获取索引i的前驱结点
int EarClip::getPrev(std::vector<int> &indices, int i)
{
    if(i < 0 || i >= indices.size())
    {
        std::cerr << i << " exceeds the scope of [0, " << indices.size()-1 << "] in getPrev" << endl;
        return -1;
    }

    return i == 0 ? indices.size()-1 : i-1;
}

// 获取索引i的后继结点
int EarClip::getNext(std::vector<int> &indices, int i)
{
    if(i < 0 || i >= indices.size())
    {
        std::cerr << i << " exceeds the scope of [0, " << indices.size()-1 << "] in getNext" << endl;
        return -1;
    }

    return i == indices.size()-1 ? 0 : i+1;
}

// 删除vector中的第i个元素
void EarClip::del_vertex(std::vector<int> &indices, int i)
{
    if(i < 0 || i >= indices.size()){
        std::cerr << i << " exceeds the scope of [0, " << indices.size()-1 << "] in del_vertex" << endl;
        return;
    }

    int j = 0;
    for(auto iter = indices.begin(); iter != indices.end(); ++iter){
        if(j == i){
            indices.erase(iter);
            return;
        }else{
            ++j;
        }
    }
}

// 某点是否在三角形内部
// 目标顶点的索引: p_index
// 三角形顶点索引: tri_v1, tri_v2, tri_v3
bool EarClip::isPointInTriangle(int p_index, int tri_v1, int tri_v2, int tri_v3)
{
    // 首先构造向量
    Eigen::Vector3d v_12, v_13, v_p1, v_p2, v_p3;
    v_12 << vertices[tri_v2]->point.x - vertices[tri_v1]->point.x,
            vertices[tri_v2]->point.y - vertices[tri_v1]->point.y,
            vertices[tri_v2]->point.z - vertices[tri_v1]->point.z;
    v_13 << vertices[tri_v3]->point.x - vertices[tri_v1]->point.x,
            vertices[tri_v3]->point.y - vertices[tri_v1]->point.y,
            vertices[tri_v3]->point.z - vertices[tri_v1]->point.z;

    v_p1 << vertices[tri_v1]->point.x - vertices[p_index]->point.x,
            vertices[tri_v1]->point.y - vertices[p_index]->point.y,
            vertices[tri_v1]->point.z - vertices[p_index]->point.z;
    v_p2 << vertices[tri_v2]->point.x - vertices[p_index]->point.x,
            vertices[tri_v2]->point.y - vertices[p_index]->point.y,
            vertices[tri_v2]->point.z - vertices[p_index]->point.z;
    v_p3 << vertices[tri_v3]->point.x - vertices[p_index]->point.x,
            vertices[tri_v3]->point.y - vertices[p_index]->point.y,
            vertices[tri_v3]->point.z - vertices[p_index]->point.z;

    double eps = 0.000001;
    double S = v_12.cross(v_13).norm() *0.5;
    double S1 = v_p1.cross(v_p2).norm() *0.5;
    double S2 = v_p1.cross(v_p3).norm() *0.5;
    double S3 = v_p2.cross(v_p3).norm() *0.5;
    return fabs(S1+S2+S3-S) <= eps;
}

// 获取多边形顶点的内角
// vertices: 顶点序列，逆时针排列
// normal: 多边形法向量，指向物体外部
// index: 目标顶点索引
// prev: index的前驱
// next: index的后继
double EarClip::getInnerAngle(vector<Vertex *> &vertices, Eigen::Vector3f &normal, int index, int prev, int next)
{
    // 首先排除冗余点情况
    // 三角形T：<index, prev, next>
    // T的最小边长不能小于阈值，如此便可排除向量为0的情况
    float T_min_side_length = 0.001;    // 1毫米
    float d1 = distP2P(vertices[index]->point, vertices[prev]->point);
    float d2 = distP2P(vertices[index]->point, vertices[next]->point);
    float d3 = distP2P(vertices[prev]->point, vertices[next]->point);
    float d_min = d1;
    if(d_min > d2) d_min = d2;
    if(d_min > d3) d_min = d3;
    if(d_min < T_min_side_length){
        std::cerr << "有两个顶点的距离过近，为" << d_min << "m； 大于阈值T_min_side_length:" << d_min << "m" << endl;
        return -1;
    }

    // 接下来，检测三点共线的特殊情况
    Eigen::Vector3d v1, v2, n;
    v1 << vertices[prev]->point.x - vertices[index]->point.x,
          vertices[prev]->point.y - vertices[index]->point.y,
          vertices[prev]->point.z - vertices[index]->point.z;
    if(v1.norm() > 1) v1.normalize();   // 只有向量长度过长，才考虑将其归一化
    v2 << vertices[next]->point.x - vertices[index]->point.x,
          vertices[next]->point.y - vertices[index]->point.y,
          vertices[next]->point.z - vertices[index]->point.z;
    if(v2.norm() > 1) v2.normalize();
    n << normal[0], normal[1], normal[2];

    double eps = 0.000001;
    // 叉乘的范数为0，则三点共线；v1可能与v2成180°，也可能成0°。
    Eigen::Vector3d v1_cross_v2 = v1.cross(v2);
    if(v1_cross_v2.norm() <= eps){
        v1.normalize();
        v2.normalize();
        double cos_angle = v1.dot(v2);
        if(cos_angle < 0) {
            return PI;
        } else{
            /*std::cerr << "出现了尖点或尖缝， 余弦值：" << cos_angle << endl;
            std::cerr << "current point: " << index << vertices->points[index] << endl;
            std::cerr << "prev point: " << prev << vertices->points[prev] << endl;
            std::cerr << "next point: " << next << vertices->points[next] << endl;*/

            return 0;
        }
    }

    v1_cross_v2.normalize();
    // 叉乘与法向量反向，内角＜PI
    if(v1_cross_v2.dot(n) < 0){
        return acos(v1.dot(v2));
    }else{
        // 叉乘与法向量同向，内角＞PI
        return 2.0*PI - acos(v1.dot(v2));
    }
}

// 点-点距离函数
float EarClip::distP2P(pcl::PointXYZ &p1, pcl::PointXYZ &p2)
{
    float dist = (p1.x-p2.x)*(p1.x-p2.x) +
                 (p1.y-p2.y)*(p1.y-p2.y) +
                 (p1.z-p2.z)*(p1.z-p2.z);
    return pow(dist, 0.5f);
}





































