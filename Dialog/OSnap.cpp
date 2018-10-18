#include "osnap.h"

OSnap::OSnap()
{
    m_poly = NULL;
	this->delta = 0.01;
    this->max_k_count = 50;
	f1_weight = f2_weight = f3_1_weight = f3_2_weight = f8_weight = f4_weight = f5_weight = f6_weight = f7_weight = 1;
}

// 重置数据结构
void OSnap::reset()
{
	this->m_polygons = NULL;
    this->m_poly = NULL;

    this->normal_vertices.clear();
    this->p2e_vertices.clear();
    this->p2p_vertices.clear();
	normal_vertices_indices.clear();
	p2p_vertices_indices.clear();
	p2e_vertices_indices.clear();

    m_osnap_poly.c << 0, 0, 0;
    m_osnap_poly.c_ << 0, 0, 0;
    m_osnap_poly.b0 << 0, 0, 0;
    m_osnap_poly.b1 << 0, 0, 0;
    m_osnap_poly.b2 << 0, 0, 0;
    m_osnap_poly.coeff << 0, 0, 0, 0;
    m_osnap_poly.coeff_origin << 0, 0, 0, 0;
    m_osnap_poly.vertices.clear();
    m_osnap_poly.vertices_origin.clear();
    m_osnap_poly.coors.clear();

	J.resize(0, 0);
	beta.resize(0);

    //this->delta = 0.01;
}

// 根据输入的多边形对象m_poly设置目标多边形对象m_osnap_poly
void OSnap::setOSnapPoly()
{
	int vertices_size = m_poly->getSize();
	this->m_osnap_poly.vertices.resize(vertices_size);
	this->m_osnap_poly.vertices_origin.resize(vertices_size);

    // 多边形的原始顶点
    Vertex *pv;
    for(int i = 0; i < m_poly->getSize(); ++i){
        pv = i == 0 ? m_poly->start_point : pv->next_point;
		m_osnap_poly.vertices_origin[i] << pv->point.x, pv->point.y, pv->point.z;
		m_osnap_poly.vertices[i] = m_osnap_poly.vertices_origin[i];
    }
    // 多边形方程系数
    m_osnap_poly.coeff_origin << m_poly->coeff[0],
                                 m_poly->coeff[1],
                                 m_poly->coeff[2],
                                 m_poly->coeff[3];
    m_osnap_poly.coeff = m_osnap_poly.coeff_origin;
    // 设置移动坐标系
    setBases_init(m_osnap_poly);
    // 设置多边形顶点在移动坐标系中的二维坐标
    set2DCoordinates(m_osnap_poly);
}

// 利用顶点构造初始的移动坐标系
void OSnap::setBases_init(OSnapPoly &poly)
{
	Eigen::Vector3d v_centroid;
	v_centroid.setZero();	// 重置每个元素为0
    for(int i = 0; i < poly.vertices.size(); ++i)
    {
        v_centroid += poly.vertices[i];
    }

    // 移动坐标系的原点=顶点集合的重心
	v_centroid /= poly.vertices.size();
    poly.b0 = v_centroid;

    Eigen::Vector3d v, e1, e2;
    // 移动坐标系的一个轴e1指向多边形的某个顶点
    for(int i = 0; i < poly.vertices.size(); ++i)
    {
        Eigen::Vector3d v_p;
		v_p = poly.vertices[i];
        v = v_p - poly.b0;
        if(v.norm() > 0.01)
        {
            e1 = 1.0 / v.norm() * v;
            break;
        }
    }

    // 构造平面的法向量
    Eigen::Vector3d normal;
    normal << poly.coeff[0], poly.coeff[1], poly.coeff[2];

    e2 = e1.cross(normal);
    e2 = 1.0 / e2.norm() * e2;

    poly.b1 = poly.b0 + e1;
    poly.b2 = poly.b0 + e2;
}

// 设置多边形顶点在移动坐标系中的二维坐标
void OSnap::set2DCoordinates(OSnapPoly &poly)
{
    // 首先清除原有的二维坐标
    poly.coors.clear();

    Eigen::Matrix<double, 3, 2> A;
    Eigen::Vector3d e1, e2;
    e1 = poly.b1 - poly.b0;
    e2 = poly.b2 - poly.b0;
    A << e1[0], e2[0],
         e1[1], e2[1],
         e1[2], e2[2];
    Eigen::Vector2d coor;
    Coor_sigma coor_sigma;
    for(int i = 0; i < poly.vertices.size(); ++i)
    {
        Eigen::Vector3d p;
		p = poly.vertices[i];
        coor = (A.transpose() * A).inverse() * A.transpose() * (p - poly.b0);
        coor_sigma.px = coor[0];
        coor_sigma.py = coor[1];
        poly.coors.push_back(coor_sigma);
    }
}

// 根据匹配关系分类顶点
void OSnap::categoriseVertices()
{
    if(!m_poly){
        cerr << "error in OSnap::categoriseVertices(): m_poly == NULL" << endl;
		return;
	}

    Vertex *pv;
    for(int i = 0; i < m_poly->getSize(); ++i){
        pv = i == 0 ? m_poly->start_point : pv->next_point;
        if(pv->neighbor_poly_id == -1){
            // 没有匹配关系
            this->normal_vertices.push_back(pv);
			normal_vertices_indices.push_back(i);
        }else if(pv->is_vertex_neighbor){
            // 点-点匹配
            this->p2p_vertices.push_back(pv);
			p2p_vertices_indices.push_back(i);
        }else{
            // 点-边匹配
            this->p2e_vertices.push_back(pv);
			p2e_vertices_indices.push_back(i);
        }		
    }
}

// 匀速螺旋运动
// 设置变量
void OSnap::setHelicalVars( Eigen::Vector3d &c,
                     Eigen::Vector3d &c_,
                     Eigen::Matrix3d &R,
                     Eigen::Vector3d &pointOfA,
                     Eigen::Vector3d &a,
                     Eigen::Vector3d &a_,
                     double &pitch,
                     double &alpha )
{
    pitch = c.dot(c_) / c.dot(c);
    a = (1.0/c.norm()) * c;
    a_ = c_ - pitch * c;
    alpha = atan(c.norm());
    pointOfA = a.cross(a_);

    double b0 = cos(alpha*0.5);
    double b1 = a[0] * sin(alpha*0.5);
    double b2 = a[1] * sin(alpha*0.5);
    double b3 = a[2] * sin(alpha*0.5);

    double m00 = b0 * b0 + b1 * b1 + b2 * b2 + b3 * b3;
    double m11 = b0 * b0 + b1 * b1 - b2 * b2 - b3 * b3;
    double m22 = b0 * b0 - b1 * b1 + b2 * b2 - b3 * b3;
    double m33 = b0 * b0 - b1 * b1 - b2 * b2 + b3 * b3;

    R << m11, 2.0*(b1*b2 + b0*b3), 2.0*(b1*b3 - b0*b2),
         2.0*(b1*b2 - b0*b3), m22, 2.0*(b2*b3 + b0*b1),
         2.0*(b1*b3 + b0*b2), 2.0*(b2*b3 - b0*b1), m33;
    R = 1.0/m00 * R;
}

// 执行移动
void OSnap::uniformHelicalMove( Eigen::Vector3d &c,
                         Eigen::Vector3d &c_,
                         Eigen::Vector3d &x,
                         Eigen::Vector3d &x_,
                         Eigen::Matrix3d &R,
                         Eigen::Vector3d &pointOfA,
                         Eigen::Vector3d &a,
                         double &pitch,
                         double &alpha )
{
    // 如果角速度不是零向量
    if(c.norm() > 0.000000000001)
        x_ = R * (x - pointOfA) + pitch * alpha * a + pointOfA;
    else
        x_ = x + c_;
}

// 更新经过螺旋运动后的移动坐标系的基
void OSnap::updateBasesAfterHelicalMove(OSnapPoly &poly)
{
    Eigen::Matrix3d R;
    Eigen::Vector3d pointOfA;
    Eigen::Vector3d a;
    Eigen::Vector3d a_;
    double pitch;
    double alpha;
    setHelicalVars( poly.c, poly.c_, R, pointOfA, a, a_, pitch, alpha);

    Eigen::Vector3d x_;
    uniformHelicalMove( poly.c, poly.c_, poly.b0, x_, R, pointOfA, a, pitch, alpha);
    poly.b0 << x_[0], x_[1], x_[2];
    uniformHelicalMove( poly.c, poly.c_, poly.b1, x_, R, pointOfA, a, pitch, alpha);
    poly.b1 << x_[0], x_[1], x_[2];
    uniformHelicalMove( poly.c, poly.c_, poly.b2, x_, R, pointOfA, a, pitch, alpha);
    poly.b2 << x_[0], x_[1], x_[2];
}

// 更新经过螺旋运动后顶点的三维坐标
// (前提：‘移动坐标系的基’+‘二维坐标’已被更新过！)
void OSnap::updateVerticesAfterHelicalMove(OSnapPoly &poly)
{
    Eigen::Vector3d p;
    for(int i = 0; i < poly.vertices.size(); ++i)
    {
        p = poly.b0 + poly.coors[i].px * (poly.b1 - poly.b0)
                    + poly.coors[i].py * (poly.b2 - poly.b0);
        poly.vertices[i] = p;
    }
}

// 根据变量向量beta更新每个顶点在移动坐标系中的二维坐标
void OSnap::updateCoors()
{
	int start_index = 6;
	int v_index;
	int beta_index;
	// 更新具有p2p匹配关系的顶点在移动坐标系中的二维坐标
	for(int i = 0; i < this->p2p_vertices_indices.size(); ++i){
		v_index = p2p_vertices_indices[i];
		beta_index = start_index + 2*i;
		m_osnap_poly.coors[v_index].px += beta[beta_index];
		m_osnap_poly.coors[v_index].py += beta[beta_index+1];
	}

	start_index += p2p_vertices.size() * 2;
	// 更新具有p2e匹配关系的顶点在移动坐标系中的二维坐标
	for(int i = 0; i < this->p2e_vertices_indices.size(); ++i){
		v_index = p2e_vertices_indices[i];
		beta_index = start_index + 2*i;
		m_osnap_poly.coors[v_index].px += beta[beta_index];
		m_osnap_poly.coors[v_index].py += beta[beta_index+1];
	}

	/*for(int i = 0; i < this->m_osnap_poly.coors.size(); ++i){
		int j = 2*i;
		this->m_osnap_poly.coors[i].px += beta[start_index + j];
		this->m_osnap_poly.coors[i].py += beta[start_index + j + 1];
	}*/
}

// 更新后顶点的位置与原来位置的距离
double OSnap::f1(int vertex_index)
{
	//                                当前顶点的位置                        原始顶点的位置
	Eigen::Vector3d v_diff = m_osnap_poly.vertices[vertex_index] - m_osnap_poly.vertices_origin[vertex_index];
    return f1_weight * v_diff.norm();
}

// 更新后顶点位置与目标顶点位置的距离
double OSnap::f2(int vertex_index_in_p2p_vertices)
{
	// 更新后的顶点
	Eigen::Vector3d v_update, v_target, v_diff;
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	Vertex *pv = this->p2p_vertices[vertex_index_in_p2p_vertices]->neighbor_vertex;
	if(!pv){
		cerr << "OSnap::f2, pv == NULL!" << endl;
		return 0;
	}
	// 目标顶点
	v_target << pv->point.x, pv->point.y, pv->point.z;

	v_diff = v_update - v_target;
    return f2_weight * v_diff.norm();
}

// dist[更新后顶点，目标顶点的前驱节点] - dist[目标顶点，目标顶点的前驱节点]
double OSnap::f3_1(int vertex_index_in_p2p_vertices)
{
	Eigen::Vector3d v_update, v_target, v_target_prev;
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	Vertex *pv = p2p_vertices[vertex_index_in_p2p_vertices]->neighbor_vertex;
	if(!pv){
		cerr << "OSnap::f3_1, pv == NULL!" << endl;
		return 0;
	}

	v_target << pv->point.x, pv->point.y, pv->point.z;
	v_target_prev << pv->prev_point->point.x,
		             pv->prev_point->point.y,
		             pv->prev_point->point.z;
	double dist1 = (v_update - v_target_prev).norm();
	double dist2 = (v_target - v_target_prev).norm();
	return f3_1_weight * abs(dist1-dist2);
}

// dist[更新后顶点，目标顶点的后继节点] - dist[目标顶点，目标顶点的后继节点]
double OSnap::f3_2(int vertex_index_in_p2p_vertices)
{
	Eigen::Vector3d v_update, v_target, v_target_next;
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	Vertex *pv = p2p_vertices[vertex_index_in_p2p_vertices]->neighbor_vertex;
	if(!pv){
		cerr << "OSnap::f3_2, pv == NULL!" << endl;
		return 0;
	}

	v_target << pv->point.x, pv->point.y, pv->point.z;
	v_target_next << pv->next_point->point.x,
		             pv->next_point->point.y,
		             pv->next_point->point.z;
	double dist1 = (v_update - v_target_next).norm();
	double dist2 = (v_target - v_target_next).norm();
	return f3_2_weight * abs(dist1-dist2);
}

// f5的对偶版本，更新后顶点位置与两个平面交线的距离
double OSnap::f8(int vertex_index_in_p2p_vertices)
{
	// 更新后的顶点
	Eigen::Vector3d v_update;
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	// 目标多边形
	POLYGON *target_poly;

	int target_poly_id = p2p_vertices[vertex_index_in_p2p_vertices]->neighbor_poly_id;
	if(target_poly_id == -1){
		cerr << "in f8, target_poly_id == -1" << endl;
	}
	target_poly = &((*m_polygons)[target_poly_id]);

	// 获取两个多边形的交线
	Eigen::Vector3d p_base, line_dir;
	if(!getIntersectLine(m_poly->coeff, target_poly->coeff, p_base, line_dir)){
		cerr << "OSnap::f8, 多边形平行!" << endl;
		return FLT_MAX;
	}

	Eigen::Vector3d v_p, delta;
	v_p = v_update;
	delta = p_base - v_p;
	double lambda = -1 * line_dir.dot(delta);

	Eigen::Vector3d p_proj;
	p_proj = p_base + lambda * line_dir;

    return f8_weight * (v_p - p_proj).norm();
}

// 更新后顶点位置到匹配边的距离
double OSnap::f4(int vertex_index_in_p2e_vertices)
{
	Eigen::Vector3d v_update, e_begin, e_end;

	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	Vertex *pv = p2e_vertices[vertex_index_in_p2e_vertices]->neighbor_edge;
	if(!pv){
		cerr << "OSnap::f4, pv == NULL!" << endl;
		return 0;
	}

	e_begin << pv->point.x,
		       pv->point.y,
			   pv->point.z;
	e_end << pv->next_point->point.x,
		     pv->next_point->point.y,
			 pv->next_point->point.z;

	// 计算v_update到line[e_begin, e_end]的距离
	Eigen::Vector3d p_base, line_dir;
	p_base = e_begin;
	line_dir = e_end-e_begin;
	line_dir.normalize();

	Eigen::Vector3d v_p, delta;
	v_p = v_update;
	delta = p_base - v_p;
	double lambda = -1 * line_dir.dot(delta);

	Eigen::Vector3d p_proj;
	p_proj = p_base + lambda * line_dir;

    return f4_weight * (v_p - p_proj).norm();
}

// 计算两个多边形的交线
bool OSnap::getIntersectLine(Eigen::Vector4f &coeff1, Eigen::Vector4f &coeff2,
                             Eigen::Vector3d &p_base, Eigen::Vector3d &line_dir)
{
	 Eigen::Vector3d n1, n2;
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
    double min = coeff1[0];
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

// 更新后顶点位置与两个平面交线的距离
double OSnap::f5(int vertex_index_in_p2e_vertices)
{
	// 更新后的顶点
	Eigen::Vector3d v_update;
	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	// 目标多边形
	POLYGON *target_poly;
	Vertex *pv = p2e_vertices[vertex_index_in_p2e_vertices]->neighbor_edge;
	if(!pv){
		cerr << "OSnap::f5, pv == NULL!" << endl;
		return 0;
	}

	int target_poly_id = floor(pv->point.data[3]+0.5f);
	target_poly = &((*m_polygons)[target_poly_id]);

	// 获取两个多边形的交线
	Eigen::Vector3d p_base, line_dir;
	if(!getIntersectLine(m_poly->coeff, target_poly->coeff, p_base, line_dir)){
		cerr << "OSnap::f5, 多边形平行!" << endl;
		return FLT_MAX;
	}

	Eigen::Vector3d v_p, delta;
	v_p = v_update;
	delta = p_base - v_p;
	double lambda = -1 * line_dir.dot(delta);

	Eigen::Vector3d p_proj;
	p_proj = p_base + lambda * line_dir;

    return f5_weight * (v_p - p_proj).norm();
}

// 更新后顶点到匹配边两个端点距离之和 - 匹配边长度
double OSnap::f6(int vertex_index_in_p2e_vertices)
{
	Eigen::Vector3d v_update, e_begin, e_end;
	// 更新后的顶点
	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];
	v_update = m_osnap_poly.vertices[v_index];

	// 匹配边
	Vertex *pv = p2e_vertices[vertex_index_in_p2e_vertices]->neighbor_edge;
	if(!pv){
		cerr << "OSnap::f4, pv == NULL!" << endl;
		return 0;
	}

	// 匹配边的两个端点
	e_begin << pv->point.x,
		       pv->point.y,
			   pv->point.z;
	e_end << pv->next_point->point.x,
		     pv->next_point->point.y,
			 pv->next_point->point.z;

	double dist1 = (v_update - e_begin).norm();
	double dist2 = (v_update - e_end).norm();
	double dist3 = (e_end - e_begin).norm();

	return f6_weight * abs(dist1 + dist2 - dist3);
}

// 更新后顶点到匹配关系所属支撑平面的距离
double OSnap::f7(int vertex_index, Vertex *pv)
{
    //构造齐次顶点
    Eigen::Vector4d homogeneous_vertex;
    homogeneous_vertex << m_osnap_poly.vertices[vertex_index][0],
                          m_osnap_poly.vertices[vertex_index][1],
                          m_osnap_poly.vertices[vertex_index][2],
                          1;

    // 获取匹配关系所属多边形的索引
    int poly_index;
    if(pv->neighbor_poly_id == -1){
        cerr << "in OSnap::f7, pv->neighbor_poly_id == -1!" << endl;
        return 0;
    }

    poly_index = pv->neighbor_poly_id;

    // 构造方程系数
    Eigen::Vector4d coeff;
    coeff << (*m_polygons)[poly_index].coeff[0],
             (*m_polygons)[poly_index].coeff[1],
             (*m_polygons)[poly_index].coeff[2],
             (*m_polygons)[poly_index].coeff[3];

    double dist = homogeneous_vertex.dot( coeff );
    return f7_weight * abs(dist);
}

// 初始化变量向量beta,在每次优化的开始使用一次
void OSnap::initBeta()
{
	int size = this->m_osnap_poly.vertices.size();
	if(size == 0){
		cerr << "in OSnap::initBeta(), m_osnap_poly.vertices->size() == 0" << endl;
		return;
	}

	// 前六个参数：角速度，线速度
	// 后面跟2*n个参数：n个具有匹配关系的顶点，其中每个顶点相对于局部坐标系的微小偏移(delta_x, delta_y)
	size = 6 + (p2p_vertices.size() + p2e_vertices.size())*2;

	beta.resize(size);
	beta.setZero();	// 将参数向量中所有的元素均初始化为0
	for(int i = 0; i < beta.size(); ++i){
		beta[i] = 0.0001;
	}
}

// 更新移动坐标系
void OSnap::updateMovingCoor()
{
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 首先从beta中提取角速度c与线速度c_
	/*m_osnap_poly.c << beta[0], beta[1], beta[2];
	m_osnap_poly.c_ << beta[3], beta[4], beta[5];*/
	m_osnap_poly.c[0] = beta[0];
	m_osnap_poly.c[1] = beta[1];
	m_osnap_poly.c[2] = beta[2];
	m_osnap_poly.c_[0] = beta[3];
	m_osnap_poly.c_[1] = beta[4];
	m_osnap_poly.c_[2] = beta[5];

	// 利用c与c_设置匀速螺旋运动的相关变量
	setHelicalVars(m_osnap_poly.c, m_osnap_poly.c_, R, pointOfA, a, a_, pitch, alpha);

	// 更新经过螺旋运动后的移动坐标系的基
	Eigen::Vector3d x_;
    uniformHelicalMove( m_osnap_poly.c, m_osnap_poly.c_, m_osnap_poly.b0, x_, R, pointOfA, a, pitch, alpha);
    m_osnap_poly.b0 = x_;
    uniformHelicalMove( m_osnap_poly.c, m_osnap_poly.c_, m_osnap_poly.b1, x_, R, pointOfA, a, pitch, alpha);
    m_osnap_poly.b1 = x_;
    uniformHelicalMove( m_osnap_poly.c, m_osnap_poly.c_, m_osnap_poly.b2, x_, R, pointOfA, a, pitch, alpha);
    m_osnap_poly.b2 = x_;	
}

// 更新顶点位置
void OSnap::updateVertices()
{
	// 更新顶点二维坐标
	updateCoors();
	// 更新顶点位置
	updateVerticesAfterHelicalMove(this->m_osnap_poly);
}

// 构造Jacob矩阵
void OSnap::cosntructJacobMat()
{
	/*
	* 没有匹配关系的顶点对应的目标函数为：f1
    * 点-点匹配的顶点对应的目标函数为：f1, f2, f3_1, f3_2, f7, f8
    * 点-边匹配的顶点对应的目标函数为：f1, f4, f5, f6, f7
	*/
	int rows = normal_vertices.size() +
               p2p_vertices.size() * 6 +
               p2e_vertices.size() * 5;
	int cols = 6 + (p2p_vertices.size() + p2e_vertices.size()) * 2;

	J.resize(rows, cols);
	J.setZero();	// 所有元素重置为0

	// 首先处理没有匹配关系的顶点
	// 求偏导数的步长
	double delta = this->delta;

	int v_index;

	/*
	* 没有匹配关系的顶点,其在移动坐标系的二维坐标保持不变，但是移动坐标系的变化会改变它们的位置，因而它们与原来位置的距离不能太远
	* 目标函数：f1
	* beta：前六个变量
	*/
	for(int i = 0; i < normal_vertices.size(); ++i){		
		v_index = normal_vertices_indices[i];

		pd_f1_normal_vertices(v_index, J, i, delta);
	}	// 没有匹配关系的顶点

	/*
	* 点-点匹配关系的顶点,移动坐标系有可能变动，其本身相对于移动坐标系的二维坐标也可能有一些轻微的变动
    * 目标函数: f1, f2, f3_1, f3_2, f7, f8共六个目标函数
	* beta: beta[0-5], beta[j],beta[j+1] 共八个变量
	*/
	int start_row_index = normal_vertices.size();	// 行索引
	int start_index = 6;	// 列索引
	for(int i = 0; i < p2p_vertices.size(); ++i)
	{
		v_index = p2p_vertices_indices[i];

        pd_f1_matching_vertices(v_index, J, start_row_index + 6*i,    start_index + 2*i, delta);
        pd_f2  (i, J,                       start_row_index + 6*i + 1, start_index + 2*i, delta);
        pd_f3_1(i, J,                       start_row_index + 6*i + 2, start_index + 2*i, delta);
        pd_f3_2(i, J,                       start_row_index + 6*i + 3, start_index + 2*i, delta);
        pd_f7(v_index, p2p_vertices[i], J,  start_row_index + 6*i + 4, start_index + 2*i, delta);
		pd_f8(i, J,                         start_row_index + 6*i + 5, start_index + 2*i, delta);
	}	// 点-点匹配

	/*
	* 点-边匹配关系的顶点,移动坐标系有可能变动，其本身相对于移动坐标系的二维坐标也可能有一些轻微的变动
    * 目标函数: f1, f4, f5, f6, f7 共五个目标函数
	* beta: beta[0-5], beta[j],beta[j+1] 共八个变量
	*/
    start_row_index = normal_vertices.size() + 6*p2p_vertices.size();	// 行索引
	start_index = 6 + 2*p2p_vertices.size();	// 列索引
	for(int i = 0; i < p2e_vertices.size(); ++i){
		v_index = p2e_vertices_indices[i];

        pd_f1_matching_vertices(v_index, J, start_row_index + 5*i, start_index + 2*i, delta);
        pd_f4(i, J,                         start_row_index + 5*i + 1, start_index + 2*i, delta);
        pd_f5(i, J,                         start_row_index + 5*i + 2, start_index + 2*i, delta);
        pd_f6(i, J,                         start_row_index + 5*i + 3, start_index + 2*i, delta);
        pd_f7(v_index, p2e_vertices[i], J,  start_row_index + 5*i + 4, start_index + 2*i, delta);
	}	// 点-边匹配
}

// 执行单步高斯牛顿迭代，从而更新beta
void OSnap::performGN()
{
	Eigen::MatrixXd A;
	Eigen::VectorXd b, x;
	int n = 6 + 2*(p2p_vertices.size() + p2e_vertices.size());
	A.resize(n, n);
	A.setZero();

	A = J.transpose()*J;

	b.resize(n);
	b.setZero();
	x.resize(n);
	x.setZero();

	Eigen::VectorXd r;
	fill_r(r);
    double sum_r = 0;
    for(int i = 0; i < r.size(); ++i){sum_r += r[i];}
    cout << "sum(r1+...+rn) = " << sum_r << endl;
	//cout << "r = \n" << r << endl;
	
	b = J.transpose() * r;

	//x = A.fullPivHouseholderQr().solve(b);
	//x = A.householderQr().solve(b);
	//conjugateGradient(A, x, b);

	// 备份beta的值
    Eigen::VectorXd beta_bck /*保存beta的初值*/, beta_opt/*保存beta的最优值*/;
	beta_bck.resize(beta.size());
	beta_bck = beta;
	beta_opt.resize(beta.size());

	int k_count = this->max_k_count;

    int k_count_valid = -1;  // 保存最优的迭代次数

	for(int i = 1; i <= k_count; ++i){
        this->conjugateGradient_with_max_k(A, x, b, i); // 共轭梯度迭代i次
        beta = beta-x;  // 更新beta
        fill_r(r);      // 更新目标函数

		double sum_r_tmp = 0;
		for(int j = 0; j < r.size(); ++j){sum_r_tmp += r[j];}

		if(sum_r_tmp < sum_r){	// sum_r 保存当前最小的目标函数值
			sum_r = sum_r_tmp;
			k_count_valid = i;
			beta_opt = beta;
		}

		// 重置beta;
		beta = beta_bck;
	}

	cout << "k_count_valid = " << k_count_valid << endl;
	cout << "sum_r = " << sum_r << endl;

	// 将beta设置为最优值
	if(k_count_valid > 0)
		beta = beta_opt;
	// 更新目标函数
    fill_r(r);

	sum_r = 0;
	for(int i = 0; i < r.size(); ++i){sum_r += r[i];}
    cout << "sum_r = " << sum_r << endl;

	//beta = beta - x;

	//beta = beta - (J.transpose()*J).inverse()*J.transpose()*r;
}

void OSnap::fill_r(Eigen::VectorXd &r)
{
	// 重置r
	int m = normal_vertices.size() + 6*p2p_vertices.size() + 5*p2e_vertices.size();
	r.resize(m);
	r.setZero();

	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	this->setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	int v_index;

	// 没有匹配关系的顶点
	for(int i = 0; i < normal_vertices.size(); ++i){
		v_index = normal_vertices_indices[i];
		m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                            + m_osnap_poly.coors[v_index].py * (b2 - b0);
		r[i] = f1(v_index);
	}

	// 点-点
	int start_row_index = normal_vertices.size();	// 行索引
	int start_index = 6;	// 列索引
	double px, py;
	for(int i = 0; i < p2p_vertices.size(); ++i)
	{
		v_index = p2p_vertices_indices[i];

		px = m_osnap_poly.coors[v_index].px + beta[start_index + 2*i];
		py = m_osnap_poly.coors[v_index].py + beta[start_index + 2*i + 1];

		m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

        r[start_row_index + 6*i] = f1(v_index);
        r[start_row_index + 6*i + 1] = f2(i);
        r[start_row_index + 6*i + 2] = f3_1(i);
        r[start_row_index + 6*i + 3] = f3_2(i);
        r[start_row_index + 6*i + 4] = f7(v_index, p2p_vertices[i]);
		r[start_row_index + 6*i + 5] = f8(i);
	}

	// 点-边
    start_row_index = normal_vertices.size() + 6*p2p_vertices.size();	// 行索引
	start_index = 6 + 2*p2p_vertices.size();	// 列索引
	for(int i = 0; i < p2e_vertices.size(); ++i){
		v_index = p2e_vertices_indices[i];

		px = m_osnap_poly.coors[v_index].px + beta[start_index + 2*i];
		py = m_osnap_poly.coors[v_index].py + beta[start_index + 2*i + 1];

		m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

        r[start_row_index + 5*i] = f1(v_index);
        r[start_row_index + 5*i + 1] = f4(i);
        r[start_row_index + 5*i + 2] = f5(i);
        r[start_row_index + 5*i + 3] = f6(i);
        r[start_row_index + 5*i + 4] = f7(v_index, p2e_vertices[i]);
	}	// 点-边匹配
}

// 共轭梯度解方程组
void OSnap::conjugateGradient(Eigen::MatrixXd &A, Eigen::VectorXd &x, Eigen::VectorXd &b)
{
	//cout << "A.diagonal() = " << A.diagonal() << endl;
	cout << "A.inverse() = \n" << A.inverse() << endl;
	cout << "A.inverse()*A = \n" << A.inverse()*A << endl;

	int n = b.size();
	x.resize(n);
	x.setZero();	// x的初值置为0

	Eigen::VectorXd r, r_, p;
	r.resize(n);
	p.resize(n);
	r_.resize(n);
	
	r = b - A*x;
	p = r;

	int k = 0;
	double alpha, beta;
	do{
		if(k >= 10000){
			cout << "k == 10000, in conjugateGradient()" << endl;
			break;
		}

		alpha = r.dot(r) / (p.transpose()*A*p);
		x = x + alpha * p;
		r_ = r - alpha * (A * p);

        /*if(r.norm() <= this->max_r_norm)	// 阈值不是设置得越小越好，阈值越小，迭代次数越多，会引入更多的误差
			break;*/

		beta = r_.dot(r_) / r.dot(r);
		p = r_ + beta*p;
		r = r_;
		++k;
	}while(true);

	cout << "k = " << k << endl;
}

void OSnap::conjugateGradient_with_max_k(Eigen::MatrixXd &A, Eigen::VectorXd &x, Eigen::VectorXd &b, int k_count)
{
	int n = b.size();
	x.resize(n);
	x.setZero();	// x的初值置为0

	Eigen::VectorXd r, r_, p;
	r.resize(n);
	p.resize(n);
	r_.resize(n);
	
	r = b - A*x;
	p = r;

	double alpha, beta;

	for(int k = 0; k < k_count; ++k){
		alpha = r.dot(r) / (p.transpose()*A*p);
		x = x + alpha * p;
		r_ = r - alpha * (A * p);

		beta = r_.dot(r_) / r.dot(r);
		p = r_ + beta*p;
		r = r_;
	}
}

// 最终的优化执行函数
void OSnap::optimise()
{

}

// 目标函数的偏导数
// 针对没有匹配关系的顶点：normal_vertices
void OSnap::pd_f1_normal_vertices(int v_index, Eigen::MatrixXd &J, int row, double delta)
{
	// 计算根据beta前六个元素更新的当前f1值
	// 首先计算没有更新时的f1值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	this->setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 更新m_osnap_poly.vertices[v_index]
	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	// 计算f1(beta)
	double f1_val = this->f1(v_index);	// 没有改变任何beta变量的f1值
	double f1_delta_val;				// 改变一个beta变量的f1值
	double pd;							// 偏导数

	// 下面计算f1(beta)针对beta前六个变量的偏导数
	// beta[0]///////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	m_osnap_poly.vertices[v_index] = b0 + m_osnap_poly.coors[v_index].px * (b1 - b0)
                                        + m_osnap_poly.coors[v_index].py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 5) = pd;
}

// 针对有匹配关系的顶点: p2p_vertices 和 p2e_vertices；不仅要更新坐标系，而且要更新顶点二维坐标
void OSnap::pd_f1_matching_vertices(int v_index, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 计算根据beta前六个元素更新的当前f1值以及beta[beta_index],beta[beta_index+1]两个元素更新后的f1值
	// 首先计算没有更新时的f1值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f1(beta)
	double f1_val = this->f1(v_index);	// 没有改变任何beta变量的f1值
	double f1_delta_val;				// 改变一个beta变量的f1值
	double pd;							// 偏导数

	// beta[0]//////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 0) = pd;

	// beta[1]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 1) = pd;

	// beta[2]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 2) = pd;

	// beta[3]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 3) = pd;

	// beta[4]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 4) = pd;

	// beta[5]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index + 1]//////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f1_delta_val = f1(v_index);
	pd = (f1_delta_val - f1_val) / delta;
	J(row, beta_index + 1) = pd;
}

// 针对点-点
// f2的偏导数
void OSnap::pd_f2(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];

	// 首先计算没有更新时的f2值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f2(beta)
	double f2_val = f2(vertex_index_in_p2p_vertices);	// 没有改变任何beta变量的f2值
	double f2_delta_val;								// 改变一个beta变量的f2值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f2_delta_val = f2(vertex_index_in_p2p_vertices);
	pd = (f2_delta_val - f2_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-点
// f3_1的偏导数
void OSnap::pd_f3_1(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];

	// 首先计算没有更新时的f3_1值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f3_1(beta)
	double f3_1_val = f3_1(vertex_index_in_p2p_vertices);	// 没有改变任何beta变量的f3_1值
	double f3_1_delta_val;								// 改变一个beta变量的f3_1值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_1_delta_val = f3_1(vertex_index_in_p2p_vertices);
	pd = (f3_1_delta_val - f3_1_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-点
// f3_2的偏导数
void OSnap::pd_f3_2(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];

	// 首先计算没有更新时的f3_2值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f3_2(beta)
	double f3_2_val = f3_2(vertex_index_in_p2p_vertices);	// 没有改变任何beta变量的f3_2值
	double f3_2_delta_val;								// 改变一个beta变量的f3_2值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f3_2_delta_val = f3_2(vertex_index_in_p2p_vertices);
	pd = (f3_2_delta_val - f3_2_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-点
// f8的偏导数
void OSnap::pd_f8(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2p_vertices_indices[vertex_index_in_p2p_vertices];

	// 首先计算没有更新时的f8值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f8(beta)
	double f8_val = f8(vertex_index_in_p2p_vertices);	// 没有改变任何beta变量的f8值
	double f8_delta_val;								// 改变一个beta变量的f8值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f8_delta_val = f8(vertex_index_in_p2p_vertices);
	pd = (f8_delta_val - f8_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-边
// f4的偏导数
void OSnap::pd_f4(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];

	// 首先计算没有更新时的f4值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f4(beta)
	double f4_val = f4(vertex_index_in_p2e_vertices);	// 没有改变任何beta变量的f4值
	double f4_delta_val;								// 改变一个beta变量的f4值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f4_delta_val = f4(vertex_index_in_p2e_vertices);
	pd = (f4_delta_val - f4_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-边
// f5的偏导数
void OSnap::pd_f5(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];

	// 首先计算没有更新时的f5值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f5(beta)
	double f5_val = f5(vertex_index_in_p2e_vertices);	// 没有改变任何beta变量的f5值
	double f5_delta_val;								// 改变一个beta变量的f5值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f5_delta_val = f5(vertex_index_in_p2e_vertices);
	pd = (f5_delta_val - f5_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对点-边
// f6的偏导数
void OSnap::pd_f6(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
	// 顶点在vertices中的索引v_index
	int v_index = p2e_vertices_indices[vertex_index_in_p2e_vertices];

	// 首先计算没有更新时的f6值
	// 匀速螺旋运动的相关变量
	Eigen::Vector3d c;
	Eigen::Vector3d c_;
	Eigen::Matrix3d R;
	Eigen::Vector3d pointOfA;
	Eigen::Vector3d a;
	Eigen::Vector3d a_;
	double pitch;
	double alpha;

	// 移动坐标系的基
	// 原点：b0
	// 两个方向向量：b1-b0, b2-b0
	Eigen::Vector3d b0, b1, b2;

	// 利用beta设置角速度与线速度
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	// 利用c与c_设置其它运动变量
	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	// 对坐标系的基执行匀速螺旋运动
	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	// 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
	double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	// 更新m_osnap_poly.vertices[v_index]
	// b0, b1, b2, px, py 均已被更新
	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	// 计算f6(beta)
	double f6_val = f6(vertex_index_in_p2e_vertices);	// 没有改变任何beta变量的f6值
	double f6_delta_val;								// 改变一个beta变量的f6值
	double pd;											// 偏导数

	// beta[0]///////////////////////////////////////////////////////////////////////
	c << beta[0] + delta, beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 0) = pd;

	// beta[1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1] + delta, beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 1) = pd;

	// beta[2]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2] + delta;
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 2) = pd;

	// beta[3]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3] + delta, beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 3) = pd;

	// beta[4]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4] + delta, beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 4) = pd;

	// beta[5]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5] + delta;

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, 5) = pd;

	// beta[beta_index]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, beta_index) = pd;

	// beta[beta_index+1]///////////////////////////////////////////////////////////////////////
	c << beta[0], beta[1], beta[2];
	c_ << beta[3], beta[4], beta[5];

	setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

	uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
	uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

	px = m_osnap_poly.coors[v_index].px + beta[beta_index];
	py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

	m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

	f6_delta_val = f6(vertex_index_in_p2e_vertices);
	pd = (f6_delta_val - f6_val) / delta;
	J(row, beta_index+1) = pd;
}

// 针对有匹配关系的顶点
void OSnap::pd_f7(int v_index, Vertex *pv, Eigen::MatrixXd &J, int row, int beta_index, double delta)
{
    // 首先计算没有更新时的f7值
    // 匀速螺旋运动的相关变量
    Eigen::Vector3d c;
    Eigen::Vector3d c_;
    Eigen::Matrix3d R;
    Eigen::Vector3d pointOfA;
    Eigen::Vector3d a;
    Eigen::Vector3d a_;
    double pitch;
    double alpha;

    // 移动坐标系的基
    // 原点：b0
    // 两个方向向量：b1-b0, b2-b0
    Eigen::Vector3d b0, b1, b2;

    // 利用beta设置角速度与线速度
    c << beta[0], beta[1], beta[2];
    c_ << beta[3], beta[4], beta[5];

    // 利用c与c_设置其它运动变量
    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    // 对坐标系的基执行匀速螺旋运动
    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    // 利用beta[beta_index],beta[beta_index+1]设置更新后的二维坐标
    double px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    double py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    // 更新m_osnap_poly.vertices[v_index]
    // b0, b1, b2, px, py 均已被更新
    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    // 计算f7(beta)
    double f7_val = f7(v_index, pv);	// 没有改变任何beta变量的f7值
    double f7_delta_val;								// 改变一个beta变量的f7值
    double pd;											// 偏导数

    // beta[0]///////////////////////////////////////////////////////////////////////
    c << beta[0] + delta, beta[1], beta[2];
    c_ << beta[3], beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 0) = pd;

    // beta[1]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1] + delta, beta[2];
    c_ << beta[3], beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 1) = pd;

    // beta[2]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2] + delta;
    c_ << beta[3], beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 2) = pd;

    // beta[3]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2];
    c_ << beta[3] + delta, beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 3) = pd;

    // beta[4]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2];
    c_ << beta[3], beta[4] + delta, beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 4) = pd;

    // beta[5]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2];
    c_ << beta[3], beta[4], beta[5] + delta;

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, 5) = pd;

    // beta[beta_index]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2];
    c_ << beta[3], beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index] + delta;
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1];

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, beta_index) = pd;

    // beta[beta_index+1]///////////////////////////////////////////////////////////////////////
    c << beta[0], beta[1], beta[2];
    c_ << beta[3], beta[4], beta[5];

    setHelicalVars(c, c_, R, pointOfA, a, a_, pitch, alpha);

    uniformHelicalMove(c, c_, m_osnap_poly.b0, b0, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b1, b1, R, pointOfA, a, pitch, alpha);
    uniformHelicalMove(c, c_, m_osnap_poly.b2, b2, R, pointOfA, a, pitch, alpha);

    px = m_osnap_poly.coors[v_index].px + beta[beta_index];
    py = m_osnap_poly.coors[v_index].py + beta[beta_index+1] + delta;

    m_osnap_poly.vertices[v_index] = b0 + px * (b1 - b0) + py * (b2 - b0);

    f7_delta_val = f7(v_index, pv);
    pd = (f7_delta_val - f7_val) / delta;
    J(row, beta_index+1) = pd;
}









