#ifndef OSNAP_H
#define OSNAP_H

#include "snappingvertices.h"

struct Coor_sigma	// 移动坐标系的二维坐标
{
    double px;
    double py;
};

struct OSnapPoly
{
    Eigen::Vector3d c;		// 角速度向量
    Eigen::Vector3d c_;		// 平移速度向量
    Eigen::Vector3d b0;		// 移动坐标系的原点
    Eigen::Vector3d b1;		// 移动坐标系基e1的端点
    Eigen::Vector3d b2;		// 移动坐标系基e2的端点
    Eigen::Vector4d coeff;			// 当前平面方程的系数
    Eigen::Vector4d coeff_origin;	// 原始平面方程的系数
	vector<Eigen::Vector3d> vertices;			// 平面多边形的当前顶点
    vector<Eigen::Vector3d> vertices_origin;	// 平面多边形的原始顶点
    std::vector<Coor_sigma> coors;							// 当前平面多边形中每个顶点在移动坐标系的二维坐标
};

class OSnap
{
public:
    OSnap();

	vector<POLYGON, Eigen::aligned_allocator<POLYGON>>* m_polygons;	// 输入的多边形对象集合
    POLYGON *m_poly;    // 输入的多边形对象
    OSnapPoly m_osnap_poly;     // 目标多边形对象

    vector<Vertex*> normal_vertices;    // 没有匹配关系的顶点集合
    vector<Vertex*> p2p_vertices;       // 具有点-点匹配关系的顶点集合
    vector<Vertex*> p2e_vertices;       // 具有点-边匹配关系的顶点集合

	// 在原始点云中的索引
	vector<int> normal_vertices_indices;
	vector<int> p2p_vertices_indices;
	vector<int> p2e_vertices_indices;

    Eigen::MatrixXd J;					// Jacob矩阵
    Eigen::VectorXd beta;				// 变量向量

    double delta;   // 利用数值方法求偏导的步长
    //double max_r_norm;  // 共轭梯度解方程,向量r的最大范数
	int max_k_count;	// 共轭梯度解方程的最大迭代次数

    // 重置数据结构
    void reset();

    // 根据输入的多边形对象m_poly设置目标多边形对象m_osnap_poly
    void setOSnapPoly();
    // 利用顶点构造初始的移动坐标系
    void setBases_init(OSnapPoly &poly);
    // 设置多边形顶点在移动坐标系中的二维坐标
    void set2DCoordinates(OSnapPoly &poly);

    // 根据匹配关系分类顶点
    void categoriseVertices();
	
    // 匀速螺旋运动
    // 设置变量
    void setHelicalVars( Eigen::Vector3d &c,
                         Eigen::Vector3d &c_,
                         Eigen::Matrix3d &R,
                         Eigen::Vector3d &pointOfA,
                         Eigen::Vector3d &a,
                         Eigen::Vector3d &a_,
                         double &pitch,
                         double &alpha );
    // 执行移动
    void uniformHelicalMove( Eigen::Vector3d &c,
                             Eigen::Vector3d &c_,
                             Eigen::Vector3d &x,
                             Eigen::Vector3d &x_,
                             Eigen::Matrix3d &R,
                             Eigen::Vector3d &pointOfA,
                             Eigen::Vector3d &a,
                             double &pitch,
                             double &alpha );
    // 更新经过螺旋运动后的移动坐标系的基
    void updateBasesAfterHelicalMove(OSnapPoly &poly);
    // 更新经过螺旋运动后顶点的三维坐标
    // (前提：‘移动坐标系的基’+‘二维坐标’已被更新过！)
    void updateVerticesAfterHelicalMove(OSnapPoly &poly);

	// 根据变量向量beta更新每个顶点在移动坐标系中的二维坐标
	void updateCoors();

	// 计算两个多边形的交线
	bool getIntersectLine(Eigen::Vector4f &coeff1, Eigen::Vector4f &coeff2,
                          Eigen::Vector3d &p_base, Eigen::Vector3d &line_dir);

    // 目标函数
	// 三种顶点共用
    // 更新后顶点的位置与原来位置的距离
	double f1(int vertex_index);

	// 点-点
	// 更新后顶点位置与目标顶点位置的距离
	double f2(int vertex_index_in_p2p_vertices);
	// dist[更新后顶点，目标顶点的前驱节点] - dist[目标顶点，目标顶点的前驱节点]
	double f3_1(int vertex_index_in_p2p_vertices);
	// dist[更新后顶点，目标顶点的后继节点] - dist[目标顶点，目标顶点的后继节点]
	double f3_2(int vertex_index_in_p2p_vertices);
	// f5的对偶版本，更新后顶点位置与两个平面交线的距离
	double f8(int vertex_index_in_p2p_vertices);

	// 点-边
	// 更新后顶点位置到匹配边的距离
	double f4(int vertex_index_in_p2e_vertices);
	// 更新后顶点位置与两个平面交线的距离
	double f5(int vertex_index_in_p2e_vertices);
	// 更新后顶点到匹配边两个端点距离之和 - 匹配边长度
	double f6(int vertex_index_in_p2e_vertices);

	// 有匹配关系的顶点
    // 更新后顶点到匹配关系所属支撑平面的距离
    double f7(int vertex_index, Vertex *pv);

	// 目标函数的权值
	double f1_weight;
	double f2_weight;
	double f3_1_weight;
	double f3_2_weight;
	double f8_weight;
	double f4_weight;
	double f5_weight;
	double f6_weight;
	double f7_weight;

	// 目标函数的偏导数
	// 针对没有匹配关系的顶点：normal_vertices
	void pd_f1_normal_vertices(int v_index, Eigen::MatrixXd &J, int row, double delta);

	// 针对有匹配关系的顶点: p2p_vertices 和 p2e_vertices
	void pd_f1_matching_vertices(int v_index, Eigen::MatrixXd &J, int row, int beta_index, double delta);

	// 针对点对点
	// f2的偏导数
	void pd_f2(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);
	// f3_1的偏导数
	void pd_f3_1(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);
	// f3_2的偏导数
	void pd_f3_2(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);
	// f8的偏导数
	void pd_f8(int vertex_index_in_p2p_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);

	// 针对点对边
	// f4的偏导数
	void pd_f4(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);
	// f5的偏导数
	void pd_f5(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);
	// f6的偏导数
	void pd_f6(int vertex_index_in_p2e_vertices, Eigen::MatrixXd &J, int row, int beta_index, double delta);

    // 针对有匹配关系的顶点
    // f7的偏导数
    void pd_f7(int v_index, Vertex *pv, Eigen::MatrixXd &J, int row, int beta_index, double delta);

	// 优化	
	// 初始化变量向量beta,在每次优化的开始使用一次
	void initBeta();
	// 构造Jacob矩阵
	void cosntructJacobMat();
	// 执行单步高斯牛顿迭代，从而更新beta
	void fill_r(Eigen::VectorXd &r);
	void conjugateGradient(Eigen::MatrixXd &A, Eigen::VectorXd &x, Eigen::VectorXd &b);
	void conjugateGradient_with_max_k(Eigen::MatrixXd &A, Eigen::VectorXd &x, Eigen::VectorXd &b, int k_count);
	void performGN();
	// 更新移动坐标系
	void updateMovingCoor();
	// 更新顶点位置
	void updateVertices();

	// 最终的优化执行函数
	void optimise();
};

#endif // OSNAP_H
