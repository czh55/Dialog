#ifndef HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H
#define HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H

#include <QDialog>
#include "HoleFilling.h"

struct Vertex_info
{
    int poly_id;    // 多边形编号
    Vertex2* v;     // 顶点位置
};

namespace Ui {
class HoleFillingInteractiveOperationDialog;
}

class HoleFillingInteractiveOperationDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit HoleFillingInteractiveOperationDialog(QWidget *parent = 0);
    ~HoleFillingInteractiveOperationDialog();

    void * m_pclviewer;
    char * m_CMD;           // 命令行缓冲器

    int selected_poly_id;   // 选中的多边形id
    Vertex2 * m_vertex1;    // 选中的顶点1
    Vertex2 * m_vertex2;    // 选中的顶点2

    vector<Vertex2*> line_segs;  // 顶点1与顶点2构成的线段序列

    void updateInfo();

    ////////////////////////////////////////////////////////
    // 构造三角化点集
    vector<Vertex2*> m_points_set_for_tri;  // 最初利用模式(1,2,3)选择的一个点集，它框定了进行三角化需要检测的点与线段的范围
    // 选择孔洞中的一个顶点，从而识别出孔洞边界
    Vertex_info m_vi_one_vertex;
    Eigen::Vector3d m_camera_pos;
    // 模式1：仅利用一个顶点生成孔洞边界点集
    bool genPointsSetWithOneVertex();
	// Input: v: 当前顶点；poly：v的隶属多边形
	// Output: next_v: 下一个顶点；next_poly：next_v的隶属多边形
	// 应用该方法的前提是需要保证目标多边形的v的前驱节点务必在边界上！
    bool findNextPoint(Vertex2* v, POLYGON2* poly, Vertex2* &next_v, POLYGON2* &next_poly);
	// 点向平面的投影
	void getProjPointOnPlane(Eigen::Vector3d &v_pos, Eigen::Vector4d &coeff, Eigen::Vector3d &proj_pos);
	// e2顺时针旋转到e1转动的角度
	double getRotationAngle(Eigen::Vector3d &e1, Eigen::Vector3d &e2, Eigen::Vector3d &normal);
    // 模式2： 利用一个起点与一个终点生成孔洞边界点集
    Vertex_info m_v_begin_mode2;
    Vertex_info m_v_end_mode2;
    bool genPointsSetWithTwoVertices();
    // 模式3：利用四个顶点生成孔洞边界点集，适用于连接两个分离多边形的情况
    Vertex_info m_v1_mode3;
    Vertex_info m_v2_mode3;
    Vertex_info m_v3_mode3;
    Vertex_info m_v4_mode3;
    bool genPointsSetWithFourVertices();

	// 三角化
	// 获取点到线段的距离
	double getDistP2LS(Vertex2* l_begin, Vertex2* l_end, Vertex2* v);
	// 获取指定线段到某个点集中每个顶点的距离，并按照从小到大的顺序进行排列
	void getDistInfo(Vertex2* l_begin, Vertex2 *l_end, vector<Vertex2*> &points, vector<int> &indices, vector<float> &dists);
	// 三角形合法性检查
	bool isValidTriangle(Vertex2* l_begin, Vertex2* l_end, Vertex2* v, 
	                     vector<Vertex2*> &m_points_set_for_tri, int poly_start_index,
	                     Eigen::Vector3d &normal);
	// 判断顶点是否在线段上
	bool isVertexOnLingSeg(Vertex2* l_begin, Vertex2* l_end, Vertex2* v);
	// 判断顶点是否在三角形内部
	bool isVertexInTriangle(Vertex2* l_begin, Vertex2* l_end, Vertex2* v, Vertex2* test_point);
	// 两个线段是否平行
	bool isTwoLineSegsParal(Vertex2* l1_begin, Vertex2* l1_end, Vertex2* l2_begin, Vertex2* l2_end);
    // 两个线段是否相交(首位相连不算相交)
	bool isTwoLineSegsIntersect(Vertex2* l1_begin, Vertex2* l1_end, Vertex2* l2_begin, Vertex2* l2_end);
	// 利用三个顶点构造新的三角形多边形
	void genNewTriPly(Vertex2* v0, Vertex2* v1, Vertex2* v2, Eigen::Vector3d &normal);

    // 水平方向，垂直方向延伸模型边界
    // 投影的线段序列是否合法
    // 有问题的线段索引: problem_i, problem_j
    bool isProjLineSegsValid(vector<Vertex2*> &line_segs, int &problem_i, int &problem_j);
    // 两个线段的方向向量是否平行(i,j均为线段的起点)
    bool isBothLineSegsParal(vector<Vertex2*> &line_segs, int i, int j);
    // 两条线段是否部分重叠
    bool isBothLineSegsOverlap(vector<Vertex2*> &line_segs, int i, int j);
    // 线段i与线段i+1是否具有合法的位置关系
    bool isAdjacentLineSegsHasValidPosRelation(vector<Vertex2*> &line_segs, int i);
	//调用本函数时默认模型内部已经不存在孔洞了
	// input：v，模型边界上的顶点
	// output: next_v, 模型边界上的顶点
	bool findNextPointVertexOnModelBorder(Vertex2* v, Vertex2* &next_v);
	// 判断顶点v是否具有单向的连接关系
	bool isVertexHasSingleConnRelation(Vertex2*v);
    // 执行函数
    void performExtention(int extend_type);
    // 模型封闭的数据结构
    vector<vector<Vertex2*>> m_four_vertices;
    Vertex2* m_bottom_left_vertex;
    Vertex2* m_bottom_right_vertex;
    double m_min_y;
    double m_min_x;
    double m_max_x;
    double m_min_z;
    int plane_type; // 0:底面；1：左侧立面;2:右侧立面；3：后侧立面

	// 视图更新函数
	void updateViewer();// 全局更新

private slots:
    void on_pushButton_select_poly_clicked();

    void on_pushButton_simplify_cancel_poly_selection_clicked();

    void on_pushButton_select_p1_clicked();

    void on_pushButton_select_p2_clicked();

    void on_pushButton_gener_line_seg_clicked();

    void on_pushButton_del_lineseg_clicked();

    void on_pushButton_switch_ling_seg_clicked();

    void on_pushButton_simplify_line_segs_clicked();

    void on_pushButton_find_prev_point_clicked();

    void on_pushButton_find_next_point_clicked();

    void on_pushButton_search_angle2angle_vertices_clicked();

    void on_pushButton_select_one_vertex_clicked();

    void on_pushButton_get_camera_pos_clicked();

    void on_pushButton_gen_points_set_for_tri_clicked();

    void on_pushButton_test_findNextPoint_clicked();

    void on_pushButton_select_two_vertices_clicked();

    void on_pushButton_sel_v1_for_mode2_clicked();

    void on_pushButton_sel_v2_for_mode2_clicked();

    void on_pushButton_select_four_vertices_clicked();

    void on_pushButton_sel_v1_for_mode3_clicked();

    void on_pushButton_sel_v2_for_mode3_clicked();

    void on_pushButton_sel_v3_for_mode3_clicked();

    void on_pushButton_sel_v4_for_mode3_clicked();

    void on_pushButton_test_clicked();

    void on_pushButton_set_camera_pos_clicked();

    void on_pushButton_get_poly_num_clicked();

    void on_pushButton_get_ver_num_clicked();

    void on_pushButton_cancel_poly_wireframe_display_clicked();

    void on_pushButton_display_poly_wireframe_clicked();

    void on_pushButton_gen_proj_line_segs_clicked();

    void on_pushButton_extend_to_bottom_clicked();

    void on_pushButton_extend_to_left_clicked();

    void on_pushButton_extend_to_right_clicked();

    void on_pushButton_extend_to_backward_clicked();

    //void on_pushButton_close_the_model_clicked();

    //void on_pushButton_gen_the_two_border_vertices_clicked();

    //void on_pushButton_gen_bottom_plane_clicked();

    //void on_pushButton_gen_left_plane_clicked();

    //void on_pushButton_gen_right_plane_clicked();

    //void on_pushButton_gen_backward_plane_clicked();

    void on_pushButton_clicked();

    void on_pushButton_del_poly_clicked();

    void on_pushButton_CDT_Hole_Filling_clicked();

    void on_pushButton_dist_V2V_clicked();

    void on_pushButton_test_ply_op_clicked();

    void on_pushButton_insert_vertices_clicked();

    void on_pushButton_gen_two_key_vertices_clicked();

    void on_pushButton_gen_bottom_plane_clicked();

    void on_pushButton__insert_interior_vertices_clicked();

    void on_pushButton_display_ply_normal_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_3_clicked();

private:
    Ui::HoleFillingInteractiveOperationDialog *ui;    
};

#endif // HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H
