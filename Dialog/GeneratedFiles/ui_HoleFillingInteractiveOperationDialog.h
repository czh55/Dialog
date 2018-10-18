/********************************************************************************
** Form generated from reading UI file 'HoleFillingInteractiveOperationDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H
#define UI_HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_HoleFillingInteractiveOperationDialog
{
public:
    QPushButton *pushButton_simplify_line_segs;
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton_select_poly;
    QLabel *label_poly_id;
    QLineEdit *lineEdit_poly_id;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *pushButton_select_p1;
    QLabel *label_p1;
    QLineEdit *lineEdit_p1;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *pushButton_select_p2;
    QLabel *label_p2;
    QLineEdit *lineEdit_p2;
    QPushButton *pushButton_simplify_cancel_poly_selection;
    QPushButton *pushButton_gener_line_seg;
    QPushButton *pushButton_switch_ling_seg;
    QPushButton *pushButton_del_lineseg;
    QPushButton *pushButton_find_prev_point;
    QPushButton *pushButton_find_next_point;
    QPushButton *pushButton_search_angle2angle_vertices;
    QFrame *line;
    QFrame *line_3;
    QPushButton *pushButton_gen_points_set_for_tri;
    QWidget *layoutWidget1;
    QHBoxLayout *horizontalLayout_9;
    QFrame *line_4;
    QLabel *label;
    QFrame *line_2;
    QWidget *layoutWidget2;
    QHBoxLayout *horizontalLayout_10;
    QPushButton *pushButton_get_camera_pos;
    QLineEdit *lineEdit_camera_pos;
    QPushButton *pushButton_test_findNextPoint;
    QFrame *line_7;
    QFrame *line_8;
    QWidget *layoutWidget3;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *pushButton_select_one_vertex;
    QRadioButton *radioButton_select_one_vertex;
    QLineEdit *lineEdit_one_vertex;
    QFrame *line_5;
    QHBoxLayout *horizontalLayout_5;
    QPushButton *pushButton_select_two_vertices;
    QRadioButton *radioButton_select_two_vertices;
    QPushButton *pushButton_sel_v1_for_mode2;
    QLineEdit *lineEdit_two_vertices_1;
    QLabel *label_2;
    QPushButton *pushButton_sel_v2_for_mode2;
    QLineEdit *lineEdit_two_vertices_2;
    QFrame *line_6;
    QHBoxLayout *horizontalLayout_8;
    QPushButton *pushButton_select_four_vertices;
    QRadioButton *radioButton_select_four_vertices;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_6;
    QPushButton *pushButton_sel_v1_for_mode3;
    QLineEdit *lineEdit_three_vertices_1;
    QLabel *label_3;
    QPushButton *pushButton_sel_v2_for_mode3;
    QLineEdit *lineEdit_three_vertices_2;
    QHBoxLayout *horizontalLayout_7;
    QPushButton *pushButton_sel_v3_for_mode3;
    QLineEdit *lineEdit_three_vertices_3;
    QLabel *label_4;
    QPushButton *pushButton_sel_v4_for_mode3;
    QLineEdit *lineEdit_three_vertices_4;
    QPushButton *pushButton_test;
    QPushButton *pushButton_set_camera_pos;
    QPushButton *pushButton_gen_proj_line_segs;
    QPushButton *pushButton_cancel_poly_wireframe_display;
    QPushButton *pushButton_display_poly_wireframe;
    QPushButton *pushButton_get_ver_num;
    QPushButton *pushButton_get_poly_num;
    QPushButton *pushButton_extend_to_bottom;
    QPushButton *pushButton_extend_to_backward;
    QPushButton *pushButton_extend_to_left;
    QPushButton *pushButton_extend_to_right;
    QPushButton *pushButton;
    QPushButton *pushButton_del_poly;
    QPushButton *pushButton_CDT_Hole_Filling;
    QLabel *label_5;
    QPushButton *pushButton_dist_V2V;
    QPushButton *pushButton_test_ply_op;
    QFrame *line_9;
    QPushButton *pushButton_insert_vertices;
    QFrame *line_10;
    QPushButton *pushButton_gen_two_key_vertices;
    QPushButton *pushButton_gen_bottom_plane;
    QPushButton *pushButton_4;
    QPushButton *pushButton_5;
    QPushButton *pushButton_6;
    QPushButton *pushButton__insert_interior_vertices;
    QPushButton *pushButton_display_ply_normal;
    QPushButton *pushButton_2;
    QPushButton *pushButton_3;

    void setupUi(QDialog *HoleFillingInteractiveOperationDialog)
    {
        if (HoleFillingInteractiveOperationDialog->objectName().isEmpty())
            HoleFillingInteractiveOperationDialog->setObjectName(QStringLiteral("HoleFillingInteractiveOperationDialog"));
        HoleFillingInteractiveOperationDialog->resize(746, 662);
        HoleFillingInteractiveOperationDialog->setSizeGripEnabled(false);
        HoleFillingInteractiveOperationDialog->setModal(false);
        pushButton_simplify_line_segs = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_simplify_line_segs->setObjectName(QStringLiteral("pushButton_simplify_line_segs"));
        pushButton_simplify_line_segs->setGeometry(QRect(210, 110, 75, 41));
        layoutWidget = new QWidget(HoleFillingInteractiveOperationDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(93, 10, 501, 89));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        pushButton_select_poly = new QPushButton(layoutWidget);
        pushButton_select_poly->setObjectName(QStringLiteral("pushButton_select_poly"));

        horizontalLayout->addWidget(pushButton_select_poly);

        label_poly_id = new QLabel(layoutWidget);
        label_poly_id->setObjectName(QStringLiteral("label_poly_id"));

        horizontalLayout->addWidget(label_poly_id);

        lineEdit_poly_id = new QLineEdit(layoutWidget);
        lineEdit_poly_id->setObjectName(QStringLiteral("lineEdit_poly_id"));
        lineEdit_poly_id->setEnabled(false);

        horizontalLayout->addWidget(lineEdit_poly_id);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        pushButton_select_p1 = new QPushButton(layoutWidget);
        pushButton_select_p1->setObjectName(QStringLiteral("pushButton_select_p1"));

        horizontalLayout_2->addWidget(pushButton_select_p1);

        label_p1 = new QLabel(layoutWidget);
        label_p1->setObjectName(QStringLiteral("label_p1"));

        horizontalLayout_2->addWidget(label_p1);

        lineEdit_p1 = new QLineEdit(layoutWidget);
        lineEdit_p1->setObjectName(QStringLiteral("lineEdit_p1"));
        lineEdit_p1->setEnabled(false);

        horizontalLayout_2->addWidget(lineEdit_p1);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        pushButton_select_p2 = new QPushButton(layoutWidget);
        pushButton_select_p2->setObjectName(QStringLiteral("pushButton_select_p2"));

        horizontalLayout_3->addWidget(pushButton_select_p2);

        label_p2 = new QLabel(layoutWidget);
        label_p2->setObjectName(QStringLiteral("label_p2"));

        horizontalLayout_3->addWidget(label_p2);

        lineEdit_p2 = new QLineEdit(layoutWidget);
        lineEdit_p2->setObjectName(QStringLiteral("lineEdit_p2"));
        lineEdit_p2->setEnabled(false);

        horizontalLayout_3->addWidget(lineEdit_p2);


        verticalLayout->addLayout(horizontalLayout_3);

        pushButton_simplify_cancel_poly_selection = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_simplify_cancel_poly_selection->setObjectName(QStringLiteral("pushButton_simplify_cancel_poly_selection"));
        pushButton_simplify_cancel_poly_selection->setGeometry(QRect(10, 60, 81, 41));
        pushButton_gener_line_seg = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_gener_line_seg->setObjectName(QStringLiteral("pushButton_gener_line_seg"));
        pushButton_gener_line_seg->setGeometry(QRect(10, 110, 75, 41));
        pushButton_switch_ling_seg = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_switch_ling_seg->setObjectName(QStringLiteral("pushButton_switch_ling_seg"));
        pushButton_switch_ling_seg->setGeometry(QRect(110, 110, 75, 41));
        pushButton_del_lineseg = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_del_lineseg->setObjectName(QStringLiteral("pushButton_del_lineseg"));
        pushButton_del_lineseg->setGeometry(QRect(320, 110, 75, 41));
        pushButton_find_prev_point = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_find_prev_point->setObjectName(QStringLiteral("pushButton_find_prev_point"));
        pushButton_find_prev_point->setGeometry(QRect(420, 110, 75, 41));
        pushButton_find_next_point = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_find_next_point->setObjectName(QStringLiteral("pushButton_find_next_point"));
        pushButton_find_next_point->setGeometry(QRect(520, 110, 75, 41));
        pushButton_search_angle2angle_vertices = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_search_angle2angle_vertices->setObjectName(QStringLiteral("pushButton_search_angle2angle_vertices"));
        pushButton_search_angle2angle_vertices->setGeometry(QRect(0, 160, 111, 41));
        line = new QFrame(HoleFillingInteractiveOperationDialog);
        line->setObjectName(QStringLiteral("line"));
        line->setGeometry(QRect(0, 150, 731, 16));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        line_3 = new QFrame(HoleFillingInteractiveOperationDialog);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setGeometry(QRect(0, 410, 731, 16));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);
        pushButton_gen_points_set_for_tri = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_gen_points_set_for_tri->setObjectName(QStringLiteral("pushButton_gen_points_set_for_tri"));
        pushButton_gen_points_set_for_tri->setGeometry(QRect(0, 430, 111, 41));
        layoutWidget1 = new QWidget(HoleFillingInteractiveOperationDialog);
        layoutWidget1->setObjectName(QStringLiteral("layoutWidget1"));
        layoutWidget1->setGeometry(QRect(0, 200, 731, 20));
        horizontalLayout_9 = new QHBoxLayout(layoutWidget1);
        horizontalLayout_9->setObjectName(QStringLiteral("horizontalLayout_9"));
        horizontalLayout_9->setContentsMargins(0, 0, 0, 0);
        line_4 = new QFrame(layoutWidget1);
        line_4->setObjectName(QStringLiteral("line_4"));
        line_4->setFrameShape(QFrame::HLine);
        line_4->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line_4);

        label = new QLabel(layoutWidget1);
        label->setObjectName(QStringLiteral("label"));
        label->setAlignment(Qt::AlignCenter);

        horizontalLayout_9->addWidget(label);

        line_2 = new QFrame(layoutWidget1);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);

        horizontalLayout_9->addWidget(line_2);

        layoutWidget2 = new QWidget(HoleFillingInteractiveOperationDialog);
        layoutWidget2->setObjectName(QStringLiteral("layoutWidget2"));
        layoutWidget2->setGeometry(QRect(0, 380, 381, 25));
        horizontalLayout_10 = new QHBoxLayout(layoutWidget2);
        horizontalLayout_10->setObjectName(QStringLiteral("horizontalLayout_10"));
        horizontalLayout_10->setContentsMargins(0, 0, 0, 0);
        pushButton_get_camera_pos = new QPushButton(layoutWidget2);
        pushButton_get_camera_pos->setObjectName(QStringLiteral("pushButton_get_camera_pos"));

        horizontalLayout_10->addWidget(pushButton_get_camera_pos);

        lineEdit_camera_pos = new QLineEdit(layoutWidget2);
        lineEdit_camera_pos->setObjectName(QStringLiteral("lineEdit_camera_pos"));

        horizontalLayout_10->addWidget(lineEdit_camera_pos);

        pushButton_test_findNextPoint = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_test_findNextPoint->setObjectName(QStringLiteral("pushButton_test_findNextPoint"));
        pushButton_test_findNextPoint->setGeometry(QRect(540, 430, 101, 41));
        line_7 = new QFrame(HoleFillingInteractiveOperationDialog);
        line_7->setObjectName(QStringLiteral("line_7"));
        line_7->setGeometry(QRect(0, 360, 731, 16));
        line_7->setFrameShape(QFrame::HLine);
        line_7->setFrameShadow(QFrame::Sunken);
        line_8 = new QFrame(HoleFillingInteractiveOperationDialog);
        line_8->setObjectName(QStringLiteral("line_8"));
        line_8->setGeometry(QRect(0, 540, 731, 16));
        line_8->setFrameShape(QFrame::HLine);
        line_8->setFrameShadow(QFrame::Sunken);
        layoutWidget3 = new QWidget(HoleFillingInteractiveOperationDialog);
        layoutWidget3->setObjectName(QStringLiteral("layoutWidget3"));
        layoutWidget3->setGeometry(QRect(0, 223, 731, 142));
        verticalLayout_3 = new QVBoxLayout(layoutWidget3);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        pushButton_select_one_vertex = new QPushButton(layoutWidget3);
        pushButton_select_one_vertex->setObjectName(QStringLiteral("pushButton_select_one_vertex"));

        horizontalLayout_4->addWidget(pushButton_select_one_vertex);

        radioButton_select_one_vertex = new QRadioButton(layoutWidget3);
        radioButton_select_one_vertex->setObjectName(QStringLiteral("radioButton_select_one_vertex"));

        horizontalLayout_4->addWidget(radioButton_select_one_vertex);

        lineEdit_one_vertex = new QLineEdit(layoutWidget3);
        lineEdit_one_vertex->setObjectName(QStringLiteral("lineEdit_one_vertex"));
        lineEdit_one_vertex->setEnabled(true);
        lineEdit_one_vertex->setReadOnly(true);

        horizontalLayout_4->addWidget(lineEdit_one_vertex);


        verticalLayout_3->addLayout(horizontalLayout_4);

        line_5 = new QFrame(layoutWidget3);
        line_5->setObjectName(QStringLiteral("line_5"));
        line_5->setFrameShape(QFrame::HLine);
        line_5->setFrameShadow(QFrame::Sunken);

        verticalLayout_3->addWidget(line_5);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        pushButton_select_two_vertices = new QPushButton(layoutWidget3);
        pushButton_select_two_vertices->setObjectName(QStringLiteral("pushButton_select_two_vertices"));

        horizontalLayout_5->addWidget(pushButton_select_two_vertices);

        radioButton_select_two_vertices = new QRadioButton(layoutWidget3);
        radioButton_select_two_vertices->setObjectName(QStringLiteral("radioButton_select_two_vertices"));

        horizontalLayout_5->addWidget(radioButton_select_two_vertices);

        pushButton_sel_v1_for_mode2 = new QPushButton(layoutWidget3);
        pushButton_sel_v1_for_mode2->setObjectName(QStringLiteral("pushButton_sel_v1_for_mode2"));

        horizontalLayout_5->addWidget(pushButton_sel_v1_for_mode2);

        lineEdit_two_vertices_1 = new QLineEdit(layoutWidget3);
        lineEdit_two_vertices_1->setObjectName(QStringLiteral("lineEdit_two_vertices_1"));
        lineEdit_two_vertices_1->setReadOnly(true);

        horizontalLayout_5->addWidget(lineEdit_two_vertices_1);

        label_2 = new QLabel(layoutWidget3);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_5->addWidget(label_2);

        pushButton_sel_v2_for_mode2 = new QPushButton(layoutWidget3);
        pushButton_sel_v2_for_mode2->setObjectName(QStringLiteral("pushButton_sel_v2_for_mode2"));

        horizontalLayout_5->addWidget(pushButton_sel_v2_for_mode2);

        lineEdit_two_vertices_2 = new QLineEdit(layoutWidget3);
        lineEdit_two_vertices_2->setObjectName(QStringLiteral("lineEdit_two_vertices_2"));
        lineEdit_two_vertices_2->setReadOnly(true);

        horizontalLayout_5->addWidget(lineEdit_two_vertices_2);


        verticalLayout_3->addLayout(horizontalLayout_5);

        line_6 = new QFrame(layoutWidget3);
        line_6->setObjectName(QStringLiteral("line_6"));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);

        verticalLayout_3->addWidget(line_6);

        horizontalLayout_8 = new QHBoxLayout();
        horizontalLayout_8->setObjectName(QStringLiteral("horizontalLayout_8"));
        pushButton_select_four_vertices = new QPushButton(layoutWidget3);
        pushButton_select_four_vertices->setObjectName(QStringLiteral("pushButton_select_four_vertices"));

        horizontalLayout_8->addWidget(pushButton_select_four_vertices);

        radioButton_select_four_vertices = new QRadioButton(layoutWidget3);
        radioButton_select_four_vertices->setObjectName(QStringLiteral("radioButton_select_four_vertices"));

        horizontalLayout_8->addWidget(radioButton_select_four_vertices);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        pushButton_sel_v1_for_mode3 = new QPushButton(layoutWidget3);
        pushButton_sel_v1_for_mode3->setObjectName(QStringLiteral("pushButton_sel_v1_for_mode3"));

        horizontalLayout_6->addWidget(pushButton_sel_v1_for_mode3);

        lineEdit_three_vertices_1 = new QLineEdit(layoutWidget3);
        lineEdit_three_vertices_1->setObjectName(QStringLiteral("lineEdit_three_vertices_1"));
        lineEdit_three_vertices_1->setReadOnly(true);

        horizontalLayout_6->addWidget(lineEdit_three_vertices_1);

        label_3 = new QLabel(layoutWidget3);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_6->addWidget(label_3);

        pushButton_sel_v2_for_mode3 = new QPushButton(layoutWidget3);
        pushButton_sel_v2_for_mode3->setObjectName(QStringLiteral("pushButton_sel_v2_for_mode3"));

        horizontalLayout_6->addWidget(pushButton_sel_v2_for_mode3);

        lineEdit_three_vertices_2 = new QLineEdit(layoutWidget3);
        lineEdit_three_vertices_2->setObjectName(QStringLiteral("lineEdit_three_vertices_2"));
        lineEdit_three_vertices_2->setReadOnly(true);

        horizontalLayout_6->addWidget(lineEdit_three_vertices_2);


        verticalLayout_2->addLayout(horizontalLayout_6);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        pushButton_sel_v3_for_mode3 = new QPushButton(layoutWidget3);
        pushButton_sel_v3_for_mode3->setObjectName(QStringLiteral("pushButton_sel_v3_for_mode3"));

        horizontalLayout_7->addWidget(pushButton_sel_v3_for_mode3);

        lineEdit_three_vertices_3 = new QLineEdit(layoutWidget3);
        lineEdit_three_vertices_3->setObjectName(QStringLiteral("lineEdit_three_vertices_3"));
        lineEdit_three_vertices_3->setReadOnly(true);

        horizontalLayout_7->addWidget(lineEdit_three_vertices_3);

        label_4 = new QLabel(layoutWidget3);
        label_4->setObjectName(QStringLiteral("label_4"));

        horizontalLayout_7->addWidget(label_4);

        pushButton_sel_v4_for_mode3 = new QPushButton(layoutWidget3);
        pushButton_sel_v4_for_mode3->setObjectName(QStringLiteral("pushButton_sel_v4_for_mode3"));

        horizontalLayout_7->addWidget(pushButton_sel_v4_for_mode3);

        lineEdit_three_vertices_4 = new QLineEdit(layoutWidget3);
        lineEdit_three_vertices_4->setObjectName(QStringLiteral("lineEdit_three_vertices_4"));
        lineEdit_three_vertices_4->setReadOnly(true);

        horizontalLayout_7->addWidget(lineEdit_three_vertices_4);


        verticalLayout_2->addLayout(horizontalLayout_7);


        horizontalLayout_8->addLayout(verticalLayout_2);


        verticalLayout_3->addLayout(horizontalLayout_8);

        pushButton_test = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_test->setObjectName(QStringLiteral("pushButton_test"));
        pushButton_test->setGeometry(QRect(450, 430, 81, 41));
        pushButton_set_camera_pos = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_set_camera_pos->setObjectName(QStringLiteral("pushButton_set_camera_pos"));
        pushButton_set_camera_pos->setGeometry(QRect(400, 380, 121, 23));
        pushButton_gen_proj_line_segs = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_gen_proj_line_segs->setObjectName(QStringLiteral("pushButton_gen_proj_line_segs"));
        pushButton_gen_proj_line_segs->setGeometry(QRect(0, 550, 111, 41));
        pushButton_cancel_poly_wireframe_display = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_cancel_poly_wireframe_display->setObjectName(QStringLiteral("pushButton_cancel_poly_wireframe_display"));
        pushButton_cancel_poly_wireframe_display->setGeometry(QRect(140, 160, 121, 41));
        pushButton_display_poly_wireframe = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_display_poly_wireframe->setObjectName(QStringLiteral("pushButton_display_poly_wireframe"));
        pushButton_display_poly_wireframe->setGeometry(QRect(290, 160, 121, 41));
        pushButton_get_ver_num = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_get_ver_num->setObjectName(QStringLiteral("pushButton_get_ver_num"));
        pushButton_get_ver_num->setGeometry(QRect(440, 160, 61, 41));
        pushButton_get_poly_num = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_get_poly_num->setObjectName(QStringLiteral("pushButton_get_poly_num"));
        pushButton_get_poly_num->setGeometry(QRect(530, 160, 81, 41));
        pushButton_extend_to_bottom = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_extend_to_bottom->setObjectName(QStringLiteral("pushButton_extend_to_bottom"));
        pushButton_extend_to_bottom->setGeometry(QRect(220, 620, 81, 31));
        pushButton_extend_to_backward = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_extend_to_backward->setObjectName(QStringLiteral("pushButton_extend_to_backward"));
        pushButton_extend_to_backward->setGeometry(QRect(220, 560, 81, 31));
        pushButton_extend_to_left = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_extend_to_left->setObjectName(QStringLiteral("pushButton_extend_to_left"));
        pushButton_extend_to_left->setGeometry(QRect(140, 590, 75, 31));
        pushButton_extend_to_right = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_extend_to_right->setObjectName(QStringLiteral("pushButton_extend_to_right"));
        pushButton_extend_to_right->setGeometry(QRect(300, 590, 75, 31));
        pushButton = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setGeometry(QRect(10, 10, 81, 41));
        pushButton_del_poly = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_del_poly->setObjectName(QStringLiteral("pushButton_del_poly"));
        pushButton_del_poly->setGeometry(QRect(610, 110, 111, 41));
        pushButton_CDT_Hole_Filling = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_CDT_Hole_Filling->setObjectName(QStringLiteral("pushButton_CDT_Hole_Filling"));
        pushButton_CDT_Hole_Filling->setGeometry(QRect(180, 430, 151, 41));
        label_5 = new QLabel(HoleFillingInteractiveOperationDialog);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(120, 442, 54, 20));
        pushButton_dist_V2V = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_dist_V2V->setObjectName(QStringLiteral("pushButton_dist_V2V"));
        pushButton_dist_V2V->setGeometry(QRect(620, 160, 101, 41));
        pushButton_test_ply_op = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_test_ply_op->setObjectName(QStringLiteral("pushButton_test_ply_op"));
        pushButton_test_ply_op->setGeometry(QRect(340, 430, 101, 41));
        line_9 = new QFrame(HoleFillingInteractiveOperationDialog);
        line_9->setObjectName(QStringLiteral("line_9"));
        line_9->setGeometry(QRect(0, 470, 731, 16));
        line_9->setFrameShape(QFrame::HLine);
        line_9->setFrameShadow(QFrame::Sunken);
        pushButton_insert_vertices = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_insert_vertices->setObjectName(QStringLiteral("pushButton_insert_vertices"));
        pushButton_insert_vertices->setGeometry(QRect(0, 490, 151, 41));
        line_10 = new QFrame(HoleFillingInteractiveOperationDialog);
        line_10->setObjectName(QStringLiteral("line_10"));
        line_10->setGeometry(QRect(380, 550, 16, 111));
        line_10->setFrameShape(QFrame::VLine);
        line_10->setFrameShadow(QFrame::Sunken);
        pushButton_gen_two_key_vertices = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_gen_two_key_vertices->setObjectName(QStringLiteral("pushButton_gen_two_key_vertices"));
        pushButton_gen_two_key_vertices->setGeometry(QRect(400, 550, 61, 101));
        pushButton_gen_bottom_plane = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_gen_bottom_plane->setObjectName(QStringLiteral("pushButton_gen_bottom_plane"));
        pushButton_gen_bottom_plane->setGeometry(QRect(470, 550, 75, 51));
        pushButton_4 = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_4->setObjectName(QStringLiteral("pushButton_4"));
        pushButton_4->setGeometry(QRect(550, 550, 75, 51));
        pushButton_5 = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_5->setObjectName(QStringLiteral("pushButton_5"));
        pushButton_5->setGeometry(QRect(470, 600, 75, 51));
        pushButton_6 = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_6->setObjectName(QStringLiteral("pushButton_6"));
        pushButton_6->setGeometry(QRect(550, 600, 75, 51));
        pushButton__insert_interior_vertices = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton__insert_interior_vertices->setObjectName(QStringLiteral("pushButton__insert_interior_vertices"));
        pushButton__insert_interior_vertices->setGeometry(QRect(180, 490, 171, 41));
        pushButton_display_ply_normal = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_display_ply_normal->setObjectName(QStringLiteral("pushButton_display_ply_normal"));
        pushButton_display_ply_normal->setGeometry(QRect(644, 430, 101, 41));
        pushButton_2 = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));
        pushButton_2->setGeometry(QRect(0, 600, 111, 41));
        pushButton_3 = new QPushButton(HoleFillingInteractiveOperationDialog);
        pushButton_3->setObjectName(QStringLiteral("pushButton_3"));
        pushButton_3->setGeometry(QRect(640, 550, 91, 101));

        retranslateUi(HoleFillingInteractiveOperationDialog);

        QMetaObject::connectSlotsByName(HoleFillingInteractiveOperationDialog);
    } // setupUi

    void retranslateUi(QDialog *HoleFillingInteractiveOperationDialog)
    {
        HoleFillingInteractiveOperationDialog->setWindowTitle(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\251\272\346\264\236\344\277\256\345\244\215\346\216\247\345\210\266\351\235\242\347\211\210", Q_NULLPTR));
        pushButton_simplify_line_segs->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\272\277\346\256\265\345\214\226\347\256\200", Q_NULLPTR));
        pushButton_select_poly->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\200\211\346\213\251\345\244\232\350\276\271\345\275\242", Q_NULLPTR));
        label_poly_id->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\244\232\350\276\271\345\275\242\347\274\226\345\217\267\357\274\232", Q_NULLPTR));
        pushButton_select_p1->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\200\211\346\213\251\351\241\266\347\202\2711", Q_NULLPTR));
        label_p1->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2711\357\274\232", Q_NULLPTR));
        pushButton_select_p2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\200\211\346\213\251\351\241\266\347\202\2712", Q_NULLPTR));
        label_p2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2712\357\274\232", Q_NULLPTR));
        pushButton_simplify_cancel_poly_selection->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\217\226\346\266\210\345\244\232\350\276\271\345\275\242 ", Q_NULLPTR));
        pushButton_gener_line_seg->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\347\272\277\346\256\265", Q_NULLPTR));
        pushButton_switch_ling_seg->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\210\207\346\215\242\347\272\277\346\256\265", Q_NULLPTR));
        pushButton_del_lineseg->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\217\226\346\266\210\347\272\277\346\256\265", Q_NULLPTR));
        pushButton_find_prev_point->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\211\215\344\270\200\344\270\252", Q_NULLPTR));
        pushButton_find_next_point->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\220\216\344\270\200\344\270\252", Q_NULLPTR));
        pushButton_search_angle2angle_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\220\234\347\264\242\350\247\222\345\257\271\350\247\222\351\241\266\347\202\271", Q_NULLPTR));
        pushButton_gen_points_set_for_tri->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\236\204\351\200\240\344\270\211\350\247\222\345\214\226\347\202\271\351\233\206", Q_NULLPTR));
        label->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\236\204\351\200\240\344\270\211\350\247\222\345\214\226\347\202\271\351\233\206", Q_NULLPTR));
        pushButton_get_camera_pos->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\350\216\267\345\217\226\345\275\223\345\211\215camera\344\275\215\347\275\256", Q_NULLPTR));
        pushButton_test_findNextPoint->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\265\213\350\257\225\344\270\213\344\270\200\344\270\252\351\241\266\347\202\271", Q_NULLPTR));
        pushButton_select_one_vertex->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\250\241\345\274\217\344\270\200", Q_NULLPTR));
        radioButton_select_one_vertex->setText(QString());
        pushButton_select_two_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\250\241\345\274\217\344\272\214", Q_NULLPTR));
        radioButton_select_two_vertices->setText(QString());
        pushButton_sel_v1_for_mode2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2711", Q_NULLPTR));
        label_2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "------->", Q_NULLPTR));
        pushButton_sel_v2_for_mode2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2712", Q_NULLPTR));
        pushButton_select_four_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\250\241\345\274\217\344\270\211", Q_NULLPTR));
        radioButton_select_four_vertices->setText(QString());
        pushButton_sel_v1_for_mode3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2711", Q_NULLPTR));
        label_3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "------->", Q_NULLPTR));
        pushButton_sel_v2_for_mode3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2712", Q_NULLPTR));
        pushButton_sel_v3_for_mode3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2713", Q_NULLPTR));
        label_4->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "------->", Q_NULLPTR));
        pushButton_sel_v4_for_mode3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\2714", Q_NULLPTR));
        pushButton_test->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\265\213\350\257\225\345\261\200\351\203\250\346\212\225\345\275\261", Q_NULLPTR));
        pushButton_set_camera_pos->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\350\256\276\347\275\256camera\344\275\215\347\275\256", Q_NULLPTR));
        pushButton_gen_proj_line_segs->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\346\212\225\345\275\261\347\272\277\346\256\265", Q_NULLPTR));
        pushButton_cancel_poly_wireframe_display->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\217\226\346\266\210\346\230\276\347\244\272\345\244\232\350\276\271\345\275\242\347\272\277\346\241\206", Q_NULLPTR));
        pushButton_display_poly_wireframe->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\230\276\347\244\272\345\244\232\350\276\271\345\275\242\347\272\277\346\241\206", Q_NULLPTR));
        pushButton_get_ver_num->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\271\346\225\260\351\207\217", Q_NULLPTR));
        pushButton_get_poly_num->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\244\232\350\276\271\345\275\242\346\225\260\351\207\217", Q_NULLPTR));
        pushButton_extend_to_bottom->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\220\221\344\270\213\345\273\266\344\274\270", Q_NULLPTR));
        pushButton_extend_to_backward->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\220\221\345\220\216\345\273\266\344\274\270", Q_NULLPTR));
        pushButton_extend_to_left->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\220\221\345\267\246\345\273\266\344\274\270", Q_NULLPTR));
        pushButton_extend_to_right->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\220\221\345\217\263\345\273\266\344\274\270", Q_NULLPTR));
        pushButton->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "CDT \344\270\211\350\247\222\345\214\226", Q_NULLPTR));
        pushButton_del_poly->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\210\240\351\231\244\345\244\232\350\276\271\345\275\242", Q_NULLPTR));
        pushButton_CDT_Hole_Filling->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\237\272\344\272\216CDT\347\232\204\344\270\211\350\247\222\345\214\226\345\255\224\346\264\236\344\277\256\345\244\215", Q_NULLPTR));
        label_5->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", ">>>>>>>>", Q_NULLPTR));
        pushButton_dist_V2V->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\351\241\266\347\202\271\350\267\235\347\246\273", Q_NULLPTR));
        pushButton_test_ply_op->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\265\213\350\257\225\345\244\232\350\276\271\345\275\242\346\223\215\344\275\234", Q_NULLPTR));
        pushButton_insert_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\234\250\345\244\232\350\276\271\345\275\242\350\276\271\344\270\212\346\267\273\345\212\240\351\241\266\347\202\271", Q_NULLPTR));
        pushButton_gen_two_key_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\344\270\244\344\270\252\n"
"\345\205\263\351\224\256\351\241\266\347\202\271\n"
"\346\205\216\351\207\215\347\202\271\345\207\273", Q_NULLPTR));
        pushButton_gen_bottom_plane->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\345\272\225\351\235\242", Q_NULLPTR));
        pushButton_4->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\345\267\246\344\276\247\n"
"\347\253\213\351\235\242", Q_NULLPTR));
        pushButton_5->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\345\217\263\344\276\247\n"
"\347\253\213\351\235\242", Q_NULLPTR));
        pushButton_6->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\347\224\237\346\210\220\345\220\216\344\276\247\n"
"\347\253\213\351\235\242", Q_NULLPTR));
        pushButton__insert_interior_vertices->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\345\234\250\345\244\232\350\276\271\345\275\242\345\206\205\351\203\250\346\267\273\345\212\240\351\241\266\347\202\271", Q_NULLPTR));
        pushButton_display_ply_normal->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\230\276\347\244\272\346\263\225\345\220\221\351\207\217", Q_NULLPTR));
        pushButton_2->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\250\241\345\236\213\347\251\272\351\227\264\350\214\203\345\233\264", Q_NULLPTR));
        pushButton_3->setText(QApplication::translate("HoleFillingInteractiveOperationDialog", "\346\250\241\345\236\213\345\260\201\351\227\255\346\200\247\n"
"\346\243\200\346\265\213", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class HoleFillingInteractiveOperationDialog: public Ui_HoleFillingInteractiveOperationDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HOLEFILLINGINTERACTIVEOPERATIONDIALOG_H
