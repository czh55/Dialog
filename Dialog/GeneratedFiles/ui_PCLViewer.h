/********************************************************************************
** Form generated from reading UI file 'PCLViewer.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PCLVIEWER_H
#define UI_PCLVIEWER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_PCLViewer
{
public:
    QAction *OpenAction;
    QAction *setBgColorAction;
    QAction *setCoordinateAxesAction;
    QAction *removePointCloudAction;
    QAction *verticesNormalEstimationAction;
    QAction *clearScreenAction;
    QAction *ObservePolyVerNormalAction;
    QAction *simplifyVerticesSizeAction;
    QAction *RGtoDetectLineSegsAction;
    QAction *RGlineSegDisplayAction;
    QAction *doSimplifyVerticesSizeAction;
    QAction *removeAllPointCloudsAction;
    QAction *removeAllShapesAction;
    QAction *removeAllInitialVerticesAction;
    QAction *showNeighborRegionAction;
    QAction *mergeVerticesParamsAction;
    QAction *showMatchingAction;
    QAction *mergeMatchingRelationsAction;
    QAction *showMatchingClustersAction;
    QAction *mergeVertexEdgeMatchingAction;
    QAction *viewSpecificPolyAction;
    QAction *action_4;
    QAction *dealWithPolyDegenerationAction;
    QAction *performSnappingVerticesAction;
    QAction *goOnPerformSnappingVerticesAction;
    QAction *observePolyIDAction;
    QAction *performOptimisePositionAction;
    QAction *action;
    QAction *action_2;
    QAction *action_5;
    QAction *action_6;
    QAction *action_3;
    QAction *action_7;
    QAction *action_8;
    QAction *action_9;
    QAction *action_dda;
    QAction *action_11;
    QAction *action_10;
    QAction *action_12;
    QAction *action_13;
    QAction *action_14;
    QAction *action_16;
    QAction *action_1;
    QAction *action_17;
    QAction *action_18;
    QAction *action_19;
    QAction *action_DDA;
    QAction *action_20;
    QAction *action_21;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    QVTKWidget *qvtkWidget;
    QMenuBar *menuBar;
    QMenu *menu;
    QMenu *menu_2;
    QMenu *menu_3;
    QMenu *menu_5;
    QMenu *menu_6;
    QMenu *menu_4;
    QMenu *menu_7;
    QMenu *menu_8;
    QMenu *menu_9;
    QMenu *menu_10;
    QMenu *menu_11;
    QMenu *menu_12;
    QMenu *menu_13;

    void setupUi(QMainWindow *PCLViewer)
    {
        if (PCLViewer->objectName().isEmpty())
            PCLViewer->setObjectName(QStringLiteral("PCLViewer"));
        PCLViewer->resize(989, 498);
        PCLViewer->setMinimumSize(QSize(0, 0));
        PCLViewer->setMaximumSize(QSize(5000, 5000));
        OpenAction = new QAction(PCLViewer);
        OpenAction->setObjectName(QStringLiteral("OpenAction"));
        setBgColorAction = new QAction(PCLViewer);
        setBgColorAction->setObjectName(QStringLiteral("setBgColorAction"));
        setCoordinateAxesAction = new QAction(PCLViewer);
        setCoordinateAxesAction->setObjectName(QStringLiteral("setCoordinateAxesAction"));
        removePointCloudAction = new QAction(PCLViewer);
        removePointCloudAction->setObjectName(QStringLiteral("removePointCloudAction"));
        verticesNormalEstimationAction = new QAction(PCLViewer);
        verticesNormalEstimationAction->setObjectName(QStringLiteral("verticesNormalEstimationAction"));
        clearScreenAction = new QAction(PCLViewer);
        clearScreenAction->setObjectName(QStringLiteral("clearScreenAction"));
        ObservePolyVerNormalAction = new QAction(PCLViewer);
        ObservePolyVerNormalAction->setObjectName(QStringLiteral("ObservePolyVerNormalAction"));
        simplifyVerticesSizeAction = new QAction(PCLViewer);
        simplifyVerticesSizeAction->setObjectName(QStringLiteral("simplifyVerticesSizeAction"));
        RGtoDetectLineSegsAction = new QAction(PCLViewer);
        RGtoDetectLineSegsAction->setObjectName(QStringLiteral("RGtoDetectLineSegsAction"));
        RGlineSegDisplayAction = new QAction(PCLViewer);
        RGlineSegDisplayAction->setObjectName(QStringLiteral("RGlineSegDisplayAction"));
        doSimplifyVerticesSizeAction = new QAction(PCLViewer);
        doSimplifyVerticesSizeAction->setObjectName(QStringLiteral("doSimplifyVerticesSizeAction"));
        removeAllPointCloudsAction = new QAction(PCLViewer);
        removeAllPointCloudsAction->setObjectName(QStringLiteral("removeAllPointCloudsAction"));
        removeAllShapesAction = new QAction(PCLViewer);
        removeAllShapesAction->setObjectName(QStringLiteral("removeAllShapesAction"));
        removeAllInitialVerticesAction = new QAction(PCLViewer);
        removeAllInitialVerticesAction->setObjectName(QStringLiteral("removeAllInitialVerticesAction"));
        showNeighborRegionAction = new QAction(PCLViewer);
        showNeighborRegionAction->setObjectName(QStringLiteral("showNeighborRegionAction"));
        mergeVerticesParamsAction = new QAction(PCLViewer);
        mergeVerticesParamsAction->setObjectName(QStringLiteral("mergeVerticesParamsAction"));
        showMatchingAction = new QAction(PCLViewer);
        showMatchingAction->setObjectName(QStringLiteral("showMatchingAction"));
        mergeMatchingRelationsAction = new QAction(PCLViewer);
        mergeMatchingRelationsAction->setObjectName(QStringLiteral("mergeMatchingRelationsAction"));
        showMatchingClustersAction = new QAction(PCLViewer);
        showMatchingClustersAction->setObjectName(QStringLiteral("showMatchingClustersAction"));
        mergeVertexEdgeMatchingAction = new QAction(PCLViewer);
        mergeVertexEdgeMatchingAction->setObjectName(QStringLiteral("mergeVertexEdgeMatchingAction"));
        viewSpecificPolyAction = new QAction(PCLViewer);
        viewSpecificPolyAction->setObjectName(QStringLiteral("viewSpecificPolyAction"));
        action_4 = new QAction(PCLViewer);
        action_4->setObjectName(QStringLiteral("action_4"));
        dealWithPolyDegenerationAction = new QAction(PCLViewer);
        dealWithPolyDegenerationAction->setObjectName(QStringLiteral("dealWithPolyDegenerationAction"));
        performSnappingVerticesAction = new QAction(PCLViewer);
        performSnappingVerticesAction->setObjectName(QStringLiteral("performSnappingVerticesAction"));
        goOnPerformSnappingVerticesAction = new QAction(PCLViewer);
        goOnPerformSnappingVerticesAction->setObjectName(QStringLiteral("goOnPerformSnappingVerticesAction"));
        observePolyIDAction = new QAction(PCLViewer);
        observePolyIDAction->setObjectName(QStringLiteral("observePolyIDAction"));
        performOptimisePositionAction = new QAction(PCLViewer);
        performOptimisePositionAction->setObjectName(QStringLiteral("performOptimisePositionAction"));
        action = new QAction(PCLViewer);
        action->setObjectName(QStringLiteral("action"));
        action_2 = new QAction(PCLViewer);
        action_2->setObjectName(QStringLiteral("action_2"));
        action_5 = new QAction(PCLViewer);
        action_5->setObjectName(QStringLiteral("action_5"));
        action_6 = new QAction(PCLViewer);
        action_6->setObjectName(QStringLiteral("action_6"));
        action_3 = new QAction(PCLViewer);
        action_3->setObjectName(QStringLiteral("action_3"));
        action_7 = new QAction(PCLViewer);
        action_7->setObjectName(QStringLiteral("action_7"));
        action_8 = new QAction(PCLViewer);
        action_8->setObjectName(QStringLiteral("action_8"));
        action_9 = new QAction(PCLViewer);
        action_9->setObjectName(QStringLiteral("action_9"));
        action_dda = new QAction(PCLViewer);
        action_dda->setObjectName(QStringLiteral("action_dda"));
        action_11 = new QAction(PCLViewer);
        action_11->setObjectName(QStringLiteral("action_11"));
        action_10 = new QAction(PCLViewer);
        action_10->setObjectName(QStringLiteral("action_10"));
        action_12 = new QAction(PCLViewer);
        action_12->setObjectName(QStringLiteral("action_12"));
        action_13 = new QAction(PCLViewer);
        action_13->setObjectName(QStringLiteral("action_13"));
        action_14 = new QAction(PCLViewer);
        action_14->setObjectName(QStringLiteral("action_14"));
        action_16 = new QAction(PCLViewer);
        action_16->setObjectName(QStringLiteral("action_16"));
        action_1 = new QAction(PCLViewer);
        action_1->setObjectName(QStringLiteral("action_1"));
        action_17 = new QAction(PCLViewer);
        action_17->setObjectName(QStringLiteral("action_17"));
        action_18 = new QAction(PCLViewer);
        action_18->setObjectName(QStringLiteral("action_18"));
        action_19 = new QAction(PCLViewer);
        action_19->setObjectName(QStringLiteral("action_19"));
        action_DDA = new QAction(PCLViewer);
        action_DDA->setObjectName(QStringLiteral("action_DDA"));
        action_20 = new QAction(PCLViewer);
        action_20->setObjectName(QStringLiteral("action_20"));
        action_21 = new QAction(PCLViewer);
        action_21->setObjectName(QStringLiteral("action_21"));
        centralwidget = new QWidget(PCLViewer);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        horizontalLayout = new QHBoxLayout(centralwidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        qvtkWidget = new QVTKWidget(centralwidget);
        qvtkWidget->setObjectName(QStringLiteral("qvtkWidget"));

        horizontalLayout->addWidget(qvtkWidget);

        PCLViewer->setCentralWidget(centralwidget);
        menuBar = new QMenuBar(PCLViewer);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 989, 23));
        menu = new QMenu(menuBar);
        menu->setObjectName(QStringLiteral("menu"));
        menu_2 = new QMenu(menuBar);
        menu_2->setObjectName(QStringLiteral("menu_2"));
        menu_3 = new QMenu(menuBar);
        menu_3->setObjectName(QStringLiteral("menu_3"));
        menu_5 = new QMenu(menu_3);
        menu_5->setObjectName(QStringLiteral("menu_5"));
        menu_6 = new QMenu(menu_3);
        menu_6->setObjectName(QStringLiteral("menu_6"));
        menu_4 = new QMenu(menuBar);
        menu_4->setObjectName(QStringLiteral("menu_4"));
        menu_7 = new QMenu(menuBar);
        menu_7->setObjectName(QStringLiteral("menu_7"));
        menu_8 = new QMenu(menu_7);
        menu_8->setObjectName(QStringLiteral("menu_8"));
        menu_9 = new QMenu(menu_7);
        menu_9->setObjectName(QStringLiteral("menu_9"));
        menu_10 = new QMenu(menuBar);
        menu_10->setObjectName(QStringLiteral("menu_10"));
        menu_11 = new QMenu(menuBar);
        menu_11->setObjectName(QStringLiteral("menu_11"));
        menu_12 = new QMenu(menuBar);
        menu_12->setObjectName(QStringLiteral("menu_12"));
        menu_13 = new QMenu(menu_12);
        menu_13->setObjectName(QStringLiteral("menu_13"));
        PCLViewer->setMenuBar(menuBar);

        menuBar->addAction(menu->menuAction());
        menuBar->addAction(menu_2->menuAction());
        menuBar->addAction(menu_4->menuAction());
        menuBar->addAction(menu_3->menuAction());
        menuBar->addAction(menu_7->menuAction());
        menuBar->addAction(menu_10->menuAction());
        menuBar->addAction(menu_11->menuAction());
        menuBar->addAction(menu_12->menuAction());
        menu->addAction(OpenAction);
        menu->addAction(action_9);
        menu_2->addAction(setBgColorAction);
        menu_2->addAction(setCoordinateAxesAction);
        menu_2->addAction(removePointCloudAction);
        menu_2->addAction(clearScreenAction);
        menu_2->addAction(removeAllPointCloudsAction);
        menu_2->addAction(removeAllShapesAction);
        menu_3->addAction(menu_5->menuAction());
        menu_3->addAction(menu_6->menuAction());
        menu_3->addSeparator();
        menu_3->addAction(doSimplifyVerticesSizeAction);
        menu_3->addAction(removeAllInitialVerticesAction);
        menu_5->addAction(verticesNormalEstimationAction);
        menu_5->addAction(ObservePolyVerNormalAction);
        menu_6->addAction(RGtoDetectLineSegsAction);
        menu_6->addAction(RGlineSegDisplayAction);
        menu_4->addAction(simplifyVerticesSizeAction);
        menu_4->addAction(mergeVerticesParamsAction);
        menu_7->addAction(showNeighborRegionAction);
        menu_7->addAction(showMatchingAction);
        menu_7->addAction(showMatchingClustersAction);
        menu_7->addAction(mergeMatchingRelationsAction);
        menu_7->addAction(mergeVertexEdgeMatchingAction);
        menu_7->addAction(menu_8->menuAction());
        menu_7->addSeparator();
        menu_7->addAction(performSnappingVerticesAction);
        menu_7->addAction(goOnPerformSnappingVerticesAction);
        menu_7->addSeparator();
        menu_7->addAction(menu_9->menuAction());
        menu_7->addSeparator();
        menu_7->addAction(action_5);
        menu_7->addSeparator();
        menu_7->addAction(action_11);
        menu_7->addSeparator();
        menu_7->addAction(action_6);
        menu_7->addSeparator();
        menu_7->addAction(action_20);
        menu_8->addAction(viewSpecificPolyAction);
        menu_8->addAction(action_4);
        menu_8->addAction(dealWithPolyDegenerationAction);
        menu_9->addAction(observePolyIDAction);
        menu_9->addAction(performOptimisePositionAction);
        menu_9->addAction(action_10);
        menu_10->addAction(action);
        menu_10->addAction(action_2);
        menu_11->addAction(action_3);
        menu_11->addAction(action_7);
        menu_11->addAction(action_8);
        menu_11->addAction(action_dda);
        menu_12->addAction(action_12);
        menu_12->addAction(action_14);
        menu_12->addAction(action_13);
        menu_12->addAction(menu_13->menuAction());
        menu_12->addAction(action_19);
        menu_12->addAction(action_DDA);
        menu_12->addSeparator();
        menu_12->addAction(action_21);
        menu_13->addAction(action_16);
        menu_13->addAction(action_1);
        menu_13->addAction(action_17);
        menu_13->addAction(action_18);

        retranslateUi(PCLViewer);

        QMetaObject::connectSlotsByName(PCLViewer);
    } // setupUi

    void retranslateUi(QMainWindow *PCLViewer)
    {
        PCLViewer->setWindowTitle(QApplication::translate("PCLViewer", "Rockmass3DReconstruction", Q_NULLPTR));
        OpenAction->setText(QApplication::translate("PCLViewer", "\346\211\223\345\274\200\345\244\232\350\276\271\345\275\242\346\225\260\346\215\256", Q_NULLPTR));
#ifndef QT_NO_TOOLTIP
        OpenAction->setToolTip(QApplication::translate("PCLViewer", "\346\211\223\345\274\200\345\244\232\350\276\271\345\275\242\346\226\207\344\273\266", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        setBgColorAction->setText(QApplication::translate("PCLViewer", "\350\203\214\346\231\257\350\211\262", Q_NULLPTR));
        setCoordinateAxesAction->setText(QApplication::translate("PCLViewer", "\345\235\220\346\240\207\350\275\264", Q_NULLPTR));
        removePointCloudAction->setText(QApplication::translate("PCLViewer", "\347\247\273\351\231\244\346\214\207\345\256\232\347\202\271\344\272\221\346\210\226\345\275\242\347\212\266", Q_NULLPTR));
        verticesNormalEstimationAction->setText(QApplication::translate("PCLViewer", "\346\263\225\345\220\221\351\207\217\344\274\260\350\256\241", Q_NULLPTR));
        clearScreenAction->setText(QApplication::translate("PCLViewer", "\346\270\205\345\261\217", Q_NULLPTR));
        ObservePolyVerNormalAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\347\211\271\345\256\232\345\244\232\350\276\271\345\275\242\347\232\204\351\241\266\347\202\271\346\263\225\345\220\221\351\207\217", Q_NULLPTR));
        simplifyVerticesSizeAction->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\214\226\347\256\200", Q_NULLPTR));
        RGtoDetectLineSegsAction->setText(QApplication::translate("PCLViewer", "\345\214\272\345\237\237\347\224\237\351\225\277", Q_NULLPTR));
        RGlineSegDisplayAction->setText(QApplication::translate("PCLViewer", "\346\230\276\347\244\272", Q_NULLPTR));
        doSimplifyVerticesSizeAction->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\214\226\347\256\200(\344\270\200\346\255\245\346\211\247\350\241\214)", Q_NULLPTR));
        removeAllPointCloudsAction->setText(QApplication::translate("PCLViewer", "\347\247\273\351\231\244\346\211\200\346\234\211\347\202\271\344\272\221", Q_NULLPTR));
        removeAllShapesAction->setText(QApplication::translate("PCLViewer", "\347\247\273\351\231\244\346\211\200\346\234\211\345\275\242\347\212\266", Q_NULLPTR));
        removeAllInitialVerticesAction->setText(QApplication::translate("PCLViewer", "\347\247\273\351\231\244\345\210\235\345\247\213\351\241\266\347\202\271\351\233\206\345\220\210", Q_NULLPTR));
        showNeighborRegionAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\351\241\266\347\202\271\347\232\204\351\202\273\345\237\237\345\215\212\345\276\204", Q_NULLPTR));
        mergeVerticesParamsAction->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\220\210\345\271\266", Q_NULLPTR));
        showMatchingAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\345\214\271\351\205\215\345\205\263\347\263\273(\347\202\271-\347\202\271\357\274\214\347\202\271-\350\276\271)", Q_NULLPTR));
        mergeMatchingRelationsAction->setText(QApplication::translate("PCLViewer", "\345\220\210\345\271\266\345\214\271\351\205\215\345\205\263\347\263\273(\347\202\271-\347\202\271)", Q_NULLPTR));
        showMatchingClustersAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\345\214\271\351\205\215\347\232\204\351\241\266\347\202\271\347\260\207", Q_NULLPTR));
        mergeVertexEdgeMatchingAction->setText(QApplication::translate("PCLViewer", "\345\220\210\345\271\266\345\214\271\351\205\215\345\205\263\347\263\273(\347\202\271-\350\276\271\357\274\211", Q_NULLPTR));
        viewSpecificPolyAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\347\211\271\345\256\232\345\244\232\350\276\271\345\275\242", Q_NULLPTR));
        action_4->setText(QApplication::translate("PCLViewer", "\350\275\256\347\225\252\346\237\245\347\234\213\347\211\271\345\256\232\345\244\232\350\276\271\345\275\242", Q_NULLPTR));
        dealWithPolyDegenerationAction->setText(QApplication::translate("PCLViewer", "\345\244\204\347\220\206\345\244\232\350\276\271\345\275\242\351\200\200\345\214\226\346\203\205\345\206\265", Q_NULLPTR));
        performSnappingVerticesAction->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\220\210\345\271\266(\344\270\200\346\255\245\346\211\247\350\241\214)", Q_NULLPTR));
        goOnPerformSnappingVerticesAction->setText(QApplication::translate("PCLViewer", "\347\273\247\347\273\255\346\211\247\350\241\214", Q_NULLPTR));
        observePolyIDAction->setText(QApplication::translate("PCLViewer", "\346\237\245\347\234\213\345\244\232\350\276\271\345\275\242id", Q_NULLPTR));
        performOptimisePositionAction->setText(QApplication::translate("PCLViewer", "\346\211\247\350\241\214\344\275\215\347\275\256\344\274\230\345\214\226", Q_NULLPTR));
        action->setText(QApplication::translate("PCLViewer", "\347\224\237\346\210\220\346\265\213\350\257\225\345\244\232\350\276\271\345\275\242", Q_NULLPTR));
        action_2->setText(QApplication::translate("PCLViewer", "\346\211\247\350\241\214\344\274\230\345\214\226", Q_NULLPTR));
        action_5->setText(QApplication::translate("PCLViewer", "\345\244\232\350\276\271\345\275\242\344\270\211\350\247\222\345\210\206\345\211\262(\346\230\276\347\244\272\345\271\263\351\235\242\346\225\210\346\236\234)", Q_NULLPTR));
        action_6->setText(QApplication::translate("PCLViewer", "\345\234\250\346\234\211\350\247\222\345\257\271\350\276\271\345\205\263\347\263\273\347\232\204\347\272\277\346\256\265\344\270\212\345\242\236\345\212\240\351\241\266\347\202\271\357\274\210\345\255\224\346\264\236\344\277\256\345\244\215\344\271\213\345\211\215\345\277\205\347\202\271\357\274\211", Q_NULLPTR));
        action_3->setText(QApplication::translate("PCLViewer", "\346\236\204\351\200\240\345\205\250\345\261\200\347\202\271\351\233\206", Q_NULLPTR));
        action_7->setText(QApplication::translate("PCLViewer", "\350\277\233\345\205\245\344\272\244\344\272\222\346\250\241\345\274\217", Q_NULLPTR));
        action_8->setText(QApplication::translate("PCLViewer", "\344\277\235\345\255\230\346\250\241\345\236\213\346\225\260\346\215\256", Q_NULLPTR));
        action_9->setText(QApplication::translate("PCLViewer", "\345\262\251\344\275\223\346\250\241\345\236\213\346\225\260\346\215\256", Q_NULLPTR));
        action_dda->setText(QApplication::translate("PCLViewer", "\344\277\235\345\255\230\346\210\220dda\346\225\260\346\215\256", Q_NULLPTR));
        action_11->setText(QApplication::translate("PCLViewer", "\346\243\200\346\237\245\346\250\241\345\236\213\345\220\210\346\263\225\346\200\247", Q_NULLPTR));
        action_10->setText(QApplication::translate("PCLViewer", "\346\211\247\350\241\214\345\244\232\350\276\271\345\275\242\345\210\240\351\231\244", Q_NULLPTR));
        action_12->setText(QApplication::translate("PCLViewer", "\345\212\240\350\275\275\344\274\230\345\212\277\347\273\223\346\236\204\351\235\242", Q_NULLPTR));
        action_13->setText(QApplication::translate("PCLViewer", "\350\257\206\345\210\253\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
        action_14->setText(QApplication::translate("PCLViewer", "\345\212\240\350\275\275\345\271\263\351\235\242\347\202\271\344\272\221", Q_NULLPTR));
        action_16->setText(QApplication::translate("PCLViewer", "\345\207\206\345\244\207", Q_NULLPTR));
        action_1->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\2711", Q_NULLPTR));
        action_17->setText(QApplication::translate("PCLViewer", "\351\241\266\347\202\2712", Q_NULLPTR));
        action_18->setText(QApplication::translate("PCLViewer", "\347\224\237\346\210\220\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
        action_19->setText(QApplication::translate("PCLViewer", "\346\230\276\347\244\272\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
        action_DDA->setText(QApplication::translate("PCLViewer", "\344\277\235\345\255\230\345\262\251\344\275\223\346\250\241\345\236\213DDA\346\225\260\346\215\256", Q_NULLPTR));
        action_20->setText(QApplication::translate("PCLViewer", "\345\261\225\347\244\272\347\262\276\345\215\216\345\214\271\351\205\215\345\205\263\347\263\273", Q_NULLPTR));
        action_21->setText(QApplication::translate("PCLViewer", "\346\211\213\345\212\250\346\214\207\345\256\232\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
        menu->setTitle(QApplication::translate("PCLViewer", "\346\226\207\344\273\266", Q_NULLPTR));
        menu_2->setTitle(QApplication::translate("PCLViewer", "\346\230\276\347\244\272", Q_NULLPTR));
        menu_3->setTitle(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\214\226\347\256\200", Q_NULLPTR));
        menu_5->setTitle(QApplication::translate("PCLViewer", "\346\263\225\345\220\221\351\207\217\344\274\260\350\256\241", Q_NULLPTR));
        menu_6->setTitle(QApplication::translate("PCLViewer", "\345\214\272\345\237\237\347\224\237\351\225\277\346\217\220\345\217\226\350\277\236\347\273\255\347\272\277\346\256\265", Q_NULLPTR));
        menu_4->setTitle(QApplication::translate("PCLViewer", "\345\217\202\346\225\260\350\256\276\347\275\256", Q_NULLPTR));
        menu_7->setTitle(QApplication::translate("PCLViewer", "\351\241\266\347\202\271\345\220\210\345\271\266", Q_NULLPTR));
        menu_8->setTitle(QApplication::translate("PCLViewer", "\345\244\204\347\220\206\345\244\232\350\276\271\345\275\242\351\200\200\345\214\226", Q_NULLPTR));
        menu_9->setTitle(QApplication::translate("PCLViewer", "\344\274\230\345\214\226\347\211\271\345\256\232\345\244\232\350\276\271\345\275\242\347\232\204\344\275\215\347\275\256", Q_NULLPTR));
        menu_10->setTitle(QApplication::translate("PCLViewer", "\346\265\213\350\257\225\344\274\230\345\214\226", Q_NULLPTR));
        menu_11->setTitle(QApplication::translate("PCLViewer", "\345\255\224\346\264\236\344\277\256\345\244\215", Q_NULLPTR));
        menu_12->setTitle(QApplication::translate("PCLViewer", "\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
        menu_13->setTitle(QApplication::translate("PCLViewer", "\344\272\244\344\272\222\345\274\217\346\217\222\345\205\245\344\270\215\350\277\236\347\273\255\351\235\242", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PCLViewer: public Ui_PCLViewer {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PCLVIEWER_H
