/********************************************************************************
** Form generated from reading UI file 'HiDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_HIDIALOG_H
#define UI_HIDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_HiDialogClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *HiDialogClass)
    {
        if (HiDialogClass->objectName().isEmpty())
            HiDialogClass->setObjectName(QStringLiteral("HiDialogClass"));
        HiDialogClass->resize(600, 400);
        menuBar = new QMenuBar(HiDialogClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        HiDialogClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(HiDialogClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        HiDialogClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(HiDialogClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        HiDialogClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(HiDialogClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        HiDialogClass->setStatusBar(statusBar);

        retranslateUi(HiDialogClass);

        QMetaObject::connectSlotsByName(HiDialogClass);
    } // setupUi

    void retranslateUi(QMainWindow *HiDialogClass)
    {
        HiDialogClass->setWindowTitle(QApplication::translate("HiDialogClass", "HiDialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class HiDialogClass: public Ui_HiDialogClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_HIDIALOG_H
