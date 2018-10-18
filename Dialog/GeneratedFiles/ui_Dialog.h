/********************************************************************************
** Form generated from reading UI file 'Dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_H
#define UI_DIALOG_H

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

class Ui_DialogClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *DialogClass)
    {
        if (DialogClass->objectName().isEmpty())
            DialogClass->setObjectName(QStringLiteral("DialogClass"));
        DialogClass->resize(600, 400);
        menuBar = new QMenuBar(DialogClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        DialogClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(DialogClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        DialogClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(DialogClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        DialogClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(DialogClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        DialogClass->setStatusBar(statusBar);

        retranslateUi(DialogClass);

        QMetaObject::connectSlotsByName(DialogClass);
    } // setupUi

    void retranslateUi(QMainWindow *DialogClass)
    {
        DialogClass->setWindowTitle(QApplication::translate("DialogClass", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class DialogClass: public Ui_DialogClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_H
