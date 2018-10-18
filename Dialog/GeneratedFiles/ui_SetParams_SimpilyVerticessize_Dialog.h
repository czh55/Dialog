/********************************************************************************
** Form generated from reading UI file 'SetParams_SimpilyVerticessize_Dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SETPARAMS_SIMPILYVERTICESSIZE_DIALOG_H
#define UI_SETPARAMS_SIMPILYVERTICESSIZE_DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QTableView>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SetParams_SimpilyVerticesSize_Dialog
{
public:
    QWidget *widget;
    QVBoxLayout *verticalLayout;
    QTableView *tableView;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *SetParams_SimpilyVerticesSize_Dialog)
    {
        if (SetParams_SimpilyVerticesSize_Dialog->objectName().isEmpty())
            SetParams_SimpilyVerticesSize_Dialog->setObjectName(QStringLiteral("SetParams_SimpilyVerticesSize_Dialog"));
        SetParams_SimpilyVerticesSize_Dialog->resize(718, 428);
        widget = new QWidget(SetParams_SimpilyVerticesSize_Dialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(10, 10, 701, 411));
        verticalLayout = new QVBoxLayout(widget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        tableView = new QTableView(widget);
        tableView->setObjectName(QStringLiteral("tableView"));

        verticalLayout->addWidget(tableView);

        buttonBox = new QDialogButtonBox(widget);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(SetParams_SimpilyVerticesSize_Dialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), SetParams_SimpilyVerticesSize_Dialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), SetParams_SimpilyVerticesSize_Dialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(SetParams_SimpilyVerticesSize_Dialog);
    } // setupUi

    void retranslateUi(QDialog *SetParams_SimpilyVerticesSize_Dialog)
    {
        SetParams_SimpilyVerticesSize_Dialog->setWindowTitle(QApplication::translate("SetParams_SimpilyVerticesSize_Dialog", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class SetParams_SimpilyVerticesSize_Dialog: public Ui_SetParams_SimpilyVerticesSize_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SETPARAMS_SIMPILYVERTICESSIZE_DIALOG_H
