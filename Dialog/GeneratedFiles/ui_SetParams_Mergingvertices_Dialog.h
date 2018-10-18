/********************************************************************************
** Form generated from reading UI file 'SetParams_Mergingvertices_Dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SETPARAMS_MERGINGVERTICES_DIALOG_H
#define UI_SETPARAMS_MERGINGVERTICES_DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QTableView>

QT_BEGIN_NAMESPACE

class Ui_SetParams_MergingVertices_Dialog
{
public:
    QTableView *tableView;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *SetParams_MergingVertices_Dialog)
    {
        if (SetParams_MergingVertices_Dialog->objectName().isEmpty())
            SetParams_MergingVertices_Dialog->setObjectName(QStringLiteral("SetParams_MergingVertices_Dialog"));
        SetParams_MergingVertices_Dialog->resize(667, 425);
        tableView = new QTableView(SetParams_MergingVertices_Dialog);
        tableView->setObjectName(QStringLiteral("tableView"));
        tableView->setGeometry(QRect(1, 1, 661, 381));
        buttonBox = new QDialogButtonBox(SetParams_MergingVertices_Dialog);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setGeometry(QRect(500, 390, 156, 23));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        retranslateUi(SetParams_MergingVertices_Dialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), SetParams_MergingVertices_Dialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), SetParams_MergingVertices_Dialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(SetParams_MergingVertices_Dialog);
    } // setupUi

    void retranslateUi(QDialog *SetParams_MergingVertices_Dialog)
    {
        SetParams_MergingVertices_Dialog->setWindowTitle(QApplication::translate("SetParams_MergingVertices_Dialog", "Dialog", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class SetParams_MergingVertices_Dialog: public Ui_SetParams_MergingVertices_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SETPARAMS_MERGINGVERTICES_DIALOG_H
