/********************************************************************************
** Form generated from reading UI file 'ControlOptimizationSnappingDialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CONTROLOPTIMIZATIONSNAPPINGDIALOG_H
#define UI_CONTROLOPTIMIZATIONSNAPPINGDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ControlOptimizationSnappingDialog
{
public:
    QWidget *widget;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton;
    QLabel *label;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *ControlOptimizationSnappingDialog)
    {
        if (ControlOptimizationSnappingDialog->objectName().isEmpty())
            ControlOptimizationSnappingDialog->setObjectName(QStringLiteral("ControlOptimizationSnappingDialog"));
        ControlOptimizationSnappingDialog->resize(392, 134);
        widget = new QWidget(ControlOptimizationSnappingDialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(0, 0, 391, 131));
        verticalLayout = new QVBoxLayout(widget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        pushButton = new QPushButton(widget);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        horizontalLayout->addWidget(pushButton);

        label = new QLabel(widget);
        label->setObjectName(QStringLiteral("label"));
        QFont font;
        font.setPointSize(18);
        label->setFont(font);

        horizontalLayout->addWidget(label);


        verticalLayout->addLayout(horizontalLayout);

        buttonBox = new QDialogButtonBox(widget);
        buttonBox->setObjectName(QStringLiteral("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(ControlOptimizationSnappingDialog);
        QObject::connect(buttonBox, SIGNAL(accepted()), ControlOptimizationSnappingDialog, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), ControlOptimizationSnappingDialog, SLOT(reject()));

        QMetaObject::connectSlotsByName(ControlOptimizationSnappingDialog);
    } // setupUi

    void retranslateUi(QDialog *ControlOptimizationSnappingDialog)
    {
        ControlOptimizationSnappingDialog->setWindowTitle(QApplication::translate("ControlOptimizationSnappingDialog", "Dialog", Q_NULLPTR));
        pushButton->setText(QApplication::translate("ControlOptimizationSnappingDialog", "\346\211\247\350\241\214\344\274\230\345\214\226", Q_NULLPTR));
        label->setText(QApplication::translate("ControlOptimizationSnappingDialog", "\350\277\255\344\273\243\347\254\2540\346\254\241", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ControlOptimizationSnappingDialog: public Ui_ControlOptimizationSnappingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CONTROLOPTIMIZATIONSNAPPINGDIALOG_H
