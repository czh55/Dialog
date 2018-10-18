#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_HiDialog.h"

class HiDialog : public QMainWindow
{
	Q_OBJECT

public:
	HiDialog(QWidget *parent = Q_NULLPTR);

private:
	Ui::HiDialogClass ui;
};
