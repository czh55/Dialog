#ifndef SETPARAMS_MERGINGVERTICES_DIALOG_H
#define SETPARAMS_MERGINGVERTICES_DIALOG_H

#include <QDialog>
#include <QStandardItemModel>
#include <QString>
#include <QSettings>

namespace Ui {
class SetParams_MergingVertices_Dialog;
}

class SetParams_MergingVertices_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit SetParams_MergingVertices_Dialog(QWidget *parent = 0);
    ~SetParams_MergingVertices_Dialog();

private slots:
    void on_buttonBox_accepted();

private:
    Ui::SetParams_MergingVertices_Dialog *ui;
    QStandardItemModel *model;
    QSettings *configIniRead;

    void initTable();
    void setTableItem(int row, QString str1, QString str2);
};

#endif // SETPARAMS_MERGINGVERTICES_DIALOG_H
