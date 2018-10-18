#ifndef CONTROLOPTIMIZATIONSNAPPINGDIALOG_H
#define CONTROLOPTIMIZATIONSNAPPINGDIALOG_H

#include <QDialog>

namespace Ui {
class ControlOptimizationSnappingDialog;
}

class ControlOptimizationSnappingDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit ControlOptimizationSnappingDialog(QWidget *parent = 0);
    ~ControlOptimizationSnappingDialog();
    void resetCount();
    void setPolyID(int id);
    void *pclviewer;
    
private slots:
    void on_pushButton_clicked();

    void on_buttonBox_accepted();

private:
    Ui::ControlOptimizationSnappingDialog *ui;
    int m_count;
    int poly_id;
};

#endif // CONTROLOPTIMIZATIONSNAPPINGDIALOG_H
