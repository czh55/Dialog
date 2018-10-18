//以下这句话解决这个错误：https://stackoverflow.com/questions/18642155/no-override-found-for-vtkpolydatamapper
#define vtkRenderingCore_AUTOINIT 2(vtkRenderingOpenGL2, vtkInteractionStyle)

#include <QApplication>
#include "pclviewer.h"
//#include "MainWindow.h"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	/*MainWindow w;
	w.show();*/

	PCLViewer v;
	v.show();
	v.init();
	return a.exec();
}
