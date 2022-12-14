#ifndef MDIWINDOW_H
#define MDIWINDOW_H

#include <Standard_WarningsDisable.hxx>
#include <QMainWindow>
#include <Standard_WarningsRestore.hxx>

#include "CommonSample.h"

class DocumentCommon;
class View;

class COMMONSAMPLE_EXPORT MDIWindow: public QMainWindow
{
    Q_OBJECT

public:
  MDIWindow( DocumentCommon* aDocument, QWidget* parent, Qt::WindowFlags wflags );
  MDIWindow( View* aView, DocumentCommon* aDocument, QWidget* parent, Qt::WindowFlags wflags );
  ~MDIWindow() override;

	DocumentCommon*            getDocument();
	void                       fitAll();
  QSize              sizeHint() const override;

signals:
  void                       selectionChanged();
  void                       message(const QString&, int );
	void                       sendCloseView(MDIWindow* theView);

public slots:
  void                       closeEvent(QCloseEvent* e) override;
  void                       onWindowActivated ();
  void                       dump();

protected:
  void                       createViewActions();
  void                       createRaytraceActions();

protected:
  DocumentCommon*            myDocument;
  View*                      myView;
};

#endif