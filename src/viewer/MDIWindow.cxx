#include "MDIWindow.h"

#include "View.h"

MDIWindow::MDIWindow(View* aView,
                     DocumentCommon* aDocument, 
                     QWidget* parent, 
                     Qt::WindowFlags wflags )
: QMainWindow( parent, wflags )
{
  myView = aView;
  myDocument = aDocument;
}

MDIWindow::MDIWindow( DocumentCommon* aDocument, QWidget* parent, Qt::WindowFlags wflags)
: QMainWindow( parent, wflags )
{
  auto *vb = new QFrame( this );

  auto *layout = new QVBoxLayout( vb );
  layout->setMargin( 0 );

  vb->setFrameStyle( QFrame::StyledPanel | QFrame::Sunken );

  setCentralWidget( vb );

  myDocument = aDocument;
  myView = new View (myDocument->getContext(), vb);
  layout->addWidget (myView);

  connect( myView, SIGNAL( selectionChanged() ),
           this,   SIGNAL( selectionChanged() ) );

  createViewActions();
  createRaytraceActions();

  resize( sizeHint() );

  setFocusPolicy( Qt::StrongFocus );
}

MDIWindow::~MDIWindow()
= default;

DocumentCommon* MDIWindow::getDocument()
{
  return myDocument;
}

void MDIWindow::closeEvent(QCloseEvent* )
{
  emit sendCloseView(this);
}

void MDIWindow::fitAll()
{
  myView->fitAll();
}

void MDIWindow::createViewActions()
{
  // populate a tool bar with some actions
  QToolBar* aToolBar = addToolBar( tr( "View Operations" ) );
  
  QList<QAction*>* aList = myView->getViewActions();
  aToolBar->addActions( *aList );

  aToolBar->toggleViewAction()->setVisible(false);
  aList->at(View::ViewHlrOffId)->setChecked( true );
}

void MDIWindow::createRaytraceActions()
{
  // populate a tool bar with some actions
  QToolBar* aToolBar = addToolBar( tr( "Ray-tracing Options" ) );
  
  QList<QAction*>* aList = myView->getRaytraceActions();
  aToolBar->addActions( *aList );

  aToolBar->toggleViewAction()->setVisible (true);
  aList->at (View::ToolRaytracingId)->setChecked (false);
  aList->at (View::ToolShadowsId)->setChecked (true);
  aList->at (View::ToolReflectionsId)->setChecked (false);
  aList->at (View::ToolAntialiasingId)->setChecked (false);
}

void MDIWindow::onWindowActivated ()
{
  getDocument()->getApplication()->onSelectionChanged();
}

void MDIWindow::dump()
{
  QString filter = "Images Files (*.bmp *.ppm *.png *.jpg *.tiff *.tga *.gif *.exr)";
  QFileDialog fd ( nullptr );
  fd.setModal( true );
  fd.setNameFilter ( filter );
  fd.setWindowTitle( QObject::tr("INF_APP_EXPORT") );
  fd.setFileMode( QFileDialog::AnyFile );
  int ret = fd.exec(); 

  /* update the desktop after the dialog is closed */
  qApp->processEvents();

  QStringList fileNames;
  fileNames = fd.selectedFiles();

  QString file ( (ret == QDialog::Accepted && !fileNames.isEmpty() )? fileNames[0] : QString::null);
  if ( !file.isNull() )
  {
    // QApplication::setOverrideCursor( Qt::WaitCursor );
    if ( !QFileInfo( file ).completeSuffix().length() )
      file += QString( ".bmp" );

    const TCollection_AsciiString anUtf8Path (file.toUtf8().data());

    bool res = myView->dump( anUtf8Path.ToCString() );
    QApplication::restoreOverrideCursor();                
    if ( !res )
    {
      QWidgetList list = qApp->allWidgets();
      QWidget* mainWidget = nullptr;
      for(auto i : list)
      {
        if( qobject_cast<ApplicationCommonWindow*>( i ) )
        mainWidget = qobject_cast<ApplicationCommonWindow*>( i );
      }

      QMessageBox::information ( mainWidget, QObject::tr("TIT_ERROR"), QObject::tr("INF_ERROR"), QObject::tr("BTN_OK"),
                                 QString::null, QString::null, 0, 0 );
      qApp->processEvents();
    }
  }
}

QSize MDIWindow::sizeHint() const
{
  return { 450, 300 };
}
