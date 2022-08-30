#include "DocumentCommon.h"

#include "ApplicationCommon.h"
#include "Transparency.h"
#include "Material.h"

#include <QColorDialog>

#include <Aspect_DisplayConnection.hxx>
#include <OpenGl_GraphicDriver.hxx>

#if !defined(_WIN32) && !defined(__WIN32__) && (!defined(__APPLE__) || defined(MACOSX_USE_GLX))

#endif

// =======================================================================
// function : Viewer
// purpose  :
// =======================================================================
Handle(V3d_Viewer) DocumentCommon::Viewer(const Standard_ExtString,
                                          const Standard_CString,
                                          const Standard_Real theViewSize,
                                          const V3d_TypeOfOrientation theViewProj,
                                          const Standard_Boolean theComputedMode,
                                          const Standard_Boolean theDefaultComputedMode) {
    static Handle(OpenGl_GraphicDriver) aGraphicDriver;

    if (aGraphicDriver.IsNull()) {
        Handle(Aspect_DisplayConnection) aDisplayConnection;
#if !defined(_WIN32) && !defined(__WIN32__) && (!defined(__APPLE__) || defined(MACOSX_USE_GLX))
        aDisplayConnection = new Aspect_DisplayConnection(OSD_Environment("DISPLAY").Value());
#endif
        aGraphicDriver = new OpenGl_GraphicDriver(aDisplayConnection);
    }

    Handle(V3d_Viewer) aViewer = new V3d_Viewer(aGraphicDriver);

    aViewer->SetDefaultViewSize(theViewSize);
    aViewer->SetDefaultViewProj(theViewProj);
    aViewer->SetComputedMode(theComputedMode);
    aViewer->SetDefaultComputedMode(theDefaultComputedMode);
//    Quantity_Color clr =aViewer->DefaultBackgroundColor();
//    std::cout << "aVIEWER "<<clr.Red()<<"\n";
//    aViewer->SetDefaultBackgroundColor(clr);
    return aViewer;
}

DocumentCommon::DocumentCommon(const int theIndex, ApplicationCommonWindow *app)
        : QObject(app),
          myApp(app),
          myIndex(theIndex),
          myNbViews(0) {
    TCollection_ExtendedString a3DName("Visu3D");

    myViewer = Viewer(a3DName.ToExtString(), "", 1000.0, V3d_XposYnegZpos, Standard_True, Standard_True);

    myViewer->SetDefaultLights();
    myViewer->SetLightOn();

    myContext = new AIS_InteractiveContext(myViewer);

    firstStart = true;

    str2clr["IfcWindow"] = std::array<double, 3>{1, 0, 0};
    str2clr["IfcWall"] = std::array<double, 3>{.8, .8, .8};
    str2clr["IfcSite"] = std::array<double, 3>{.75, .8, .65};
    str2clr["IfcSlab"] = std::array<double, 3>{.4, .4, .4};
    str2clr["IfcWallStandardCase"] = std::array<double, 3>{.9, .9, .9};
    str2clr["IfcWall"] = std::array<double, 3>{.9, .9, .9};
    str2clr["IfcWindow"] = std::array<double, 3>{.75, .8, .75};
    str2clr["IfcDoor"] = std::array<double, 3>{.55, .3, .15};
    str2clr["IfcBeam"] = std::array<double, 3>{.75, .7, .7};
    str2clr["IfcRailing"] = std::array<double, 3>{.65, .6, .6};
    str2clr["IfcMember"] = std::array<double, 3>{.65, .6, .6};
    str2clr["IfcPlate"] = std::array<double, 3>{.8, .8, .8};
    str2clr["IfcSpace"] = std::array<double, 3>{.9, .0, .1};
    str2clr["IfcFurniture"] = std::array<double, 3>{.7, .4, .3};
    str2clr["IfcTank"] = std::array<double, 3>{1.0, 1.0, 1.0};
    str2clr["IfcFlowFitting"] = std::array<double, 3>{.0, .0, .5};
    str2clr["IfcRelSpaceBoundary"] = std::array<double, 3>{.5, .8, .9};
    str2clr["RED"] = std::array<double, 3>{1, 0, 0};
    str2clr["BLUE"] = std::array<double, 3>{0, 0, 1};
    str2clr["DARKBLUE"] = std::array<double, 3>{0, 0, 0.5};
    str2clr["GREEN"] = std::array<double, 3>{0, 1, 0};
    str2clr["GREY"] = std::array<double, 3>{0.8, 0.8, 0.8};
    str2clr["MIDGREY"] = std::array<double, 3>{0.5, 0.5, 0.5};
    str2clr["DARKGREY"] = std::array<double, 3>{0.3, 0.3, 0.3};
    str2clr["DARKERGREY"] = std::array<double, 3>{0.1, 0.1, 0.1};
    str2clr["YELLOW"] = std::array<double, 3>{1, 1, 0};
    str2clr["BLACK"] = std::array<double, 3>{0, 0, 0};
    str2clr["LIGHTGREEN"] = std::array<double, 3>{0.56, 0.93, 0.56};
    str2clr["DARKGREEN"] = std::array<double, 3>{0.95, 0.97, 0.91};
    str2clr["LIGHTRED"] = std::array<double, 3>{1, 0.52, 0.52};
    str2clr["DARKRED"] = std::array<double, 3>{0.5, 0, 0};
    str2clr["LIGHTBROWN"] = std::array<double, 3>{1, 0.89, 0.77};
    str2clr["BROWN"] = std::array<double, 3>{0.55, 0.27, 0.07};
    str2clr["ORANGE"] = std::array<double, 3>{1, 0.56, 0.0};
    str2clr["CYAN"] = std::array<double, 3>{0, 1, 1};
    str2clr["PURPLE"] = std::array<double, 3>{0.31, 0.15, 0.54};
    str2clr["EGG"] = std::array<double, 3>{0.81, 0.71, 0.4};
    str2clr["OPENSTUDIORED"] = std::array<double, 3>{0.4, 0.2, 0.2};
}

DocumentCommon::~DocumentCommon() {
}

ApplicationCommonWindow *DocumentCommon::getApplication() {
    return myApp;
}

MDIWindow *DocumentCommon::createNewMDIWindow() {
    QMdiArea *ws = myApp->getWorkspace();
    return new MDIWindow(this, ws, 0);
}

void DocumentCommon::onCreateNewView() {
    QMdiArea *ws = myApp->getWorkspace();
    MDIWindow *w = createNewMDIWindow();

    if (!w)
        return;

    ws->addSubWindow(w);
    myViews.append(w);

    connect(w, SIGNAL(selectionChanged()),
            this, SIGNAL(selectionChanged()));
    connect(w, SIGNAL(message( const QString&, int )),
            myApp->statusBar(), SLOT(showMessage( const QString&, int )));
    connect(w, SIGNAL(sendCloseView(MDIWindow * )),
            this, SLOT(onCloseView(MDIWindow * )));

    QString aName;
    w->setWindowTitle(aName.sprintf("Fenster Nummer %d:%d", myIndex, ++myNbViews));
    w->setWindowState(Qt::WindowMaximized);

    QString dir = ApplicationCommonWindow::getResourceDir() + QString("/");

    w->setWindowIcon(QPixmap(dir + QObject::tr("ICON_DOC")));

    if (ws->subWindowList().isEmpty()) {
        // Due to strange Qt4.2.3 feature the child window icon is not drawn
        // in the main menu if showMaximized() is called for a non-visible child window
        // Therefore calling show() first...
        w->show();
        w->showMaximized();
    } else
        w->show();

    w->setFocus();

    getApplication()->onSelectionChanged();
}

void DocumentCommon::onCloseView(MDIWindow *theView) {
    removeView(theView);
    if (countOfWindow() == 0)
            emit sendCloseDocument(this);
}

void DocumentCommon::removeView(MDIWindow *theView) {
    if (myViews.count(theView)) {
        myViews.removeAll(theView);
        delete theView;
    }
}

void DocumentCommon::removeViews() {
    while (myViews.count()) {
        removeView(myViews.first());
    }
}

int DocumentCommon::countOfWindow() {
    return myViews.count();
}

Handle(AIS_InteractiveContext) DocumentCommon::getContext() {
    return myContext;
}

void DocumentCommon::fitAll() {
    QList<MDIWindow *>::iterator i;
    for (i = myViews.begin(); i != myViews.end(); i++)
        (*i)->fitAll();
}

void DocumentCommon::onWireframe() {
    // QApplication::setOverrideCursor(Qt::WaitCursor);
    for (myContext->InitSelected(); myContext->MoreSelected(); myContext->NextSelected())
        myContext->SetDisplayMode(myContext->SelectedInteractive(), 0, false);
    myContext->UpdateCurrentViewer();
    getApplication()->onSelectionChanged();
    QApplication::restoreOverrideCursor();
}

void DocumentCommon::onShading() {
    // QApplication::setOverrideCursor(Qt::WaitCursor);
    for (myContext->InitSelected(); myContext->MoreSelected(); myContext->NextSelected())
        myContext->SetDisplayMode(myContext->SelectedInteractive(), 1, false);
    myContext->UpdateCurrentViewer();
    getApplication()->onSelectionChanged();
    QApplication::restoreOverrideCursor();
}

void DocumentCommon::onColor() {
    QColor aColor;
    myContext->InitSelected();
    Handle(AIS_InteractiveObject) Current = myContext->SelectedInteractive();
    if (Current->HasColor()) {
        Quantity_Color aShapeColor;
        myContext->Color(Current, aShapeColor);
        aColor.setRgb((Standard_Integer) (aShapeColor.Red() * 255), (Standard_Integer) (aShapeColor.Green() * 255),
                      (Standard_Integer) (aShapeColor.Blue() * 255));
    } else
        aColor.setRgb(255, 255, 255);

    QColor aRetColor = QColorDialog::getColor(aColor);
    if (aRetColor.isValid()) {
        Quantity_Color color(aRetColor.red() / 255., aRetColor.green() / 255.,
                             aRetColor.blue() / 255., Quantity_TOC_RGB);
        for (; myContext->MoreSelected(); myContext->NextSelected()) {
            myContext->SetColor(myContext->SelectedInteractive(), color, Standard_False);
            myContext->SelectedInteractive()->Attributes()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
        }
        myContext->UpdateCurrentViewer();
    }
}

void DocumentCommon::onMaterial(int theMaterial) {
    for (myContext->InitSelected(); myContext->MoreSelected(); myContext->NextSelected())
        myContext->SetMaterial(myContext->SelectedInteractive(), (Graphic3d_NameOfMaterial) theMaterial, Standard_False);
    myContext->UpdateCurrentViewer();
}

void DocumentCommon::onMaterial() {
    DialogMaterial *m = new DialogMaterial();
    connect(m, SIGNAL(sendMaterialChanged(int)), this, SLOT(onMaterial(int)));
    m->exec();
}

void DocumentCommon::onTransparency(int theTrans) {
    for (myContext->InitSelected(); myContext->MoreSelected(); myContext->NextSelected())
        myContext->SetTransparency(myContext->SelectedInteractive(), ((Standard_Real) theTrans) / 10.0, Standard_False);
    myContext->UpdateCurrentViewer();
}

void DocumentCommon::onTransparency() {
    DialogTransparency *aDialog = new DialogTransparency();
    connect(aDialog, SIGNAL(sendTransparencyChanged(int)), this, SLOT(onTransparency(int)));
    aDialog->exec();
}

void DocumentCommon::onDelete() {
    myContext->EraseSelected(Standard_False);
    myContext->ClearSelected(Standard_False);
    myContext->UpdateCurrentViewer();
    getApplication()->onSelectionChanged();
}
