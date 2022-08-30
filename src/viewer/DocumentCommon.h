#ifndef DOCUMENTCOMMON_H
#define DOCUMENTCOMMON_H

#include "MDIWindow.h"

#include <Standard_WarningsDisable.hxx>
#include <QObject>
#include <QList>
#include <Standard_WarningsRestore.hxx>

#include <AIS_InteractiveContext.hxx>
#include <V3d_Viewer.hxx>
#include <unordered_map>

class ApplicationCommonWindow;

class COMMONSAMPLE_EXPORT DocumentCommon : public QObject {
Q_OBJECT

public:
    DocumentCommon(int, ApplicationCommonWindow *);

    ~DocumentCommon() override;

    ApplicationCommonWindow *getApplication();

    Handle(AIS_InteractiveContext) getContext();

    void removeView(MDIWindow *);

    void removeViews();

    int countOfWindow();

    void fitAll();

    bool firstStart;
    std::unordered_map<std::string, std::array<double, 3>> str2clr;

protected:
    virtual MDIWindow *createNewMDIWindow();

signals:

    void selectionChanged();

    void sendCloseDocument(DocumentCommon *);

public slots:

    virtual void onCloseView(MDIWindow *);

    virtual void onCreateNewView();

    virtual void onMaterial();

    virtual void onMaterial(int);

    virtual void onDelete();

    void onWireframe();

    void onShading();

    void onColor();

    void onTransparency();

    void onTransparency(int);

private:
    Handle(V3d_Viewer) Viewer(Standard_ExtString theName,
                              Standard_CString theDomain,
                              Standard_Real theViewSize,
                              V3d_TypeOfOrientation theViewProj,
                              Standard_Boolean theComputedMode,
                              Standard_Boolean theDefaultComputedMode);

protected:
    ApplicationCommonWindow *myApp;
    QList<MDIWindow *> myViews;
    Handle(V3d_Viewer) myViewer;
    Handle(AIS_InteractiveContext) myContext;
    int myIndex;
    int myNbViews;
};

#endif


