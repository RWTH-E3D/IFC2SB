#ifndef VIEW_H
#define VIEW_H

#include "ViewerIncludes.h"


class TopoDS_Shape;

class QRubberBand;

//class COMMONSAMPLE_EXPORT View: public QWidget
class View : public QWidget {
Q_OBJECT
protected:
    enum CurrentAction3d {
        CurAction3d_Nothing, CurAction3d_DynamicZooming,
        CurAction3d_WindowZooming, CurAction3d_DynamicPanning,
        CurAction3d_GlobalPanning, CurAction3d_DynamicRotation
    };

public:
    enum ViewAction {
        ViewFitAllId, ViewFitAreaId, ViewZoomId, ViewPanId, ViewGlobalPanId,
        ViewFrontId, ViewBackId, ViewTopId, ViewBottomId, ViewLeftId, ViewRightId,
        ViewAxoId, ViewRotationId, ViewResetId, ViewHlrOffId, ViewHlrOnId, ViewHideTriedron,
        ViewTakeSnapshot, ViewClip, BGColorId
    };
    enum RaytraceAction {
        ToolRaytracingId, ToolShadowsId, ToolReflectionsId, ToolAntialiasingId
    };

    View(Handle(AIS_InteractiveContext) theContext, QWidget *parent);

    ~View() override;

    virtual void init();

    bool dump(Standard_CString theFile);

    QList<QAction *> *getViewActions();

    QList<QAction *> *getRaytraceActions();

    void noActiveActions();

    void EnableRaytracing();

    void DisableRaytracing();

    void SetRaytracedShadows(bool theState);

    void SetRaytracedReflections(bool theState);

    void SetRaytracedAntialiasing(bool theState);

    Standard_EXPORT

    Standard_EXPORT

    Standard_EXPORT

    QPaintEngine *paintEngine() const override;

signals:

    void selectionChanged();

public slots:

    void fitAll();

    void fitArea();

    void zoom();

    void pan();

    void globalPan();

    void front();

    void back();

    void top();

    void bottom();

    void left();

    void right();

    void axo();

    void rotation();

    void reset();

    void hlrOn();

    void hlrOff();

    void hideTriedron();

    void showTriedron();

    void takeSnapshot();

    void clip();

    void updateToggled(bool);

    void onBackground();

    void onEnvironmentMap();

    void onRaytraceAction();

protected:
    void paintEvent(QPaintEvent *) override;

    void resizeEvent(QResizeEvent *) override;

    void mousePressEvent(QMouseEvent *) override;

    void mouseReleaseEvent(QMouseEvent *) override;

    void mouseMoveEvent(QMouseEvent *) override;

    virtual void addItemInPopup(QMenu *);

    Handle(AIS_InteractiveContext) &getContext();

    void activateCursor(CurrentAction3d);

    void Popup(int x, int y);

    virtual void onLButtonDown(int nFlags, QPoint point);

    virtual void onMButtonDown(int nFlags, QPoint point);

    virtual void onRButtonDown(int nFlags, QPoint point);

    virtual void onLButtonUp(Qt::MouseButtons nFlags, QPoint point);

    virtual void onMButtonUp(Qt::MouseButtons nFlags, QPoint point);

    virtual void onRButtonUp(Qt::MouseButtons nFlags, QPoint point);

    virtual void onMouseMove(Qt::MouseButtons nFlags, QPoint point);

    void wheelEvent(QWheelEvent *event) override;

private:
    static void initCursors();

    void initViewActions();

    void initRaytraceActions();

    void DragEvent(int x, int y, int TheState);

    void InputEvent(int x, int y);

    void MoveEvent(int x, int y);

    void MultiMoveEvent(int x, int y);

    void MultiDragEvent(int x, int y, int TheState);

    void MultiInputEvent(int x, int y);

    void DrawRectangle(int MinX, int MinY,
                       int MaxX, int MaxY, bool Draw);

private:
    bool myIsRaytracing;
    bool myIsShadowsEnabled;
    bool myIsReflectionsEnabled;
    bool myIsAntialiasingEnabled;

    // set when a rect is used for selection or magnify
    Handle(V3d_View) myView;
    Handle(AIS_InteractiveContext) myContext;
    CurrentAction3d myCurrentMode;
    Standard_Integer myXmin;
    Standard_Integer myYmin;
    Standard_Integer myXmax;
    Standard_Integer myYmax;
    Standard_Real myCurZoom;
    Standard_Boolean myHlrModeIsOn;
    QList<QAction *> *myViewActions;
    QList<QAction *> *myRaytraceActions;
    QMenu *myBackMenu;
    QRubberBand *myRectBand; //!< selection rectangle rubber band
};

#endif