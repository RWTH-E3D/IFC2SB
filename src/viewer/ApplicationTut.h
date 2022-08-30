#ifndef APPLICATIONTUT_H
#define APPLICATIONTUT_H

#include "ViewerIncludes.h"


class ApplicationTut : public ApplicationCommonWindow {
Q_OBJECT

public:

    ApplicationTut(std::list<viewerHelper::DisplayShapes> shapes);

    ApplicationTut(std::list<viewerHelper::DisplayShapes_SB> shapes_spaces);

    ~ApplicationTut();

    std::list<viewerHelper::DisplayShapes> shapes;
    std::list<viewerHelper::DisplayShapes_SB> shapes_SB;
    std::list<viewerHelper::DisplayShapes_SB *> selected_shapes_SB;

    static QString getTutResourceDir();

    virtual void updateFileActions();

public slots:

    void onMakeBottleAction();

    void printInfo();

    void showSelected();

    void selectAll();

    void deselectAll();

    void changeTransparency();

    void changeTransparency(int theTrans);

    void changeColor();

    void changeEdgeColor();

    void findClicked();

    void findUnclicked();

private:
    void createMakeBottleOperation();

    void createMakeBottleOperation_SB();

    void update_selected_shapes();

    std::set<std::string> categories;
    QList<QCheckBox *> checkBoxes;

    QGroupBox *groupBox;
    QLineEdit *lineEdit;

private:
    QToolBar *myMakeBottleBar;
};

#endif