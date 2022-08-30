#ifndef DOCUMENTTUT_H
#define DOCUMENTTUT_H

#include "ViewerIncludes.h"


class SurfConstruction;

class DocumentTut : public DocumentCommon {
Q_OBJECT

public:
    DocumentTut(int, ApplicationCommonWindow *);

    ~DocumentTut() override;

    void onMakeBottle(const std::list<viewerHelper::DisplayShapes>& shapes);

    void create_AIS_of_SB(std::list<viewerHelper::DisplayShapes_SB> &shapes, std::set<std::string> selected_checkBoxes);

    void show_and_hide(std::list<viewerHelper::DisplayShapes_SB> &shapes, std::set<std::string>& selected_checkBoxes);

private:

    static TopoDS_Shape construct_vectors(gp_Pnt pnt, gp_Vec vec, double scale);
};

#endif
