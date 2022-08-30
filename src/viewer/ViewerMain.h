#ifndef ERIC_VIEWERMAIN_H
#define ERIC_VIEWERMAIN_H

#include "ViewerIncludes.h"
#include "ApplicationTut.h"

class ViewerMain {

public:
    static void start_viewer(std::list<viewerHelper::DisplayShapes> shapes);

    static void start_viewer(std::list<viewerHelper::DisplayShapes_SB> shapes);
};


#endif
