#include "ViewerMain.h"

#include <utility>

void ViewerMain::start_viewer(std::list<viewerHelper::DisplayShapes> shapes) {

    int argc_changed = 1;
    char *argv_changed[] = {(char *) ("SB Viewer ErF")};

    //printf("argc_changed: %d\n", argc_changed);
    //for (int i = 0; i < argc_changed; i++) {
    //    printf("argv[%d]: %s\n", i, argv_changed[i]);
    //}

    TCollection_AsciiString aPlugindsDirName = OSD_Environment("QTDIR").Value();
    if (!aPlugindsDirName.IsEmpty())
        QApplication::addLibraryPath(QString(aPlugindsDirName.ToCString()) + "/plugins");

    QApplication a(argc_changed, argv_changed);

    QString resDir = ApplicationCommonWindow::getResourceDir();
    QString resTutDir = ApplicationTut::getTutResourceDir();

    QTranslator strTrans(nullptr);
    Standard_Boolean isOK = strTrans.load("Common-string", resDir);
    if (isOK)
        QApplication::installTranslator(&strTrans);

    QTranslator iconTrans(nullptr);
    isOK = iconTrans.load("Common-icon", resDir);
    if (isOK)
        QApplication::installTranslator(&iconTrans);

    QTranslator strTutTrans(nullptr);
    isOK = strTutTrans.load("Tutorial-string", resTutDir);
    if (isOK)
        QApplication::installTranslator(&strTutTrans);

    QTranslator iconTutTrans(nullptr);
    isOK = iconTutTrans.load("Tutorial-icon", resTutDir);
    if (isOK)
        QApplication::installTranslator(&iconTutTrans);

    QObject::connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

    auto *mw = new ApplicationTut(std::move(shapes));
    mw->setWindowTitle(QObject::tr("Space Boundary Viewer"));
    //mw->resize(3800,3000);
    mw->setWindowState(Qt::WindowMaximized);
    QString aResName(resDir + QString("/") + QObject::tr("ICON_SAMPLE"));
    mw->setWindowIcon(QPixmap(aResName));

    mw->show();

    QApplication::exec();
}

void ViewerMain::start_viewer(std::list<viewerHelper::DisplayShapes_SB> shapes) {

    int argc_changed = 1;
    char *argv_changed[] = {(char *) ("SB Viewer")};

    TCollection_AsciiString aPlugindsDirName = OSD_Environment("QTDIR").Value();
    if (!aPlugindsDirName.IsEmpty())
        QApplication::addLibraryPath(QString(aPlugindsDirName.ToCString()) + "/plugins");

    QApplication a(argc_changed, argv_changed);

    QString resDir = ApplicationCommonWindow::getResourceDir();
    QString resTutDir = ApplicationTut::getTutResourceDir();

    QTranslator strTrans(nullptr);
    Standard_Boolean isOK = strTrans.load("Common-string", resDir);
    if (isOK)
        QApplication::installTranslator(&strTrans);

    QTranslator iconTrans(nullptr);
    isOK = iconTrans.load("Common-icon", resDir);
    if (isOK)
        QApplication::installTranslator(&iconTrans);

    QTranslator strTutTrans(nullptr);
    isOK = strTutTrans.load("Tutorial-string", resTutDir);
    if (isOK)
        QApplication::installTranslator(&strTutTrans);

    QTranslator iconTutTrans(nullptr);
    isOK = iconTutTrans.load("Tutorial-icon", resTutDir);
    if (isOK)
        QApplication::installTranslator(&iconTutTrans);

    QObject::connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));

    auto *mw = new ApplicationTut(std::move(shapes));
    mw->setWindowTitle(QObject::tr("Space Boundary Viewer"));
    mw->setWindowState(Qt::WindowMaximized);
    QString aResName(resDir + QString("/") + QObject::tr("ICON_SAMPLE"));
    mw->setWindowIcon(QPixmap(aResName));

    mw->show();

    QApplication::exec();
}