#include "ApplicationTut.h"
#include "Transparency.h"
#include <QColorDialog>

ApplicationTut::ApplicationTut(std::list<viewerHelper::DisplayShapes> _shapes) : ApplicationCommonWindow() {
    shapes = std::move(_shapes);
    createMakeBottleOperation();
}

ApplicationTut::ApplicationTut(std::list<viewerHelper::DisplayShapes_SB> _shapes_SB) : ApplicationCommonWindow() {
    shapes_SB = std::move(_shapes_SB);

    categories.insert("Normals");
    for (const auto &SB: shapes_SB)
        for (const auto &c: SB.categories)
            categories.insert(c);

    createMakeBottleOperation_SB();
}

ApplicationTut::~ApplicationTut() {}

void ApplicationTut::createMakeBottleOperation() {

    QPixmap MakeBottleIcon;
    QString dir = getTutResourceDir() + QString("/");
    MakeBottleIcon = QPixmap(dir + QObject::tr("ICON_MAKE_BOTTLE"));

    auto *MakeBottleAction = new QAction(MakeBottleIcon, QObject::tr("Show Shapes"), this);
    MakeBottleAction->setToolTip(QObject::tr("Show Shapes!"));
    MakeBottleAction->setStatusTip(QObject::tr("Show Shapes!"));
    MakeBottleAction->setShortcut(QObject::tr("CTRL+M"));

    connect(MakeBottleAction, SIGNAL(triggered()), this, SLOT(onMakeBottleAction()));

    myMakeBottleBar = addToolBar(tr("Mach mal eine Bottle"));
    insertToolBar(getCasCadeBar(), myMakeBottleBar);
    myMakeBottleBar->addAction(MakeBottleAction);
    myMakeBottleBar->hide();
}

void ApplicationTut::createMakeBottleOperation_SB() {

    QPixmap MakeBottleIcon;
    QString dir = getTutResourceDir() + QString("/");
    MakeBottleIcon = QPixmap(dir + QObject::tr("ICON_MAKE_BOTTLE"));

    myMakeBottleBar = addToolBar(tr("Mach mal eine Bottle"));
    addToolBar(Qt::LeftToolBarArea, myMakeBottleBar);
    /* toolbar break to prevent minimized widgets */
    insertToolBarBreak(myMakeBottleBar);

    /* custom group box */
    groupBox = new QGroupBox(tr(" Categories"));

    auto vbox = new QVBoxLayout();
    vbox->setSpacing(0);
    vbox->setContentsMargins(0, 0, 0, 0);
    vbox->setAlignment(Qt::AlignTop);

    auto *scrollArea = new QScrollArea();
    scrollArea->setStyleSheet("QScrollArea {background-color:gray;}");
    scrollArea->setWidgetResizable(true);

    auto *deselectButton = new QPushButton(tr("&Deselect"));
    connect(deselectButton, SIGNAL(clicked()), this, SLOT(findUnclicked()));

    auto *findButton = new QPushButton(tr("&Find"));
    connect(findButton, SIGNAL(clicked()), this, SLOT(findClicked()));

    auto *btn = new QPushButton(tr("Show selected"), this);
    connect(btn, SIGNAL(clicked()), this, SLOT(showSelected()));

    auto *btn1 = new QPushButton(tr("Check all"), this);
    connect(btn1, SIGNAL(clicked()), this, SLOT(selectAll()));

    auto *btn2 = new QPushButton(tr("Decheck all"), this);
    connect(btn2, SIGNAL(clicked()), this, SLOT(deselectAll()));

    auto *btn3 = new QPushButton(tr("Transparency"), this);
    connect(btn3, SIGNAL(clicked()), this, SLOT(changeTransparency()));

    auto *btn5 = new QPushButton(tr("Color"), this);
    connect(btn5, SIGNAL(clicked()), this, SLOT(changeColor()));

    auto *btn6 = new QPushButton(tr("Edge Color"), this);
    connect(btn6, SIGNAL(clicked()), this, SLOT(changeEdgeColor()));

    auto *btn4 = new QPushButton(tr("Print info of mouse-selected"), this);
    connect(btn4, SIGNAL(clicked()), this, SLOT(printInfo()));

    myMakeBottleBar->addWidget(btn);
    myMakeBottleBar->addWidget(btn1);
    myMakeBottleBar->addWidget(btn2);
    myMakeBottleBar->addWidget(btn3);
    myMakeBottleBar->addWidget(btn4);
    myMakeBottleBar->addWidget(btn5);
    myMakeBottleBar->addWidget(btn6);

    //auto *findLabel = new QLabel(tr("Enter key:"));
    //myMakeBottleBar->addWidget(findLabel);
    myMakeBottleBar->addWidget(findButton);
    myMakeBottleBar->addWidget(deselectButton);
    lineEdit = new QLineEdit;
    myMakeBottleBar->addWidget(lineEdit);

    for (const auto &it: categories) {
        QString qstr = QString::fromStdString(it);
        auto *checkBox = new QCheckBox(qstr);
        checkBox->setCheckable(true);
        //checkBox->setFixedHeight(10);

        QFont font;
        //font.setFamily(QStringLiteral("Arial"));
        font.setPointSize(10);
        checkBox->setFont(font);

        checkBoxes.push_back(checkBox);
        vbox->addWidget(checkBox);
    }

    groupBox->setLayout(vbox);
    scrollArea->setWidget(groupBox);
    myMakeBottleBar->addWidget(scrollArea);

    myMakeBottleBar->hide();
}

void ApplicationTut::updateFileActions() {
    if (getWorkspace()->subWindowList().isEmpty()) {
        if (!isDocument())
            myMakeBottleBar->show();
        else
            myMakeBottleBar->hide();
    }
    ApplicationCommonWindow::updateFileActions();
}

void ApplicationTut::onMakeBottleAction() {
    QMdiArea *ws = ApplicationCommonWindow::getWorkspace();
    auto *doc = (DocumentTut *) (qobject_cast<MDIWindow *>(ws->activeSubWindow()->widget())->getDocument());
    statusBar()->showMessage(QObject::tr("INF_MAKE_BOTTLE"), 5000);
    doc->onMakeBottle(shapes);
    statusBar()->showMessage(QObject::tr("INF_DONE"));
}

QString ApplicationTut::getTutResourceDir() {
    static QString resDir(OSD_Environment("CSF_TutorialResourcesDefaults").Value().ToCString());
    if (resDir.isEmpty())
        resDir = QString(OSD_Environment("CSF_OCCTResourcePath").Value().ToCString()) + "/samples";
    return resDir;
}

void ApplicationTut::printInfo() {

    QMdiArea *ws = ApplicationCommonWindow::getWorkspace();
    auto *doc = (DocumentTut *) (qobject_cast<MDIWindow *>(ws->activeSubWindow()->widget())->getDocument());

    doc->getContext()->InitSelected();
    Handle(AIS_InteractiveObject) Current = doc->getContext()->SelectedInteractive();

    if (Current.IsNull()) return;

    for (const auto &shape: shapes_SB) {
        if (shape.AIS == Current) {
            std::cout << "Info: " << std::endl;
            for (auto &c: shape.categories)
                std::cout << "\t" << c << std::endl;

            TopExp_Explorer Ex;
            std::list<gp_Pnt> pnts;

            std::set<Standard_Integer> hashes;
            for (Ex.Init(shape.shape, TopAbs_VERTEX); Ex.More(); Ex.Next()) {
                auto search = hashes.find(Ex.Current().HashCode(INT_MAX));
                if (search == hashes.end()) {
                    hashes.insert(Ex.Current().HashCode(INT_MAX));
                    pnts.push_back(BRep_Tool::Pnt(TopoDS::Vertex(Ex.Current())));
                }
            }

            for (auto &pnt: pnts)
                std::cerr << "\t\t" << pnt.X() << " " << pnt.Y() << " " << pnt.Z() << "\n";

            break;
        }

    }
}

void ApplicationTut::showSelected() {
    QMdiArea *ws = ApplicationCommonWindow::getWorkspace();
    auto *doc = (DocumentTut *) (qobject_cast<MDIWindow *>(ws->activeSubWindow()->widget())->getDocument());
    statusBar()->showMessage(QObject::tr("INF_MAKE_BOTTLE"), 5000);

    std::set<std::string> selected_checkBoxes;

    for (auto i: checkBoxes)
        if (i->isChecked())
            selected_checkBoxes.insert(i->text().toStdString());

    doc->create_AIS_of_SB(shapes_SB, selected_checkBoxes);
    statusBar()->showMessage(QObject::tr("INF_DONE"));
}

void ApplicationTut::selectAll() { for (auto i: groupBox->findChildren<QCheckBox *>()) i->setChecked(true); }

void ApplicationTut::deselectAll() { for (auto i: groupBox->findChildren<QCheckBox *>()) i->setChecked(false); }

void ApplicationTut::update_selected_shapes() {

    std::set<std::string> selected_checkBoxes;
    selected_shapes_SB.clear();

    for (auto i: checkBoxes)
        if (i->isChecked())
            selected_checkBoxes.insert(i->text().toStdString());

    for (auto &shape: shapes_SB)
        for (const auto &c: shape.categories)
            if (selected_checkBoxes.find(c) != selected_checkBoxes.end()) {
                selected_shapes_SB.push_back(&shape);
                break;
            }
}

void ApplicationTut::changeColor() {

    update_selected_shapes();

    QColor aColor;
    Standard_Real R1;
    Standard_Real G1;
    Standard_Real B1;
    aColor.setRgb((Standard_Integer) (R1 * 255), (Standard_Integer) (G1 * 255), (Standard_Integer) (B1 * 255));
    QColor aRetColor = QColorDialog::getColor(aColor);

    if (aRetColor.isValid()) {
        R1 = aRetColor.red() / 255.;
        G1 = aRetColor.green() / 255.;
        B1 = aRetColor.blue() / 255.;
        Quantity_Color clr(R1, G1, B1, Quantity_TOC_RGB);
        for (auto &shape: selected_shapes_SB) {
            shape->AIS->SetColor(clr);
            shape->AIS->Attributes()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
        }
    }
}

void ApplicationTut::changeEdgeColor() {

    update_selected_shapes();

    QColor aColor;
    Standard_Real R1;
    Standard_Real G1;
    Standard_Real B1;
    aColor.setRgb((Standard_Integer) (R1 * 255), (Standard_Integer) (G1 * 255), (Standard_Integer) (B1 * 255));
    QColor aRetColor = QColorDialog::getColor(aColor);

    if (aRetColor.isValid()) {
        R1 = aRetColor.red() / 255.;
        G1 = aRetColor.green() / 255.;
        B1 = aRetColor.blue() / 255.;
        Quantity_Color clr(R1, G1, B1, Quantity_TOC_RGB);
        for (auto &shape: selected_shapes_SB)
            shape->AIS->Attributes()->FaceBoundaryAspect()->SetColor(clr);
    }
}

void ApplicationTut::changeTransparency() {

    update_selected_shapes();

    auto *aDialog = new DialogTransparency();
    connect(aDialog, SIGNAL(sendTransparencyChanged(int)), this, SLOT(changeTransparency(int)));
    aDialog->exec();
}

void ApplicationTut::changeTransparency(int theTrans) {

    QMdiArea *ws = ApplicationCommonWindow::getWorkspace();
    auto *doc = (DocumentTut *) (qobject_cast<MDIWindow *>(ws->activeSubWindow()->widget())->getDocument());

    auto t = ((Standard_Real) theTrans) / 10.0;

    for (auto &shape: selected_shapes_SB)
        shape->AIS->SetTransparency(t);

    doc->getContext()->UpdateCurrentViewer();
}

void ApplicationTut::findClicked() {
    QString text = lineEdit->text();

    // do nothing if the string is empty
    if (text.isEmpty())
        return;

    for (auto checkbox: checkBoxes)
        if (checkbox->text().contains(text))
            checkbox->setChecked(true);
    //else
    //   checkbox->setChecked(false);

}

void ApplicationTut::findUnclicked() {
    QString text = lineEdit->text();

    // do nothing if the string is empty
    if (text.isEmpty()) return;

    for (auto checkbox: checkBoxes)
        if (checkbox->text().contains(text))
            checkbox->setCheckState(Qt::Unchecked);
}