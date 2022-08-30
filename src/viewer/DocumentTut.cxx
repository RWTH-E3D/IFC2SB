#include "DocumentTut.h"

TopoDS_Shape MakeBottle(Standard_Real myWidth, Standard_Real myHeight, Standard_Real myThickness);

DocumentTut::DocumentTut(const int theIndex, ApplicationCommonWindow *app) : DocumentCommon(theIndex, app) {}

DocumentTut::~DocumentTut() = default;

void DocumentTut::onMakeBottle(const std::list<viewerHelper::DisplayShapes> &shapes) {

    const Handle(Prs3d_Drawer) selectionStyle = new Prs3d_Drawer();
    selectionStyle->SetColor(Quantity_NOC_RED);
    getContext()->SetSelectionStyle(selectionStyle);
    getContext()->SetHighlightStyle(selectionStyle);
    getContext()->DefaultDrawer()->SetFaceBoundaryDraw(Standard_True);
    //Handle_Prs3d_LineAspect lineAspectBox = new Prs3d_LineAspect(Quantity_NOC_BLACK, Aspect_TOL_SOLID, 2.0);
    //getContext()->DefaultDrawer()->SetFaceBoundaryAspect(lineAspectBox);
    //getContext()->DefaultDrawer()->FaceBoundaryAspect()->SetWidth(10);
    //getContext()->DefaultDrawer()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
    //getContext()->DefaultDrawer()->FaceBoundaryAspect()->SetTypeOfLine(Aspect_TOL_SOLID);
    // turn up tesselation defaults, which are too conversative...
    auto chord_dev = getContext()->DefaultDrawer()->MaximalChordialDeviation() / 10.;
    getContext()->DefaultDrawer()->SetMaximalChordialDeviation(chord_dev);
    //getContext()->DefaultDrawer()
    //this->v3dview->SetBackgroundColor(Quantity_NOC_BLACK);

    for (auto &shape: shapes) {

        // skip null shapes to prevent crashing of viewer
        if (shape.shape.IsNull()) {
            std::cerr << "VIEWER: Shape is Null " << shape.shape.HashCode(INT_MAX) << std::endl;
            continue;
        }

        //***********************************************************************************
        // QApplication::setOverrideCursor(Qt::WaitCursor);
        TopoDS_Shape shape_to_display;

        if (shape.shape.ShapeType() == TopAbs_EDGE) {
            gp_Pnt P1 = BRep_Tool::Pnt(TopExp::FirstVertex(TopoDS::Edge(shape.shape), Standard_True));
            gp_Pnt P2 = BRep_Tool::Pnt(TopExp::LastVertex(TopoDS::Edge(shape.shape), Standard_True));
            shape_to_display = BRepPrimAPI_MakeCylinder(gp_Ax2(P1, gp_Dir(gp_Vec(P1, P2))), shape.edge_thickness, P1.Distance(P2)).Shape();
        } else if (shape.shape.ShapeType() == TopAbs_VERTEX) {
            gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(shape.shape));
            shape_to_display = BRepPrimAPI_MakeSphere(P, 0.06).Shape();
        } else
            shape_to_display = shape.shape;

        Handle_AIS_ColoredShape AISBauteil = new AIS_ColoredShape(shape_to_display);
        //Handle(AIS_Shape) AISBauteil = new AIS_Shape(it->shape);
        Graphic3d_NameOfMaterial mat = Graphic3d_NOM_GOLD;
        if (shape.shape.ShapeType() == TopAbs_VERTEX || shape.shape.ShapeType() == TopAbs_EDGE)
            mat = Graphic3d_NOM_PLASTER;
        getContext()->SetMaterial(AISBauteil, mat, Standard_False);
        //getContext()->SetTransparency(AISBauteil, it->transparency, Standard_False);
        getContext()->SetDisplayMode(AISBauteil, 1, Standard_False);
        //const Handle(AIS_InteractiveObject) &anIOAISBottle = AISBottle;
        //getContext()->SetSelected(anIOAISBottle, Standard_False);
        //***********************************************************************************


        //***********************************************************************************
        // color
        if (shape.clr_by_string) {
            if (str2clr.find(shape.clr_string) == str2clr.end())
                AISBauteil->SetColor(shape.clr);
                //getContext()->SetColor(AISBauteil, it->clr, Standard_False);
            else {
                std::array<double, 3> RGB = str2clr[shape.clr_string];
                Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                //getContext()->SetColor(AISBauteil, color, Standard_False);
                AISBauteil->SetColor(color);
            }
        } else if (shape.clr_random) {
            std::array<double, 3> RGB{(rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0};
            Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
            AISBauteil->SetColor(color);
        } else
            AISBauteil->SetColor(shape.clr);
        //***********************************************************************************


        //***********************************************************************************
        // line color
        if (str2clr.find(shape.line_clr_string) == str2clr.end())
            AISBauteil->Attributes()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
        else {
            std::array<double, 3> RGB = str2clr[shape.line_clr_string];
            Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
            AISBauteil->Attributes()->FaceBoundaryAspect()->SetColor(color);
        }
        //***********************************************************************************


        //***********************************************************************************
        AISBauteil->Attributes()->FaceBoundaryAspect()->SetWidth(8);
        AISBauteil->Attributes()->FaceBoundaryAspect()->SetTypeOfLine(Aspect_TOL_SOLID);
        AISBauteil->SetTransparency(shape.transparency);
        getContext()->Display(AISBauteil, Standard_False);
        //***********************************************************************************


        //***********************************************************************************
        // adaptor face normal
        if (shape.show_adaptor_face_normal) {
            Handle(AIS_Shape) AISFaceNormal = new AIS_Shape(construct_vectors(shape.center_of_face, gp_Vec(shape.adaptor_face_normal), 0.15));
            getContext()->SetMaterial(AISFaceNormal, Graphic3d_NOM_DEFAULT, Standard_False);
            getContext()->SetColor(AISFaceNormal, shape.clr_adaptor_face_normal, Standard_False);
            getContext()->SetTransparency(AISFaceNormal, 0, Standard_False);
            getContext()->SetDisplayMode(AISFaceNormal, 1, Standard_False);
            getContext()->Display(AISFaceNormal, Standard_False);
        }
        //***********************************************************************************


        //***********************************************************************************
        // edge orientation
        if (shape.show_edge_orientation) {
            TopExp_Explorer Ex;
            for (Ex.Init(shape.shape, TopAbs_EDGE); Ex.More(); Ex.Next()) {
                TopoDS_Edge E = TopoDS::Edge(Ex.Current());
                TopoDS_Vertex FirstVertex = TopExp::FirstVertex(TopoDS::Edge(Ex.Current()), Standard_True);
                gp_Pnt FirstPoint = BRep_Tool::Pnt(FirstVertex);
                TopoDS_Vertex SecondVertex = TopExp::LastVertex(TopoDS::Edge(Ex.Current()), Standard_True);
                gp_Pnt SecondPoint = BRep_Tool::Pnt(SecondVertex);
                Handle(AIS_Shape) AISFaceNormal = new AIS_Shape(construct_vectors(FirstPoint, gp_Vec(FirstPoint, SecondPoint), 0.1));
                getContext()->SetMaterial(AISFaceNormal, Graphic3d_NOM_DEFAULT, Standard_False);
                getContext()->SetColor(AISFaceNormal, shape.clr_edge_arrows, Standard_False);
                getContext()->SetTransparency(AISFaceNormal, 0, Standard_False);
                getContext()->SetDisplayMode(AISFaceNormal, 1, Standard_False);
                getContext()->Display(AISFaceNormal, Standard_False);
            }
        }
        //***********************************************************************************


        //***********************************************************************************
        // Vertices
        if (shape.show_vertices) {
            TopExp_Explorer Ex;
            TopoDS_ListOfShape elements;
            for (Ex.Init(shape.shape, TopAbs_VERTEX); Ex.More(); Ex.Next()) {
                gp_Pnt P = BRep_Tool::Pnt(TopoDS::Vertex(Ex.Current()));
                Handle_AIS_ColoredShape AISPoint = new AIS_ColoredShape(BRepPrimAPI_MakeSphere(P, shape.vertex_size).Shape());
                getContext()->SetMaterial(AISPoint, Graphic3d_NOM_GOLD, Standard_False);
                getContext()->SetDisplayMode(AISPoint, 1, Standard_False);
                std::array<double, 3> RGB{(rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0};
                Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                AISPoint->SetColor(color);
                getContext()->Display(AISPoint, Standard_False);
            }
        }
        //***********************************************************************************

    }

    emit selectionChanged();
    fitAll();
    QApplication::restoreOverrideCursor();
}

void DocumentTut::create_AIS_of_SB(std::list<viewerHelper::DisplayShapes_SB> &shapes, std::set<std::string> selected_checkBoxes) {

    const Handle(Prs3d_Drawer) selectionStyle = new Prs3d_Drawer();
    selectionStyle->SetColor(Quantity_NOC_RED);
    getContext()->SetSelectionStyle(selectionStyle);
    getContext()->SetHighlightStyle(selectionStyle);
    getContext()->DefaultDrawer()->SetFaceBoundaryDraw(Standard_True);

    auto chord_dev = getContext()->DefaultDrawer()->MaximalChordialDeviation() / 10.;
    getContext()->DefaultDrawer()->SetMaximalChordialDeviation(chord_dev);

    if (DocumentTut::firstStart) {

        for (auto &shape: shapes) {

            // skip null shapes to prevent crashing of viewer
            if (shape.shape.IsNull()) {
                std::cerr << "VIEWER: Shape is Null " << shape.shape.HashCode(INT_MAX) << std::endl;
                continue;
            }

            //***********************************************************************************
            // QApplication::setOverrideCursor(Qt::WaitCursor);
            Handle_AIS_ColoredShape AISBauteil = new AIS_ColoredShape(shape.shape);
            getContext()->SetMaterial(AISBauteil, Graphic3d_NOM_GOLD, Standard_False);
            getContext()->SetDisplayMode(AISBauteil, 1, Standard_False);
            //***********************************************************************************

            shape.AIS = AISBauteil;

            //***********************************************************************************
            // color
            if (shape.clr_by_string) {
                if (str2clr.find(shape.clr_string) == str2clr.end())
                    AISBauteil->SetColor(shape.clr);
                else {
                    std::array<double, 3> RGB = str2clr[shape.clr_string];
                    Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                    AISBauteil->SetColor(color);
                }
            } else if (shape.clr_random) {
                std::array<double, 3> RGB{(rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0, (rand() % 255 + 1) / 255.0};
                Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                AISBauteil->SetColor(color);
            } else {
                if (shape.categories.find("Shading") != shape.categories.end() || shape.categories.find("Shading:Site") != shape.categories.end() || shape.categories.find("Shading:Building") != shape.categories.end())
                    shape.clr_string = "PURPLE";
                else if (shape.categories.find("2a") != shape.categories.end() && shape.categories.find("INTERNAL") != shape.categories.end())
                    shape.clr_string = "LIGHTGREEN";
                else if (shape.categories.find("2a") != shape.categories.end() && shape.categories.find("EXTERNAL") != shape.categories.end())
                    shape.clr_string = "BLUE";
                else if (shape.categories.find("2a") != shape.categories.end() && shape.categories.find("EXTERNAL_EARTH") != shape.categories.end())
                    shape.clr_string = "BROWN";
                else if (shape.categories.find("2a") != shape.categories.end() && shape.categories.find("EXTERNAL_FIRE") != shape.categories.end())
                    shape.clr_string = "ORANGE";
                else if (shape.categories.find("2a") != shape.categories.end() && shape.categories.find("EXTERNAL_WATER") != shape.categories.end())
                    shape.clr_string = "CYAN";
                else if (shape.categories.find("2b") != shape.categories.end() && shape.categories.find("INTERNAL") != shape.categories.end())
                    shape.clr_string = "LIGHTRED";
                else if (shape.categories.find("2b") != shape.categories.end() && shape.categories.find("EXTERNAL") != shape.categories.end())
                    shape.clr_string = "RED";
                else if (shape.categories.find("2b") != shape.categories.end() && shape.categories.find("EXTERNAL_EARTH") != shape.categories.end())
                    shape.clr_string = "DARKRED";
                else if (shape.categories.find("1") != shape.categories.end() && shape.categories.find("INTERNAL") != shape.categories.end())
                    shape.clr_string = "LIGHTGREEN";
                else if (shape.categories.find("1") != shape.categories.end() && shape.categories.find("EXTERNAL") != shape.categories.end())
                    shape.clr_string = "BLUE";
                else if (shape.categories.find("1") != shape.categories.end() && shape.categories.find("EXTERNAL_EARTH") != shape.categories.end())
                    shape.clr_string = "BROWN";
                else if (shape.categories.find("1") != shape.categories.end() && shape.categories.find("EXTERNAL_FIRE") != shape.categories.end())
                    shape.clr_string = "ORANGE";
                else if (shape.categories.find("1") != shape.categories.end() && shape.categories.find("EXTERNAL_WATER") != shape.categories.end())
                    shape.clr_string = "CYAN";
                else if (shape.categories.find("SpacexSpaceVirtuals") != shape.categories.end())
                    shape.clr_string = "ORANGE";
                else
                    shape.clr_string = "BLACK";

                std::array<double, 3> RGB = str2clr[shape.clr_string];
                Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                AISBauteil->SetColor(color);
            }
            //***********************************************************************************


            //***********************************************************************************
            AISBauteil->Attributes()->FaceBoundaryAspect()->SetWidth(8);
            AISBauteil->Attributes()->FaceBoundaryAspect()->SetTypeOfLine(Aspect_TOL_SOLID);
            AISBauteil->SetTransparency(shape.transparency);
            getContext()->Display(AISBauteil, Standard_False);
            //***********************************************************************************


            //***********************************************************************************
            // line color
            if (str2clr.find(shape.line_clr_string) == str2clr.end())
                AISBauteil->Attributes()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
            else {
                std::array<double, 3> RGB = str2clr[shape.line_clr_string];
                Quantity_Color color(RGB[0], RGB[1], RGB[2], Quantity_TypeOfColor::Quantity_TOC_RGB);
                AISBauteil->Attributes()->FaceBoundaryAspect()->SetColor(color);
            }
            //***********************************************************************************


            //***********************************************************************************
            // adaptor face normal
            if (shape.show_adaptor_face_normal) {
                double scale = 0.03;
                Handle(AIS_Shape) AISFaceNormal = new AIS_Shape(construct_vectors(shape.center_of_face, gp_Vec(shape.adaptor_face_normal), scale));
                getContext()->SetMaterial(AISFaceNormal, Graphic3d_NOM_DEFAULT, Standard_False);
                getContext()->SetColor(AISFaceNormal, shape.clr_adaptor_face_normal, Standard_False);
                getContext()->SetTransparency(AISFaceNormal, 0, Standard_False);
                getContext()->SetDisplayMode(AISFaceNormal, 1, Standard_False);
                shape.AIS_NormalArrow = AISFaceNormal;
            }
                //***********************************************************************************



                //if (shape.categories.find("Rel_IfcExternalSpatialElement") != shape.categories.end())
                //    shape.line_clr_string = "DARKGREY";

            else if (shape.categories.find("IfcExternalSpatialElement") != shape.categories.end()) {
                shape.line_clr_string = "GREY";
                AISBauteil->Attributes()->FaceBoundaryAspect()->SetWidth(6);
            }


        }
    }

    show_and_hide(shapes, selected_checkBoxes);

    getContext()->UpdateCurrentViewer();
    emit selectionChanged();

    if (DocumentTut::firstStart) {
        DocumentTut::firstStart = false;
        fitAll();
    }

    QApplication::restoreOverrideCursor();

}

void DocumentTut::show_and_hide(std::list<viewerHelper::DisplayShapes_SB> &shapes, std::set<std::string> &selected_checkBoxes) {

    for (auto &shape: shapes) {

        if (shape.AIS.IsNull()) continue;

        bool to_display = false;

        for (const auto &c: shape.categories)
            if (selected_checkBoxes.find(c) != selected_checkBoxes.end()) {
                to_display = true;
                break;
            }

        if (to_display) {
            getContext()->Display(shape.AIS, Standard_False);

            if (selected_checkBoxes.find("Normals") != selected_checkBoxes.end()) {
                if (!shape.AIS_NormalArrow.IsNull()) getContext()->Display(shape.AIS_NormalArrow, Standard_False);
                else getContext()->Remove(shape.AIS_NormalArrow, Standard_False);
            }
        } else {
            getContext()->Remove(shape.AIS, Standard_False);
            getContext()->Remove(shape.AIS_NormalArrow, Standard_False);
        }
    }
}

TopoDS_Shape DocumentTut::construct_vectors(gp_Pnt center_pnt, gp_Vec vec, double scale) {

    // cylinder shaft
    double cylinder_radius = 0.2 * scale;
    vec.Scale(scale);
    Standard_Real cylinder_height = vec.Magnitude();
    gp_Ax2 cylinder_origin = gp_Ax2(center_pnt, gp_Dir(vec));
    TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(cylinder_origin, cylinder_radius, cylinder_height).Shape();

    // cone tip
    TopoDS_Shape cone = BRepPrimAPI_MakeCone(cylinder_origin, 1.5 * cylinder_radius, 0, 4 * cylinder_radius).Shape(); //
    gp_Trsf trns;
    trns.SetTranslation(vec);
    cone = BRepBuilderAPI_Transform(cone, trns).Shape();

    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    builder.Add(Comp, cylinder);
    builder.Add(Comp, cone);

    return Comp;

//    Handle(Graphic3d_Structure) aStructure = Handle(Graphic3d_Structure)::DownCast(getContext());
//    gp_Vec pnt_as_vec = gp_Vec(center_pnt.X(), center_pnt.Y(), center_pnt.Z());
//    gp_Vec start = pnt_as_vec + vec + vec;
//    gp_Pnt pnt_start = gp_Pnt(start.X(), start.Y(), start.Z());
//    Prs3d_Arrow arr;
//    arr.Draw(aStructure, pnt_start, gp_Dir(vec), 0.174533, vec.Magnitude());

}
