#ifndef VIEWERHELPER_H
#define VIEWERHELPER_H

#include <AIS_ColoredShape.hxx>
#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <Aspect_Window.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeSegment.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Graphic3d_NameOfMaterial.hxx>
#include <OSD_Environment.hxx>
#include <OSD_Environment.hxx>
#include <QAction>
#include <QApplication>
#include <QApplication>
#include <QDialog>
#include <QFileDialog>
#include <QLayout>
#include <QList>
#include <QMainWindow>
#include <QMdiSubWindow>
#include <QMessageBox>
#include <QObject>
#include <QPushButton>
#include <QSignalMapper>
#include <QTranslator>
#include <QWidget>
#include <Standard_TypeDef.hxx>
#include <Standard_Version.hxx>
#include <Standard_WarningsDisable.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <V3d_View.hxx>
#include <cstdlib>
#include <gp.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax2d.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <set>

namespace viewerHelper {

    struct DisplayShapes {
        TopoDS_Shape shape;
        Quantity_Color clr;
        bool clr_by_string;
        bool clr_random;
        std::string clr_string;
        double transparency;
        std::string line_clr_string;

        // normal
        gp_Dir adaptor_face_normal;
        gp_Pnt center_of_face;

        bool show_adaptor_face_normal;
        Quantity_Color clr_adaptor_face_normal;
        Quantity_Color clr_edge_arrows;

        bool show_edge_orientation;
        bool show_vertices;

        double edge_thickness; // when shape is an edge
        double vertex_size; // when show_vertices

        DisplayShapes() {
            clr = Quantity_NOC_GRAY;
            transparency = 0;
            show_adaptor_face_normal = false;
            show_edge_orientation = false;
            show_vertices = false;
            clr_by_string = false;
            clr_random = false;
            clr_adaptor_face_normal = Quantity_NOC_BLACK;
            clr_edge_arrows = Quantity_NOC_BLACK;
            line_clr_string = "BLACK";
            edge_thickness = 0.01;
            vertex_size = 0.01;
        }
    };

    /* Derived struct for SB viewer */
    struct DisplayShapes_SB : DisplayShapes {
        std::string sb_type;
        Handle_AIS_ColoredShape AIS;
        Handle_AIS_Shape AIS_NormalArrow;
        std::set<std::string> categories;

        DisplayShapes_SB() {
            clr = Quantity_NOC_GRAY;
            transparency = 0;
            show_adaptor_face_normal = false;
            show_edge_orientation = false;
            show_vertices = false;
            clr_by_string = false;
            clr_random = false;
            clr_adaptor_face_normal = Quantity_NOC_BLACK;
            clr_edge_arrows = Quantity_NOC_BLACK;
            AIS = nullptr;
            AIS_NormalArrow = nullptr;
        }
    };
}

#endif //VIEWERHELPER_H