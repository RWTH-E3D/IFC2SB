#include "SB_StlAPI_Writer.h"
#include "SB_RWStl.h"

#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>

SB_StlAPI_Writer::SB_StlAPI_Writer() : myASCIIMode(Standard_True) {}

Standard_Boolean SB_StlAPI_Writer::Write(const TopoDS_Shape &theShape, const Standard_CString theFileName, const Standard_CString solid_name, bool overwrite, const Message_ProgressRange &theProgress) const {

    Standard_Integer aNbNodes = 0;
    Standard_Integer aNbTriangles = 0;

    // calculate total number of the nodes and triangles
    for (TopExp_Explorer anExpSF(theShape, TopAbs_FACE); anExpSF.More(); anExpSF.Next()) {
        TopLoc_Location aLoc;
        Handle(Poly_Triangulation) aTriangulation = BRep_Tool::Triangulation(TopoDS::Face(anExpSF.Current()), aLoc);
        if (!aTriangulation.IsNull()) {
            aNbNodes += aTriangulation->NbNodes();
            aNbTriangles += aTriangulation->NbTriangles();
        }
    }

    // no triangulation on the shape
    if (aNbTriangles == 0)
        return Standard_False;

    // create temporary triangulation
    Handle(Poly_Triangulation) aMesh = new Poly_Triangulation(aNbNodes, aNbTriangles, Standard_False);
    Standard_Integer aNbFacesNoTri = 0; // count faces missing triangulation
    Standard_Integer aNodeOffset = 0; // fill temporary triangulation
    Standard_Integer aTriangleOffet = 0;

    for (TopExp_Explorer anExpSF(theShape, TopAbs_FACE); anExpSF.More(); anExpSF.Next()) {
        const TopoDS_Shape &aFace = anExpSF.Current();
        TopLoc_Location aLoc;
        Handle(Poly_Triangulation) aTriangulation = BRep_Tool::Triangulation(TopoDS::Face(aFace), aLoc);
        if (aTriangulation.IsNull()) {
            ++aNbFacesNoTri;
            continue;
        }

        const Poly_Array1OfTriangle &aTriangles = aTriangulation->Triangles();

        // copy nodes
        gp_Trsf aTrsf = aLoc.Transformation();

        for (Standard_Integer aNodeIter = 1; aNodeIter <= aTriangulation->NbNodes(); ++aNodeIter) {
            gp_Pnt aPnt = aTriangulation->Node(aNodeIter);
            aPnt.Transform(aTrsf);
            aMesh->ChangeNode(aNodeIter + aNodeOffset) = aPnt;
        }

        // copy triangles
        const TopAbs_Orientation anOrientation = anExpSF.Current().Orientation();

        for (Standard_Integer aTriIter = 1; aTriIter <= aTriangulation->NbTriangles(); ++aTriIter) {
            Poly_Triangle aTri = aTriangulation->Triangle(aTriIter);

            Standard_Integer anId[3];
            aTri.Get(anId[0], anId[1], anId[2]);
            if (anOrientation == TopAbs_REVERSED) {
                // Swap 1, 2.
                Standard_Integer aTmpIdx = anId[1];
                anId[1] = anId[2];
                anId[2] = aTmpIdx;
            }

            // Update nodes according to the offset.
            anId[0] += aNodeOffset;
            anId[1] += aNodeOffset;
            anId[2] += aNodeOffset;

            aTri.Set(anId[0], anId[1], anId[2]);
            aMesh->ChangeTriangle(aTriIter + aTriangleOffet) = aTri;
        }

        aNodeOffset += aTriangulation->NbNodes();
        aTriangleOffet += aTriangulation->NbTriangles();
    }

    OSD_Path aPath(theFileName);
    Standard_Boolean isDone = (myASCIIMode ? SB_RWStl::WriteAscii(aMesh, aPath, solid_name, overwrite, theProgress) : SB_RWStl::WriteBinary(aMesh, aPath, theProgress));

    if (isDone && (aNbFacesNoTri > 0)) {
        // Print warning with number of faces missing triangulation
        TCollection_AsciiString aWarningMsg = TCollection_AsciiString("Warning: ") + TCollection_AsciiString(aNbFacesNoTri) + TCollection_AsciiString((aNbFacesNoTri == 1) ? " face has" : " faces have") + TCollection_AsciiString(" been skipped due to null triangulation");
        Message::DefaultMessenger()->Send(aWarningMsg, Message_Warning);
    }

    return isDone;
}