#ifndef SB_RWSTL_H
#define SB_RWSTL_H

#include <Message_ProgressIndicator.hxx>
#include <OSD_Path.hxx>
#include <Poly_Triangulation.hxx>
#include <Standard_Macro.hxx>

class SB_RWStl {
public:

    Standard_EXPORT static Standard_Boolean WriteBinary(const Handle(Poly_Triangulation) &theMesh, const OSD_Path &thePath, const Message_ProgressRange &theProgress = Message_ProgressRange());

    Standard_EXPORT static Standard_Boolean WriteAscii(const Handle(Poly_Triangulation) &theMesh, const OSD_Path &thePath, Standard_CString solid_name, bool overwrite, const Message_ProgressRange& theProgress = Message_ProgressRange());

private:

    static Standard_Boolean writeASCII(const Handle(Poly_Triangulation) &theMesh, FILE *theFile, const Message_ProgressRange& theProgress, Standard_CString solid_name);

    static Standard_Boolean writeBinary(const Handle(Poly_Triangulation) &theMesh, FILE *theFile, const Message_ProgressRange& theProgress);

    static double rounded(double d, unsigned int n);
};

#endif //SB_RWSTL_H
