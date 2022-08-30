#ifndef SB_STLAPI_WRITER_H
#define SB_STLAPI_WRITER_H

#include <Message_ProgressScope.hxx>
#include <Message_ProgressIndicator.hxx>

class TopoDS_Shape;

class SB_StlAPI_Writer {
public:

    DEFINE_STANDARD_ALLOC

    Standard_EXPORT SB_StlAPI_Writer();

    Standard_Boolean &ASCIIMode() { return myASCIIMode; }

    Standard_EXPORT Standard_Boolean Write(const TopoDS_Shape &theShape, Standard_CString theFileName, Standard_CString solid_name = "", bool overwrite = false, const Message_ProgressRange &theProgress = Message_ProgressRange()) const;

private:
    Standard_Boolean myASCIIMode;
};

#endif //SB_STLAPI_WRITER_H
