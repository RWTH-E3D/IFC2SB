#include "SB_RWStl.h"

#include <Message_ProgressScope.hxx>
#include <OSD_File.hxx>
#include <OSD_OpenFile.hxx>

namespace {

    const Standard_Integer THE_STL_SIZEOF_FACET = 50;
    const Standard_Integer IND_THRESHOLD = 1000; // increment the indicator every 1k triangles

    //! Writing a Little Endian 32 bits integer
    inline void convertInteger(const Standard_Integer theValue,
                               Standard_Character *theResult) {
        union {
            Standard_Integer i;
            Standard_Character c[4];
        } anUnion{};
        anUnion.i = theValue;

        theResult[0] = anUnion.c[0];
        theResult[1] = anUnion.c[1];
        theResult[2] = anUnion.c[2];
        theResult[3] = anUnion.c[3];
    }

    //! Writing a Little Endian 32 bits float
    inline void convertDouble(const Standard_Real theValue,
                              Standard_Character *theResult) {
        union {
            Standard_ShortReal i;
            Standard_Character c[4];
        } anUnion{};
        anUnion.i = (Standard_ShortReal) theValue;

        theResult[0] = anUnion.c[0];
        theResult[1] = anUnion.c[1];
        theResult[2] = anUnion.c[2];
        theResult[3] = anUnion.c[3];
    }

}

Standard_Boolean SB_RWStl::WriteBinary(const Handle(Poly_Triangulation) &theMesh, const OSD_Path &thePath, const Message_ProgressRange &theProgress) {
    if (theMesh.IsNull() || theMesh->NbTriangles() <= 0)
        return Standard_False;

    TCollection_AsciiString aPath;
    thePath.SystemName(aPath);

    FILE *aFile = OSD_OpenFile(aPath, "wb");
    if (aFile == nullptr)
        return Standard_False;

    Standard_Boolean isOK = writeBinary(theMesh, aFile, theProgress);

    fclose(aFile);
    return isOK;
}

Standard_Boolean SB_RWStl::WriteAscii(const Handle(Poly_Triangulation) &theMesh, const OSD_Path &thePath, const Standard_CString solid_name, bool overwrite, const Message_ProgressRange &theProgress) {
    if (theMesh.IsNull() || theMesh->NbTriangles() <= 0)
        return Standard_False;

    TCollection_AsciiString aPath;
    thePath.SystemName(aPath);

    const char *mode;
    if (overwrite)
        mode = "w";
    else
        mode = "a";

    FILE *aFile = OSD_OpenFile(aPath, mode);
    if (aFile == nullptr)
        return Standard_False;

    Standard_Boolean isOK = writeASCII(theMesh, aFile, theProgress, solid_name);
    fclose(aFile);
    return isOK;
}

Standard_Boolean SB_RWStl::writeASCII(const Handle(Poly_Triangulation) &theMesh, FILE *theFile, const Message_ProgressRange &theProgress, const Standard_CString solid_name) {
    // note that space after 'solid' is necessary for many systems
    std::string temp = "solid ";
    temp += solid_name;
    temp += +"\n";

    Standard_CString solid_line = temp.c_str();

    fwrite(solid_line, 1, temp.size(), theFile);

    char aBuffer[512];
    memset(aBuffer, 0, sizeof(aBuffer));

    const Standard_Integer NBTriangles = theMesh->NbTriangles();
    Message_ProgressScope aPS(theProgress, "Triangles", NBTriangles);

    const Poly_Array1OfTriangle &aTriangles = theMesh->Triangles();
    Standard_Integer anElem[3] = {0, 0, 0};
    for (Standard_Integer aTriIter = 1; aTriIter <= NBTriangles; ++aTriIter) {
        const Poly_Triangle &aTriangle = aTriangles(aTriIter);
        aTriangle.Get(anElem[0], anElem[1], anElem[2]);

        const gp_Pnt aP1 = theMesh->Node(anElem[0]);
        const gp_Pnt aP2 = theMesh->Node(anElem[1]);
        const gp_Pnt aP3 = theMesh->Node(anElem[2]);

        const gp_Vec aVec1(aP1, aP2);
        const gp_Vec aVec2(aP1, aP3);
        gp_Vec aVNorm = aVec1.Crossed(aVec2);
        if (aVNorm.SquareMagnitude() > gp::Resolution())
            aVNorm.Normalize();
        else
            aVNorm.SetCoord(0.0, 0.0, 0.0);

        //  unsigned int plc = 4;

        Sprintf(aBuffer,
                " facet normal % 12e % 12e % 12e\n"
                "   outer loop\n"
                "     vertex % 12e % 12e % 12e\n"
                "     vertex % 12e % 12e % 12e\n"
                "     vertex % 12e % 12e % 12e\n"
                "   endloop\n"
                " endfacet\n",
                aVNorm.X(), aVNorm.Y(), aVNorm.Z(),
//                rounded(aP1.X(), plc), rounded(aP1.Y(), plc), rounded(aP1.Z(), plc),
//                rounded(aP2.X(), plc), rounded(aP2.Y(), plc), rounded(aP2.Z(), plc),
//                rounded(aP3.X(), plc), rounded(aP3.Y(), plc), rounded(aP3.Z(), plc));
                aP1.X(), aP1.Y(), aP1.Z(),
                aP2.X(), aP2.Y(), aP2.Z(),
                aP3.X(), aP3.Y(), aP3.Z());

        if (fprintf(theFile, "%s", aBuffer) < 0)
            return Standard_False;

        // update progress only per 1k triangles
        if ((aTriIter % IND_THRESHOLD) == 0) {
            if (!aPS.More()) return Standard_False;
            aPS.Next(IND_THRESHOLD);
        }
    }

    if (fwrite("endsolid\n", 1, 9, theFile) != 9)
        return Standard_False;

    return Standard_True;
}

double SB_RWStl::rounded(double d, unsigned int n) {
    unsigned int f = pow(10, n);
    return std::round(d * f) / f;
}

Standard_Boolean SB_RWStl::writeBinary(const Handle(Poly_Triangulation) &theMesh, FILE *theFile, const Message_ProgressRange &theProgress) {
    char aHeader[80] = "STL Exported by OpenCASCADE [www.opencascade.com]";
    if (fwrite(aHeader, 1, 80, theFile) != 80) return Standard_False;

    const Standard_Integer aNBTriangles = theMesh->NbTriangles();
    Message_ProgressScope aPS(theProgress, "Triangles", aNBTriangles);

    const Standard_Size aNbChunkTriangles = 4096;
    const Standard_Size aChunkSize = aNbChunkTriangles * THE_STL_SIZEOF_FACET;
    NCollection_Array1<Standard_Character> aData(1, aChunkSize);
    Standard_Character *aDataChunk = &aData.ChangeFirst();

    const Poly_Array1OfTriangle &aTriangles = theMesh->Triangles();

    Standard_Character aConv[4];
    convertInteger(aNBTriangles, aConv);
    if (fwrite(aConv, 1, 4, theFile) != 4) return Standard_False;

    Standard_Size aByteCount = 0;
    for (Standard_Integer aTriIter = 1; aTriIter <= aNBTriangles; ++aTriIter) {
        Standard_Integer id[3];
        const Poly_Triangle &aTriangle = aTriangles(aTriIter);
        aTriangle.Get(id[0], id[1], id[2]);

        const gp_Pnt aP1 = theMesh->Node(id[0]);
        const gp_Pnt aP2 = theMesh->Node(id[1]);
        const gp_Pnt aP3 = theMesh->Node(id[2]);

        gp_Vec aVec1(aP1, aP2);
        gp_Vec aVec2(aP1, aP3);
        gp_Vec aVNorm = aVec1.Crossed(aVec2);
        if (aVNorm.SquareMagnitude() > gp::Resolution())
            aVNorm.Normalize();
        else
            aVNorm.SetCoord(0.0, 0.0, 0.0);

        convertDouble(aVNorm.X(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aVNorm.Y(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aVNorm.Z(), &aDataChunk[aByteCount]);
        aByteCount += 4;

        convertDouble(aP1.X(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP1.Y(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP1.Z(), &aDataChunk[aByteCount]);
        aByteCount += 4;

        convertDouble(aP2.X(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP2.Y(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP2.Z(), &aDataChunk[aByteCount]);
        aByteCount += 4;

        convertDouble(aP3.X(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP3.Y(), &aDataChunk[aByteCount]);
        aByteCount += 4;
        convertDouble(aP3.Z(), &aDataChunk[aByteCount]);
        aByteCount += 4;

        aDataChunk[aByteCount] = 0;
        aByteCount += 1;
        aDataChunk[aByteCount] = 0;
        aByteCount += 1;

        // Chunk is filled. Dump it to the file.
        if (aByteCount == aChunkSize) {
            if (fwrite(aDataChunk, 1, aChunkSize, theFile) != aChunkSize)
                return Standard_False;


            aByteCount = 0;
        }

        // update progress only per 1k triangles
        if ((aTriIter % IND_THRESHOLD) == 0) {
            if (!aPS.More()) return Standard_False;
            aPS.Next(IND_THRESHOLD);
        }
    }

    // Write last part if necessary.
    if (fwrite(aDataChunk, 1, aByteCount, theFile) != aByteCount)
        return Standard_False;

    return Standard_True;
}
