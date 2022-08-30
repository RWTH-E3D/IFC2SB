// Copyright 2022 Eric Fichter
#ifndef HEADERS_H
#define HEADERS_H

// Standard
#include <memory>
#include <regex>
#include <stack>
#include <thread>
#include <unordered_map>

// CGAL
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel CGAL_Kernel;
typedef CGAL_Kernel::Point_3 CGAL_Point;
typedef CGAL_Kernel::Segment_3 CGAL_Segment;
typedef CGAL_Kernel::Triangle_3 CGAL_Triangle;

// External scripts
#include "clipper.hpp" // Clipper http://www.angusj.com/delphi/clipper.php
#include "extern/rtree.h" // RTree https://github.com/nushoin/RTree
#include <tbb/parallel_for_each.h> // TBB
#include "extern/concave_polygon.h" // Convex Decomposition

// Boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// OpenCascade
#include <BOPAlgo_CellsBuilder.hxx>
#include <BOPAlgo_MakerVolume.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Edge.hxx>
#include <BRepCheck_Face.hxx>
#include <BRepCheck_Shell.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepGProp.hxx>
#include <BRepLProp_SLProps.hxx>
#include <BRepOffsetAPI_MakeOffset.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_OBB.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GProp_GProps.hxx>
#include <GeomLProp_SLProps.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <ShapeAnalysis.hxx>
#include <ShapeAnalysis_CheckSmallFace.hxx>
#include <ShapeAnalysis_Edge.hxx>
#include <ShapeAnalysis_ShapeContents.hxx>
#include <ShapeAnalysis_ShapeTolerance.hxx>
#include <ShapeAnalysis_Shell.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ShapeAnalysis_Wire.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Version.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>

// IfcOpenShell
#include <IfcGeomTree.h>
#include "ifcgeom_schema_agnostic/Serialization.h"
#include <ifcgeom_schema_agnostic/IfcGeomIterator.h>
#include <ifcparse/utils.h>
#include <Ifc2x3.h>
#include <Ifc4.h>
#include <Ifc4x1.h>
#include <Ifc4x2.h>
#include <Ifc4x3_rc1.h>

// IFC2SB
#include "definitions.h"
#include "octree/octree_interface.h"
#include "intersector/IntersectorInterface.h"
#include "stl_export/SB_StlAPI_Writer.h"
#include "Intersector.h"
#include "ShapeChecker.h"
#include "ifc_creator.h"
#include "Clipper.h"
#include "Space.h"
#include "Graph.h"
#include "Clip.h"
#include "Kernel.h"
#include "Product.h"
#include "oFace.h"
#include "cFace.h"
#include "Topo.h"
#include "sFace.h"
#include "ShapeHealing.h"
#include "IfcCheck.h"
#include "Layer.h"
#include "Topo2D.h"
#include "graph_2D.h"

#ifdef VISUALIZATION
#include "Viewer.h"
#include "viewer/ViewerMain.h"
#endif

#endif //HEADERS_H


