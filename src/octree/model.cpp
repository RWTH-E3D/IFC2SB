/*

Load geometric model data as vertices and faces.

- Ralf-Peter Mundani, mundani@mytum.de

-----------------------------------------------------------------------

Copyright (c) 2010 Technische Universit�t M�nchen.  All rights reserved.   
  
Permission to use, copy, modify and distribute this software and its   
documentation for any purpose is hereby granted without fee, provided   
that the above copyright notice and this permission notice appear in   
all copies of this software and that you do not sell the software.   
  
THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,   
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY   
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

*/
#include "octree_includes.h"

//=============================================================================
// @ARG: none
//=============================================================================
model::model() {
    nvertex = nface = nedge = 0;
    vlist = nullptr;
    elist = nullptr;
    flist = nullptr;
    minx = maxx = miny = maxy = minz = maxz = 0;
    bbox_diffx = bbox_diffy = bbox_diffz = bbox_norm = 0;
}

//=============================================================================
// @ARG: none
//=============================================================================
model::~model() {
    if (vlist != nullptr) {
        delete[] vlist;
        vlist = nullptr;
    }
    if (elist != nullptr) {
        delete[] elist;
        elist = nullptr;
    }
    if (flist != nullptr) {
        delete[] flist;
        flist = nullptr;
    }
}

//=============================================================================
// @ARG: name of TRI file to load
//=============================================================================
void model::load_model(const char *filename, int oct_depth, bool extend_bbox_for_floodfill) {
    OFCGeometry geom;
    OFCInterfaceEAF interf_eaf;

    std::string s = filename;
    s = s.substr(s.length() - 3);
    if (s != "stl" && s != "STL") {
        std::cerr << "Filetype " << s << " is unknown. Abort." << std::endl;
        abort();
    }

    // read in geometry
    interf_eaf.inputGeometry(geom, filename);

    // set reference to maps (i.e. shortcuts)
    const std::map<int, OFCNode> &nds = geom.getConstNdsMap();
    const std::map<int, OFCTri> &trs = geom.getConstTrsMap();

    //******************************
    /*
    for (map<int, OFCNode>::const_iterator it = nds.begin(); it != nds.end(); ++it) {
        cout << it->first << " " << (it->second).id;
        for (int i = 0; i < 3; ++i)
            cout << " " << (it->second).coord[i];
        cout << endl;
    }
    for (map<int, OFCTri>::const_iterator it = trs.begin(); it != trs.end(); ++it) {
        cout << it->first << " " << (it->second).id;
        for (int i = 0; i < 3; ++i)
            cout << " " << (it->second).nd[i] << " " << ((it->second).nd[i])->id;
        cout << endl;

        list<OFCAttrib *>::const_iterator lit;
        for (lit = (it->second).attrib.begin(); lit != (it->second).attrib.end(); ++lit) {
            OFCAttrib *attr = *lit;
            cout << attr->type << " " << attr->key.length() << endl;
            string val = *((string *) attr->getAttrib());
            cout << attr->type << "\t" << "\t" << attr->key << "\t" << val << endl;
        }
    }*/
    //******************************

    /* allocate memory for vertex and face list */
    nvertex = nds.size();
    nface = trs.size();
    vlist = new Vertex[nvertex];
    flist = new Face[nface];

    std::cout << "\n... Number of vertices: " << nvertex << ", number of faces: " << nface << std::endl;

    // restore node information
    std::map<int, OFCNode>::const_iterator nit;
    for (nit = nds.begin(); nit != nds.end(); ++nit) {
        const OFCNode &nd = nit->second;

        // point coordinates
        vlist[nd.id].x = nd.coord[0];
        vlist[nd.id].y = nd.coord[1];
        vlist[nd.id].z = nd.coord[2];

        // attributes are ignored at this moment
    }

    // check if any tri attributes are present
    std::map<int, OFCTri>::const_iterator tit;
    for (tit = trs.begin(); tit != trs.end(); ++tit) {
        const OFCTri &tr = tit->second;

        // node ids for tris
        flist[tr.id].v1 = tr.nd[0]->id;
        flist[tr.id].v2 = tr.nd[1]->id;
        flist[tr.id].v3 = tr.nd[2]->id;

        flist[tr.id].nx = flist[tr.id].ny = flist[tr.id].nz = 0;

        flist[tr.id].minx = flist[tr.id].miny = flist[tr.id].minz = flist[tr.id].maxx = flist[tr.id].maxy = flist[tr.id].maxz = 0;
        flist[tr.id].A = 0;

        // add attribute "guid"
        std::list<OFCAttrib *>::const_iterator lit;
        for (lit = tr.attrib.begin(); lit != tr.attrib.end(); ++lit) {
            OFCAttrib *attr = *lit;
            std::string val = *((std::string *) attr->getAttrib());
            flist[tr.id].guid = val;
        }

    }

    /* compute face normals */
    __calc_normals();

    /* compute model's bounding box and normalise model to unit cube */
    __comp_BBOX();

    double Lx = maxx - minx;
    double Ly = maxy - miny;
    double Lz = maxz - minz;
    std::cout << "... Original geometry data bounding box\t\t: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")" << std::endl;

    if (extend_bbox_for_floodfill) {
        // For flood fill algorithm bounding box has to be extended, so that geometry does not touch border
        // condition for flood fill is that gap between boundary and geometry is bigger than minimal octnode size
        // because after extending the bounding box, the minimal size increases also a little bit
        // the following condition must therefore be fulfilled: gap >! smin + 2*smin/(2^d). Implemented is with d-1
        double sminx = Lx / pow(2, oct_depth) + 2.01 * Lx / pow(2, 2 * (oct_depth - 1));
        double sminy = Ly / pow(2, oct_depth) + 2.01 * Ly / pow(2, 2 * (oct_depth - 1));
        double sminz = Lz / pow(2, oct_depth) + 2.01 * Lz / pow(2, 2 * (oct_depth - 1));
        minx -= sminx;
        maxx += sminx;
        miny -= sminy;
        maxy += sminy;
        minz -= sminz;
        maxz += sminz;
        Lx = maxx - minx;
        Ly = maxy - miny;
        Lz = maxz - minz;
        std::cout << "... Extended original geometry data bounding box: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")" << std::endl;
    } else { // minor margin, otherwise triangles at border are ignored
        double d = 0.01;
        minx -= d;
        maxx += d;
        miny -= d;
        maxy += d;
        minz -= d;
        maxz += d;
    }

    // Minimal voxel size
    double dmin = std::max(std::max(Lx, Ly), Lz) / pow(2, oct_depth);
    std::cout << "... Minimal voxel size\t: " << dmin << std::endl;

    __normalise_model();
    __comp_BBOX();
    std::cout << "... Normalized geometry data bounding box\t: " << minx << " " << maxx << " | " << miny << " " << maxy << " | " << minz << " " << maxz << std::endl;

    /* compute edges for wireframe display */
    //__comp_edges();

    /* compute face bounding boxes */
    __comp_bbox();

}

//=============================================================================
// @ARG: name of TRI file to load
//=============================================================================
void model::load_model(const std::vector<std::tuple<double, double, double>> &vertices, const std::vector<std::tuple<int, int, int>> &faces, const std::vector<std::string> &attrs, int oct_depth, bool extend_bbox_for_floodfill) {
    OFCGeometry geom;
    OFCInterfaceEAF interf_eaf;

    // read in geometry
    interf_eaf.inputGeometry(geom, vertices, faces, attrs);

    // set reference to maps (i.e. shortcuts)
    const std::map<int, OFCNode> &nds = geom.getConstNdsMap();
    const std::map<int, OFCTri> &trs = geom.getConstTrsMap();

    //******************************
    /*
    for (map<int, OFCNode>::const_iterator it = nds.begin(); it != nds.end(); ++it) {
        cout << it->first << " " << (it->second).id;
        for (int i = 0; i < 3; ++i)
            cout << " " << (it->second).coord[i];
        cout << endl;
    }
    for (map<int, OFCTri>::const_iterator it = trs.begin(); it != trs.end(); ++it) {
        cout << it->first << " " << (it->second).id;
        for (int i = 0; i < 3; ++i)
            cout << " " << (it->second).nd[i] << " " << ((it->second).nd[i])->id;
        cout << endl;

        list<OFCAttrib *>::const_iterator lit;
        for (lit = (it->second).attrib.begin(); lit != (it->second).attrib.end(); ++lit) {
            OFCAttrib *attr = *lit;
            cout << attr->type << " " << attr->key.length() << endl;
            string val = *((string *) attr->getAttrib());
            cout << attr->type << "\t" << "\t" << attr->key << "\t" << val << endl;
        }
    }*/
    //******************************

    /* allocate memory for vertex and face list */
    nvertex = nds.size();
    nface = trs.size();
    vlist = new Vertex[nvertex];
    flist = new Face[nface];

#ifdef PRINT
    std::cout << "\t\t... Number of vertices: " << nvertex << ", number of faces: " << nface << "\n";
#endif
    // restore node information
    std::map<int, OFCNode>::const_iterator nit;
    for (nit = nds.begin(); nit != nds.end(); ++nit) {
        const OFCNode &nd = nit->second;

        // point coordinates
        vlist[nd.id].x = nd.coord[0];
        vlist[nd.id].y = nd.coord[1];
        vlist[nd.id].z = nd.coord[2];

        // attributes are ignored at this moment
    }

    // check if any tri attributes are present
    std::map<int, OFCTri>::const_iterator tit;
    for (tit = trs.begin(); tit != trs.end(); ++tit) {
        const OFCTri &tr = tit->second;

        // node ids for tris
        flist[tr.id].v1 = tr.nd[0]->id;
        flist[tr.id].v2 = tr.nd[1]->id;
        flist[tr.id].v3 = tr.nd[2]->id;

        flist[tr.id].nx = flist[tr.id].ny = flist[tr.id].nz = 0;

        flist[tr.id].minx = flist[tr.id].miny = flist[tr.id].minz = flist[tr.id].maxx = flist[tr.id].maxy = flist[tr.id].maxz = 0;
        flist[tr.id].A = 0;

        // add attribute "guid"
        std::list<OFCAttrib *>::const_iterator lit;
        for (lit = tr.attrib.begin(); lit != tr.attrib.end(); ++lit) {
            OFCAttrib *attr = *lit;
            std::string val = *((std::string *) attr->getAttrib());
            flist[tr.id].guid = val;
        }

    }

    /* compute face normals */
    __calc_normals();

    /* compute model's bounding box and normalise model to unit cube */
    __comp_BBOX();

    double Lx = maxx - minx;
    double Ly = maxy - miny;
    double Lz = maxz - minz;
#ifdef PRINT
    std::cout << "\t\t... Original geometry data bounding box\t\t: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")\n";
#endif
    if (extend_bbox_for_floodfill) {
        // For flood fill algorithm bounding box has to be extended, so that geometry does not touch border
        // condition for flood fill is that gap between boundary and geometry is bigger than minimal octnode size
        // because after extending the bounding box, the minimal size increases also a little bit
        // the following condition must therefore be fulfilled: gap >! smin + 2*smin/(2^d). Implemented is with d-1
        double sminx = Lx / pow(2, oct_depth) + 2.01 * Lx / pow(2, 2 * (oct_depth - 1));
        double sminy = Ly / pow(2, oct_depth) + 2.01 * Ly / pow(2, 2 * (oct_depth - 1));
        double sminz = Lz / pow(2, oct_depth) + 2.01 * Lz / pow(2, 2 * (oct_depth - 1));
        minx -= sminx;
        maxx += sminx;
        miny -= sminy;
        maxy += sminy;
        minz -= sminz;
        maxz += sminz;
    } else { // minor margin, otherwise triangles at border are ignored
        double d = 0.01;
        minx -= d;
        maxx += d;
        miny -= d;
        maxy += d;
        minz -= d;
        maxz += d;
    }

#ifdef PRINT
    Lx = maxx - minx;
    Ly = maxy - miny;
    Lz = maxz - minz;
    std::cout << "\t\t... Extended original geometry data bounding box: " << Lx << " (" << minx << " -> " << maxx << ") | " << Ly << " (" << miny << " -> " << maxy << ") | " << Lz << " (" << minz << " -> " << maxz << ")\n";

    // Minimal voxel size
    double dmin = std::max({Lx, Ly, Lz}) / pow(2, oct_depth);
    std::cout << "\t\t... Minimal voxel size\t: " << dmin << "\n";
#endif

    __normalise_model();
    __comp_BBOX();
#ifdef PRINT
    std::cout << "\t\t... Normalized geometry data bounding box\t: " << minx << " " << maxx << " | " << miny << " " << maxy << " | " << minz << " " << maxz << "\n";
#endif
    /* compute face bounding boxes */
    __comp_bbox();

}

//=============================================================================
// @ARG: none
//=============================================================================
void model::__calc_normals() {
    int i;                // loop counter
    double nx, ny, nz;    // face normal's components
    double norm;        // normalisation factor

    for (i = 0; i < nface; ++i) {
        __vec_cross(vlist[flist[i].v2].x - vlist[flist[i].v1].x,
                    vlist[flist[i].v2].y - vlist[flist[i].v1].y,
                    vlist[flist[i].v2].z - vlist[flist[i].v1].z,
                    vlist[flist[i].v3].x - vlist[flist[i].v2].x,
                    vlist[flist[i].v3].y - vlist[flist[i].v2].y,
                    vlist[flist[i].v3].z - vlist[flist[i].v2].z,
                    &nx, &ny, &nz);
        norm = sqrt(nx * nx + ny * ny + nz * nz);
        flist[i].nx = nx / norm;
        flist[i].ny = ny / norm;
        flist[i].nz = nz / norm;
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
void model::__normalise_model() {
    double diffx = (minx + maxx) / 2;    // translation difference in x-direction
    double diffy = (miny + maxy) / 2;    // translation difference in y-direction
    double diffz = (minz + maxz) / 2;    // translation difference in z-direction
    double norm = (maxx - minx) / 2;    // normalise factor
    int i;                            // loop counter

    if ((maxy - miny) / 2 > norm) { norm = (maxy - miny) / 2; }
    if ((maxz - minz) / 2 > norm) { norm = (maxz - minz) / 2; }
    for (i = 0; i < nvertex; ++i) {
        vlist[i].x = (vlist[i].x - diffx) / norm;
        vlist[i].y = (vlist[i].y - diffy) / norm;
        vlist[i].z = (vlist[i].z - diffz) / norm;
    }

    bbox_diffx = diffx;
    bbox_diffy = diffy;
    bbox_diffz = diffz;
    bbox_norm = norm;
}

//=============================================================================
// @ARG: vector vec1's components
// @ARG: vector vec2's components
// @ARG: face normal's components
//=============================================================================
void model::__vec_cross(const double vec1x, const double vec1y, const double vec1z,
                        const double vec2x, const double vec2y, const double vec2z,
                        double *nx, double *ny, double *nz) {
    *nx = vec1y * vec2z - vec1z * vec2y;
    *ny = vec1z * vec2x - vec1x * vec2z;
    *nz = vec1x * vec2y - vec1y * vec2x;
}

//=============================================================================
// @ARG: none
//=============================================================================
void model::__comp_edges() {
    int i;            // loop counter
    int idx(0);        // index of edge list for inserting next element

    /* allocate memory for edge list and initialise it*/
    nedge = nface * 3 / 2;
    elist = new Edge[nedge];
    for (i = 0; i < nedge; ++i) {
        elist[i].v1 = MODEL_ERROR;
        elist[i].v2 = MODEL_ERROR;
    }

    /* insert edges in edge list */
    for (i = 0; i < nface; ++i) {
        if (__check_edge(flist[i].v1, flist[i].v2) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v1;
            elist[idx].v2 = flist[i].v2;
            ++idx;
        }
        if (__check_edge(flist[i].v2, flist[i].v3) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v2;
            elist[idx].v2 = flist[i].v3;
            ++idx;
        }
        if (__check_edge(flist[i].v3, flist[i].v1) != MODEL_ERROR) {
            elist[idx].v1 = flist[i].v3;
            elist[idx].v2 = flist[i].v1;
            ++idx;
        }
    }
}

//=============================================================================
// @ARG: vertex indices delimiting edge
//=============================================================================
int model::__check_edge(const int v1, const int v2) {
    int idx(0);        // index of edge list for checking given vertices

    while (idx < nedge) {
        if (elist[idx].v1 == MODEL_ERROR) { return idx; }
        else if (elist[idx].v1 == v1 || elist[idx].v2 == v1) {
            if (elist[idx].v1 == v2 || elist[idx].v2 == v2) { return MODEL_ERROR; }
        }
        ++idx;
    }
    return MODEL_ERROR;
}

//=============================================================================
// @ARG: none
//=============================================================================
void model::__comp_BBOX() {
    int i;        // loop counter

    minx = vlist[0].x;
    miny = vlist[0].y;
    minz = vlist[0].z;
    maxx = vlist[0].x;
    maxy = vlist[0].y;
    maxz = vlist[0].z;
    for (i = 1; i < nvertex; ++i) {
        if (vlist[i].x < minx) { minx = vlist[i].x; } else if (vlist[i].x > maxx) { maxx = vlist[i].x; }
        if (vlist[i].y < miny) { miny = vlist[i].y; } else if (vlist[i].y > maxy) { maxy = vlist[i].y; }
        if (vlist[i].z < minz) { minz = vlist[i].z; } else if (vlist[i].z > maxz) { maxz = vlist[i].z; }
    }
}

//=============================================================================
// @ARG: none
//=============================================================================
void model::__comp_bbox() {
    int i;        // loop counter

    for (i = 0; i < nface; ++i) {
        flist[i].minx = vlist[flist[i].v1].x;
        flist[i].miny = vlist[flist[i].v1].y;
        flist[i].minz = vlist[flist[i].v1].z;
        flist[i].maxx = vlist[flist[i].v1].x;
        flist[i].maxy = vlist[flist[i].v1].y;
        flist[i].maxz = vlist[flist[i].v1].z;
        __check_vertex(i, vlist[flist[i].v2].x, vlist[flist[i].v2].y, vlist[flist[i].v2].z);
        __check_vertex(i, vlist[flist[i].v3].x, vlist[flist[i].v3].y, vlist[flist[i].v3].z);
    }
}

//=============================================================================
// @ARG: face index
// @ARG: vertex coordinates to be processed
//=============================================================================
void model::__check_vertex(const int idx, const double x, const double y, const double z) {
    if (x < flist[idx].minx) { flist[idx].minx = x; } else if (x > flist[idx].maxx) { flist[idx].maxx = x; }
    if (y < flist[idx].miny) { flist[idx].miny = y; } else if (y > flist[idx].maxy) { flist[idx].maxy = y; }
    if (z < flist[idx].minz) { flist[idx].minz = z; } else if (z > flist[idx].maxz) { flist[idx].maxz = z; }
}