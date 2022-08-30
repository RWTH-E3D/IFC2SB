// https://raw.githubusercontent.com/mjjq/ConvexDecomposition/master/src/ConcavePolygon.h

#ifndef CONCAVE_POLY_H
#define CONCAVE_POLY_H

#include <vector>
#include <cmath>
#include <map>
#include <iostream>

namespace cxdsb {

    struct Vec2 {
        double x;
        double y;

        static double length(Vec2 const &v) {
            return sqrt(v.x * v.x + v.y * v.y);
        }

        static Vec2 norm(Vec2 const &v) {
            if (length(v) < 1e-30)
                return {0.0f, 0.0f};
        }

        static double dot(Vec2 const &v1, Vec2 const &v2) {
            return v1.x * v2.x + v1.y * v2.y;
        }

        static double square(Vec2 const &v) {
            return dot(v, v);
        }

        static double cross(Vec2 const &v1, Vec2 const &v2) {
            return v1.x * v2.y - v1.y * v2.x;
        }

        Vec2 operator-(Vec2 const &v1) const {
            return {x - v1.x, y - v1.y};
        }

        Vec2 operator+(Vec2 const &v1) const {
            return {x + v1.x, y + v1.y};
        }

        Vec2 operator*(double const &f) const {
            return {f * x, f * y};
        }

        Vec2 operator/(double const &f) const {
            return {x / f, y / f};
        }

        static double getSignedArea(Vec2 const &v1,
                                   Vec2 const &v2) {
            return (v2.x - v1.x) * (v2.y + v1.y);
        }
    };

    struct Vertex {
        Vec2 position;

        Vertex() {}

        Vertex(Vec2 const &_position) : position{_position} {}

        static double getHandedness(Vertex const &v1,
                                   Vertex const &v2,
                                   Vertex const &v3) {
            Vec2 edge1 = v2.position - v1.position;
            Vec2 edge2 = v3.position - v2.position;

            return Vec2::cross(edge1, edge2);
        }
    };

    struct SliceVertex : public Vertex {
        SliceVertex() {}

        SliceVertex(Vec2 const &_position) : Vertex{_position} {}

        int index;
        double distanceToSlice;
    };

    struct LineSegment {
        Vec2 startPos;
        Vec2 finalPos;

        LineSegment() {}

        LineSegment(Vec2 const &_startPos,
                    Vec2 const &_finalPos) : startPos{_startPos},
                                             finalPos{_finalPos} {}

        Vec2 direction() const {
            return finalPos - startPos;
        }

        Vec2 normalisedDirection() const {
            return Vec2::norm(finalPos - startPos);
        }

        LineSegment operator+(LineSegment const &ls) const {
            Vec2 newStartPos = (startPos + ls.startPos) / 2.0f;
            Vec2 newFinalPos = (finalPos + ls.finalPos) / 2.0f;

            return {newStartPos, newFinalPos};
        }

        static std::pair<bool, Vec2> intersects(LineSegment s1, LineSegment s2) {
            const double TOLERANCE = 1e-7;

            Vec2 p1 = s1.startPos;
            Vec2 p2 = s2.startPos;
            Vec2 d1 = s1.direction();
            Vec2 d2 = s2.direction();

            if (std::abs(Vec2::cross(d1, d2)) < 1e-30)
                return {false, {0.0f, 0.0f}};

            double t1 = Vec2::cross(p2 - p1, d2) / Vec2::cross(d1, d2);

            if ((t1 < (0.0f - TOLERANCE)) || (t1 > (1.0f + TOLERANCE)))
                return {false, {0.0f, 0.0f}};

            Vec2 pIntersect = p1 + d1 * t1;

            double t2 = Vec2::dot(pIntersect - p2,
                                 s2.finalPos - p2);

            if (t2 < (0.0f - TOLERANCE) || t2 / Vec2::square(s2.finalPos - p2) >= 1.0f - TOLERANCE)
                return {false, {0.0f, 0.0f}};

            return {true, pIntersect};
        }
    };

    class ConcavePolygon {
        typedef std::vector<Vertex> VertexArray;
        typedef std::vector<ConcavePolygon> PolygonArray;

        VertexArray vertices;
        PolygonArray subPolygons;


        static int mod(int x, int m) {
            int r = x % m;
            return r < 0 ? r + m : r;
        }

        static void flipPolygon(VertexArray &_verts) {
            int iMax = _verts.size() / 2;

            if (_verts.size() % 2 != 0)
                iMax += 1;

            for (int i = 1; i < iMax; ++i) {
                std::swap(_verts[i], _verts[_verts.size() - i]);
            }
        }

        static bool checkIfRightHanded(VertexArray &_verts) {
            if (_verts.size() < 3)
                return false;

            double signedArea = 0.0f;

            for (unsigned int i = 0; i < _verts.size(); ++i) {
                signedArea += Vec2::getSignedArea(_verts[i].position,
                                                  _verts[mod(i + 1, _verts.size())].position);
            }

            if (signedArea < 0.0f)
                return true;

            return false;
        }


        static bool isVertexInCone(LineSegment const &ls1,
                            LineSegment const &ls2,
                            Vec2 const &origin,
                            Vertex const &vert) {
            Vec2 relativePos = vert.position - origin;

            double ls1Product = Vec2::cross(relativePos, ls1.direction());
            double ls2Product = Vec2::cross(relativePos, ls2.direction());

            if (ls1Product < 0.0f && ls2Product > 0.0f)
                return true;

            return false;
        }

        typedef std::vector<int> IntArray;

        static IntArray findVerticesInCone(LineSegment const &ls1,
                                    LineSegment const &ls2,
                                    Vec2 const &origin,
                                    VertexArray const &inputVerts) {
            IntArray result;

            for (unsigned int i = 0; i < inputVerts.size(); ++i) {
                if (isVertexInCone(ls1, ls2, origin, inputVerts[i]))
                    result.push_back(i);
            }
            return result;
        }

        static bool checkVisibility(Vec2 const &originalPosition,
                             Vertex const &vert,
                             VertexArray const &polygonVertices) {
            LineSegment ls(originalPosition, vert.position);
            VertexIntMap intersectingVerts = verticesAlongLineSegment(ls, polygonVertices);

            //std::cout << intersectingVerts.size() << " intverts\n";

            if (intersectingVerts.size() > 3)
                return false;

            return true;
        }


        static int getBestVertexToConnect(IntArray const &indices,
                                   VertexArray const &polygonVertices,
                                   Vec2 const &origin) {
            if (indices.size() == 1) {
                if (checkVisibility(origin, polygonVertices[indices[0]], polygonVertices))
                    return indices[0];
            } else if (indices.size() > 1) {
                for (int index : indices) {
                    int vertSize = polygonVertices.size();

                    Vertex prevVert = polygonVertices[mod(index - 1, vertSize)].position;
                    Vertex currVert = polygonVertices[index].position;
                    Vertex nextVert = polygonVertices[mod(index + 1, vertSize)].position;

                    LineSegment ls1(prevVert.position, currVert.position);
                    LineSegment ls2(nextVert.position, currVert.position);

                    if ((Vertex::getHandedness(prevVert, currVert, nextVert) < 0.0f) &&
                        isVertexInCone(ls1, ls2, polygonVertices[index].position, Vertex(origin)) &&
                        checkVisibility(origin, polygonVertices[index], polygonVertices))
                        return index;
                }

                for (int index : indices) {
                    int vertSize = polygonVertices.size();

                    Vertex prevVert = polygonVertices[mod(index - 1, vertSize)].position;
                    Vertex currVert = polygonVertices[index].position;
                    Vertex nextVert = polygonVertices[mod(index + 1, vertSize)].position;

                    LineSegment ls1(prevVert.position, currVert.position);
                    LineSegment ls2(nextVert.position, currVert.position);

                    if ((Vertex::getHandedness(prevVert, currVert, nextVert) < 0.0f) &&
                        checkVisibility(origin, polygonVertices[index], polygonVertices))
                        return index;
                }


                double minDistance = 1e+15;
                int closest = indices[0];
                for (int index : indices) {
                    double currDistance = Vec2::square(polygonVertices[index].position - origin);
                    if (currDistance < minDistance) {
                        minDistance = currDistance;
                        closest = index;
                    }
                }

                return closest;
            }


            return -1;
        }

        void convexDecomp(VertexArray const &_vertices, unsigned int &count) {
            if (!subPolygons.empty()) {
                return;
            }

            int reflexIndex = findFirstReflexVertex(_vertices);
            if (reflexIndex == -1)
                return;

            Vec2 prevVertPos = _vertices[mod(reflexIndex - 1, _vertices.size())].position;
            Vec2 currVertPos = _vertices[reflexIndex].position;
            Vec2 nextVertPos = _vertices[mod(reflexIndex + 1, _vertices.size())].position;

            LineSegment ls1(prevVertPos, currVertPos);
            LineSegment ls2(nextVertPos, currVertPos);

            IntArray vertsInCone = findVerticesInCone(ls1, ls2, currVertPos, _vertices);

            int bestVert = -1;

            if (!vertsInCone.empty()) {
                bestVert = getBestVertexToConnect(vertsInCone, _vertices, currVertPos);
                if (bestVert != -1) {
                    LineSegment newLine(currVertPos, _vertices[bestVert].position);

                    slicePolygon(newLine);
                }
            }
            if (vertsInCone.empty() || bestVert == -1) {
                LineSegment newLine(currVertPos, (ls1.direction() + ls2.direction()) * 1e+10);
                slicePolygon(newLine);
            }

            for (auto & subPolygon : subPolygons) {
                subPolygon.convexDecomp(count);
            }
        }

        static int findFirstReflexVertex(VertexArray const &_vertices) {
            for (unsigned int i = 0; i < _vertices.size(); ++i) {
                double handedness = Vertex::getHandedness(_vertices[mod(i - 1, _vertices.size())],
                                                         _vertices[i],
                                                         _vertices[mod(i + 1, _vertices.size())]);
                if (handedness < 0.0f)
                    return i;
            }

            return -1;
        }

        void flipPolygon() {
            flipPolygon(vertices);
        }

        typedef std::map<int, Vertex> VertexIntMap;
        typedef std::pair<int, Vertex> VertexIntPair;


        static VertexIntMap cullByDistance(VertexIntMap const &input,
                                    Vec2 const &origin,
                                    int const &maxVertsToKeep) {
            if (maxVertsToKeep >= (int) input.size())
                return input;

            std::vector<SliceVertex> sliceVertices;

            for (const auto & it : input) {
                SliceVertex vert(it.second.position);
                vert.index = it.first;
                vert.distanceToSlice = Vec2::square(it.second.position - origin);

                sliceVertices.push_back(vert);
            }

            for (unsigned int i = 1; i < sliceVertices.size(); ++i)
                for (unsigned int j = i; j > 0 && sliceVertices[j].distanceToSlice < sliceVertices[j - 1].distanceToSlice; --j)
                    std::swap(sliceVertices[j], sliceVertices[j - 1]);

            sliceVertices.erase(sliceVertices.begin() + maxVertsToKeep, sliceVertices.end());

            for (unsigned int i = 1; i < sliceVertices.size(); ++i)
                for (unsigned int j = i; j > 0 && sliceVertices[j].index < sliceVertices[j - 1].index; --j)
                    std::swap(sliceVertices[j], sliceVertices[j - 1]);

            VertexIntMap result;
            for (auto & sliceVertice : sliceVertices) {
                result.insert({sliceVertice.index, Vertex(sliceVertice.position)});
            }

            return result;
        }

        static VertexIntMap verticesAlongLineSegment(LineSegment const &segment,
                                              VertexArray const &_vertices) {
            VertexIntMap result;

            LineSegment tempSegment;

            for (unsigned int i = 0; i < _vertices.size(); ++i) {
                tempSegment.startPos = _vertices[i].position;
                tempSegment.finalPos = _vertices[mod(i + 1, _vertices.size())].position;

                std::pair<bool, Vec2> intersectionResult = LineSegment::intersects(segment, tempSegment);

                if (intersectionResult.first) {
                    result.insert({i, Vertex(intersectionResult.second)});
                }
            }

            return result;
        }

    public:
        ConcavePolygon(VertexArray const &_vertices) : vertices{_vertices} {
            if (vertices.size() > 2)
                if (!checkIfRightHanded())
                    flipPolygon();
        }

        ConcavePolygon() {}

        bool checkIfRightHanded() {
            return checkIfRightHanded(vertices);
        }

        void slicePolygon(int vertex1, int vertex2) {
            if (vertex1 == vertex2 ||
                vertex2 == vertex1 + 1 ||
                vertex2 == vertex1 - 1)
                return;

            if (vertex1 > vertex2)
                std::swap(vertex1, vertex2);

            VertexArray returnVerts;
            VertexArray newVerts;
            for (int i = 0; i < (int) vertices.size(); ++i) {
                if (i == vertex1 || i == vertex2) {
                    returnVerts.push_back(vertices[i]);
                    newVerts.push_back(vertices[i]);
                } else if (i > vertex1 && i < vertex2)
                    returnVerts.push_back(vertices[i]);
                else
                    newVerts.push_back(vertices[i]);
            }

            subPolygons.push_back(ConcavePolygon(returnVerts));
            subPolygons.push_back(ConcavePolygon(newVerts));
        }

        void slicePolygon(LineSegment segment) {
            if (!subPolygons.empty()) {
                subPolygons[0].slicePolygon(segment);
                subPolygons[1].slicePolygon(segment);
                return;
            }

            const double TOLERANCE = 1e-7;

            VertexIntMap slicedVertices = verticesAlongLineSegment(segment, vertices);
            slicedVertices = cullByDistance(slicedVertices, segment.startPos, 2);

            if (slicedVertices.size() < 2)
                return;

            VertexArray leftVerts;
            VertexArray rightVerts;

            for (int i = 0; i < (int) vertices.size(); ++i) {
                Vec2 relativePosition = vertices[i].position - segment.startPos;

                auto it = slicedVertices.begin();

                double perpDistance = std::abs(Vec2::cross(relativePosition, segment.direction()));

                if (perpDistance > TOLERANCE ||
                    (perpDistance <= TOLERANCE && (slicedVertices.find(i) == slicedVertices.end()))
                        ) {
                    //std::cout << relCrossProd << ", i: " << i << "\n";
                    if ((i > it->first) && (i <= (++it)->first)) {
                        leftVerts.push_back(vertices[i]);
                        //std::cout << i << " leftVertAdded\n";
                    } else {
                        rightVerts.push_back(vertices[i]);
                        //std::cout << i << " rightVertAdded\n";
                    }

                }

                if (slicedVertices.find(i) != slicedVertices.end()) {
                    rightVerts.push_back(slicedVertices[i]);
                    leftVerts.push_back(slicedVertices[i]);
                }
            }

            subPolygons.push_back(ConcavePolygon(leftVerts));
            subPolygons.push_back(ConcavePolygon(rightVerts));
        }

        void convexDecomp(unsigned int &count) {

            count++;

            if (count > 300) {
                std::cerr << "[Warning] Decomposition stopped at recursion 300" << std::endl; // Often because of collinear points or to many splits
                return;
            }

            if (vertices.size() > 3)
                convexDecomp(vertices, count);
        }

        VertexArray const &getVertices() const {
            return vertices;
        }

        ConcavePolygon const &getSubPolygon(int subPolyIndex) const {
            if (!subPolygons.empty() && subPolyIndex < (int) subPolygons.size())
                return subPolygons[subPolyIndex];

            return *this;
        }

        int getNumberSubPolys() const {
            return subPolygons.size();
        }

        void returnLowestLevelPolys(std::vector<ConcavePolygon> &returnArr) {
            if (!subPolygons.empty()) {
                subPolygons[0].returnLowestLevelPolys(returnArr);
                subPolygons[1].returnLowestLevelPolys(returnArr);
            } else
                returnArr.push_back(*this);
        }

        void reset() {
            if (!subPolygons.empty()) {
                subPolygons[0].reset();
                subPolygons[1].reset();
                subPolygons.clear();
            }
        }

        Vec2 getPoint(unsigned int index) const {
            if (index >= 0 && index < vertices.size())
                return vertices[index].position;

            return {0.0f, 0.0f};
        }

        int getPointCount() const {
            return vertices.size();
        }
    };

}

#endif // CONCAVE_POLY_H

