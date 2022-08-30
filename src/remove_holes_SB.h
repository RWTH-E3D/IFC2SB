// Copyright 2022 Eric Fichter
#ifndef REMOVE_HOLES_BY_TRIANGULATION_H
#define REMOVE_HOLES_BY_TRIANGULATION_H

#include <algorithm>
#include <array>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "extern/triangle.h"


namespace ifc2sb_triangle_hole_remove {

    typedef std::array<double, 2> Point;
    typedef std::list<Point> PointChain;
    typedef std::vector<std::pair<PointChain, Point>> PointChainsWithSpots;
    typedef std::pair<PointChain, Point> PointChainWithSpot;

    struct Triangle {

        int id, v1, v2, v3, region;
        std::pair<int, int> s1f, s1b, s2f, s2b, s3f, s3b;
        bool visited, tested;
        std::list<Triangle *> nbs;

        Triangle(int _id, int _v1, int _v2, int _v3, int _region) : id(_id), v1(_v1), v2(_v2), v3(_v3), region(_region) {
            s1f = std::make_pair(v1, v2);
            s1b = std::make_pair(v2, v1);
            s2f = std::make_pair(v2, v3);
            s2b = std::make_pair(v3, v2);
            s3f = std::make_pair(v3, v1);
            s3b = std::make_pair(v1, v3);
            visited = false;
            tested = false;
            nbs.clear();
        }
    };

    void report(struct triangulateio *io, int reporttriangles) {

        int i, j;

        for (i = 0; i < io->numberofpoints; i++) {
            printf("Point %4d:", i);
            for (j = 0; j < 2; j++) {
                printf("  %.6g", io->pointlist[i * 2 + j]);
            }
            if (io->numberofpointattributes > 0) {
                printf("   attributes");
            }
            for (j = 0; j < io->numberofpointattributes; j++) {
                printf("  %.6g",
                       io->pointattributelist[i * io->numberofpointattributes + j]);
            }
            printf("\n");
        }
        printf("\n");

        if (reporttriangles) {

            for (i = 0; i < io->numberoftriangles; i++) {
                printf("Triangle %4d points:", i);
                for (j = 0; j < io->numberofcorners; j++) {
                    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
                }
                if (io->numberoftriangleattributes > 0) {
                    printf("   attributes");
                }
                for (j = 0; j < io->numberoftriangleattributes; j++) {
                    printf("  %.6g", io->triangleattributelist[i * io->numberoftriangleattributes + j]);
                }
                printf("   neighbors:");
                for (j = 0; j < 3; j++) {
                    printf("  %4d", io->neighborlist[i * 3 + j]);
                }
                printf("\n");

            }
            printf("\n");
        }

        for (i = 0; i < io->numberofsegments; i++) {
            printf("Segment %4d points:", i);
            for (j = 0; j < 2; j++) {
                printf("  %4d", io->segmentlist[i * 2 + j]);
            }
            printf("\n");
        }
        printf("\n");

        for (i = 0; i < io->numberofholes; i++) {
            fprintf(stderr, "Hole %4d:", i);
            for (j = 0; j < 2; j++) {
                fprintf(stderr, "  %.6g", io->holelist[i * 2 + j]);
            }
            fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }

    void add_marking_points(struct triangulateio &in, const PointChainWithSpot &wire, const PointChainsWithSpots &holes, const PointChainsWithSpots &children) {

        int h = 0;

        // add hole points
        for (const auto &hole: holes) {
            in.holelist[h] = hole.second[0];
            in.holelist[h + 1] = hole.second[1];
            h += 2;
        }

        // add region point of wire
        in.regionlist[0] = wire.second[0];
        in.regionlist[1] = wire.second[1];
        in.regionlist[2] = 99999999;
        in.regionlist[3] = 99;

        int r = 4;

        // add region points
        // 0 is opening, 1 is normal
        for (unsigned int i = 0; i < children.size(); i++) {
            auto &child = children[i];
            in.regionlist[r] = child.second[0];
            in.regionlist[r + 1] = child.second[1];
            in.regionlist[r + 2] = i; // Regional attribute
            in.regionlist[r + 3] = 99; // Area constraint that will not be used.
            r += 4;
        }

    }

    void add_polygon(struct triangulateio &in, const PointChain &L, int &c, int &e) {

        // add points
        for (const auto &P: L) {
            in.pointlist[c] = P[0];
            in.pointlist[c + 1] = P[1];
            c += 2;
        }

        // add segments. ATTENTION indices start with 1 instead of 0
        int t = e / 2;
        for (int i = 1; i < L.size(); i++) {
            in.segmentlist[e] = t + i;
            in.segmentlist[e + 1] = t + i + 1;
            e += 2;
        }
        in.segmentlist[e] = t + static_cast<int>(L.size());
        in.segmentlist[e + 1] = t + 1;
        e += 2;
    }

    void calc_array_sizes(struct triangulateio &in, const PointChainWithSpot &wire, const PointChainsWithSpots &holes, const PointChainsWithSpots &children) {

        int c = static_cast<int>(wire.first.size());
        for (auto &v: holes) c += static_cast<int>(v.first.size());
        for (auto &v: children) c += static_cast<int>(v.first.size());

        in.numberofpoints = c;
        in.numberofsegments = c;
        in.numberofholes = static_cast<int>(holes.size());
        in.numberofregions = static_cast<int>(children.size() + 1);

        in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
        in.holelist = (double *) malloc(in.numberofholes * 2 * sizeof(double));
        in.regionlist = (double *) malloc(in.numberofregions * 4 * sizeof(double));
    }

    void free_memory(struct triangulateio &in, struct triangulateio &out) {
        free(in.pointlist);
        free(in.segmentlist);
        free(in.holelist);
        free(in.regionlist);

        free(out.pointlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        free(out.neighborlist);
    }

    void allocate_memory(struct triangulateio &in, struct triangulateio &out) {
        in.pointlist = nullptr;
        in.pointattributelist = nullptr;
        in.pointmarkerlist = nullptr;
        in.numberofpoints = 0;
        in.numberofpointattributes = 0;
        in.trianglelist = nullptr;
        in.triangleattributelist = nullptr;
        in.trianglearealist = nullptr;
        in.neighborlist = nullptr;
        in.numberoftriangles = 0;
        in.numberofcorners = 0;
        in.numberoftriangleattributes = 0;
        in.segmentlist = nullptr;
        in.segmentmarkerlist = nullptr;
        in.numberofsegments = 0;
        in.holelist = nullptr;
        in.numberofholes = 0;
        in.regionlist = nullptr;
        in.numberofregions = 0;
        in.edgelist = nullptr;
        in.edgemarkerlist = nullptr;
        in.normlist = nullptr;
        in.numberofedges = 0;

        out.pointlist = nullptr;
        out.pointattributelist = nullptr;
        out.pointmarkerlist = nullptr;
        out.trianglelist = nullptr;
        out.triangleattributelist = nullptr;
        out.trianglearealist = nullptr;
        out.neighborlist = nullptr;
        out.segmentlist = nullptr;
        out.segmentmarkerlist = nullptr;
        out.holelist = nullptr;
        out.regionlist = nullptr;
        out.edgelist = nullptr;
        out.edgemarkerlist = nullptr;
        out.normlist = nullptr;

        /*
        out.pointlist = (double *) nullptr;
        out.trianglelist = (int *) nullptr;
        out.triangleattributelist = (double *) nullptr;
        out.segmentlist = (int *) nullptr;
        out.segmentmarkerlist = (int *) nullptr;
        out.neighborlist = (int *) nullptr;
        */

        /*
        out.numberofedges = 0;
        out.numberofpoints = 0;
        out.numberofpointattributes = 0;
        out.numberofholes = 0;
        out.numberofregions = 0;
        out.numberofsegments = 0;
        out.numberoftriangles = 0;
        out.numberofcorners = 0;
        out.numberoftriangleattributes = 0;
        */
    }

    std::list<Triangle> get_triangles(struct triangulateio &out) {

        std::list<Triangle> TRIANGLES;
        std::map<int, Triangle *> M;

        for (int i = 0; i < out.numberoftriangles; i++) {
            int pos = i * 3;
            TRIANGLES.emplace_back(i + 1, out.trianglelist[pos], out.trianglelist[pos + 1], out.trianglelist[pos + 2], out.triangleattributelist[i * out.numberoftriangleattributes]);
            M[TRIANGLES.back().id] = &TRIANGLES.back();
        }

        // nbs
        for (int i = 0; i < out.numberoftriangles; i++) {

            auto &T = M[i + 1];

            for (int j = 0; j < 3; j++) {
                int v = out.neighborlist[i * 3 + j];
                if (v > 0) T->nbs.push_back(M[v]);
            }
        }

        return TRIANGLES;
    }

    Triangle *get_current_triangle(std::list<Triangle> &TRIANGLES) {

        for (auto &t: TRIANGLES)
            if (!t.visited)
                return &t;
        return nullptr;
    }

    bool get_loop_current(const std::map<std::pair<int, int>, std::list<bool>> &SINGLE_SEGMENTS, std::pair<int, int> &l) {

        for (auto &s: SINGLE_SEGMENTS)
            if (!s.second.front()) {
                l = s.first;
                return true;
            }
        return false;
    }

    bool get_loop_next(const std::map<std::pair<int, int>, std::list<bool>> &SINGLE_SEGMENTS, const std::pair<int, int> &l, std::pair<int, int> &key) {

        for (auto &SEG: SINGLE_SEGMENTS)
            if (l.second == SEG.first.first && !SEG.second.front()) {
                key = SEG.first;
                return true;
            }
        return false;
    }

    std::map<std::pair<int, int>, std::list<bool>> get_cluster_segments(const std::set<Triangle *> &TRIANGLES) {

        std::map<std::pair<int, int>, std::list<bool>> SEGMENTS;

        for (auto &t: TRIANGLES) {
            if (SEGMENTS.find(t->s1b) != SEGMENTS.end()) SEGMENTS[t->s1b].push_back(false);
            else SEGMENTS[t->s1f].push_back(false);

            if (SEGMENTS.find(t->s2b) != SEGMENTS.end()) SEGMENTS[t->s2b].push_back(false);
            else SEGMENTS[t->s2f].push_back(false);

            if (SEGMENTS.find(t->s3b) != SEGMENTS.end()) SEGMENTS[t->s3b].push_back(false);
            else SEGMENTS[t->s3f].push_back(false);
        }

        return SEGMENTS;
    }

    std::map<std::pair<int, int>, std::list<bool>> get_all_border_segments(const std::map<std::pair<int, int>, std::list<bool>> &ALL_SEGMENTS) {

        std::map<std::pair<int, int>, std::list<bool>> SINGLE_SEGMENTS;
        for (auto &m: ALL_SEGMENTS)
            if (m.second.size() == 1)
                SINGLE_SEGMENTS[m.first] = m.second;
        return SINGLE_SEGMENTS;
    }

    unsigned int get_number_of_loops(std::map<std::pair<int, int>, std::list<bool>> &SINGLE_SEGMENTS, bool &gibraltar) {

        unsigned int n_loops = 0;
        gibraltar = false;

        while (true) {

            std::pair<int, int> l;
            if (!get_loop_current(SINGLE_SEGMENTS, l)) break;

            SINGLE_SEGMENTS[l] = {true};
            n_loops += 1;

            //std::list<int> gibraltar1, gibraltar2;

            while (true) {

                std::pair<int, int> n;
                if (!get_loop_next(SINGLE_SEGMENTS, l, n)) break;

                SINGLE_SEGMENTS[n] = {true};
                //gibraltar1.push_back(n.first);
                //gibraltar2.push_back(n.second);

                l = n;
            }

/*            if (std::adjacent_find(gibraltar1.begin(), gibraltar1.end()) != gibraltar1.end()) {
                gibraltar = true;
                return 0;
            }
            if (std::adjacent_find(gibraltar2.begin(), gibraltar2.end()) != gibraltar2.end()) {
                gibraltar = true;
                return 0;
            }*/
        }

        return n_loops;
    }

    bool is_ok(const std::set<Triangle *> &Cluster, Triangle *next, bool &gibraltar) {

        std::set<Triangle *> TEMP = Cluster;
        TEMP.insert(next);

        auto ALL_SEGMENTS = get_cluster_segments(TEMP);
        auto SINGLE_SEGMENTS = get_all_border_segments(ALL_SEGMENTS);

        // vertex should only occur twice in loop
        std::map<int, std::list<bool>> M;
        for (auto &SEGMENT: SINGLE_SEGMENTS) {
            M[SEGMENT.first.first].push_back(true);
            M[SEGMENT.first.second].push_back(true);
            if (M[SEGMENT.first.first].size() > 2 || M[SEGMENT.first.second].size() > 2) return false;
        }

        unsigned int n_loops = get_number_of_loops(SINGLE_SEGMENTS, gibraltar);

        if (n_loops != 1 || gibraltar)
            return false;
        else
            return true;
    }

    bool find_ok_scenario_search_all_nbs(const std::map<int, std::set<Triangle *>> &all_nbs, std::set<Triangle *> &Cluster, Triangle *&next) {

        for (auto &m: all_nbs)
            for (auto &nb: m.second) {
                if (nb->visited || nb->tested) continue;

                bool gibraltar;

                if (is_ok(Cluster, nb, gibraltar)) {
                    next = nb;
                    return true;
                } else nb->tested = true;
            }

        return false;
    }

    std::vector<std::set<Triangle *>> get_clusters(std::list<Triangle> &TRIANGLES) {

        std::vector<std::set<Triangle *>> CLUSTERS;

        while (true) {

            auto current = get_current_triangle(TRIANGLES);
            if (current == nullptr) break;

            CLUSTERS.emplace_back();

            CLUSTERS.back().insert(current);
            current->visited = true;

            for (auto &Tri: TRIANGLES)
                Tri.tested = false;

            std::map<int, std::set<Triangle *>> all_nbs;

            while (true) {

                for (auto &nb: current->nbs)
                    if (!nb->visited and !nb->tested)
                        all_nbs[nb->region].insert(nb);

                Triangle *next = nullptr;
                if (!find_ok_scenario_search_all_nbs(all_nbs, CLUSTERS.back(), next))
                    break;

                CLUSTERS.back().insert(next);
                next->visited = true;
                current = next;
            }
        }

        return CLUSTERS;
    }

    std::list<int> get_loop(std::map<std::pair<int, int>, std::list<bool>> &SINGLE_SEGMENTS) {

        std::pair<int, int> l;
        get_loop_current(SINGLE_SEGMENTS, l);
        SINGLE_SEGMENTS[l] = {true};

        std::list<int> gibraltar1, gibraltar2;

        while (true) {

            std::pair<int, int> n;
            if (!get_loop_next(SINGLE_SEGMENTS, l, n)) break;

            SINGLE_SEGMENTS[n] = {true};
            gibraltar1.push_back(n.first);
            gibraltar2.push_back(n.second);

            l = n;
        }

        gibraltar1.push_back(gibraltar2.back());
        gibraltar1.push_back(gibraltar1.front());

        return gibraltar1;
    }

    std::vector<std::list<int>> get_loops(const std::vector<std::set<Triangle *>> &CLUSTERS) {

        std::vector<std::list<int>> loops;

        for (auto &C: CLUSTERS) {

            auto ALL_SEGMENTS = get_cluster_segments(C);
            auto SINGLE_SEGMENTS = get_all_border_segments(ALL_SEGMENTS);

            std::list<int> l = get_loop(SINGLE_SEGMENTS);
            loops.push_back(l);

        }
        return loops;
    }

    bool orientation(const PointChain &t) {

        // true for counterclockwise winding
        // https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order

        std::vector<Point> l(t.begin(), t.end());
        double v = 0;

        for (unsigned int i = 0; i < l.size() - 1; i++)
            v += (l[i + 1][0] - l[i][0]) * (l[i + 1][1] + l[i][1]);
        v += (l[0][0] - l.back()[0]) * (l[0][1] + l.back()[1]);

        return v < 0;
    }

    std::list<std::pair<PointChain, std::set<unsigned int>>> get_result(const std::vector<std::set<Triangle *>> &CLUSTERS, const std::vector<std::list<int>> &LOOPS, struct triangulateio &out) {

        std::list<std::pair<PointChain, std::set<unsigned int>>> RESULT;

        for (unsigned int m = 0; m < CLUSTERS.size(); m++) {
            auto &C = CLUSTERS[m];
            auto &L = LOOPS[m];

            PointChain chain;
            std::set<unsigned int> regions;

            for (auto &T: C)
                if (T->region < 99999998) regions.insert(T->region);

            for (auto &i: L)
                chain.push_back({out.pointlist[(i - 1) * 2], out.pointlist[(i - 1) * 2 + 1]});

            RESULT.emplace_back(chain, regions);
        }

        return RESULT;
    }

    std::list<PointChain> simplify_result(const std::list<std::pair<PointChain, std::set<unsigned int>>> &res) {

        std::list<PointChain> RESULT;

        for (const auto &r: res) {

            PointChain chain;
            for (const auto &P_2D: r.first)
                chain.push_back({P_2D[0], P_2D[1]});

            RESULT.push_back(chain);
        }

        return RESULT;
    }

    std::list<PointChain> run(const PointChainWithSpot &wire, const PointChainsWithSpots &holes, const PointChainsWithSpots &children) {

        /*
         * The whole procedure is not really time efficient at the moment.
         * For the test buildings it does the job.
         * Most time-consuming are the "ok checks" before adding new triangle to a cluster.
         * That's because the segments of a cluster are not stored, but are calculated everytime a new triangle is tested.
         * Also, the check for loops can simplified.
         * IF this is done, quality mesh 'q' should be used instead of 'S0'
         */

        struct triangulateio in, out;

        allocate_memory(in, out);
        calc_array_sizes(in, wire, holes, children);

        int c(0), e(0);
        add_polygon(in, wire.first, c, e);

        for (const auto &hole: holes)
            add_polygon(in, hole.first, c, e);

        for (const auto &child: children)
            add_polygon(in, child.first, c, e);

        add_marking_points(in, wire, holes, children);

        std::string s = "pAnIQ"; // 'I' cannot be used when 'q' is used. std::string s = "pS0AnQ";
        char *cstr = &s[0];

        try { triangulate(cstr, &in, &out, nullptr); }
        catch (...) {
            std::cerr << "[Warning] Triangulation failed!" << std::endl;
            std::list<PointChain> RESULT2;
            return RESULT2;
        }

        auto TRIANGLES = get_triangles(out);
        auto CLUSTERS = get_clusters(TRIANGLES);
        auto LOOPS = get_loops(CLUSTERS);
        auto RESULT = get_result(CLUSTERS, LOOPS, out);
        auto RESULT2 = simplify_result(RESULT);

/*        for (auto &T: TRIANGLES) {
            std::cerr << "TRI: " << T.id << "\t" << T.v1 << "\t" << T.v2 << "\t" << T.v3 << "\tR: " << T.region << "\tNBS: " << T.nbs.size() << "\t";
            for (auto &nb: T.nbs)
                std::cerr << "\t" << nb->id;
            std::cerr << std::endl;
        }*/
/*        for (auto &S: SEGMENTS)
            for (auto &T: S.second)
                std::cerr << "SEG: " << S.first.first << "\t" << S.first.second << "\t" << T.id << std::endl;*/

/*        for (unsigned int m = 0; m < CLUSTERS.size(); m++) {
            auto &C = CLUSTERS[m];
            auto &L = LOOPS[m];
            std::cerr << "CLUSTER: " << std::endl;
            for (auto &T: C)
                std::cerr << "\t" << T->id << "(" << T->region << ")";
            std::cerr << std::endl;
            for (auto &i: L) {
                std::cerr << "\t" << i << "(";
                for (int j = 0; j < 2; j++) {
                    std::cerr << " " << out.pointlist[(i - 1) * 2 + j] << " ";
                }
                std::cerr << ")";
            }
            std::cerr << std::endl;
        }*/

/*        for (auto &C: CLUSTERS) {
            std::cerr << "CLUSTER: " << std::endl;
            for (auto &T: C)
                std::cerr << "\t" << T->id << "(" << T->region << ")";
            std::cerr << std::endl;
        }

        for (auto &L: LOOPS) {
            std::cerr << "LOOP: " << std::endl;
            for (auto &i: L) {
                std::cerr << "\t" << i << "(";
                for (int j = 0; j < 2; j++) {
                    std::cerr << " " << out.pointlist[(i - 1) * 2 + j] << " ";
                }
                std::cerr << ")";
            }
            std::cerr << std::endl;
        }*/

/*        std::cerr << "IN" << std::endl;
        report(&in, 0);
        std::cerr << "OUT" << std::endl;
        report(&out, 1);*/

        free_memory(in, out);

        return RESULT2;
    }

}

#endif //REMOVE_HOLES_BY_TRIANGULATION_H