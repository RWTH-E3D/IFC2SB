#ifndef AKENINE_INCLUDES_H
#define AKENINE_INCLUDES_H

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#include <cmath>
#include <cstdio>

namespace akine_moeller_isct {

#define X_octree 0
#define Y_octree 1
#define Z_octree 2

#define CROSS_octree(dest, v1, v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT_octree(v1, v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB_octree(dest, v1, v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

#define FINDMINMAX_octree(x0, x1, x2, min, max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

    inline int planeBoxOverlap(const float normal[3], const float vert[3], const float maxbox[3])    // -NJMP-
    {
        int q;
        float vmin[3], vmax[3], v;
        for (q = X_octree; q <= Z_octree; q++) {
            v = vert[q];                    // -NJMP-
            if (normal[q] > 0.0f) {
                vmin[q] = -maxbox[q] - v;    // -NJMP-
                vmax[q] = maxbox[q] - v;    // -NJMP-
            } else {
                vmin[q] = maxbox[q] - v;    // -NJMP-
                vmax[q] = -maxbox[q] - v;    // -NJMP-
            }
        }
        if (DOT_octree(normal, vmin) > 0.0f) return 0;    // -NJMP-
        if (DOT_octree(normal, vmax) >= 0.0f) return 1;    // -NJMP-

        return 0;
    }


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)               \
    p0 = a*v0[Y_octree] - b*v0[Z_octree];                       \
    p2 = a*v2[Y_octree] - b*v2[Z_octree];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[Y_octree] + fb * boxhalfsize[Z_octree];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)               \
    p0 = a*v0[Y_octree] - b*v0[Z_octree];                       \
    p1 = a*v1[Y_octree] - b*v1[Z_octree];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[Y_octree] + fb * boxhalfsize[Z_octree];   \
    if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)               \
    p0 = -a*v0[X_octree] + b*v0[Z_octree];                   \
    p2 = -a*v2[X_octree] + b*v2[Z_octree];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[X_octree] + fb * boxhalfsize[Z_octree];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)               \
    p0 = -a*v0[X_octree] + b*v0[Z_octree];                   \
    p1 = -a*v1[X_octree] + b*v1[Z_octree];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X_octree] + fb * boxhalfsize[Z_octree];   \
    if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)               \
    p1 = a*v1[X_octree] - b*v1[Y_octree];                       \
    p2 = a*v2[X_octree] - b*v2[Y_octree];                       \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
    rad = fa * boxhalfsize[X_octree] + fb * boxhalfsize[Y_octree];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)               \
    p0 = a*v0[X_octree] - b*v0[Y_octree];                   \
    p1 = a*v1[X_octree] - b*v1[Y_octree];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X_octree] + fb * boxhalfsize[Y_octree];   \
    if(min>rad || max<-rad) return 0;

    inline int triBoxOverlap(const float boxcenter[3], float boxhalfsize[3], float triverts[3][3]) {

        /*    use separating axis theorem to test overlap between triangle and box */
        /*    need to test for overlap in these directions: */
        /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
        /*       we do not even need to test these) */
        /*    2) normal of the triangle */
        /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
        /*       this gives 3x3=9 more tests */
        float v0[3], v1[3], v2[3];
//   float axis[3];
        float min, max, p0, p1, p2, rad, fex, fey, fez;        // -NJMP- "d" local variable removed
        float normal[3], e0[3], e1[3], e2[3];

        /* This is the fastest branch on Sun */
        /* move everything so that the boxcenter is in (0,0,0) */
        SUB_octree(v0, triverts[0], boxcenter);
        SUB_octree(v1, triverts[1], boxcenter);
        SUB_octree(v2, triverts[2], boxcenter);

        /* compute triangle edges */
        SUB_octree(e0, v1, v0);      /* tri edge 0 */
        SUB_octree(e1, v2, v1);      /* tri edge 1 */
        SUB_octree(e2, v0, v2);      /* tri edge 2 */

        /* Bullet 3:  */
        /*  test the 9 tests first (this was faster) */
        fex = fabsf(e0[X_octree]);
        fey = fabsf(e0[Y_octree]);
        fez = fabsf(e0[Z_octree]);
        AXISTEST_X01(e0[Z_octree], e0[Y_octree], fez, fey);
        AXISTEST_Y02(e0[Z_octree], e0[X_octree], fez, fex);
        AXISTEST_Z12(e0[Y_octree], e0[X_octree], fey, fex);

        fex = fabsf(e1[X_octree]);
        fey = fabsf(e1[Y_octree]);
        fez = fabsf(e1[Z_octree]);
        AXISTEST_X01(e1[Z_octree], e1[Y_octree], fez, fey);
        AXISTEST_Y02(e1[Z_octree], e1[X_octree], fez, fex);
        AXISTEST_Z0(e1[Y_octree], e1[X_octree], fey, fex);

        fex = fabsf(e2[X_octree]);
        fey = fabsf(e2[Y_octree]);
        fez = fabsf(e2[Z_octree]);
        AXISTEST_X2(e2[Z_octree], e2[Y_octree], fez, fey);
        AXISTEST_Y1(e2[Z_octree], e2[X_octree], fez, fex);
        AXISTEST_Z12(e2[Y_octree], e2[X_octree], fey, fex);

        /* Bullet 1: */
        /*  first test overlap in the {x,y,z}-directions */
        /*  find min, max of the triangle each direction, and test for overlap in */
        /*  that direction -- this is equivalent to testing a minimal AABB around */
        /*  the triangle against the AABB */

        /* test in X-direction */
        FINDMINMAX_octree(v0[X_octree], v1[X_octree], v2[X_octree], min, max);
        if (min > boxhalfsize[X_octree] || max < -boxhalfsize[X_octree]) return 0;

        /* test in Y-direction */
        FINDMINMAX_octree(v0[Y_octree], v1[Y_octree], v2[Y_octree], min, max);
        if (min > boxhalfsize[Y_octree] || max < -boxhalfsize[Y_octree]) return 0;

        /* test in Z-direction */
        FINDMINMAX_octree(v0[Z_octree], v1[Z_octree], v2[Z_octree], min, max);
        if (min > boxhalfsize[Z_octree] || max < -boxhalfsize[Z_octree]) return 0;

        /* Bullet 2: */
        /*  test if the box intersects the plane of the triangle */
        /*  compute plane equation of triangle: normal*x+d=0 */
        CROSS_octree(normal, e0, e1);
        // -NJMP- (line removed here)
        if (!planeBoxOverlap(normal, v0, boxhalfsize)) return 0;    // -NJMP-

        return 1;   /* box and triangle overlaps */
    }
}

#endif //AKENINE_INCLUDES_H