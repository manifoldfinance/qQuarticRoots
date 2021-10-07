#include <stdio.h>
#include <stdlib.h>

#include "quartic_real.h"
#include "k.h"

K quarticRoots(K ka, K kb, K kc, K kd, K ke) {
    double poly[5];
    poly[4] = ka->f;
    poly[3] = kb->f;
    poly[2] = kc->f;
    poly[1] = kd->f;
    poly[0] = ke->f;
    double sols[4];
    const int num_sols = solve_real_poly(4, poly, sols);
    K retz = ktn(KF,4);
    kF(retz)[0]=sols[0];
    kF(retz)[1]=sols[1];
    kF(retz)[2]=sols[2];
    kF(retz)[3]=sols[3];
    return(retz);
}