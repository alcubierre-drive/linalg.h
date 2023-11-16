#define LINALG_EXPORT static inline
#include "linalg.h"
LINALG_IMPLEMENT

#include <stdio.h>

int main() {
    claV2cd_t vec = {1.,2.};
    claM2cd_t mat = claM2cd_id();
    claV2cd_t vp = claM2cd_matvec( mat, vec );
    printf("vp = %.2f %.2f\n", creal(vp.v[0]), creal(vp.v[1]));
}
