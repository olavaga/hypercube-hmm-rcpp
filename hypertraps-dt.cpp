#include <Rcpp.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "hypertraps-dt.c"

#ifdef __cplusplus
}
#endif

// [[Rcpp::export]]
int hypertrapsrcpp(char[] fn, int seed, int lengthindex, int kernelindex, int losses) {
    return hypertrapsc(fn, seed, lengthindex, kernelindex, losses);
}
