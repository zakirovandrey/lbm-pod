#ifdef SEMI_LAGRANGIAN
#include "streaming-pod.cu"
#else
#include "streaming-usual.cu"
struct MomentsMatrix{};
#endif
