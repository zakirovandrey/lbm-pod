#define DIM 3
#define QN 64
#define QN_IN_DIM 512   // QN*2^DIM

const int ec1=1;
const int ec2=4;
const int ec1d=ec1*ec1;
const int ec2d=ec2*ec2;

#include "REFERENCE_TEMP_AND_WEIGHTS-SET-4.h"

//#define NON_ISOTHERMAL_RELAXATION

constexpr const int3 _e[QN] = {
  {+ec1, +ec1, +ec1 }, {+ec1, +ec1, -ec1}, {+ec1, +ec1, +ec2}, {+ec1, +ec1, -ec2}, // 30--34
  {+ec1, -ec1, +ec1 }, {+ec1, -ec1, -ec1}, {+ec1, -ec1, +ec2}, {+ec1, -ec1, -ec2}, // 35--39
  {+ec1, +ec2, +ec1 }, {+ec1, +ec2, -ec1}, {+ec1, +ec2, +ec2}, {+ec1, +ec2, -ec2}, // 40--44
  {+ec1, -ec2, +ec1 }, {+ec1, -ec2, -ec1}, {+ec1, -ec2, +ec2}, {+ec1, -ec2, -ec2}, // 45--49
  {-ec1, +ec1, +ec1 }, {-ec1, +ec1, -ec1}, {-ec1, +ec1, +ec2}, {-ec1, +ec1, -ec2}, // 55--59
  {-ec1, -ec1, +ec1 }, {-ec1, -ec1, -ec1}, {-ec1, -ec1, +ec2}, {-ec1, -ec1, -ec2}, // 60--64
  {-ec1, +ec2, +ec1 }, {-ec1, +ec2, -ec1}, {-ec1, +ec2, +ec2}, {-ec1, +ec2, -ec2}, // 65--69
  {-ec1, -ec2, +ec1 }, {-ec1, -ec2, -ec1}, {-ec1, -ec2, +ec2}, {-ec1, -ec2, -ec2}, // 70--74
  {+ec2, +ec1, +ec1 }, {+ec2, +ec1, -ec1}, {+ec2, +ec1, +ec2}, {+ec2, +ec1, -ec2}, // 80--84
  {+ec2, -ec1, +ec1 }, {+ec2, -ec1, -ec1}, {+ec2, -ec1, +ec2}, {+ec2, -ec1, -ec2}, // 85--89
  {+ec2, +ec2, +ec1 }, {+ec2, +ec2, -ec1}, {+ec2, +ec2, +ec2}, {+ec2, +ec2, -ec2}, // 90--94
  {+ec2, -ec2, +ec1 }, {+ec2, -ec2, -ec1}, {+ec2, -ec2, +ec2}, {+ec2, -ec2, -ec2}, // 95--99
  {-ec2, +ec1, +ec1 }, {-ec2, +ec1, -ec1}, {-ec2, +ec1, +ec2}, {-ec2, +ec1, -ec2}, // 105--109
  {-ec2, -ec1, +ec1 }, {-ec2, -ec1, -ec1}, {-ec2, -ec1, +ec2}, {-ec2, -ec1, -ec2}, // 110--114
  {-ec2, +ec2, +ec1 }, {-ec2, +ec2, -ec1}, {-ec2, +ec2, +ec2}, {-ec2, +ec2, -ec2}, // 115--119
  {-ec2, -ec2, +ec1 }, {-ec2, -ec2, -ec1}, {-ec2, -ec2, +ec2}, {-ec2, -ec2, -ec2}, // 120--124
};

#define WEIGHTS_MANUALLY
constexpr ftype _w[QN] = {
  W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
};


