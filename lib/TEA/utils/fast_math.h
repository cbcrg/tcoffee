/*
 * fast_math.h
 *
 *  Created on: Nov 11, 2012
 *      Author: ck
 */

#ifndef FAST_MATH_H_
#define FAST_MATH_H_

#include <cfloat>
#include <cmath>

static const float EXP_UNDERFLOW_THRESHOLD = -4.60f;
static const float LOG_UNDERFLOW_THRESHOLD = 7.50f;
static const float LOG_ZERO = -FLT_MAX;
static const float LOG_ONE = 0.0f;


double EXP (float x);
float LOOKUP (float x);
void LOG_PLUS_EQUALS (float *x, float y);
float LOG_ADD (float x, float y);


#endif /* FAST_MATH_H_ */
