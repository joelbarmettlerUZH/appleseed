
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2013 Francois Beaune, Jupiter Jazz Limited
// Copyright (c) 2014-2015 Francois Beaune, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#ifndef APPLESEED_FOUNDATION_MATH_MIS_H
#define APPLESEED_FOUNDATION_MATH_MIS_H

// appleseed.foundation headers.
#include "foundation/math/minmax.h"

// Standard headers.
#include <algorithm>
#include <cassert>
#include <cmath>

namespace foundation
{

//
// Multiple Importance Sampling (MIS) heuristics.
//
// Reference:
//
//   Robust Monte Carlo Methods For Light Transport Simulation
//   http://graphics.stanford.edu/papers/veach_thesis/thesis.pdf
//

// Balance heuristic.
template <typename T> T mis_balance(const T q1, const T q2);
template <typename T> T mis_balance(const T q1, const T q2, const T q3);

// Power heuristic. beta is >= 0.
template <typename T> T mis_power(const T q1, const T q2, const T beta);
template <typename T> T mis_power(const T q1, const T q2, const T q3, const T beta);

// Power heuristic with beta = 2.
template <typename T> T mis_power2(const T q1, const T q2);
template <typename T> T mis_power2(const T q1, const T q2, const T q3);

// Cutoff heuristic. The cutoff threshold alpha is in [0,1].
template <typename T> T mis_cutoff(const T q1, const T q2, const T alpha);
template <typename T> T mis_cutoff(const T q1, const T q2, const T q3, const T alpha);

// Maximum heuristic.
template <typename T> T mis_maximum(const T q1, const T q2);
template <typename T> T mis_maximum(const T q1, const T q2, const T q3);


//
// Implementation.
//

template <typename T>
inline T mis_balance(const T q1, const T q2)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q1 + q2 > T(0.0));

    return q1 / (q1 + q2);
}

template <typename T>
inline T mis_balance(const T q1, const T q2, const T q3)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q3 >= T(0.0));
    assert(q1 + q2 + q3 > T(0.0));

    return q1 / (q1 + q2 + q3);
}

template <typename T>
inline T mis_power(const T q1, const T q2, const T beta)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q1 + q2 > T(0.0));

    assert(beta >= T(0.0));

    const T q1_pow = std::pow(q1, beta);
    const T q2_pow = std::pow(q2, beta);

    return q1_pow / (q1_pow + q2_pow);
}

template <typename T>
inline T mis_power(const T q1, const T q2, const T q3, const T beta)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q3 >= T(0.0));
    assert(q1 + q2 + q3 > T(0.0));

    assert(beta >= T(0.0));

    const T q1_pow = std::pow(q1, beta);
    const T q2_pow = std::pow(q2, beta);
    const T q3_pow = std::pow(q3, beta);

    return q1_pow / (q1_pow + q2_pow + q3_pow);
}

template <typename T>
inline T mis_power2(const T q1, const T q2)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q1 + q2 > T(0.0));

    const T q1_pow = q1 * q1;
    const T q2_pow = q2 * q2;

    return q1_pow / (q1_pow + q2_pow);
}

template <typename T>
inline T mis_power2(const T q1, const T q2, const T q3)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q3 >= T(0.0));
    assert(q1 + q2 + q3 > T(0.0));

    const T q1_pow = q1 * q1;
    const T q2_pow = q2 * q2;
    const T q3_pow = q3 * q3;

    return q1_pow / (q1_pow + q2_pow + q3_pow);
}

template <typename T>
inline T mis_cutoff(const T q1, const T q2, const T alpha)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q1 + q2 > T(0.0));

    assert(alpha >= T(0.0));
    assert(alpha <= T(1.0));

    const T cutoff = std::max(q1, q2) * alpha;

    if (q1 < cutoff)
         return T(0.0);
    else if (q2 < cutoff)
         return T(1.0);
    else return q1 / (q1 + q2);
}

template <typename T>
inline T mis_cutoff(const T q1, const T q2, const T q3, const T alpha)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q3 >= T(0.0));
    assert(q1 + q2 + q3 > T(0.0));

    assert(alpha >= T(0.0));
    assert(alpha <= T(1.0));

    const T cutoff = max(q1, q2, q3) * alpha;

    if (q1 < cutoff)
         return T(0.0);

    T den = q1;
    if (q2 >= cutoff) den += q2;
    if (q3 >= cutoff) den += q3;

    return q1 / den;
}

template <typename T>
inline T mis_maximum(const T q1, const T q2)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q1 + q2 > T(0.0));

    return q1 >= q2 ? T(1.0) : T(0.0);
}

template <typename T>
inline T mis_maximum(const T q1, const T q2, const T q3)
{
    assert(q1 >= T(0.0));
    assert(q2 >= T(0.0));
    assert(q3 >= T(0.0));
    assert(q1 + q2 + q3 > T(0.0));

    return q1 >= q2 && q1 >= q3 ? T(1.0) : T(0.0);
}

}       // namespace foundation

#endif  // !APPLESEED_FOUNDATION_MATH_MIS_H
