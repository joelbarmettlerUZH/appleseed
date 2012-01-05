
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2011 Francois Beaune
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

#ifndef APPLESEED_FOUNDATION_UTILITY_BENCHMARK_BENCHMARKSUITEREPOSITORY_H
#define APPLESEED_FOUNDATION_UTILITY_BENCHMARK_BENCHMARKSUITEREPOSITORY_H

// appleseed.foundation headers.
#include "foundation/core/concepts/singleton.h"

// Forward declarations.
namespace foundation    { class BenchmarkResult; }
namespace foundation    { class BenchmarkSuite; }
namespace foundation    { class IFilter; }

//
// On Windows, define FOUNDATIONDLL to __declspec(dllexport) when building the DLL
// and to __declspec(dllimport) when building an application using the DLL.
// Other platforms don't use this export mechanism and the symbol FOUNDATIONDLL is
// defined to evaluate to nothing.
//

#ifndef FOUNDATIONDLL
#ifdef _WIN32
#ifdef APPLESEED_FOUNDATION_EXPORTS
#define FOUNDATIONDLL __declspec(dllexport)
#else
#define FOUNDATIONDLL __declspec(dllimport)
#endif
#else
#define FOUNDATIONDLL
#endif
#endif

namespace foundation
{

//
// The (unique) benchmark suite repository, as a collection of benchmark suites.
//

class FOUNDATIONDLL BenchmarkSuiteRepository
  : public Singleton<BenchmarkSuiteRepository>
{
  public:
    // Register a benchmark suite.
    void register_suite(BenchmarkSuite* suite);

    // Run all the registered benchmark suites.
    void run(BenchmarkResult& result) const;

    // Run those benchmark suites whose name pass a given filter.
    void run(
        const IFilter&      filter,
        BenchmarkResult&    result) const;

  private:
    friend class Singleton<BenchmarkSuiteRepository>;

    struct Impl;
    Impl* impl;

    // Constructor.
    BenchmarkSuiteRepository();

    // Destructor.
    ~BenchmarkSuiteRepository();
};

}       // namespace foundation

#endif  // !APPLESEED_FOUNDATION_UTILITY_BENCHMARK_BENCHMARKSUITEREPOSITORY_H
