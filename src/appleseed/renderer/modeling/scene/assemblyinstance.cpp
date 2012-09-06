
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2012 Francois Beaune, Jupiter Jazz Limited
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

// Interface header.
#include "assemblyinstance.h"

// appleseed.renderer headers.
#include "renderer/modeling/scene/assembly.h"
#include "renderer/modeling/scene/objectinstance.h"
#include "renderer/utility/bbox.h"

// Standard headers.
#include <string>

using namespace foundation;
using namespace std;

namespace renderer
{

//
// AssemblyInstance class implementation.
//

struct AssemblyInstance::Impl
{
    string m_assembly_name;
};

namespace
{
    const UniqueID g_class_uid = new_guid();
}

AssemblyInstance::AssemblyInstance(
    const char*         name,
    const ParamArray&   params,
    const char*         assembly_name)
  : Entity(g_class_uid, params)
  , impl(new Impl())
{
    set_name(name);

    impl->m_assembly_name = assembly_name;

    m_assembly = 0;
}

AssemblyInstance::~AssemblyInstance()
{
    delete impl;
}

void AssemblyInstance::release()
{
    delete this;
}

const char* AssemblyInstance::get_assembly_name() const
{
    return impl->m_assembly_name.c_str();
}

Assembly* AssemblyInstance::find_assembly() const
{
    const Entity* parent = get_parent();

    while (parent)
    {
        Assembly* assembly =
            static_cast<const Assembly*>(parent)
                ->assemblies().get_by_name(impl->m_assembly_name.c_str());

        if (assembly)
            return assembly;

        parent = parent->get_parent();
    }

    return 0;
}

GAABB3 AssemblyInstance::compute_local_bbox() const
{
    // In many places, we need the parent-space bounding box of an assembly instance
    // before input binding is performed, i.e. before the instantiated assembly is
    // bound to the instance. Therefore we manually look the assembly up through the
    // assembly hierarchy instead of simply using m_assembly.

    const Assembly* assembly = find_assembly();

    if (assembly == 0)
        return GAABB3::invalid();

    const ObjectInstanceContainer& object_instances = assembly->object_instances();

    GAABB3 bbox =
        renderer::compute_parent_bbox<GAABB3>(
            object_instances.begin(),
            object_instances.end());

    const AssemblyInstanceContainer& assembly_instances = assembly->assembly_instances();

    bbox.insert(
        renderer::compute_parent_bbox<GAABB3>(
            assembly_instances.begin(),
            assembly_instances.end()));

    return bbox;
}

GAABB3 AssemblyInstance::compute_parent_bbox() const
{
    return m_transform_sequence.to_parent(compute_local_bbox());
}

void AssemblyInstance::unbind_assembly()
{
    m_assembly = 0;
}

void AssemblyInstance::bind_assembly(const AssemblyContainer& assemblies)
{
    if (m_assembly == 0)
        m_assembly = assemblies.get_by_name(impl->m_assembly_name.c_str());
}

void AssemblyInstance::check_assembly() const
{
    if (m_assembly == 0)
        throw ExceptionUnknownEntity(impl->m_assembly_name.c_str());
}

bool AssemblyInstance::on_frame_begin(const Project& project)
{
    if (!m_transform_sequence.prepare())
    {
        RENDERER_LOG_ERROR("assembly instance \"%s\" has one or more invalid transforms.", get_name());
        return false;
    }

    return true;
}

void AssemblyInstance::on_frame_end(const Project& project)
{
}


//
// AssemblyInstanceFactory class implementation.
//

auto_release_ptr<AssemblyInstance> AssemblyInstanceFactory::create(
    const char*         name,
    const ParamArray&   params,
    const char*         assembly_name)
{
    return
        auto_release_ptr<AssemblyInstance>(
            new AssemblyInstance(
                name,
                params,
                assembly_name));
}

}   // namespace renderer
