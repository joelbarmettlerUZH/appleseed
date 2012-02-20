
//
// This source file is part of appleseed.
// Visit http://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2010-2012 Francois Beaune, Jupiter Jazz
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
#include "assemblytree.h"

// appleseed.renderer headers.
#include "renderer/kernel/shading/shadingpoint.h"
#include "renderer/modeling/object/iregion.h"
#include "renderer/modeling/object/object.h"
#include "renderer/modeling/object/regionkit.h"
#include "renderer/modeling/scene/assembly.h"
#include "renderer/modeling/scene/assemblyinstance.h"
#include "renderer/utility/bbox.h"

// appleseed.foundation headers.
#include "foundation/math/intersection.h"
#include "foundation/math/permutation.h"
#include "foundation/platform/system.h"
#include "foundation/platform/timer.h"
#include "foundation/utility/string.h"

// Standard headers.
#include <algorithm>
#include <utility>

using namespace foundation;
using namespace std;

namespace renderer
{

//
// Specializations.
//

typedef bvh::SAHPartitioner<
    AssemblyTree
> AssemblyTreePartitioner;

typedef bvh::Builder<
    AssemblyTree,
    AssemblyTreePartitioner
> AssemblyTreeBuilder;

typedef bvh::TreeStatistics<
    AssemblyTree,
    AssemblyTreeBuilder
> AssemblyTreeStatistics;


//
// AssemblyTree class implementation.
//

AssemblyTree::AssemblyTree(const Scene& scene)
  : TreeType(AlignedAllocator<void>(System::get_l1_data_cache_line_size()))
  , m_scene(scene)
{
    update();
}

AssemblyTree::~AssemblyTree()
{
    // Log a progress message.
    RENDERER_LOG_INFO("deleting assembly tree...");

    // Delete region trees.
    for (each<RegionTreeContainer> i = m_region_trees; i; ++i)
        delete i->second;
    m_region_trees.clear();

    // Delete triangle trees.
    for (each<TriangleTreeContainer> i = m_triangle_trees; i; ++i)
        delete i->second;
    m_triangle_trees.clear();
}

void AssemblyTree::update()
{
    rebuild_assembly_tree();
    update_child_trees();
}

size_t AssemblyTree::get_memory_size() const
{
    return
          TreeType::get_memory_size()
        - sizeof(*static_cast<const TreeType*>(this))
        + sizeof(*this)
        + m_assembly_instances.capacity() * sizeof(UniqueID);
}

namespace
{
    // Collect all assemblies from a given scene.
    void collect_assemblies(const Scene& scene, vector<UniqueID>& assemblies)
    {
        assert(assemblies.empty());

        assemblies.reserve(scene.assembly_instances().size());

        for (const_each<AssemblyInstanceContainer> i = scene.assembly_instances(); i; ++i)
            assemblies.push_back(i->get_assembly_uid());

        sort(assemblies.begin(), assemblies.end());

        const vector<UniqueID>::iterator new_end =
            unique(assemblies.begin(), assemblies.end());

        assemblies.erase(new_end, assemblies.end());
    }

    // Collect all regions of all objects from a given assembly.
    void collect_regions(const Assembly& assembly, RegionInfoVector& regions)
    {
        assert(regions.empty());

        const ObjectInstanceContainer& object_instances = assembly.object_instances();
        const size_t object_instance_count = object_instances.size();

        // Collect all regions of all object instances of this assembly.
        for (size_t obj_inst_index = 0; obj_inst_index < object_instance_count; ++obj_inst_index)
        {
            // Retrieve the object instance and its transformation.
            const ObjectInstance* object_instance = object_instances.get_by_index(obj_inst_index);
            assert(object_instance);
            const Transformd& transform = object_instance->get_transform();

            // Retrieve the object.
            Object& object = object_instance->get_object();

            // Retrieve the region kit of the object.
            Access<RegionKit> region_kit(&object.get_region_kit());

            // Collect all regions of the object.
            for (size_t region_index = 0; region_index < region_kit->size(); ++region_index)
            {
                // Retrieve the region.
                const IRegion* region = (*region_kit)[region_index];

                // Compute the assembly space bounding box of the region.
                const GAABB3 region_bbox =
                    transform.transform_to_parent(region->get_local_bbox());

                regions.push_back(
                    RegionInfo(
                        obj_inst_index,
                        region_index,
                        region_bbox));
            }
        }
    }
}

Lazy<TriangleTree>* AssemblyTree::create_triangle_tree(const Assembly& assembly) const
{
    // Compute the assembly space bounding box of the assembly.
    const GAABB3 assembly_bbox =
        get_parent_bbox<GAABB3>(
            assembly.object_instances().begin(),
            assembly.object_instances().end());

    RegionInfoVector regions;
    collect_regions(assembly, regions);

    auto_ptr<ILazyFactory<TriangleTree> > triangle_tree_factory(
        new TriangleTreeFactory(
            TriangleTree::Arguments(
                assembly.get_uid(),
                assembly_bbox,
                assembly,
                regions)));

    return new Lazy<TriangleTree>(triangle_tree_factory);
}

Lazy<RegionTree>* AssemblyTree::create_region_tree(const Assembly& assembly) const
{
    auto_ptr<ILazyFactory<RegionTree> > region_tree_factory(
        new RegionTreeFactory(
            RegionTree::Arguments(
                assembly.get_uid(),
                assembly)));

    return new Lazy<RegionTree>(region_tree_factory);
}

void AssemblyTree::rebuild_assembly_tree()
{
    // Clear the current tree.
    clear();
    m_assembly_instances.clear();

    // Insert all assembly instances of the scene into the tree.
    for (const_each<AssemblyInstanceContainer> i = m_scene.assembly_instances(); i; ++i)
    {
        // Retrieve the assembly instance.
        const AssemblyInstance& assembly_instance = *i;

        // Retrieve the assembly.
        const Assembly& assembly = assembly_instance.get_assembly();

        // Skip empty assemblies.
        if (assembly.object_instances().empty())
            continue;

        // Insert the bounding box of the assembly instance into the tree's root leaf.
        insert(AABB3d(assembly_instance.compute_parent_bbox()));

        // Store the assembly instance.
        m_assembly_instances.push_back(assembly_instance.get_uid());
    }

    // Log a progress message.
    RENDERER_LOG_INFO(
        "building assembly tree (%s %s)...",
        pretty_int(size()).c_str(),
        plural(size(), "assembly instance").c_str());

    // Build the assembly tree.
    AssemblyTreePartitioner partitioner(1);
    AssemblyTreeBuilder builder;
    builder.build<DefaultWallclockTimer>(*this, partitioner);

    if (!m_assembly_instances.empty())
    {
        const vector<size_t>& ordering = partitioner.get_item_ordering();

        assert(m_assembly_instances.size() == ordering.size());

        // Reorder the assembly instances according to the tree ordering.
        vector<UniqueID> temp_assembly_instances(ordering.size());
        small_item_reorder(
            &m_assembly_instances[0],
            &temp_assembly_instances[0],
            &ordering[0],
            ordering.size());

        // Store assembly instances in the tree leaves whenever possible.
        store_assembly_instances_in_leaves();
    }

    // Collect and print assembly tree statistics.
    AssemblyTreeStatistics tree_stats(*this, builder);
    RENDERER_LOG_DEBUG("assembly tree statistics:");
    tree_stats.print(global_logger());
}

void AssemblyTree::store_assembly_instances_in_leaves()
{
    size_t leaf_count = 0;
    size_t fat_leaf_count = 0;

    const size_t node_count = m_nodes.size();

    for (size_t i = 0; i < node_count; ++i)
    {
        NodeType& node = m_nodes[i];

        if (node.is_leaf())
        {
            ++leaf_count;

            const size_t item_count = node.get_item_count();

            if (item_count <= NodeType::MaxUserDataSize / sizeof(UniqueID))
            {
                ++fat_leaf_count;

                const size_t item_begin = node.get_item_index();
                UniqueID* user_data = &node.get_user_data<UniqueID>();

                for (size_t j = 0; j < item_count; ++j)
                    user_data[j] = m_assembly_instances[item_begin + j];
            }
        }
    }

    RENDERER_LOG_DEBUG(
        "fat assembly tree leaves: %s",
        pretty_percent(fat_leaf_count, leaf_count).c_str());
}

void AssemblyTree::update_child_trees()
{
    // Collect all assemblies in the scene.
    vector<UniqueID> assemblies;
    collect_assemblies(m_scene, assemblies);

    // Create or update the child tree of each assembly.
    for (const_each<vector<UniqueID> > i = assemblies; i; ++i)
    {
        // Retrieve the assembly.
        const UniqueID assembly_uid = *i;
        const Assembly& assembly = *m_scene.assemblies().get_by_uid(assembly_uid);

        // Retrieve the current version ID of the assembly.
        const VersionID current_version_id = assembly.get_version_id();

        // Retrieve the stored version ID of the assembly.
        const AssemblyVersionMap::const_iterator stored_version_it =
            m_assembly_versions.find(assembly_uid);

        if (stored_version_it == m_assembly_versions.end())
        {
            // No tree for this assembly yet, create one.
            if (assembly.is_flushable())
            {
                m_region_trees.insert(
                    make_pair(assembly_uid, create_region_tree(assembly)));
            }
            else
            {
                m_triangle_trees.insert(
                    make_pair(assembly_uid, create_triangle_tree(assembly)));
            }
        }
        else if (current_version_id != stored_version_it->second)
        {
            // The tree corresponding to this assembly is out-of-date.
            if (assembly.is_flushable())
            {
                const RegionTreeContainer::iterator region_tree_it =
                    m_region_trees.find(assembly_uid);
                delete region_tree_it->second;
                region_tree_it->second = create_region_tree(assembly);
            }
            else
            {
                const TriangleTreeContainer::iterator triangle_tree_it =
                    m_triangle_trees.find(assembly_uid);
                delete triangle_tree_it->second;
                triangle_tree_it->second = create_triangle_tree(assembly);
            }
        }

        // Update the stored version ID of the assembly.
        m_assembly_versions[assembly_uid] = current_version_id;
    }
}


//
// Utility function to transform a ray to the space of an assembly instance.
//

namespace
{
    void transform_ray_to_assembly_instance_space(
        const AssemblyInstance&         assembly_instance,
        const ShadingPoint*             parent_shading_point,
        const ShadingRay::RayType&      input_ray,
        ShadingRay::RayType&            output_ray)
    {
        // Retrieve the transformation of the assembly instance.
        const Transformd& transform = assembly_instance.get_transform();

        // Transform the ray direction.
        output_ray.m_dir = transform.transform_vector_to_local(input_ray.m_dir);

        if (parent_shading_point &&
            parent_shading_point->get_assembly_instance_uid() == assembly_instance.get_uid())
        {
            // The caller provided the previous intersection, and we are about
            // to intersect the assembly instance that contains the previous
            // intersection. Use the properly offset intersection point as the
            // origin of the child ray.
            output_ray.m_org = parent_shading_point->get_offset_point(output_ray.m_dir);
        }
        else
        {
            // The caller didn't provide the previous intersection, or we are
            // about to intersect an assembly instance that does not contain
            // the previous intersection: proceed normally.
            output_ray.m_org = transform.transform_point_to_local(input_ray.m_org);
        }
    }
}


//
// AssemblyLeafVisitor class implementation.
//

bool AssemblyLeafVisitor::visit(
    const AssemblyTree::AABBVector&     bboxes,
    const AssemblyTree::NodeType&       node,
    const ShadingRay::RayType&          ray,
    const ShadingRay::RayInfoType&      ray_info,
    double&                             distance
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
    , bvh::TraversalStatistics&         stats
#endif
    )
{
    // Retrieve the assembly instances for this leaf.
    const size_t assembly_instance_index = node.get_item_index();
    const size_t assembly_instance_count = node.get_item_count();
    const UniqueID* assembly_instance_uids =
        assembly_instance_count <= AssemblyTree::NodeType::MaxUserDataSize / sizeof(UniqueID)
            ? &node.get_user_data<UniqueID>()                           // items are stored in the leaf node
            : &m_tree.m_assembly_instances[assembly_instance_index];    // items are stored in the tree

    for (size_t i = 0; i < assembly_instance_count; ++i)
    {
        // Retrieve the assembly instance.
        const UniqueID assembly_instance_uid = assembly_instance_uids[i];
        const AssemblyInstance& assembly_instance =
            *m_tree.m_scene.assembly_instances().get_by_uid(assembly_instance_uid);

        // Transform the ray to assembly instance space.
        ShadingPoint local_shading_point;
        transform_ray_to_assembly_instance_space(
            assembly_instance,
            m_parent_shading_point,
            ray,
            local_shading_point.m_ray);
        local_shading_point.m_ray.m_tmin = ray.m_tmin;
        local_shading_point.m_ray.m_tmax = ray.m_tmax;
        local_shading_point.m_ray.m_time = m_shading_point.m_ray.m_time;
        local_shading_point.m_ray.m_flags = m_shading_point.m_ray.m_flags;
        const RayInfo3d local_ray_info(local_shading_point.m_ray);

        FOUNDATION_BVH_TRAVERSAL_STATS(stats.m_intersected_items.insert(1));

        if (assembly_instance.get_assembly().is_flushable())
        {
            // Retrieve the region tree of this assembly.
            const RegionTree& region_tree =
                *m_region_tree_cache.access(
                    assembly_instance.get_assembly_uid(),
                    m_tree.m_region_trees);

            // Check the intersection between the ray and the region tree.
            RegionLeafVisitor visitor(
                local_shading_point,
                m_triangle_tree_cache
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
                , m_triangle_tree_stats
#endif
                );
            RegionLeafIntersector intersector;
            intersector.intersect(
                region_tree,
                local_shading_point.m_ray,
                local_ray_info,
                visitor);
        }
        else
        {
            // Retrieve the triangle tree of this assembly.
            const TriangleTree* triangle_tree =
                m_triangle_tree_cache.access(
                    assembly_instance.get_assembly_uid(),
                    m_tree.m_triangle_trees);

            if (triangle_tree)
            {
                // Check the intersection between the ray and the triangle tree.
                TriangleTreeIntersector intersector;
                TriangleLeafVisitor visitor(*triangle_tree, local_shading_point);
                intersector.intersect(
                    *triangle_tree,
                    local_shading_point.m_ray,
                    local_ray_info,
                    visitor
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
                    , m_triangle_tree_stats
#endif
                    );
                visitor.read_hit_triangle_data();
            }
        }

        // Keep track of the closest hit.
        if (local_shading_point.m_hit && local_shading_point.m_ray.m_tmax < m_shading_point.m_ray.m_tmax)
        {
            m_shading_point.m_ray.m_tmax = local_shading_point.m_ray.m_tmax;
            m_shading_point.m_hit = true;
            m_shading_point.m_bary = local_shading_point.m_bary;
            m_shading_point.m_asm_instance_uid = assembly_instance_uid;
            m_shading_point.m_object_instance_index = local_shading_point.m_object_instance_index;
            m_shading_point.m_region_index = local_shading_point.m_region_index;
            m_shading_point.m_triangle_index = local_shading_point.m_triangle_index;
            m_shading_point.m_triangle_support_plane = local_shading_point.m_triangle_support_plane;
        }
    }

    // Continue traversal.
    distance = m_shading_point.m_ray.m_tmax;
    return true;
}


//
// AssemblyLeafProbeVisitor class implementation.
//

bool AssemblyLeafProbeVisitor::visit(
    const AssemblyTree::AABBVector&     bboxes,
    const AssemblyTree::NodeType&       node,
    const ShadingRay::RayType&          ray,
    const ShadingRay::RayInfoType&      ray_info,
    double&                             distance
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
    , bvh::TraversalStatistics&         stats
#endif
    )
{
    // Retrieve the assembly instances for this leaf.
    const size_t assembly_instance_count = node.get_item_count();
    const UniqueID* assembly_instance_uids =
        assembly_instance_count <= AssemblyTree::NodeType::MaxUserDataSize / sizeof(UniqueID)
            ? &node.get_user_data<UniqueID>()                           // items are stored in the leaf node
            : &m_tree.m_assembly_instances[node.get_item_index()];      // items are stored in the tree

    for (size_t i = 0; i < assembly_instance_count; ++i)
    {
        // Retrieve the assembly instance.
        const UniqueID assembly_instance_uid = assembly_instance_uids[i];
        const AssemblyInstance& assembly_instance =
            *m_tree.m_scene.assembly_instances().get_by_uid(assembly_instance_uid);

        // Transform the ray to assembly instance space.
        ShadingRay::RayType local_ray;
        transform_ray_to_assembly_instance_space(
            assembly_instance,
            m_parent_shading_point,
            ray,
            local_ray);
        local_ray.m_tmin = ray.m_tmin;
        local_ray.m_tmax = ray.m_tmax;
        const RayInfo3d local_ray_info(local_ray);

        FOUNDATION_BVH_TRAVERSAL_STATS(stats.m_intersected_items.insert(1));

        if (assembly_instance.get_assembly().is_flushable())
        {
            // Retrieve the region tree of this assembly.
            const RegionTree& region_tree =
                *m_region_tree_cache.access(
                    assembly_instance.get_assembly_uid(),
                    m_tree.m_region_trees);

            // Check the intersection between the ray and the region tree.
            RegionLeafProbeVisitor visitor(
                m_triangle_tree_cache
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
                , m_triangle_tree_stats
#endif
                );
            RegionLeafProbeIntersector intersector;
            intersector.intersect(
                region_tree,
                local_ray,
                local_ray_info,
                visitor);
        
            // Terminate traversal if there was a hit.
            if (visitor.hit())
            {
                m_hit = true;
                return false;
            }
        }
        else
        {
            // Retrieve the triangle tree of this leaf.
            const TriangleTree* triangle_tree =
                m_triangle_tree_cache.access(
                    assembly_instance.get_assembly_uid(),
                    m_tree.m_triangle_trees);

            if (triangle_tree)
            {
                // Check the intersection between the ray and the triangle tree.
                TriangleTreeProbeIntersector intersector;
                TriangleLeafProbeVisitor visitor(*triangle_tree);
                intersector.intersect(
                    *triangle_tree,
                    local_ray,
                    local_ray_info,
                    visitor
#ifdef FOUNDATION_BVH_ENABLE_TRAVERSAL_STATS
                    , m_triangle_tree_stats
#endif
                    );

                // Terminate traversal if there was a hit.
                if (visitor.hit())
                {
                    m_hit = true;
                    return false;
                }
            }
        }
    }

    // Continue traversal.
    distance = ray.m_tmax;
    return true;
}

}   // namespace renderer
