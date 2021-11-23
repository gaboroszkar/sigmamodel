#ifndef CONFIGURATION_CLUSTER_H
#define CONFIGURATION_CLUSTER_H

#include "configuration.h"
#include "point.h"
#include <list>

template <typename Value_type, typename Connect_functor>
class Configuration_cluster : public virtual Configuration<Value_type>
{

protected:
    static std::set<std::shared_ptr<Point<Value_type>>> build_cluster(const std::shared_ptr<Point<Value_type>>& starting_point, Connect_functor& connect_functor)
    {
        std::set<std::shared_ptr<Point<Value_type>>> cluster;
        build_cluster(starting_point, connect_functor, cluster);
        return cluster;
    }

private:
#ifdef TAIL_RECURSIVE_CLUSTER_BUILDING
    static void build_cluster(std::shared_ptr<Point<Value_type>> point, Connect_functor& connect_functor, std::set<std::shared_ptr<Point<Value_type>>>& cluster)
    {
        std::list<std::pair<std::shared_ptr<Point<Value_type>>, std::shared_ptr<Point<Value_type>>>> new_edges;
        new_edges.push_back(std::pair<std::shared_ptr<Point<Value_type>>, std::shared_ptr<Point<Value_type>>>(nullptr, point));
        build_cluster(new_edges, connect_functor, cluster);
    }
    static void build_cluster(const std::list<std::pair<std::shared_ptr<Point<Value_type>>, std::shared_ptr<Point<Value_type>>>>& edges, Connect_functor& connect_functor, std::set<std::shared_ptr<Point<Value_type>>>& cluster)
    {
        std::list<std::pair<std::shared_ptr<Point<Value_type>>, std::shared_ptr<Point<Value_type>>>> new_edges;
        for (auto edge : edges)
            if (edge.first == nullptr || connect_functor(edge.first->value(), edge.second->value()))
            {
                cluster.insert(edge.second);
                for (auto neighbor : edge.second->neighbors())
                    if (cluster.count(neighbor) == 0)
                        new_edges.push_back(std::pair<std::shared_ptr<Point<Value_type>>, std::shared_ptr<Point<Value_type>>>(edge.second, neighbor));
            }

        if (new_edges.empty())
            return;
        return build_cluster(new_edges, connect_functor, cluster);
    }
#else
    static void build_cluster(const std::shared_ptr<Point<Value_type>>& point, Connect_functor& connect_functor, std::set<std::shared_ptr<Point<Value_type>>>& cluster)
    {
        cluster.insert(point);
        for (auto neighbor : point->neighbors())
            if (cluster.count(neighbor) == 0)
                if (connect_functor(point->value(), neighbor->value()))
                    build_cluster(neighbor, connect_functor, cluster);
    }
#endif
};

#endif
