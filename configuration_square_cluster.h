#ifndef CONFIGURATION_SQUARE_CLUSTER_H
#define CONFIGURATION_SQUARE_CLUSTER_H

#include "configuration_square.h"
#include "configuration_cluster.h"
#include "point.h"
#include <random>

template <typename Value_type, typename Connect_functor>
class Configuration_square_cluster : public virtual Configuration_square<Value_type>, public Configuration_cluster<Value_type, Connect_functor>
{

public:
    Configuration_square_cluster() : Configuration_square<Value_type>() {}
    Configuration_square_cluster(int linear_size) : Configuration_square<Value_type>(linear_size) {}
};

#endif
