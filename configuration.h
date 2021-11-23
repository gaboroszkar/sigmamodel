#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "point.h"
#include <memory>
#include <istream>
#include <ostream>

template <typename Value_type>
class Configuration
{

public:
    virtual void write(std::ostream& os) const = 0;
    virtual void read(std::istream& is) = 0;
    virtual void update() = 0;

protected:
    inline static void set_neighbor(const std::shared_ptr<Point<Value_type>>& point_0, const std::shared_ptr<Point<Value_type>>& point_1)
    {
        point_0->set_neighbor(point_1);
        point_1->set_neighbor(point_0);
    }
};

#endif
