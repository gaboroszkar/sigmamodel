#ifndef POINT_H
#define POINT_H

#include <set>
#include <memory>
#include "serialization.h"

template <typename Value_type>
class Configuration;

template <typename Value_type>
class Point
{

    friend class Configuration<Value_type>;

public:
    Point() {}
    Point(const Value_type& value) : m_value{value} {}

    inline Value_type& value()
    {
        return m_value;
    }

    inline const std::set<std::shared_ptr<Point>>& neighbors() const
    {
        return m_neighbors;
    }

    inline void write(std::ostream& os) const
    {
        ::write(m_value, os);
    }

    inline void read(std::istream& is)
    {
        ::read(m_value, is);
    }

private:
    inline virtual void set_neighbor(const std::shared_ptr<Point>& neighbor)
    {
        m_neighbors.insert(neighbor);
    }

    Value_type m_value;
    std::set<std::shared_ptr<Point>> m_neighbors;
};

#endif
