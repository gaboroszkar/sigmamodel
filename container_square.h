#ifndef CONTAINER_SQUARE_H
#define CONTAINER_SQUARE_H

#include <vector>

class Container_square_not_initialized : public std::exception
{
public:
    const char* what () const throw ()
    {
        return "Not initialized";
    }
};

template <typename Value_type>
class Container_square
{

public:
    Container_square() : m_linear_size{0} {}
    Container_square(int linear_size) : m_values(linear_size * linear_size), m_linear_size{linear_size} {}

    inline Value_type& value(int x, int y)
    {
        if (linear_size() == 0)
            throw Container_square_not_initialized();
        return m_values.at(index(x, y));
    }

    inline const Value_type& value(int x, int y) const
    {
        if (linear_size() == 0)
            throw Container_square_not_initialized();
        return m_values.at(index(x, y));
    }

    inline int linear_size() const
    {
        return m_linear_size;
    }

protected:
    void clear_and_resize(int linear_size)
    {
        m_values.clear();
        m_values.resize(linear_size * linear_size);
        m_linear_size = linear_size;
    }

private:
    inline int normalize(int x) const
    {
        if (x < 0)
            x += linear_size() * (1 - (x / linear_size()));
        return x % linear_size();
    }

    inline int index(int x, int y) const
    {
        return normalize(x) * linear_size() + normalize(y);
    }

    std::vector<Value_type> m_values;
    int m_linear_size;
};

#endif
