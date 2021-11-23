#ifndef CONFIGURATION_SQUARE_H
#define CONFIGURATION_SQUARE_H

#include "container_square.h"
#include "point.h"
#include "configuration.h"

template <typename Value_type>
class Configuration_square : public Container_square<std::shared_ptr<Point<Value_type>>>, public virtual Configuration<Value_type>
{

public:
    Configuration_square() : Container_square<std::shared_ptr<Point<Value_type>>>() {}

    Configuration_square(int linear_size) : Container_square<std::shared_ptr<Point<Value_type>>>(linear_size)
    {
        initialize();
    }

    virtual void write(std::ostream& os) const
    {
        ::write(this->linear_size(), os);
        for (int x = 0; x != this->linear_size(); ++x)
            for (int y = 0; y != this->linear_size(); ++y)
                ::write(this->value(x, y)->value(), os);
    }

    virtual void read(std::istream& is)
    {
        int linear_size;
        ::read(linear_size, is);
        this->clear_and_resize(linear_size);
        initialize();
        for (int x = 0; x != this->linear_size(); ++x)
            for (int y = 0; y != this->linear_size(); ++y)
                ::read(this->value(x, y)->value(), is);
    }

private:
    void initialize()
    {
        for (int x = 0; x != this->linear_size(); ++x)
            for (int y = 0; y != this->linear_size(); ++y)
                this->value(x, y) = std::make_shared<Point<Value_type>>();

        for (int x = 0; x != this->linear_size(); ++x)
        {
            for (int y = 0; y != this->linear_size(); ++y)
            {
                this->set_neighbor(this->value(x, y), this->value(x + 1, y));
                this->set_neighbor(this->value(x, y), this->value(x, y + 1));
            }
        }
    }
};

#endif
