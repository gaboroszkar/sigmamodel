#ifndef VECTOR_H
#define VECTOR_H

#include <array>
#include <cmath>
#include "serialization.h"

template <typename Value_type, int dimension>
class Vector
{

public:
    Vector() {}
    Vector(const std::array<Value_type, dimension>& values) : m_values(values) {}

    inline Value_type operator[](int i) const
    {
        return m_values[i];
    }

    inline Value_type& operator[](int i)
    {
        return m_values[i];
    }

    void write(std::ostream& os) const
    {
        for (int i = 0; i != dimension; ++i)
            ::write(m_values[i], os);
    }

    void read(std::istream& is)
    {
        for (int i = 0; i != dimension; ++i)
            ::read(m_values[i], is);
    }

protected:
    std::array<Value_type, dimension> m_values;
};

template <typename Value_type, int dimension>
Vector<Value_type, dimension> operator *(Value_type scalar, Vector<Value_type, dimension> vector)
{
    for (int i = 0; i != dimension; ++i)
        vector[i] = scalar * vector[i];
    return vector;
}

template <typename Value_type, int dimension>
Vector<Value_type, dimension> operator /(Vector<Value_type, dimension> vector, Value_type scalar)
{
    for (int i = 0; i != dimension; ++i)
        vector[i] = vector[i] / scalar;
    return vector;
}

template <typename Value_type, int dimension>
Value_type operator *(const Vector<Value_type, dimension>& vector_0, const Vector<Value_type, dimension>& vector_1)
{
    Value_type scalar_result = 0.0;
    for (int i = 0; i != dimension; ++i)
        scalar_result = scalar_result + vector_0[i] * vector_1[i];
    return scalar_result;
}

template <typename Value_type, int dimension>
Vector<Value_type, dimension> operator +(Vector<Value_type, dimension> vector_0, const Vector<Value_type, dimension>& vector_1)
{
    for (int i = 0; i != dimension; ++i)
        vector_0[i] = vector_0[i] + vector_1[i];
    return vector_0;
}

template <typename Value_type, int dimension>
Vector<Value_type, dimension> operator -(Vector<Value_type, dimension> vector_0, const Vector<Value_type, dimension>& vector_1)
{
    for (int i = 0; i != dimension; ++i)
        vector_0[i] = vector_0[i] - vector_1[i];
    return vector_0;
}

template <typename Value_type>
inline Vector<Value_type, 3> vector_cross(const Vector<Value_type, 3>& vector_0, const Vector<Value_type, 3>& vector_1)
{
    return Vector<Value_type, 3>({vector_0[1] * vector_1[2] - vector_0[2] * vector_1[1],
            vector_0[2] * vector_1[0] - vector_0[0] * vector_1[2],
            vector_0[0] * vector_1[1] - vector_0[1] * vector_1[0]});
}

template <typename Value_type, int dimension>
inline Value_type vector_length(const Vector<Value_type, dimension>& vector)
{
    return std::sqrt(vector * vector);
}

template <typename Value_type, int dimension>
inline Vector<Value_type, dimension> vector_normalize(const Vector<Value_type, dimension>& vector)
{
    return vector / vector_length(vector);
}

#endif
