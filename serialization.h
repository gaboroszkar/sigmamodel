#ifndef SERIALIZATION
#define SERIALIZATION

#include <ostream>
#include <istream>

template <typename Value_type>
class Has_write
{

private:
    Has_write();
    template <typename T> static void sfinae(void (T::*)(std::ostream&) const);
    template <typename T> static auto test(void*) -> decltype(sfinae(&T::write), std::true_type());
    template <typename T> static std::false_type test(...);

public:
    typedef decltype(test<Value_type>(0)) Has;
};

template <typename Value_type, typename Has>
class Write
{

public:
    static void write(const Value_type& value, std::ostream& os)
    {
        os.write(reinterpret_cast<const char*>(&value), sizeof(Value_type));
    }

private:
    Write();
};

template <typename Value_type>
class Write<Value_type, std::true_type>
{

public:
    static void write(const Value_type& value, std::ostream& os)
    {
        value.write(os);
    }

private:
    Write();
};

template <typename Value_type>
void write(const Value_type& value, std::ostream& os)
{
    Write<Value_type, typename Has_write<Value_type>::Has>::write(value, os);
}

template <typename Value_type>
class Has_read
{

private:
    Has_read();
    template <typename T> static void sfinae(void (T::*)(std::istream&));
    template <typename T> static auto test(void*) -> decltype(sfinae(&T::read), std::true_type());
    template <typename T> static std::false_type test(...);

public:
    typedef decltype(test<Value_type>(0)) Has;
};

template <typename Value_type, typename Has>
class Read
{

public:
    static void read(Value_type& value, std::istream& is)
    {
        is.read(reinterpret_cast<char*>(&value), sizeof(Value_type));
    }

private:
    Read();
};

template <typename Value_type>
class Read<Value_type, std::true_type>
{

public:
    static void read(Value_type& value, std::istream& is)
    {
        value.read(is);
    }

private:
    Read();
};

template <typename Value_type>
void read(Value_type& value, std::istream& is)
{
    Read<Value_type, typename Has_read<Value_type>::Has>::read(value, is);
}

#endif
