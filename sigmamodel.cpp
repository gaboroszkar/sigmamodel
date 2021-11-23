#include "configuration_square_cluster.h"
#include "vector.h"
#include <random>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cstring>

constexpr double pi = 3.14159265358979323846;
constexpr double epsilon = 0.00000001;

class Vector3_default : public Vector<double, 3>
{

public:
    Vector3_default() : Vector<double, 3>({0.0, 0.0, 1.0}) {}

    Vector3_default(const Vector<double, 3>& v) : Vector<double, 3>(v) {}
    Vector3_default(const std::array<double, 3>& v) : Vector<double, 3>(v) {}
};

template <typename Random_generator>
Vector3_default random_vector_in_one_sphere(Random_generator& random_generator)
{
    std::uniform_real_distribution<> distribution(-1.0, 1.0);
    Vector3_default v({distribution(random_generator), distribution(random_generator), distribution(random_generator)});
    if (vector_length(v) > 1.0)
        return random_vector_in_one_sphere(random_generator);
    return v;
}

template <typename Random_generator>
Vector3_default random_normalized_vector(Random_generator& random_generator)
{
    Vector3_default v = random_vector_in_one_sphere(random_generator);
    if (vector_length(v) < epsilon)
        return random_normalized_vector(random_generator);
    return vector_normalize(v);
}

template <typename Random_generator>
std::array<Vector3_default, 3> random_normalized_perpendicular_vectors(Random_generator& random_generator)
{
    Vector3_default r_0 = random_normalized_vector(random_generator);

    Vector3_default r_n = random_normalized_vector(random_generator);
    Vector3_default r_1 = vector_normalize(vector_cross(r_0, r_n));

    Vector3_default r_2 = vector_normalize(vector_cross(r_0, r_1));
    std::uniform_int_distribution<> distribution(0, 1);
    if (distribution(random_generator))
        r_2 = -1.0 * r_2;

    return {r_0, r_1, r_2};
}

template <typename Random_generator>
class Connect
{

public:
    Connect(double beta, Vector3_default r, Random_generator& random_generator) : m_beta{beta}, m_r{r}, m_random_generator(random_generator) {}

    bool operator()(Vector3_default point_0, Vector3_default point_1)
    {
        double x = (point_0 * m_r) * (point_1 * m_r);
        std::uniform_real_distribution<> distribution(0.0, 1.0);
        if (x < 0.0)
            return false;
        return distribution(m_random_generator) < 1.0 - std::exp(-2.0 * m_beta * x);
    }

private:
    const double m_beta;
    Vector3_default m_r;
    Random_generator& m_random_generator;
};

class Configuration_square_parameters : public virtual Configuration_square<Vector3_default>
{
public:
    Configuration_square_parameters() : Configuration_square<Vector3_default>() {}
    Configuration_square_parameters(int linear_size) : Configuration_square<Vector3_default>(linear_size) {}

    virtual double beta() const = 0;
    virtual Vector3_default external_field() const = 0;
    virtual std::string algorithm_string() const = 0;
};

template <typename Random_generator>
class Configuration_metropolis : public Configuration_square_parameters
{

public:
    Configuration_metropolis(Random_generator& random_generator)
        : Configuration_square<Vector3_default>(),
          Configuration_square_parameters(),
          m_random_generator(random_generator) {}
    Configuration_metropolis(int linear_size, double beta, double external_field_length, Random_generator& random_generator)
        : Configuration_square<Vector3_default>(linear_size),
          Configuration_square_parameters(linear_size),
          m_beta(beta),
          m_external_field(Vector3_default({0.0, 0.0, external_field_length})),
          m_random_generator(random_generator) {}

    void write(std::ostream& os) const
    {
        Configuration_square<Vector3_default>::write(os);
        ::write(m_beta, os);
        ::write(m_external_field, os);
    }

    virtual void read(std::istream& is)
    {
        Configuration_square<Vector3_default>::read(is);
        ::read(m_beta, is);
        ::read(m_external_field, is);
    }

    void update()
    {
        std::uniform_int_distribution<> distribution_linear_size(0, this->linear_size() - 1);
        std::shared_ptr<Point<Vector3_default>> p = this->value(distribution_linear_size(m_random_generator), distribution_linear_size(m_random_generator));

        Vector3_default r_0 = p->value();
        double local_action_0 = get_local_action(*p);

        Vector3_default r_1 = random_normalized_vector(m_random_generator);
        p->value() = r_1;
        double local_action_1 = get_local_action(*p);

        std::uniform_real_distribution<> distribution_1(0.0, 1.0);
        double ds = local_action_1 - local_action_0;
        if (std::min(1.0, std::exp(-ds)) < distribution_1(m_random_generator))
            p->value() = r_0;
    }

    double beta() const
    {
        return m_beta;
    }

    Vector3_default external_field() const
    {
        return m_external_field;
    }

    std::string algorithm_string() const
    {
        return "metropolis";
    }

private:
    double get_local_action(Point<Vector3_default>& p)
    {
        double local_action = 0.0;
        for (auto& neighbor : p.neighbors())
            local_action -= m_beta * p.value() * neighbor->value();
        local_action -= p.value() * m_external_field;
        return local_action;
    }

    double m_beta;
    Vector3_default m_external_field;
    Random_generator& m_random_generator;
};

template <typename Random_generator>
class Configuration_wolff : public Configuration_square_cluster<Vector3_default, Connect<Random_generator>>, public Configuration_square_parameters
{

public:
    Configuration_wolff(Random_generator& random_generator)
        : Configuration_square<Vector3_default>(),
          Configuration_square_parameters(),
          Configuration_square_cluster<Vector3_default, Connect<Random_generator>>(), m_random_generator(random_generator) {}
    Configuration_wolff(int linear_size, double beta, double external_field_length, Random_generator& random_generator)
        : Configuration_square<Vector3_default>(linear_size),
          Configuration_square_parameters(linear_size),
          Configuration_square_cluster<Vector3_default, Connect<Random_generator>>(linear_size),
          m_beta(beta),
          m_external_field(Vector3_default({0.0, 0.0, external_field_length})),
          m_random_generator(random_generator) {}

    void write(std::ostream& os) const
    {
        Configuration_square<Vector3_default>::write(os);
        ::write(m_beta, os);
        ::write(m_external_field, os);
    }

    virtual void read(std::istream& is)
    {
        Configuration_square<Vector3_default>::read(is);
        ::read(m_beta, is);
        ::read(m_external_field, is);
    }

    void update()
    {
        std::array<Vector3_default, 3> r = random_normalized_perpendicular_vectors(m_random_generator);
        update_single_component(r[0]);
        update_single_component(r[1]);
        update_single_component(r[2]);
    }

    double beta() const
    {
        return m_beta;
    }

    Vector3_default external_field() const
    {
        return m_external_field;
    }

    std::string algorithm_string() const
    {
        return "wolff";
    }

private:
    void update_single_component(Vector3_default r)
    {
        Connect<Random_generator> connect(m_beta, r, m_random_generator);
        std::uniform_int_distribution<> distribution_linear_size(0, this->linear_size() - 1);
        std::shared_ptr<Point<Vector3_default>> starting_point = this->value(distribution_linear_size(m_random_generator), distribution_linear_size(m_random_generator));
        std::set<std::shared_ptr<Point<Vector3_default>>> cluster = this->build_cluster(starting_point, connect);

        double cluster_magnetization_r = 0.0;
        for (const auto& point : cluster)
            cluster_magnetization_r += point->value() * r;

        double flip_probability_nominator = std::exp(-2.0 * (m_external_field * r) * cluster_magnetization_r);
        double flip_probability = flip_probability_nominator / (flip_probability_nominator + 1.0);
        std::uniform_real_distribution<> distribution(0.0, 1.0);
        if (distribution(m_random_generator) < flip_probability)
            for (auto& point : cluster)
                point->value() = vector_normalize(point->value() - (2.0 * (r * point->value())) * r);
    }

    double m_beta;
    Vector3_default m_external_field;
    Random_generator& m_random_generator;
};

template <typename Random_generator>
class Configuration_shadow : public Configuration_square_cluster<Vector3_default, Connect<Random_generator>>, public Configuration_square_parameters
{

public:
    Configuration_shadow(Random_generator& random_generator)
        : Configuration_square<Vector3_default>(),
          Configuration_square_parameters(),
          Configuration_square_cluster<Vector3_default, Connect<Random_generator>>(),
          m_random_generator(random_generator) {}
    Configuration_shadow(int linear_size, double beta, double external_field_length, Random_generator& random_generator)
        : Configuration_square<Vector3_default>(linear_size),
          Configuration_square_parameters(linear_size),
          Configuration_square_cluster<Vector3_default, Connect<Random_generator>>(linear_size),
          m_beta(beta),
          m_external_point(std::make_shared<Point<Vector3_default>>()),
          m_random_generator(random_generator)
    {
        m_external_point->value() = Vector3_default({0.0, 0.0, external_field_length / beta});
        set_neighbors_external_point();
    }

    void write(std::ostream& os) const
    {
        Configuration_square<Vector3_default>::write(os);
        ::write(m_beta, os);
        ::write(m_external_point->value(), os);
    }

    virtual void read(std::istream& is)
    {
        Configuration_square<Vector3_default>::read(is);
        ::read(m_beta, is);
        Vector3_default external_point_value;
        ::read(external_point_value, is);
        m_external_point = std::make_shared<Point<Vector3_default>>(external_point_value);
        set_neighbors_external_point();

    }

    void update()
    {
        std::array<Vector3_default, 3> r = random_normalized_perpendicular_vectors(m_random_generator);
        update_single_component(r[0]);
        update_single_component(r[1]);
        update_single_component(r[2]);
    }

    double beta() const
    {
        return m_beta;
    }

    Vector3_default external_field() const
    {
        return m_external_point->value();
    }

    std::string algorithm_string() const
    {
        return "shadow";
    }

private:
    void set_neighbors_external_point()
    {
        for (int x = 0; x != this->linear_size(); ++x)
            for (int y = 0; y != this->linear_size(); ++y)
                this->set_neighbor(this->value(x, y), m_external_point);
    }

    void update_single_component(Vector3_default r)
    {
        Connect<Random_generator> connect(m_beta, r, m_random_generator);
        std::uniform_int_distribution<> distribution_linear_size(0, this->linear_size() - 1);
        std::shared_ptr<Point<Vector3_default>> starting_point = this->value(distribution_linear_size(m_random_generator), distribution_linear_size(m_random_generator));
        std::set<std::shared_ptr<Point<Vector3_default>>> cluster = this->build_cluster(starting_point, connect);

        if (cluster.count(m_external_point) == 0)
            for (auto& point : cluster)
                point->value() = vector_normalize(point->value() - (2.0 * (r * point->value())) * r);
    }

    double m_beta;
    std::shared_ptr<Point<Vector3_default>> m_external_point;
    Random_generator& m_random_generator;
};

Vector3_default magnetization(Configuration_square_parameters& configuration)
{
    Vector3_default m({0.0, 0.0, 0.0});
    for (int x = 0; x != configuration.linear_size(); ++x)
        for (int y = 0; y != configuration.linear_size(); ++y)
            m = m + configuration.value(x, y)->value();

    return m / static_cast<double>(configuration.linear_size() * configuration.linear_size());
}

std::vector<std::array<double, 2>> get_two_point_function(const Configuration_square_parameters& configuration)
{
    bool external_field_non_null = vector_length(configuration.external_field()) > epsilon;
    Vector3_default external_field_direction;
    if (external_field_non_null)
        external_field_direction = vector_normalize(configuration.external_field());

    std::vector<std::array<double, 2>> two_point_function(configuration.linear_size());
    for (int tau = 0; tau != configuration.linear_size(); ++tau)
    {
        two_point_function[tau] = std::array<double, 2>({0.0, 0.0});
        for (int x = 0; x != configuration.linear_size(); ++x)
        {
            for (int y = 0; y != configuration.linear_size(); ++y)
            {
                double all = configuration.value(x, tau)->value() * configuration.value(y, 0)->value();
                if (external_field_non_null)
                {
                    double parallel = (configuration.value(x, tau)->value() * external_field_direction) * (configuration.value(y, 0)->value() * external_field_direction);
                    two_point_function[tau][0] += parallel;
                    two_point_function[tau][1] += all - parallel;
                }
                else
                {
                    two_point_function[tau][0] += all / 3.0;
                    two_point_function[tau][1] += all / 1.5;
                }
            }
        }
        two_point_function[tau][0] /= configuration.linear_size() * configuration.linear_size();
        two_point_function[tau][1] /= configuration.linear_size() * configuration.linear_size();
    }

    return two_point_function;
}

double get_moment_0(const std::vector<std::array<double, 2>>& two_point_function, int i)
{
    double moment_0 = 0.0;
    for (int x = 0; x != two_point_function.size(); ++x)
        moment_0 += two_point_function[x][i];
    return moment_0 / two_point_function.size();
}

double get_moment_2(const std::vector<std::array<double, 2>>& two_point_function, int i)
{
    double moment_2 = 0.0;
    for (int x = 0; x != two_point_function.size(); ++x)
        moment_2 += std::pow(2.0 * std::sin(pi * x / two_point_function.size()), 2) * two_point_function[x][i];
    return (moment_2 / std::pow(2.0 * pi, 2)) / two_point_function.size();
}

enum class Algorithm
{
    METROPOLIS,
    WOLFF,
    SHADOW,
    INVALID
};

Algorithm parse_algorithm(char* arg)
{
    if (std::strcmp(arg, "metropolis") == 0)
    {
        return Algorithm::METROPOLIS;
    }
    else if (std::strcmp(arg, "wolff") == 0)
    {
        return Algorithm::WOLFF;
    }
    else if (std::strcmp(arg, "shadow") == 0)
    {
        return Algorithm::SHADOW;
    }
    else
    {
        return Algorithm::INVALID;
    }
}

int main(int argc, char* argv[])
{
    if (argc < 5 || 8 < argc)
    {
        std::cout << "Usage: \n";
        std::cout << "    " << argv[0] << " sweeps correlation_length algorithm linear_size beta external_field_length [configuration_output_filename]" << '\n';
        std::cout << "    " << argv[0] << " sweeps correlation_length algorithm configuration_input_filename [configuration_output_filename]" << '\n';
        return 1;
    }

    int sweeps = std::stoi(argv[1]);
    int correlation_length = std::stoi(argv[2]);
    Algorithm algorithm = parse_algorithm(argv[3]);

    std::random_device random_device;
    using random_generator_type = std::ranlux48;
    random_generator_type random_generator(random_device());

    std::shared_ptr<Configuration_square_parameters> configuration;

    const char* configuration_input_filename = nullptr;
    if (argc == 5 || argc == 6)
        configuration_input_filename = argv[4];
    if (configuration_input_filename)
    {
        switch (algorithm)
        {
            case Algorithm::METROPOLIS:
                configuration = std::make_shared<Configuration_metropolis<random_generator_type>>(random_generator);
                break;
            case Algorithm::WOLFF:
                configuration = std::make_shared<Configuration_wolff<random_generator_type>>(random_generator);
                break;
            case Algorithm::SHADOW:
                configuration = std::make_shared<Configuration_shadow<random_generator_type>>(random_generator);
                break;
        }

        std::ifstream is;
        is.open(configuration_input_filename, std::ios::in | std::ios::binary);
        if (is.peek() == std::ifstream::traits_type::eof())
        {
            std::cout << "Empty input file '" << configuration_input_filename << "'\n";
            return 1;
        }
        read(*configuration, is);
        is.close();
    }

    if (argc == 7 || argc == 8)
    {
        int linear_size = std::stoi(argv[4]);
        double beta = std::stod(argv[5]);
        double external_field_length = std::stod(argv[6]);

        switch (algorithm)
        {
            case Algorithm::METROPOLIS:
                configuration = std::make_shared<Configuration_metropolis<random_generator_type>>(linear_size, beta, external_field_length, random_generator);
                break;
            case Algorithm::WOLFF:
                configuration = std::make_shared<Configuration_wolff<random_generator_type>>(linear_size, beta, external_field_length, random_generator);
                break;
            case Algorithm::SHADOW:
                configuration = std::make_shared<Configuration_shadow<random_generator_type>>(linear_size, beta, external_field_length, random_generator);
                break;
        }
    }

    const char* configuration_output_filename = nullptr;
    if (argc == 6 || argc == 8)
        configuration_output_filename = argv[argc - 1];

    std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);

    std::cout << "# " << std::ctime(&now_time) << '\n';
    std::cout << "# sweeps " << sweeps << '\n';
    std::cout << "# correlation_length " << static_cast<int>(correlation_length) << '\n';
    std::cout << "# algorithm_string " << configuration->algorithm_string() << '\n';
    std::cout << "# algorithm " << static_cast<int>(algorithm) << '\n';
    std::cout << "# linear_size " << configuration->linear_size() << '\n';
    std::cout << "# beta " << configuration->beta() << '\n';
    std::cout << "# external_field_length " << vector_length(configuration->external_field()) << '\n';
    std::cout << "# configuration_output_filename ";
    if (configuration_output_filename != nullptr)
        std::cout << '"' <<  configuration_output_filename << '"';
    else
        std::cout << "none";
    std::cout << '\n';
    std::cout << "# configuration_input_filename ";
    if (configuration_input_filename != nullptr)
        std::cout << '"' <<  configuration_input_filename << '"';
    else
        std::cout << "none";
    std::cout << '\n';


    std::chrono::high_resolution_clock::time_point t_0 = std::chrono::high_resolution_clock::now();

    std::cout << "# magnetization moment_0_parallel moment_0_perpendicular moment_2_parallel moment_2_perpendicular\n";
    for (int i = 0; i != sweeps; ++i)
    {
        for (int j = 0; j != correlation_length; ++j)
            configuration->update();

        std::vector<std::array<double, 2>> two_point_function;
        two_point_function = get_two_point_function(*configuration);

        double moment_0_parallel = get_moment_0(two_point_function, 0);
        double moment_0_perpendicular = get_moment_0(two_point_function, 1);
        double moment_2_parallel = get_moment_2(two_point_function, 0);
        double moment_2_perpendicular = get_moment_2(two_point_function, 1);

        double current_magnetization = vector_length(magnetization(*configuration));

        std::cout << current_magnetization << ' ' << moment_0_parallel << ' ' << moment_0_perpendicular << ' ' << moment_2_parallel << ' ' << moment_2_perpendicular << '\n';
    }

    std::chrono::high_resolution_clock::time_point t_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = t_1 - t_0;

    std::cout << "# timespan " << time_span.count() << '\n';

    if (configuration_output_filename)
    {
        std::ofstream os;
        os.open(configuration_output_filename, std::ios::out | std::ios::binary);
        write(*configuration, os);
        os.close();
    }

    return 0;
}
