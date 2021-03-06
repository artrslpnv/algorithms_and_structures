#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double E = 0.0000001;

class Vector {
public:
    Vector() = default;

    ~Vector() {};

    explicit Vector(double x, double y, double z) noexcept : _x(x), _y(y), _z(z) {}

    Vector(const Vector &vector2) {
        _x = vector2._x;
        _y = vector2._y;
        _z = vector2._z;
    }

    Vector(Vector &&vector) noexcept {
        _x = vector._x;
        _y = vector._y;
        _z = vector._z;
    }

    Vector &operator=(const Vector &vector2) {
        _x = vector2._x;
        _y = vector2._y;
        _z = vector2._z;
        return *this;
    }


    friend double scalar_multiply(const Vector &vector1, const Vector &vector2) {
        return vector1._x * vector2._x + vector1._y * vector2._y + vector1._z * vector2._z;
    };

    friend Vector vector_multiply(const Vector &vector1, const Vector &vector2) {
        Vector result;
        result._x = vector1._y * vector2._z - vector1._z * vector2._y;
        result._y = vector1._z * vector2._x - vector2._z * vector1._x;
        result._z = vector1._x * vector2._y - vector2._x * vector1._y;
        return result;
    }

    friend Vector operator*(double constant, const Vector &vector) {
        Vector result;
        result._x = constant * vector._x;
        result._y = constant * vector._y;
        result._z = constant * vector._z;
        return result;

    };

    friend Vector operator/(const Vector &vector, double constant) {
        Vector result;
        result._x = vector._x / constant;
        result._y = vector._y / constant;
        result._z = vector._z / constant;
        return result;
    };

    Vector &operator*=(double constant) {
        _x *= constant;
        _y *= constant;
        _z *= constant;
        return *this;
    }

    friend Vector operator+(const Vector &vector1, const Vector &vector2) {
        Vector result;
        result._x = vector1._x + vector2._x;
        result._y = vector1._y + vector2._y;
        result._z = vector1._z + vector2._z;
        return result;
    }

    double len() const {
        return sqrt(_x * _x + _y * _y + _z * _z);
    }

    friend Vector operator-(const Vector &vector1, const Vector &vector2) {
        Vector result;
        result._x = vector1._x - vector2._x;
        result._y = vector1._y - vector2._y;
        result._z = vector1._z - vector2._z;
        return result;
    }

    friend double find_cos_of_angle_between(const Vector &vector1, const Vector &vector2) {
        double scalar = scalar_multiply(vector1, vector2);
        double len_of_first = vector1.len();
        double len_of_second = vector2.len();
        return scalar / (len_of_first * len_of_second);
    }

    friend double find_sin_of_angle_between(const Vector &vector1, const Vector &vector2) {
        Vector vector3 = vector_multiply(vector1, vector2);
        return vector3.len() / (vector1.len() * vector2.len());
    }

    double _x;
    double _y;
    double _z;
};

class Segment {
public:
    friend class Vector;

    Segment(double x, double y, double z, double x1, double y1, double z1) {
        from = Vector(x, y, z);
        to = Vector(x1, y1, z1);
    }

    double len() {
        return sqrt((to._x - from._x) * (to._x - from._x) + (to._y - from._y) * (to._y - from._y) +
                    (to._z - from._z) * (to._z - from._z));
    }

    Vector from;
    Vector to;
};

double distance_between_two_dots(Vector dot1, Vector dot2) {
    return (dot1 - dot2).len();
}

template<typename T>
class UnimodalFunc;

template<typename T>
double ternar_search(Segment segment1, UnimodalFunc<T> function);

template<typename T>
class UnimodalFunc {
public:
    explicit UnimodalFunc(T &t) : segment_or_vector(t) {};

    double operator()(Vector vector) {
        return distance(vector, segment_or_vector);
    }

private:
    T segment_or_vector;

    double distance(Vector vector, Vector vector1) {
        return distance_between_two_dots(vector, vector1);
    }

    double distance(Vector vector, Segment segment) {
        UnimodalFunc<Vector> func(vector);
        return ternar_search(segment, func);
    }
};

template<class T>
double ternar_search(Segment segment1, UnimodalFunc<T> function) {
    Vector l1 = segment1.from;
    Vector r1 = segment1.to;
    Vector m1 = l1 + (r1 - l1) / 3;
    Vector m2 = r1 - (r1 - l1) / 3;
    double distance = 0;
    while ((m1 - m2).len() > E) {

        double distance_m1 = function(m1);
        double distance_m2 = function(m2);

        if (distance_m1 > distance_m2) {
            l1 = m1;
            distance = distance_m2;
        } else {
            r1 = m2;
            distance = distance_m1;
        }
        m1 = l1 + (r1 - l1) / 3;
        m2 = r1 - (r1 - l1) / 3;
    }
    return distance;
}

int main() {
    double x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    std::cin >> x >> y >> z >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
    Segment segment = Segment(x, y, z, x1, y1, z1);
    Segment segment1 = Segment(x2, y2, z2, x3, y3, z3);
    UnimodalFunc<Segment> unimodalFunc(segment1);
    std::cout << std::setprecision(10) << ternar_search(segment, unimodalFunc);
    return 0;
}

