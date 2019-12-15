#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

const double E = 0.00000001;

class Vector {
public:
    Vector() = default;

    ~Vector() {};

    explicit Vector(double x, double y, double z) noexcept : _x(x), _y(y), _z(z) {}

    explicit Vector(double x, double y) : _z(0), _x(x), _y(y) {}

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

    friend bool is_less_polar_angle(Vector v1, Vector v2) {
        if (v1._x > 0 && v2._x < 0) {
            return true;
        }
        if (v1._x < 0 && v2._x > 0) {
            return false;
        }

        double tg1 = v1._y / v1._x;
        double tg2 = v2._y / v2._x;

        if (fabs(v1._x) < E) {
            if (v1._y > 0) {
                return v2._x < 0;
            } else {
                return false;
            }
        }
        if (fabs(v2._x) < E) {
            if (v2._y > 0) {
                return v1._x > 0;
            } else {
                return fabs(tg1 - tg2) > E;
            }
        }

        if (fabs(tg1 - tg2) < E) {
            return false;
        }
        return tg1 < tg2;

    };

    double _x = 0;
    double _y = 0;
    double _z = 0;
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

class Polygon {
public:
    Polygon() = default;

    explicit Polygon(std::vector <Vector> &_polygon_array) : polygon_array(_polygon_array) {}

    int find_extremal_dot_index() const {
        Vector extremal_dot = polygon_array[0];
        int index = 0;
        for (int i = 1; i < polygon_array.size(); ++i) {
            if ((polygon_array[i]._x - extremal_dot._x > E) && polygon_array[i]._x < extremal_dot._x) {
                extremal_dot = polygon_array[i];
                index = i;
            } else if (polygon_array[i]._x - extremal_dot._x < E) {
                if (polygon_array[i]._y < extremal_dot._y) {
                    extremal_dot = polygon_array[i];
                    index = i;
                }
            }
        }
        return index;
    };

    friend std::vector <Vector> Minkovskiy_sum(Polygon &polygon1, Polygon &polygon2) {
        std::vector <Vector> normalized_polygon1 = get_neccessary_permutation(polygon1);
        std::vector <Vector> normalized_polygon2 = get_neccessary_permutation(polygon2);
        normalized_polygon1.push_back(normalized_polygon1[0]);
        normalized_polygon2.push_back(normalized_polygon2[0]);
        std::vector <Vector> result;
        int i = 0, j = 0;
        while (i < normalized_polygon1.size() - 1 && j < normalized_polygon2.size() - 1) {
            result.push_back(normalized_polygon1[i] + normalized_polygon2[j]);
            if (is_less_polar_angle(normalized_polygon1[i + 1] - normalized_polygon1[i],
                                    normalized_polygon2[j + 1] - normalized_polygon2[j])) {
                ++i;
            } else if (is_less_polar_angle(normalized_polygon2[j + 1] - normalized_polygon2[j],
                                           normalized_polygon1[i + 1] - normalized_polygon1[i])) {
                ++j;
            } else {
                ++i;
                ++j;
            }
        }
        return result;
    }

    friend bool Is_zero_zero_in_polygon(Polygon &polygon) {
        size_t i = 0;
        std::vector<double> crossed;

        for (; i < polygon.polygon_array.size(); ++i) {
            size_t next = (i != polygon.polygon_array.size() - 1 ? i + 1 : 0);
            if (fabs(polygon.polygon_array[i]._y) < E) {
                crossed.push_back(polygon.polygon_array[i]._x);
            }
            if (polygon.polygon_array[i]._y * polygon.polygon_array[next]._y < 0) {
                crossed.push_back((polygon.polygon_array[i]._x + polygon.polygon_array[next]._x) / 2);
            }
        }
        if (crossed.size() == 0) {
            return false;
        }
        double max_x_coordinate = crossed[0];
        double min_x_coordinate = crossed[0];
        for (double d: crossed) {
            if (d - max_x_coordinate > E) {
                max_x_coordinate = d;
            }
            if (min_x_coordinate - d > E) {
                min_x_coordinate = d;
            }
        }
        if (max_x_coordinate * min_x_coordinate < 0 || fabs(max_x_coordinate * min_x_coordinate) < E * E) {
            return true;
        } else { return false; }

    }

    friend std::vector <Vector> get_neccessary_permutation(const Polygon &polygon) {
        std::vector <Vector> r;
        size_t min = polygon.find_extremal_dot_index();
        for (size_t i = 0; i < polygon.polygon_array.size(); ++i) {
            r.push_back(polygon.polygon_array[(i + min + 1) % polygon.polygon_array.size()]);
        }
        std::reverse(r.begin(), r.end());
        return r;
    }

    std::vector <Vector> polygon_array;
};

int main() {
    int n, m;
    std::cin >> n;
    std::vector <Vector> first(n);
    for (size_t i = 0; i < n; ++i) {
        double x, y;
        std::cin >> x >> y;
        first[i] = Vector(x, y);
    }
    std::cin >> m;
    std::vector <Vector> second(m);
    for (size_t i = 0; i < m; ++i) {
        double x, y;
        std::cin >> x >> y;
        x = -x;
        y = -y;
        second[i] = Vector(x, y);
    }
    Polygon polygon1 = Polygon(first);
    Polygon polygon2 = Polygon(second);
    std::vector <Vector> a = Minkovskiy_sum(polygon1, polygon2);
    Polygon polygon(a);
    std::cout << (Is_zero_zero_in_polygon(polygon) ? "YES" : "NO");

}

