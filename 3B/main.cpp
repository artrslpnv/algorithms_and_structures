#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <queue>
#include <unordered_set>

const double E = 0.0001;


struct Point {
    double x = 0;
    double y = 0;
    double z = 0;
    int index = -1;

    Point() = default;

    Point(double _x, double _y, double _z, int _index) : x(_x), y(_y), z(_z), index(_index) {}

    friend bool operator<(const Point &p1, const Point &p2) {
        if (p1.z != p2.z) {
            return p1.z < p2.z;
        } else if (p1.y != p2.y) {
            return p1.y < p2.y;
        } else {
            return p1.x < p2.x;
        }
    }

    friend bool operator==(const Point &p1, const Point &p2) {
        return p1.x == p2.x && p1.y == p2.y && p1.z == p2.z;
    }
};

class Vector {
public:
    Vector() = default;

    Vector(Point p) {
        _x = p.x;
        _y = p.y;
        _z = p.z;
    }

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

double scalar_multiply(const Vector &vector1, const Vector &vector2) {
    return vector1._x * vector2._x + vector1._y * vector2._y + vector1._z * vector2._z;
};

Vector vector_multiply(const Vector &vector1, const Vector &vector2) {
    Vector result;
    result._x = vector1._y * vector2._z - vector1._z * vector2._y;
    result._y = vector1._z * vector2._x - vector2._z * vector1._x;
    result._z = vector1._x * vector2._y - vector2._x * vector1._y;
    return result;
}

double find_cos_of_angle_between(const Vector &vector1, const Vector &vector2) {
    double scalar = scalar_multiply(vector1, vector2);
    double len_of_first = vector1.len();
    double len_of_second = vector2.len();
    return scalar / (len_of_first * len_of_second);
}

double find_sin_of_angle_between(const Vector &vector1, const Vector &vector2) {
    Vector vector3 = vector_multiply(vector1, vector2);
    return vector3.len() / (vector1.len() * vector2.len());
}

struct Edge {
    Edge(Point _p1, Point _p2) {

        p1 = _p1, p2 = _p2;

    }

    Point p1;
    Point p2;

    friend bool operator==(const Edge &face1, const Edge &face2) {
        return (face1.p1 == face2.p1 && face1.p2 == face2.p2) || (face1.p2 == face2.p1 && face1.p1 == face2.p2);
    }
};

struct Face {
    Point p1;
    Point p2;
    Point p3;

    Face(Point _p1, Point _p2, Point _p3) : p1(_p1), p2(_p2), p3(_p3) {
        Order();
    }

    void Order() {
        if (p2.index < p1.index && p2.index < p3.index) {
            std::swap(p1, p2);
            std::swap(p2, p3);
        } else if (p3.index < p1.index && p3.index < p2.index) {
            std::swap(p2, p3);
            std::swap(p1, p2);
        }
    }

    Vector find_normal() {
        Vector p1p2 = Vector(p2) - Vector(p1);
        Vector p2p3 = Vector(p3) - Vector(p2);
        return vector_multiply(p1p2, p2p3);
    }

    friend bool operator<(const Face &face1, const Face &face2) {
        if (face1.p1.index < face2.p1.index) {
            return true;
        } else if (face1.p1.index == face2.p1.index) {
            if (face1.p2.index < face2.p2.index) {
                return true;
            } else if (face1.p2.index == face2.p2.index) {
                return face1.p3.index < face2.p3.index;
            } else { return false; }
        } else { return false; }
    }

    friend bool operator==(const Face &face1, const Face &face2) {
        return face1.p1 == face2.p1 && face1.p2 == face2.p2 && face1.p3 == face2.p3;
    }
};


int find_the_first_hull_dot(const std::vector<Point> &points) {
    int result = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i] < points[result]) {
            result = i;
        }
    }
    return result;
}

struct Face_hasher {
    size_t operator()(const Face &face) const noexcept {
        return (face.p1.x + face.p1.y + face.p2.z) * 100;
    }
};

struct Edge_hasher {
    size_t operator()(const Edge &edge) const noexcept {
        return (edge.p1.x + edge.p2.y) * 100;
    }
};

int find_the_second_dot(const std::vector<Point> &points, Point first_dot) {
    Vector v1(1, 0, 0);
    Vector n_Oxy(0, 0, 1);
    double max_cos = -1.1;
    int index = -1;
    for (int i = 0; i < points.size(); ++i) {
        if (i != first_dot.index) {
            Vector v2 = Vector(points[i]) - Vector(first_dot);
            Vector v3 = vector_multiply(v1, v2);
            double cos = find_cos_of_angle_between(v3, n_Oxy);
            if (cos > max_cos) {
                max_cos = cos;
                index = i;
            }
        }
    }
    return index;
}

int find_next_dot(const std::vector<Point> &points, Point p1, Point p2, Face last_face) {
    Vector v1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    double max_cos = -1.1;
    int index = -1;
    for (int i = 0; i < points.size(); ++i) {
        if (i != p1.index && i != p2.index && i != last_face.p1.index && i != last_face.p2.index &&
            i != last_face.p3.index) {
            Vector v2 = Vector(points[i]) - Vector(p2);
            Vector v3 = vector_multiply(v1, v2);
            double cos = find_cos_of_angle_between(v3, last_face.find_normal());
            if (cos > max_cos) {
                max_cos = cos;
                index = i;
            }
        }
    }
    return index;
}

std::vector<Face> Build3dHull(std::vector<Point> &Points) {
    Point first_dot = Points[find_the_first_hull_dot(Points)];
    std::unordered_set<Face, Face_hasher> processed_faces;
    std::unordered_set<Edge, Edge_hasher> processed_edges;
    std::queue<std::pair<Edge, Face>> queue;
    std::vector<Face> hull;
    Point second_point = Points[find_the_second_dot(Points, first_dot)];
    Point Helping_dot(first_dot.x - 0.01, first_dot.y, first_dot.z, -1);
    Face last_face(Helping_dot, first_dot, second_point);
    queue.push(std::make_pair(Edge(first_dot, second_point), last_face));
    while (!queue.empty()) {
        std::pair<Edge, Face> e = queue.front();
        queue.pop();
        if (processed_edges.count(e.first) == 0) {
            Point p3 = Points[find_next_dot(Points, e.first.p1, e.first.p2, e.second)];
            Face new_face(e.first.p1, e.first.p2, p3);
            if (processed_faces.count(new_face) == 0) {
                hull.push_back(new_face);
                processed_faces.insert(new_face);
            }
            processed_edges.insert(e.first);
            queue.push(std::make_pair(Edge(e.first.p2, p3), new_face));
            queue.push(std::make_pair(Edge(p3, e.first.p1), new_face));
        }
    }
    return hull;
}

void rotate(Point &p, double angle) {
    double new_x = p.x * cos(angle) + p.z * sin(angle);
    double new_z = -p.x * sin(angle) + p.z * cos(angle);
    p.x = new_x;
    p.z = new_z;

    new_z = p.z * cos(angle) + p.y * sin(angle);
    double new_y = -p.z * sin(angle) + p.y * cos(angle);
    p.z = new_z;
    p.y = new_y;

    new_x = p.x * cos(angle) + p.y * sin(angle);
    new_y = -p.x * sin(angle) + p.y * cos(angle);
    p.x = new_x;
    p.y = new_y;
}

int main() {
    int n;
    std::cin >> n;
    for (size_t i = 0; i < n; i++) {
        int m;
        std::vector<Point> points;

        std::cin >> m;
        for (int j = 0; j < m; j++) {
            int x, y, z;
            std::cin >> x >> y >> z;
            Point p(x, y, z, j);
            rotate(p, E);
            points.push_back(p);
        }

        std::vector<Face> hull = Build3dHull(points);
        //std::sort(hull.begin(), hull.end());
        std::cout << hull.size() << "\n";
        for (Face &f : hull) {
            std::cout << 3 << " " << f.p1.index << " " << f.p2.index << " " << f.p3.index << "\n";
        }
    }
    return 0;
}