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


    friend double scalar_multiply(const Vector &vector1, const Vector &vector2) {
        return vector1._x * vector2._x + vector1._y * vector2._y + vector1._z * vector2._z;
    };


    double _x = 0;
    double _y = 0;
    double _z = 0;
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
        return (face1.p1 == face2.p1 && face1.p2 == face2.p2);
    }
};

struct Face {
    Point p1;
    Point p2;
    Point p3;

    Face() : p1(), p2(), p3() {}

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

    void NormalizeFace(std::vector<Point> &points, int first, int second, int third) {
        Face tmp_face(points[first], points[second], points[third]);
        first = tmp_face.p1.index;
        second = tmp_face.p2.index;
        third = tmp_face.p3.index;

        Vector b1(Vector(points[second]) - Vector(points[first]));
        Vector b2(Vector(points[third]) - Vector(points[second]));
        Vector normal_vector = vector_multiply(b1, b2);
        normal_vector = normal_vector / normal_vector.len();
        int one_different_point_from_points_index = -1;
        int size = points.size();
        for (int i = 0; i < size; ++i) {
            if (i != first && i != second && i != third) {
                one_different_point_from_points_index = i;
                break;
            }
        }

        bool counter = true;
        if (scalar_multiply(normal_vector,
                            Vector(points[one_different_point_from_points_index]) - Vector(points[first])) > -E) {
            normal_vector = -1 * normal_vector;
            counter = false;
        }

        if (!counter) {
            p1 = points[first];
            p2 = points[third];
            p3 = points[second];
        } else {
            p1 = points[first];
            p2 = points[second];
            p3 = points[third];
        }
    }

    friend bool operator==(const Face &face1, const Face &face2) {
        return face1.p1 == face2.p1 && face1.p2 == face2.p2 && face1.p3 == face2.p3;
    }

    friend bool operator!=(const Face &face1, const Face &face2) {
        return !(face1.p1 == face2.p1 && face1.p2 == face2.p2 && face1.p3 == face2.p3);
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
    double max_cos = -1.1;
    int index = -1;
    for (int i = 0; i < points.size(); ++i) {
        if (i != first_dot.index) {
            Vector v2 = Vector(points[i]) - Vector(first_dot);
            double cos = sqrt(v2._x * v2._x + v2._y * v2._y) / v2.len();
            if (cos > max_cos) {
                max_cos = cos;
                index = i;
            }
        }
    }
    return index;
}

int find_point(const std::vector<Point> &points, int first, int second, int third = -1) {
    int size = points.size();
    int new_point = -1;
    Vector edge_vector(Vector(points[second]) - Vector(points[first]));
    for (int i = 0; i < size; ++i) {
        if (i != first && i != second && i != third) {
            if (new_point == -1) {
                new_point = i;
                continue;
            }
            Vector v1(Vector(points[new_point]) - Vector(points[first]));
            Vector v2(Vector(points[i]) - Vector(points[first]));

            Vector normal_vector = vector_multiply(v1, v2);
            if (scalar_multiply(normal_vector, edge_vector) > -E) {
                new_point = i;
            }
        }
    }
    return new_point;
}

Face InsertFace(std::unordered_set<Edge, Edge_hasher> &ProcessedEdges, const std::vector<Point> &points, Edge edge,
                int third) {
    std::swap(edge.p1, edge.p2);
    Face result;
    if (ProcessedEdges.count(edge) == 0) {
        int point_to_add_index = find_point(points, edge.p1.index, edge.p2.index, third);
        result = Face(edge.p1, edge.p2, points[point_to_add_index]);
    }
    return result;
}

std::vector<Face> Build3dHull(std::vector<Point> &points) {
    std::vector<Face> Hull;
    int first = find_the_first_hull_dot(points);
    int second = find_the_second_dot(points, points[first]);
    int third = find_point(points, first, second);
    Face first_face(points[first], points[second], points[third]);
    first_face.NormalizeFace(points, first, second, third);
    std::queue<Face> queue;
    queue.push(first_face);
    std::unordered_set<Edge, Edge_hasher> ProcessedEdges;
    ProcessedEdges.insert(Edge(first_face.p1, first_face.p2));
    ProcessedEdges.insert(Edge(first_face.p2, first_face.p3));
    ProcessedEdges.insert(Edge(first_face.p3, first_face.p1));
    while (!queue.empty()) {
        Face cur_face = queue.front();
        queue.pop();
        Hull.push_back(cur_face);
        Edge edge1(cur_face.p1, cur_face.p2);
        Edge edge2(cur_face.p2, cur_face.p3);
        Edge edge3(cur_face.p3, cur_face.p1);
        Face face_to_check;
        Face face_to_add = InsertFace(ProcessedEdges, points, edge1, edge2.p2.index);
        if (face_to_add != face_to_check) {
            queue.push(face_to_add);
            ProcessedEdges.insert(Edge(face_to_add.p1, face_to_add.p2));
            ProcessedEdges.insert(Edge(face_to_add.p2, face_to_add.p3));
            ProcessedEdges.insert(Edge(face_to_add.p3, face_to_add.p1));
        }
        face_to_add = InsertFace(ProcessedEdges, points, edge2, edge3.p2.index);
        if (face_to_add != face_to_check) {
            ProcessedEdges.insert(Edge(face_to_add.p1, face_to_add.p2));
            ProcessedEdges.insert(Edge(face_to_add.p2, face_to_add.p3));
            ProcessedEdges.insert(Edge(face_to_add.p3, face_to_add.p1));
            queue.push(face_to_add);
        }
        face_to_add = InsertFace(ProcessedEdges, points, edge3, edge1.p2.index);
        if (face_to_add != face_to_check) { ;
            ProcessedEdges.insert(Edge(face_to_add.p1, face_to_add.p2));
            ProcessedEdges.insert(Edge(face_to_add.p2, face_to_add.p3));
            ProcessedEdges.insert(Edge(face_to_add.p3, face_to_add.p1));
            queue.push(face_to_add);
        }
    }
    return Hull;
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
        std::sort(hull.begin(), hull.end());
        std::cout << hull.size() << "\n";
        for (Face &f : hull) {
            std::cout << 3 << " " << f.p1.index << " " << f.p2.index << " " << f.p3.index << "\n";
        }
    }
    return 0;
}