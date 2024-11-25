
#include <iostream>
#include <vector>

struct Point {
    double x, y, z;
};

std::ostream& operator<<(std::ostream& os, const Point& p) {
    return os << "[ " << p.x << ", " << p.y << ", " << p.z << " ]";
}

void FromFrame3dd(double *t1, double *t2, double *t3, double *t4, double *t5, double *t6,
        double *t7, double *t8, double *t9, const Point& p1, const Point& p2, double p) {
    const double L = std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2) + std::pow(p2.z - p1.z, 2));
    const double Cx = (p2.x - p1.x) / L;
    const double Cy = (p2.y - p1.y) / L;
    const double Cz = (p2.z - p1.z) / L;

    *t1 = *t2 = *t3 = *t4 = *t5 = *t6 = *t7 = *t8 = *t9 = 0.0;

    const double Cp = std::cos(p);
    const double Sp = std::sin(p);

    if (std::abs(Cz) == 1.0) {
        *t3 =  Cz;
        *t4 = -Cz*Sp;
        *t5 =  Cp;
        *t7 = -Cz*Cp;
        *t8 = -Sp;
    } else {

        const double den = std::sqrt ( 1.0 - Cz*Cz );

        *t1 = Cx;
        *t2 = Cy;
        *t3 = Cz;

        *t4 = (-Cx*Cz*Sp - Cy*Cp)/den;
        *t5 = (-Cy*Cz*Sp + Cx*Cp)/den;
        *t6 = Sp*den;

        *t7 = (-Cx*Cz*Cp + Cy*Sp)/den;
        *t8 = (-Cy*Cz*Cp - Cx*Sp)/den;
        *t9 = Cp*den;
    }
}

std::vector<std::vector<double>> CoorTransform(const Point& p1, const Point& p2, double roll) {
    // Coefficients of the final transformation matrix:
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    double t6;
    double t7;
    double t8;
    double t9;

    FromFrame3dd(&t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8, &t9, p1, p2, roll);

    std::vector<std::vector<double>> a{};

    a.resize(3);
    for (std::size_t i = 0; i < 3; ++i) {
        a[i].resize(3);
    }

    a[0][0] = t1;
    a[0][1] = t2;
    a[0][2] = t3;
    a[1][0] = t4;
    a[1][1] = t5;
    a[1][2] = t6;
    a[2][0] = t7;
    a[2][1] = t8;
    a[2][2] = t9;

    return a;
}

void Print(const std::vector<std::vector<double>>& a) {
    for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
            std::cout << a[i][j] << "  ";
        }
        std::cout << '\n';
    }
}

void Transform(const std::vector<std::vector<double>>& a, const Point& p) {
    Point pp;

    pp.x = a[0][0]*p.x + a[0][1]*p.y + a[0][2]*p.z;
    pp.y = a[1][0]*p.x + a[1][1]*p.y + a[1][2]*p.z;
    pp.z = a[2][0]*p.x + a[2][1]*p.y + a[2][2]*p.z;

    std::cout << "Transformed " << p << " to " << pp << '\n';
}

int main(int, char**) {
    const Point p1{0, 0, 0};
    const Point p2{0, 0, 5};
    const double roll = M_PI/2;

    auto a = CoorTransform(p1, p2, roll);

    Transform(a, Point{1, 0, 0});
}
