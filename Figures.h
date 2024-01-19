//
// Created by motherlode on 20.01.23.
//

#ifndef UNTITLED4_FIGURES_H
#define UNTITLED4_FIGURES_H
#include "System.h"
typedef struct Ellipse{
    Vector center;
    double a;
    double b;
    double c;
    Matrix R;
    bool dot_in_ellipse(const Vector & dot_) {
        //std::cout<<'h'<<std::endl;
        Vector dot;
        dot = dot_ + center * (-1);
        dot = invert(R) * dot;
        if(dot.x * dot.x / (a * a) + dot.y * dot.y / (b * b)  + dot.z * dot.z / (c * c)  <= 1 ) return true;
        return false;
    }
    //x = Rx'' + C
    void print_dots() {
        std::cout<<center;
        Vector ax = R * Vector(a,0,0) + center;
        Vector ax1 = R * Vector(-a,0,0) + center;
        Vector bx = R * Vector(0,b,0) + center;
        Vector bx1 = R * Vector(0,-b,0) + center;
        Vector cx = R * Vector(0,0,c) + center;
        Vector cx1 = R * Vector(0,0,-c) + center;
        std::cout<<ax<<ax1<<bx<<bx1<<cx<<cx1<<std::endl;
    }
} Ellipse;

typedef struct Plane{
    Vector n = Vector(1, 1, 1);
    double D = 0;
    //x + z = 0
}Plane;


Vector get_contact(Ellipse & el, Plane & pl) {
    Vector dot = el.center;
    Vector p;
    double A = pl.n.x;
    double B = pl.n.y;
    double C = pl.n.z;
    double D = pl.D;
    double t = -(A * dot.x + B * dot.y + C * dot.z + D) / (A*A + B*B + C*C);
    p.x = dot.x + t * A;
    p.y = dot.y + t * B;
    p.z = dot.z + t * C;
    return p;
}

bool contact(Ellipse & el, Plane & pl) {
    if(el.center.x * pl.n.x + el.center.y * pl.n.y + el.center.z * pl.n.z + pl.D <= 0 ) return true;
    if(get_contact(el, pl).length() > max(max(el.a, el.b), el.c)) return false;
    return true;
    //return (el.dot_in_ellipse(get_contact(el, pl)));
}

void solve(Ellipse * ellipse, SystemState * x_rb, Plane plane,const Constants * data) {
    Vector p = get_contact(*ellipse, plane); /* world-space vertex location */
    Vector n = x_rb->C + p * (-1);
    n = n * (1/n.length());
            //plane.n * (1/(plane.n.length())); /* outwards pointing normal of face */
    Vector pt_velocity = x_rb ->v + vec_mult(x_rb->omega, p * (-1) + x_rb->C ); // p˙−a (t0)
    double v_rel = n * pt_velocity;//v-rel
    double epsilon = 1;
    Vector ra = p * (-1) + x_rb->C ; /* ra */
    double numerator = -(1 + epsilon) * v_rel;
    double term1 = 1 / data->mass,
            term2 = 0,
            term3 = n * vec_mult((x_rb->Iinv * vec_mult(ra , n)) , ra),
            term4 = 0;
    double j = numerator / (term1 + term2 + term3 + term4);
    Vector force = n * (j);
    x_rb->p = x_rb->p + force;
    x_rb->L = x_rb->L + vec_mult(ra, force);
    x_rb->v = x_rb->p * (1 / data->mass);
    x_rb->omega = x_rb->Iinv * x_rb->L;
}


#endif //UNTITLED4_FIGURES_H
