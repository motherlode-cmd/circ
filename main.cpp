#include <iostream>
#include <stdio.h>
#include "math.h"
#include <functional>
#include "System.h"
#include "Draw.h"
#include "Figures.h"

Matrix Star(Vector a){
    Matrix m;
    double matr[3][3] = {
            {0, -a.z, a.y},
            {a.z, 0, -a.x},
            {-a.y, a.x, 0}
    };
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            m.matrix[i][j] = matr[i][j];
        }
    }
    return m;
}

void f(const SystemState & X, SystemState & fx, const Constants * data)
{
    //Vector v = X.p / data->mass;
    //Matrix Iinv = X.R * data->Ibodyinv * transpose(X.R);
    //Vector omega = Iinv * X.L;
    fx.C = X.v;
    fx.R = Star(X.omega) * X.R;
    fx.p = X.force;
    fx.L = X.torque;
    //std::cout<<X.force.x<<' '<<X.force.y<<' '<<X.force.z<<' ';
}

double a_ralston[4][4] = {
        {0,		0, 	 0, 	0}, //
        {2.0/5, 0, 	 0, 	0},
        {(-2889 + 1428 * sqrt(5))/1024, (3875 - 1620 * sqrt(5))/1024,0, 0},
        {(-3365 + 2094 * sqrt(5))/6040, (-975 - 3046 * sqrt(5))/2552, (467040+203968*sqrt(5))/240845, 0}
};

double b_ralston[4] = {
        (263 + 24 * sqrt(5))/1812,
        (125-1000*sqrt(5))/3828,
        (1024*(3346+1623*sqrt(5)))/5924787,
        (30-4*sqrt(5))/123
};

void ralston(SystemState & x, std::function <void (const SystemState & , SystemState & , const Constants *)> f, double h, const Constants * data) {
    SystemState * k = new SystemState [4];

    for(int i = 0; i < 4; i++) {
        SystemState xt(x);
        for(int j = 0; j < 4; j++){
            xt = xt + k[j] * a_ralston[i][j];
        }
        xt.force = x.force;
        xt.torque = x.torque;
        xt.calculate_derived(data);
        f(xt,k[i], data);
    }
    for(int i = 0; i < 4; i++) {
        //std::cout<<k[i].p.x<<' '<<k[i].p.y<<' '<<k[i].p.z<<' ';
        x = x + k[i] * h * b_ralston[i];
    }
    delete [] k;
    //std::cout<<x.p.x<<' '<<x.p.y<<' '<<x.p.z<<' ';
    //std::cout<<x.p.x<<' '<<x.p.y<<' '<<x.p.z<<std::endl;
}


void start(SystemState & x, double h, double t, const Constants * data)
{
    Ellipse ellipse = {x.C, 1

                       , 2, 1, x.R};

    Plane plane = {Vector(-1,0,1), 1};
    x.force = Vector(0,0,-data->mass * 9.8);
    x.torque = Vector(0,0,0);
    x.L = Vector(0, 0,0.5);
    for(; t < 10; t+=h) {
        ralston( x, f ,h, data);
        process_GS(x.R, x.L);
        t+=h;
        ellipse.R = x.R;
        ellipse.center = x.C;
        std::cout<<x.C;
        ellipse.print_dots();
        //std::cout<<std::endl;
        if(contact(ellipse, plane)) {
            solve(&ellipse, &x, plane, data);
        }
        Drawer().drow();
    }
}
//home/motherlode/SFML-2.5.1/
//1 g++ -c main.cpp -I /home/motherlode/SFML-2.5.1/include
//2 g++ main.o -o sfml-app -lsfml-graphics -lsfml-window -lsfml-system
//3 ./sfml-app


int main()
{
    Drawer d;
    SystemState X;
    double h = 0.01;
    double t = 0;
    Constants data;
    Matrix default_matrix;
    double arr[3][3] = {{1./2, -1./2, sqrt(2)/2},
                        {sqrt(2)/2, sqrt(2)/2, 0},
                        {1./2, -1./2, -sqrt(2)/2}};
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            default_matrix.matrix[i][j] = arr[i][j];
        }
    }
    X.R = default_matrix;
    X.C = Vector(100,0,120);
    //Ellipse ellipse = {X.C, 1, 1, 1, X.R};
    //Vector c(100,100,99.5);
    //std::cout<<ellipse.dot_in_ellipse(c);
    start(X, h, t, &data);
}

