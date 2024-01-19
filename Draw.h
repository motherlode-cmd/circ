//
// Created by motherlode on 19.01.23.
//

#ifndef UNTITLED4_DRAW_H
#define UNTITLED4_DRAW_H
#include "System.h"
#include "SFML/Graphics.hpp"
typedef struct Drawer {
public:
    void drow() {
        sf::RenderWindow window(sf::VideoMode(200, 200), "SFML works!");
        sf::CircleShape shape(100.f);
        shape.setFillColor(sf::Color::Green);

        while (window.isOpen())
        {
            sf::Event event;
            while (window.pollEvent(event))
            {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();
            window.draw(shape);
            window.display();
        }

    }
}Drawer;
#endif //UNTITLED4_DRAW_H
