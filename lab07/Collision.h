#ifndef COLLISION_H
#define COLLISION_H

#include <glm/glm.hpp>
#include "build\Planet.h"
class Box;
class Sphere;
void handleBoxSphereCollision(Box& box, Sphere& sphere);
bool handlePlanetPlanetCollision(Planet& planet1, Planet& planet2);

#endif
