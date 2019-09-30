#include "Collision.h"
#include "Box.h"
#include "Sphere.h"
#include "build\Planet.h"
#include <iostream>
#include <string>
// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>


// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace glm;
/*________asteroid mass_factor_____*/
float asteroid_mass_factor=1e6f;

void handleBoxSphereCollision(Box& box, Sphere& sphere);
bool checkForBoxSphereCollision(glm::vec3& pos, const float& r,
	const float& size, glm::vec3& n);

bool handlePlanetPlanetCollision(Planet& planet1, Planet& planet2);

bool checkForPlanetPlanetCollision(glm::vec3& pos1, const float& r1,
	glm::vec3& pos2, const float& r2, glm::vec3& n1, glm::vec3& n2);

//bool collision = true;

void handleBoxSphereCollision(Box& box, Sphere& sphere) {
    vec3 n;
    if (checkForBoxSphereCollision(sphere.x, sphere.r, box.size, n)) {
        // Task 2b: define the velocity of the sphere after the collision
    }
}

bool checkForBoxSphereCollision(vec3& pos, const float& r, const float& size, vec3& n) {
    if (pos.x - r <= 0) {
        //correction
        float dis = -(pos.x - r);
        pos = pos + vec3(dis, 0, 0);

        n = vec3(-1, 0, 0);
    } else if (pos.x + r >= size) {
        //correction
        float dis = size - (pos.x + r);
        pos = pos + vec3(dis, 0, 0);

        n = vec3(1, 0, 0);
    } else if (pos.y - r <= 0) {
        //correction
        float dis = -(pos.y - r);
        pos = pos + vec3(0, dis, 0);

        n = vec3(0, -1, 0);
    } else if (pos.y + r >= size) {
        //correction
        float dis = size - (pos.y + r);
        pos = pos + vec3(0, dis, 0);

        n = vec3(0, 1, 0);
    } else if (pos.z - r <= 0) {
        //correction
        float dis = -(pos.z - r);
        pos = pos + vec3(0, 0, dis);

        n = vec3(0, 0, -1);
    } else if (pos.z + r >= size) {
        //correction
        float dis = size - (pos.z + r);
        pos = pos + vec3(0, 0, dis);

        n = vec3(0, 0, 1);
    } else {
        return false;
    }
    return true;
}


bool handlePlanetPlanetCollision(Planet& planet1, Planet& planet2) {
	vec3 n1;
	vec3 n2;
	if (checkForPlanetPlanetCollision(planet1.new_pos_moon, planet1.l, planet2.x, planet2.l, n1, n2)) {

		/*__plastic collision__*/

		vec3 vel_1_princ = planet1.v;
		vec3 vel_2_princ = planet2.v;

		//planet-1 (moon)
		
		planet1.v = (planet1.m*vel_1_princ) / (planet1.m + asteroid_mass_factor*planet2.m) + (asteroid_mass_factor*planet2.m*vel_2_princ) / (planet1.m + asteroid_mass_factor*planet2.m);
		planet1.m = planet1.m + asteroid_mass_factor*planet2.m;

		planet1.P = 0.97f*planet1.v*planet1.m;

		//planet-2 (asteroid)
		
		planet2.v = planet1.v;
		//planet2.m = planet1.m + asteroid_mass_factor* planet2.m;
		planet2.P = 0.97f*planet2.v*planet2.m;

		//collision = false;
		return true;
	}
	//return collision;
	return false;

}

bool checkForPlanetPlanetCollision(vec3& pos1, const float& r1, vec3& pos2, const float& r2, vec3& n1, vec3& n2) {
	
	
		if (abs(distance(pos2, pos1) <= r2 + r1 + 0.6)) {
			//correction
			float dis = 0.6 + r2 + r1 - (distance(pos2, pos1));
			pos1 = pos1 - normalize(pos2 - pos1)*dis / 2.f;
			pos2 = pos2 + normalize(pos2 - pos1)*dis / 2.f;
			n1 = -normalize(pos2 - pos1);
			n2 = normalize(pos2 - pos1);
			return true;
		}

		else {
			return false;
		}
	
	return false;
}