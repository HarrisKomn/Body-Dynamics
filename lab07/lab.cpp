// Include C++ headers
#include <iostream>
#include <string>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>


// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Shader loading utilities and other
#include <common/shader.h>
#include <common/util.h>
#include <common/camera.h>
#include <common/texture.h>

#include "Cube.h"
#include "Sphere.h"
#include "Box.h"
#include "MassSpringDamper.h"
#include "Collision.h"
#include "build\Planet.h"


using namespace std;
using namespace glm;

// Function prototypes
void initialize();
void createContext();
void mainLoop();
void free();
void calculateForces(int i);
void texturesAndDraw(int i, GLuint diffuseSampler , GLuint specularSampler , GLuint ambientSampler, GLuint modelMatrixLocation);
void drawAsteroidsEfe(int i, GLuint diffuseSampler, GLuint specularSampler, GLuint ambientSampler, GLuint modelMatrixLocation);
void returnToInitialConditions(int i);


#define W_WIDTH 1024
#define W_HEIGHT 768
#define TITLE "Lab 07"
#define QUESTION_A2 false
#define QUESTION_B1 false
#define QUESTION_B2 false
#define QUESTION_B3 false
#define efe false

// Global variables
GLFWwindow* window;
Camera* camera;
GLuint shaderProgram;
GLuint projectionMatrixLocation, viewMatrixLocation, modelMatrixLocation;
GLuint diffuceColorSampler, specularColorSampler,ambientColorSampler;
GLuint ambientBlack;
GLuint asterEfeDiffuse, asterEfeSpecular, asterEfeAmbient;

GLuint asteroidDiffuse, asteroidSpecular,asteroidAmbient;
GLuint efeAsteroidDiffuse, efeAsteroidSpecular;

Cube* cube;
Sphere* sphere;
Box* box;
MassSpringDamper* msd;
Planet* asteroid;
Planet* earth;
Planet* moon;
Planet * sun; 
Planet* space;
Planet* efeAsteroid;
Planet* astEfePlanet;

/*__Other planets___*/
Planet* mercury;
Planet* jupiter;
Planet* uranus;
Planet* venus;
Planet* mars;
Planet* saturn;
Planet* neptune;
Planet* pluto;

Planet* planets[14];
Planet* asteroids[2000];

GLuint triangleVAO, triangleVerticesVBO, colourVBO;
GLuint lightSun;
GLuint Colour;
GLuint l_power;
GLuint light_planet;

#define g 9.80665f
#define G float(6.674*pow(10, (-11)))
int number = 0;
int asteroid_number = 0;

void returnToInitialConditions(int i) {
	planets[i]->x = planets[i]->initial_x;
	planets[i]->v = planets[i]->initial_v;
	planets[i]->P = planets[i]->m*planets[i]->initial_v;
	planets[i]->w = planets[i]->initial_w;
	planets[i]->L = planets[i]->I*planets[i]->initial_w;

}

void drawAsteroidsEfe(int i, GLuint diffuseSampler, GLuint specularSampler, GLuint ambientSampler, GLuint modelMatrixLocation) {
	
	glUniformMatrix4fv(modelMatrixLocation, 1, GL_FALSE, &asteroids[i]->modelMatrix[0][0]);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, asteroids[i]->diffuse);
	glUniform1i(diffuseSampler, 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, asteroids[i]->specular);
	glUniform1i(specularSampler, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, asteroids[i]->ambient);
	glUniform1i(ambientSampler, 2);

	asteroids[i]->draw();

}



void texturesAndDraw(int p, GLuint diffuseSampler, GLuint specularSampler, GLuint ambientSampler,  GLuint modelLocation){


	
	glUniformMatrix4fv(modelLocation, 1, GL_FALSE, &planets[p]->modelMatrix[0][0]);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, planets[p]->diffuse);
	glUniform1i(diffuseSampler, 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, planets[p]->specular);
	glUniform1i(specularSampler, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, planets[p]->ambient);
	glUniform1i(ambientSampler, 2);

	planets[p]->draw();

}

void calculateForces(int i) {
	number = i;


	planets[i]->forcing = [&](float t, const vector<float>& y)->vector<float> {

		vector<float> f(6, 0.0f);

		float planets_planeti_x = -G*10.f*(planets[0]->m)*(planets[number]->x.x - planets[0]->x.x) / (pow(distance(planets[number]->x, planets[0]->x), 3));
		float planets_planeti_y = -G*10.f*(planets[0]->m)*(planets[number]->x.y - planets[0]->x.y) / (pow(distance(planets[number]->x, planets[0]->x), 3));
		float planets_planeti_z = -G*10.f*(planets[0]->m)*(planets[number]->x.z - planets[0]->x.z) / (pow(distance(planets[number]->x, planets[0]->x), 3));
			
		for (int j = 1; j < number; j++) {

				planets_planeti_x += -G*(planets[j]->m)*(planets[number]->x.x - planets[j]->x.x) / (pow(distance(planets[number]->x, planets[j]->x), 3));
				planets_planeti_y += -G*(planets[j]->m)*(planets[number]->x.y - planets[j]->x.y) / (pow(distance(planets[number]->x, planets[j]->x), 3));
				planets_planeti_z += -G*(planets[j]->m)*(planets[number]->x.z - planets[j]->x.z) / (pow(distance(planets[number]->x, planets[j]->x), 3));
			}

for (int j = number + 1; j < 10; j++) {

	planets_planeti_x += -G*(planets[j]->m)*(planets[number]->x.x - planets[j]->x.x) / (pow(distance(planets[number]->x, planets[j]->x), 3));
	planets_planeti_y += -G*(planets[j]->m)*(planets[number]->x.y - planets[j]->x.y) / (pow(distance(planets[number]->x, planets[j]->x), 3));
	planets_planeti_z += -G*(planets[j]->m)*(planets[number]->x.z - planets[j]->x.z) / (pow(distance(planets[number]->x, planets[j]->x), 3));
}

f[0] = (planets[number]->m)*planets_planeti_x;
f[1] = (planets[number]->m)*planets_planeti_y;
f[2] = (planets[number]->m)*planets_planeti_z;

return f;

	};

}

void createContext() {
	shaderProgram = loadShaders(
		"StandardShading.vertexshader",
		"StandardShading.fragmentshader");

	projectionMatrixLocation = glGetUniformLocation(shaderProgram, "P");
	viewMatrixLocation = glGetUniformLocation(shaderProgram, "V");
	modelMatrixLocation = glGetUniformLocation(shaderProgram, "M");
	diffuceColorSampler = glGetUniformLocation(shaderProgram, "diffuseColorSampler");
	specularColorSampler = glGetUniformLocation(shaderProgram, "specularColorSampler");
	ambientColorSampler = glGetUniformLocation(shaderProgram, "ambientColorSampler");


	Colour = glGetUniformLocation(shaderProgram, "useColour");
	lightSun= glGetUniformLocation(shaderProgram, "lightSun");
	l_power= glGetUniformLocation(shaderProgram, "power");
	light_planet= glGetUniformLocation(shaderProgram, "l_planet");

	box = new Box(0);

	  

	//space area
	space = new Planet(vec3(0, 0, 0), vec3(0, 0, 0), -0.0001f*0.000199f*vec3(0, 1, 0), 0.01/*0.00002*/, float(1.989f*pow(10, 29)), "models/sol.obj", 0., 0.f, 109 * 6378.1 *pow(10, 3));

	space->diffuse = loadSOIL("models/stars.jpg");
	space->specular = space->diffuse;
	space->ambient = space->diffuse;

	//sun
	sun = new Planet(vec3(float(2.499573*pow(10, 9)), 0, 0), vec3(0, 0, 0), -0.001f*0.000199f*vec3(0, 1, 0), 0.006/*0.00002*/, float(1.989f*pow(10, 29)), "models/earth.obj", 0., 0.f, 109 * 6378.1 *pow(10, 3));

	sun->diffuse = loadSOIL("models/2k_sun.jpg");
	sun->specular = sun->diffuse;
	sun->ambient = sun->diffuse;

	//ambientBlack = loadSOIL("models/black.jpg");

	//earth
//	if (QUESTION_B1) {
//		earth = new Planet(vec3(float(149.598023*pow(10, 9)), 0, 0), vec3(0, 0, 29771.f), (-465.1f/ float(6378.1f *pow(10, 3)))*vec3(0, 1, 0), 0.002, float(5.97f*pow(10, 24)), "models/earth.obj", 3.14f, 0.409, 6378.1f *pow(10, 3));
//	}
//	else if (QUESTION_A2) {
//		earth = new Planet(vec3(float(149.598023*pow(10, 9)), 0, 0), vec3(0, 0, 29771.f), 1000.f/**(-465.1f / float(6378.1f *pow(10, 3)))*/*vec3(0, 1, 0), 0.002, 5.f/2.f/*float(5.97f*pow(10, 24))*/, "models/earth.obj", 3.14f, 0.409,1.f/* 6378.1f *pow(10, 3)*/);
//
//		//earth = new Planet(vec3(float(149.598023*pow(10, 9)), 0, 0), vec3(0, 0, 29771.f), vec3(10, 10, 10), 0.0001, 10., "models/earth.obj", 2./*3.14/2.*/, 0.0/*vec3(sin(0.409052), cos(0.409052), 0)*/, 1. /*6.3781e6*/);
//
//	}

	earth = new Planet(vec3(float(149.598023*pow(10, 9)), 0, 0), vec3(0, 0, -29771.f), (-465.1f/ float(6378.1f *pow(10, 3)))*vec3(0, 1, 0), 0.002, float(5.97f*pow(10, 24)), "models/earth.obj", 3.14f, 0.409, 6378.1f *pow(10, 3));

	earth->diffuse = loadSOIL("models/earth_diffuse.jpg");
	earth->specular = loadSOIL("models/earthspec1k.jpg");
	earth->ambient = loadSOIL("models/night_lights.jpg");


	//moon
//	if (QUESTION_B1) {
//		moon = new Planet(vec3(float((149.598023 + 0.3626)*pow(10, 9)), 0, 0), vec3(0., 0, 1022.f/*964.f*/), (4.627f/ (1738.1f*float(pow(10, 3))))*vec3(0, 1, 0), 0.005, float(7.342f*pow(10, 22)), "models/moon.obj", 0., 0.1167f, 1738.1f*pow(10,3));
//	}
//	else if (QUESTION_A2) {
//		moon = new Planet(vec3(float((149.598023 + 0.3626)*pow(10, 9)), 0, 0), vec3(0., 0, 1022.f/*964.f*/), vec3(0, 1, 0), 0.005, float(7.342f*pow(10, 22)), "models/moon.obj", 0., 0.f, 0);
//
//		//moon = new Planet(vec3(float((149.598023 + 0.3626)*pow(10, 9)), 0, 0), vec3(0, 0, 964.f), vec3(0, 0, 0), 0.0001, float(7.342f*pow(10, 22)), "models/moon.obj", 0., 0.f, 1.737e3);
//	}
//	//moon = new Planet(vec3(float((149.598023 + 0.3626)*pow(10, 9)), 0, 0), vec3(0., 0, 1022.f/*964.f*/), vec3(0, 0, 0), 0.005, float(7.342f*pow(10, 22)), "models/moon.obj", 0., 0.f, 0);

	moon = new Planet(vec3(float((149.598023 + 0.3626)*pow(10, 9)), 0, 0), vec3(0., 0, -1022.f/*964.f*/), (-4.627f/ (1738.1f*float(pow(10, 3))))*vec3(0, 1, 0), 0.005, float(7.342f*pow(10, 22)), "models/moon.obj", 0., 0.1167f, 1738.1f*pow(10,3));


	moon->diffuse = loadSOIL("models/moon.jpg");
	moon->specular = moon->diffuse;
	moon->ambient = moon->diffuse;

	planets[0] = sun;
	planets[1] = earth;

	planets[10] = moon;
	planets[11] = space;
	//asteroid

	

	asteroid = new Planet(vec3(30.f, 0, 0), vec3(0.f, 0, 0.f), vec3(0, 0, 0), 0.003, float(7.342f*pow(10, 15)), "models/moon.obj", 0., 0.f, 1.737e3);

	asteroid->diffuse= loadSOIL("models/makemake.jpg");
	asteroid->specular = asteroid->diffuse;
	asteroid->ambient = asteroid->diffuse;
	
	
	efeAsteroid = new Planet(vec3(0, 0, 10), vec3(0.f, 0, 0.f), vec3(0, 0, 0), 0.0005, float(7.342f*pow(10, 22)), "models/asteroid_space.obj", 0., 0.f, 1.737e3);

	efeAsteroid->diffuse = loadSOIL("models/asteroid_space.jpg");
	efeAsteroid->specular = efeAsteroid->diffuse;
	efeAsteroid->ambient = efeAsteroid->diffuse;



	astEfePlanet = new Planet(vec3(0, 0, 0), vec3(0.f, 0, 0.f), vec3(0, 0, 0), 0.000001, float(7.342f*pow(10, 0)), "models/asteroid_space.obj", 0., 0.f, 1.737e0);
	asterEfeDiffuse = loadSOIL("models/grey.jpg");
	asterEfeSpecular = asterEfeDiffuse;
	asterEfeAmbient = asterEfeDiffuse;

	planets[12] = asteroid;
	planets[13] = efeAsteroid;
	
		for (int k=0; k < 1000; k++) {
		
			asteroids[k]= new Planet(vec3(float((rand() % 160) -80),  float((rand() % 160) -80), float((rand() % 160) -80)), vec3(0.f, 0, 0.f), vec3(0, 0, 0), 0.0008, float(7.342f*pow(10, 0)), "models/soft asteroid.obj", 0., 0.f, 1.737e0);
			float theta = acos(dot(normalize(vec3(asteroids[k]->x.x, 0, asteroids[k]->x.z)), vec3(1, 0, 0)));
			//float theta = acos((dot(normalize(projPlaneVert), normalize(x_vex))));

			vec3 check_normal = cross(normalize(vec3(asteroids[k]->x.x,0, asteroids[k]->x.z)), vec3(1,0,0));
	
			if (dot(vec3(0,1,0), check_normal)<0) {

				theta = -theta;
			}
			asteroids[k]->anim_theta = theta;
			asteroids[k]->anim_r = length(vec3(0, 0, 0)- vec3(asteroids[k]->x.x,0, asteroids[k]->x.z));
			asteroids[k]->diffuse = asterEfeDiffuse;
			asteroids[k]->specular = asterEfeDiffuse;
			asteroids[k]->ambient = asterEfeDiffuse;

		}
	

	

		mercury = new Planet(vec3(sun->x.x+ float((46.001200)*pow(10, 9)), 0, 0), -47362.f*vec3(0, 0, 1), (-3.026f/float(earth->radius*0.3829))*vec3(0, 1, 0), earth->l*2.*0.3829 / 10.f, float(3.3011*pow(10, 23)), "models/earth.obj", 0., 0.0006, earth->radius*0.3829);
		mercury->diffuse = loadSOIL("models/mercury/mercury_diffuse.jpg");
		mercury->specular = mercury->diffuse;
		mercury->ambient = mercury->diffuse;


		jupiter = new Planet(vec3(sun->x.x + float((740.52)*pow(10, 9)), 0, 0), -13070.f*vec3(0, 0, 1), (-12600.f / float(earth->radius*10.517))*vec3(0, 1, 0), earth->l*2.*10.f/10.f, float(1.8982*pow(10, 27)), "models/earth.obj", 0., 0.055, earth->radius*10.517);
		
		jupiter->diffuse = loadSOIL("models/jupiter/jupiter_diffuse.jpg");
		jupiter->specular = jupiter->diffuse;
		jupiter->ambient = jupiter->diffuse;
		
		
		uranus = new Planet(vec3(sun->x.x + float((2748.938461)*pow(10, 9)), 0, 0), -6800.f*vec3(0, 0, 1), (-2590.f / float(earth->radius*4.007))*vec3(0, 1, 0), earth->l*2.*4.f / 10.f, float(8.68*pow(10, 25)), "models/earth.obj", 0., 1.706, earth->radius*4.007);
		
		uranus->diffuse = loadSOIL("models/uranus/uranus_diffuse.jpg");
		uranus->specular = uranus->diffuse;
		uranus->ambient = uranus->diffuse;

		
		venus = new Planet(vec3(sun->x.x + float((107.476259)*pow(10, 9)), 0, 0), -35020.f*vec3(0, 0, 0), (-1.81f / float(earth->radius*0.9499))*vec3(0, 1, 0), earth->l*2.*0.95f / 10.f, float(4.8685*pow(10, 24)), "models/earth.obj", 0., 3.096, earth->radius*0.9499);
		
		venus->diffuse = loadSOIL("models/venus/venus_diffuse.jpg");
		venus->specular = venus->diffuse;
		venus->ambient = venus->diffuse;
		
		mars = new Planet(vec3(sun->x.x + float((206.669)*pow(10, 9)), 0, 0), -24077.f*vec3(0, 0, 1), (-241.f / float(earth->radius*0.533))*vec3(0, 1, 0), earth->l*2.*0.533f / 10.f, float(6.48*pow(10, 23)), "models/earth.obj", 0., 0.439, earth->radius*0.533);
		
		mars->diffuse = loadSOIL("models/mars/mars_diffuse.jpg");
		mars->specular = mars->diffuse;
		mars->ambient = mars->diffuse;
		
		
		saturn = new Planet(vec3(sun->x.x + float((1353.572956)*pow(10, 9)), 0, 0), -9680.f*vec3(0, 0, 1), (-9870.f / float(earth->radius*9.499))*vec3(0, 1, 0), earth->l*2.f*9.5f / 10.f, float(5.6846*pow(10, 26)), "models/earth.obj", 0.,0.467, earth->radius*9.499);
		saturn->diffuse = loadSOIL("models/saturn/saturn_diffuse.jpg");
		saturn->specular = saturn->diffuse;
		saturn->ambient = saturn->diffuse;
		
		neptune = new Planet(vec3(sun->x.x + float((4452.940833)*pow(10, 9)), 0, 0), -54300.f*vec3(0, 0, 1), (-2680.f / float(earth->radius*3.833))*vec3(0, 1, 0), earth->l*2.*3.9f / 10.f, float(1.0243*pow(10, 26)), "models/earth.obj", 0., 0.494, earth->radius*3.833);
		neptune->diffuse = loadSOIL("models/neptune/neptune_diffuse.jpg");
		neptune->specular = neptune->diffuse;
		neptune->ambient = neptune->diffuse;

		
		pluto = new Planet(vec3(sun->x.x + float((4436.824613)*pow(10, 9)), 0, 0), -46700.f*vec3(0, 0, 1), (-13.f / float(earth->radius*0.1868))*vec3(0, 1, 0), 10.f*earth->l*2.*0.19f / 10.f, float(1.305*pow(10, 22)), "models/earth.obj", 0., 2.138, earth->radius*0.1868);
		pluto->diffuse = loadSOIL("models/pluto/pluto_diffuse.jpg");
		pluto->specular = pluto->diffuse;
		pluto->ambient = pluto->diffuse;
		
		
		planets[2] = mercury;
		planets[3] = jupiter;
		planets[4] = uranus;
		planets[5] = venus;
		planets[6] = mars;
		planets[7] = saturn;
		planets[8] = neptune;
		planets[9] = pluto;
		
	


	//AXIS
	glGenVertexArrays(1, &triangleVAO);
	glBindVertexArray(triangleVAO);

	glGenBuffers(1, &triangleVerticesVBO);    
	glBindBuffer(GL_ARRAY_BUFFER, triangleVerticesVBO);
	const GLfloat tirangleVertices[] =
	{
		-500, 0, 0.0,
		500, 0, 0.0,
		0, 0, -500,
		0,0,500
	};

	glBufferData(GL_ARRAY_BUFFER, sizeof(tirangleVertices),
		&tirangleVertices[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	const GLfloat colours[] = {
		1.f,0.f,0.f,
		0.f,1.f,0.f,
		0.f,0.f,1.f,
		0.f,0.f,1.f

	};
	glGenBuffers(1, &colourVBO);    
	glBindBuffer(GL_ARRAY_BUFFER, colourVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(colours),
		&colours[0], GL_STATIC_DRAW);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(3);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}

void free() {

    glDeleteProgram(shaderProgram);
    glfwTerminate();
}

void mainLoop() {

	bool col = false;
	bool happened_col = false;
    float t = glfwGetTime();
    vec3 lightPos = vec3(0,0, 0);
    camera->position = glm::vec3(box->size / 2, box->size / 2, 20);
	bool START = false;
	float step = 500;
	int useColour = 0;
	int lightS = 0;
	float lig_pow = 0.f;
	mat4 modelMatrix = glm::mat4(1.0);
	mat4 M;
	
	/*_____________Animation_______________*/
	//Elliptical orbit: alpha,beta,gama
	float alpha = 14.9598023f/2.f ;//2*pow(10,10)
	float gama = 0.2499573f/2.f;//2*pow(10,10)
	float beta = sqrt(pow(alpha, 2) - pow(gama, 2));

	//Earth
	float ì_sun = 1.327*pow(10, 20);
	vec3 earth_pos = vec3(alpha, 0, 0);

	float earth_r = 0.0;
	float earth_theta = 0.0;
	float earth_own_angle = 0.f;

	//float earth_A = 3.14*(pow(10, 7)) / 60.0f;

	float tilted = 0.409052;//rad
	//Sun
	vec3 sun_pos = vec3(gama, 0, 0);
	
	//Moon
	float ì_earth = 3.98*pow(10, 14);
	vec3 moon_pos = vec3(0, 0, 0);

	float moon_r = 0.0;
	float moon_theta = 0.0;
	float moon_own_angle = 0.f;

	float moon_center_x = 0.f;
	float moon_center_z = 0.f;

	//elipse a,b,c

	float moon_alpha = 0.03626f;//pow(10,10)
	float moon_gama = moon_alpha-0.03565f;//pow(10,10)
	float moon_beta = sqrt(pow(moon_alpha, 2) - pow(moon_gama, 2));
	/*synchronization*/
	float Duration_earth = 60;
	float earth_A = 0;

	float P_earth = 60.f;
	float P_tilt = 26000.f;
	float P_moon = P_earth*(27.322f / 365.f);
	float theta_tilt = 0.f;

	float l_planet=0.;
	float fast_forward = 1.f;

	bool FAST_FORWARD = false;

	bool ANIMATION = false;
	bool THREE_BODY = false;
	bool COLLISION = false;
	bool N_BODY = false;
	bool RESET = true;
	
	bool showAxis = true;
	bool efeAnimation = true;

	do {
        // calculate dt
        float currentTime = glfwGetTime();
        float dt = currentTime - t;


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);


		//Fast-Forward

		if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
			fast_forward = 300.f;
		}
		if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS) {
			fast_forward = 1.f;
		}

		//efe
		if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
			efeAnimation = false;
		}
		if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
			efeAnimation = true;
		}

		if (RESET == true) {
			//Animation
			if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
				earth->update(t, dt, vec3(0, 0, 0), vec3(0, 0, 0), false, true, false, false);
				moon->update(t, dt, vec3(0, 0, 0), (earth->x)*float(pow(10, 10)), true, true, false, true);
				ANIMATION = true;
				RESET = false;
			}

			//Simulation
			//3_Body
			if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
				THREE_BODY = true;
				RESET = false;

			}

			//Collision
			if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
				THREE_BODY = true;
				COLLISION = true;
				RESET = false;

			}

			//N_body
			if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) {
				N_BODY = true;
				RESET = false;

			}
		}
			//Start-Stop
			if (ANIMATION || THREE_BODY || COLLISION || N_BODY) {
				if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS) {
					START = true;
				}

				if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
					START = false;
				}
			}

			//Reset
			if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {

				ANIMATION = false;
				THREE_BODY = false;
				COLLISION = false;
				N_BODY = false;
				START = false;
				RESET = true;
				//for (int m = 1; m < 11; m++) {
				//	returnToInitialConditions(m);
				//}
				returnToInitialConditions(1);
				returnToInitialConditions(10);
				RESET = true;
				t = 0.f;
			}
		
        // camera
        camera->update();
        mat4 projectionMatrix = camera->projectionMatrix;
        mat4 viewMatrix = camera->viewMatrix;

        glUniformMatrix4fv(viewMatrixLocation, 1, GL_FALSE, &viewMatrix[0][0]);
        glUniformMatrix4fv(projectionMatrixLocation, 1, GL_FALSE, &projectionMatrix[0][0]);

		//*/
		//AXIS
		if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
			showAxis = !showAxis;
		}
		if (showAxis) {
			useColour = 1;
			glBindVertexArray(triangleVAO);
			M = modelMatrix;
			glUniformMatrix4fv(modelMatrixLocation, 1, GL_FALSE, &M[0][0]);
			glUniform1i(Colour, useColour);
			glDrawArrays(GL_LINES, 0, 4);
			useColour = 0;
			glUniform1i(Colour, useColour);
		}
		//space	
		space->update(t, step,vec3(0,0,0), vec3(0, 0, 0),false,true,true,false);
		texturesAndDraw(11, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

		
		// sun	
		sun->update(t, step,vec3(0,0,0), vec3(0,0,0),false,true,true,false);
		lightS = 1;
		glUniform1i(lightSun, lightS);	
		glUniform1f(l_power, lig_pow);
		lig_pow += 0.01f;
		texturesAndDraw(0, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
		lightS = 0;
		glUniform1i(lightSun, lightS);
		
		l_planet = 400.f;
		glUniform1f(light_planet, l_planet);

		if (ANIMATION) {
		
			//Earth
			earth_r = ((alpha*alpha) - (gama*gama)) / (alpha - gama*cos(earth_theta));

			earth_pos.x = earth_r*cos(earth_theta);
			earth_pos.z = earth_r*sin(earth_theta);

			earth->x = earth_pos;
			earth->own_angle = earth_own_angle;

			earth->P = earth->m*vec3(0, 0, 0);
			earth->w = vec3(0, 0, 0);
			earth->L = earth->I*earth->w;

			float angle_sp = 2 * 3.14*alpha*beta / ((P_earth / (fast_forward))*pow(distance(earth_pos, sun_pos), 2));
			if (START) {

				earth->w = vec3(0, 1, 0)*float(-2 * 3.14f / (P_earth / (fast_forward * 365)));
				earth->L = earth->w;
		
				//precesion
				theta_tilt += (2.f*3.14f / ((P_earth / fast_forward)*P_tilt))*dt;

				earth_theta = earth_theta - angle_sp*(dt);
				if (earth_theta < -(300.f*3.14f * 2.f)) {
					START = false;
				}
			}

			earth->update(t, dt, vec3(0, 0, 0), vec3(0, 0, 0), false, true, false, false);
		
			mat4 extraRotate=rotate(mat4(), theta_tilt, vec3(0, 1, 0));
			earth->modelMatrix = (translate(mat4(), earth->normalized_pos)) *extraRotate*(earth->RotationalMatrix)* (earth->ScaleMatrix);		
			texturesAndDraw(1, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			//Moon
			moon_r = ((moon_alpha*moon_alpha) - (moon_gama*moon_gama)) / (moon_alpha - moon_gama*cos(moon_theta));

			moon_pos.x = moon_r*cos(moon_theta);
			moon_pos.z = moon_r*sin(moon_theta);


			moon_center_x = ((abs(distance(earth_pos, vec3(0, 0, 0))) - moon_gama) / (abs(distance(earth_pos, vec3(0, 0, 0)))))*earth_pos.x;
			moon_center_z = ((distance(earth_pos, vec3(0, 0, 0)) - moon_gama) / (distance(earth_pos, vec3(0, 0, 0))))*earth_pos.z;

			moon->x.x = moon_center_x + moon_r*cos(moon_theta);;
			moon->x.z = moon_center_z + moon_r*sin(moon_theta);
			moon->P = moon->m*vec3(0, 0, 0);
			moon->w= vec3(0, 0, 0);
			moon->L = moon->I*moon->w;
	
			float moon_sp = 2 * 3.14*moon_alpha*moon_beta / (P_moon / (fast_forward)*pow(distance(earth_pos, moon->x), 2));

			if (START) {
				moon->w = vec3(0, 1, 0)*float(-2 * 3.14f /	(P_moon / (fast_forward)));
				moon->L = earth->w;
				moon_theta = moon_theta - moon_sp*dt;
			}
			moon->update(t, dt, vec3(0, 0, 0), (earth->x)*float(pow(10, 10)), true, true, false, true);
			//moon->modelMatrix = rotate(mat4(), 0.09f, vec3(0, 0, 1))*moon->modelMatrix;

			texturesAndDraw(10, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

		}

		

		if (THREE_BODY) {
			// earth
			//*/
			earth->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			if (START) {
			

			earth->forcing = [&](float t, const vector<float>& y)->vector<float> {
				vector<float> f(6, 0.0f);
				float a_moon_earth_0 = -G*(moon->m)*(earth->x.x - moon->x.x) / (pow(distance(earth->x, moon->x), 3));
				float a_sun_earth_0 = -G*10.f*(sun->m)*(earth->x.x - sun->x.x) / (pow(distance(earth->x, sun->x), 3));

				f[0] = (earth->m)*a_moon_earth_0
					+ (earth->m)*a_sun_earth_0;

				float a_moon_earth_1 = -G*(moon->m)*(earth->x.y - moon->x.y) / (pow(distance(earth->x, moon->x), 3));
				float a_sun_earth_1 = -G*10.f*(sun->m)*(earth->x.y - sun->x.y) / (pow(distance(earth->x, sun->x), 3));

				f[1] = (earth->m)*a_moon_earth_1
					+ (earth->m)*a_sun_earth_1;

				float a_moon_earth_2 = -G*(moon->m)*(earth->x.z - moon->x.z) / (pow(distance(earth->x, moon->x), 3));
				float a_sun_earth_2 = -G*10.f*(sun->m)*(earth->x.z - sun->x.z) / (pow(distance(earth->x, sun->x), 3));

				f[2] = (earth->m)*a_moon_earth_2
					+ (earth->m)*a_sun_earth_2;

				return f;
			};

			earth->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);

			}
			texturesAndDraw(1, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
			
			//*/

			// moon
			//*/

			vec3 earth_vel = earth->P / earth->m;
			bool energies = true;
			bool coll = false;

			moon->update(t, step, earth_vel, earth->x, true, false, true, energies);
			
			if (START) {
			moon->forcing = [&](float t, const vector<float>& y)->vector<float> {
				vector<float> f(6, 0.0f);
				float a_earth_moon_0 = -G*(earth->m)*(moon->x.x - earth->x.x) / (pow(distance(moon->x, earth->x), 3));
				float a_sun_moon_0 = -G*10.f*(sun->m)*(moon->x.x - sun->x.x) / (pow(distance(moon->x, sun->x), 3));

				f[0] = (moon->m)*a_earth_moon_0
					+ (moon->m)*a_sun_moon_0;

				float a_earth_moon_1 = -G*(earth->m)*(moon->x.y - earth->x.y) / (pow(distance(moon->x, earth->x), 3));
				float a_sun_moon_1 = -G*10.f*(sun->m)*(moon->x.y - sun->x.y) / (pow(distance(moon->x, sun->x), 3));

				f[1] = (moon->m)*a_earth_moon_1
					+ (moon->m)*a_sun_moon_1;

				float a_earth_moon_2 = -G*(earth->m)*(moon->x.z - earth->x.z) / (pow(distance(moon->x, earth->x), 3));
				float a_sun_moon_2 = -G*10.f*(sun->m)*(moon->x.z - sun->x.z) / (pow(distance(moon->x, sun->x), 3));

				f[2] = (moon->m)*a_earth_moon_2
					+ (moon->m)*a_sun_moon_2;


				return f;
			};

			moon->update(t, step, earth_vel, earth->x, true, true, true, energies);

			}
			texturesAndDraw(10, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			//*/

			if (COLLISION) {
				// asteroid
				//*/
				//asteroid->P = (asteroid->m)*100000.f *normalize(moon->x - asteroid->x);
				if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
					if (col) {
						asteroid->x = vec3(10.f, 0, -20);
						asteroid->P = asteroid->m*vec3(0, 0, 0);
						asteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, false, false);

						col = false;
					}
				}
				asteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, false, false);

				if(START){
				if (!col) {

					asteroid->forcing = [&](float t, const vector<float>& y)->vector<float> {
						vector<float> f(6, 0.0f);
						float a_asteroid_moon_0 = -1000000000 * (asteroid->x.x - moon->new_pos_moon.x) / (pow(distance(asteroid->x, moon->new_pos_moon), 3));

						f[0] = a_asteroid_moon_0;



						float a_asteroid_moon_1 = -1000000000 * (asteroid->x.y - moon->new_pos_moon.y) / (pow(distance(asteroid->x, moon->new_pos_moon), 3));

						f[1] = a_asteroid_moon_1;


						float a_asteroid_moon_2 = -1000000000 * (asteroid->x.z - moon->new_pos_moon.z) / (pow(distance(asteroid->x, moon->new_pos_moon), 3));

						f[2] = a_asteroid_moon_2;
						
						return f;
					};

					asteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, false, false);

				}
				else {
					asteroid->x = moon->new_pos_moon;
					asteroid->P = vec3(0, 0, 0);
					asteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, false, false);

				}
			}
				//texturesAndDraw(12, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

				//*/
				//efe

				//*/


				if (!col) {
					efeAsteroid->x = asteroid->x + (-0.15f - asteroid->l)*normalize(asteroid->x - moon->new_pos_moon);
					efeAsteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, false, false);
					texturesAndDraw(13, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

				}
				else {

					efeAsteroid->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, false, false);

				}

				//texturesAndDraw(13, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

				//*/


				if (!col) {
					col = handlePlanetPlanetCollision(*moon, *asteroid);
				}
			}
		}
		
		if (N_BODY) {
			
			/*mercury*/
			if (START) {
				calculateForces(2);
				mercury->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				mercury->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(2, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);			

			l_planet = 200.f;
			glUniform1f(light_planet, l_planet);

			/*venus*/
			if (START) {
				calculateForces(5);
				venus->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				venus->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(5, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
			
			/*earth*/
			if (START) {
				calculateForces(1);
				earth->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				earth->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(1, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			/*moon*/
			if (START) {
				calculateForces(10);
				moon->update(t, step, earth->v, earth->x, true, true, true, true);
			}
			else {
				moon->update(t, step, earth->v, earth->x, true, false, true, true);

			}
			texturesAndDraw(10, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			/*mars*/
			if (START) {
				calculateForces(6);
				mars->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				mars->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);
			}
				texturesAndDraw(6, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			l_planet = 3500.f;
			glUniform1f(light_planet, l_planet);

			/*jupiter*/
			if (START) {
				calculateForces(3);
				jupiter->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				jupiter->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(3, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			/*saturn*/
			if (START) {
				calculateForces(7);
				saturn->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				saturn->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(7, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			l_planet = 5000.f;
			glUniform1f(light_planet, l_planet);

			/*uranus*/
			if (START) {
				calculateForces(4);
				uranus->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				uranus->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(4, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
	
			/*neptune*/
			if (START) {
				calculateForces(8);
				neptune->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				neptune->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);


			}
			texturesAndDraw(8, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);

			/*pluto*/
			if (START) {
				calculateForces(9);
				pluto->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, true, true, false);
			}
			else {
				pluto->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, true, false);

			}
			texturesAndDraw(9, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
			

		}
		
		l_planet = 500.f;
		glUniform1f(light_planet, l_planet);
		

		if (efeAnimation) {
			for (int ast = 0; ast < 1000; ast++) {
				//asteroidForces(ast);

				asteroids[ast]->x = vec3((asteroids[ast]->anim_r)*cos(asteroids[ast]->anim_theta), asteroids[ast]->x.y, (asteroids[ast]->anim_r)*sin(asteroids[ast]->anim_theta));
				asteroids[ast]->anim_theta += 0.0001;
				asteroids[ast]->update(t, step, vec3(0, 0, 0), vec3(0, 0, 0), false, false, false, false);
				drawAsteroidsEfe(ast, diffuceColorSampler, specularColorSampler, ambientColorSampler, modelMatrixLocation);
			}
		}

		t += dt;

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
        glfwWindowShouldClose(window) == 0);
}

void initialize() {
    // Initialize GLFW
    if (!glfwInit()) {
        throw runtime_error("Failed to initialize GLFW\n");
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    window = glfwCreateWindow(W_WIDTH, W_HEIGHT, TITLE, NULL, NULL);
    if (window == NULL) {
        glfwTerminate();
        throw runtime_error(string(string("Failed to open GLFW window.") +
            " If you have an Intel GPU, they are not 3.3 compatible." +
            "Try the 2.1 version.\n"));
    }
    glfwMakeContextCurrent(window);

    // Start GLEW extension handler
    glewExperimental = GL_TRUE;

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        glfwTerminate();
        throw runtime_error("Failed to initialize GLEW\n");
    }

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    // Hide the mouse and enable unlimited movement
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Set the mouse at the center of the screen
    glfwPollEvents();
    glfwSetCursorPos(window, W_WIDTH / 2, W_HEIGHT / 2);

    // Gray background color
    glClearColor(0.5f, 0.5f, 0.5f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);

    // Cull triangles which normal is not towards the camera
    //glEnable(GL_CULL_FACE);
    //glFrontFace(GL_CW);
    //glFrontFace(GL_CCW);

    // enable point size when drawing points
    glEnable(GL_PROGRAM_POINT_SIZE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// enable textures
	//glEnable(GL_TEXTURE_2D);
    // Log
    logGLParameters();

    // Create camera
    camera = new Camera(window);
}

int main(void) {
    try {
        initialize();
        createContext();
        mainLoop();
        free();
    } catch (exception& ex) {
        cout << ex.what() << endl;
        getchar();
        free();
        return -1;
    }

    return 0;
}
