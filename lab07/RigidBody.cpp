#include "RigidBody.h"
#include <common/ModelLoader.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

using namespace glm;
GLuint  textureSun;
GLuint textureSampler;
#define G float(6.674*pow(10, (-11)))
float M_earth = 5.97*pow(10, 24);
float M_sun = 1.989*pow(10, 30);
float M_moon = float(7.342f*pow(10, 22));

RigidBody::RigidBody() {

	//rigidBody = new Drawable();
	l = 1;
    m = 1;
	//extraRotation = glm::mat4(1.);
    x = v = w = P = L = vec3(0, 0, 0);

#ifdef USE_QUATERNIONS
    q = quat(0, vec3(1, 0, 0));
#else
    R = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
#endif

    I_inv = mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
}

RigidBody::~RigidBody() {
}


void RigidBody::draw(unsigned int drawable) {
	rigidBody->bind();
	rigidBody->draw();
	
	
}



std::vector<float> RigidBody::getY() {
    std::vector<float> state(STATES);
    int k = 0;

    state[k++] = x.x;
    state[k++] = x.y;
    state[k++] = x.z;

#ifdef USE_QUATERNIONS
    state[k++] = q.x;
    state[k++] = q.y;
    state[k++] = q.z;
    state[k++] = q.w;
#else
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            state[k++] = R[i][j];
        }
    }
#endif

    state[k++] = P.x;
    state[k++] = P.y;
    state[k++] = P.z;

    state[k++] = L.x;
    state[k++] = L.y;
    state[k++] = L.z;

    return state;
}

void RigidBody::setY(const std::vector<float>& y) {
    int k = 0;
    x.x = y[k++];
    x.y = y[k++];
    x.z = y[k++];

#ifdef USE_QUATERNIONS
    q.x = y[k++];
    q.y = y[k++];
    q.z = y[k++];
    q.w = y[k++];


#else
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = y[k++];
        }
    }
#endif

    P.x = y[k++];
    P.y = y[k++];
    P.z = y[k++];

    L.x = y[k++];
    L.y = y[k++];
    L.z = y[k++];

    // momentum
    v = P / m;

    // update inertia matrix
#ifdef USE_QUATERNIONS
    q = normalize(q);
    mat3 R = mat3_cast(q);
    I_inv = R * I_inv * transpose(R);
#else
    // must ensure that |R| = 1
    I_inv = R * I_inv * transpose(R);
#endif

    // angular momentum
    //w = I_inv * L;
}

std::vector<float> RigidBody::dydt(float t, const std::vector<float>& y) {
    // store initial state and override state
    std::vector<float> y0 = getY();
    setY(y);

    std::vector<float> yDot(STATES);
    int k = 0;

    //yDot = u
    yDot[k++] = v.x;
    yDot[k++] = v.y;
    yDot[k++] = v.z;

#ifdef USE_QUATERNIONS
    //q_dot = 1 / 2 * w * q;
    quat w_hat = quat(0, w);
    quat q_dot = (w_hat * q) / 2.0f;

    yDot[k++] = q_dot.x;
    yDot[k++] = q_dot.y;
    yDot[k++] = q_dot.z;
    yDot[k++] = q_dot.w;
#else
    mat3 w_hat = mat3(
        0, -w.z, w.y,
        w.z, 0, -w.x,
        -w.y, w.x, 0);
    mat3 R_dot = w_hat * R;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            yDot[k++] = R_dot[i][j];
        }
    }
#endif

    auto forces = forcing(t, y);
    //P_dot = f
    yDot[k++] = forces[0];
    yDot[k++] = forces[1];
    yDot[k++] = forces[2];

    //L_dot = tau
    yDot[k++] = forces[3];
    yDot[k++] = forces[4];
    yDot[k++] = forces[5];

    // restore initial state
    setY(y0);

    return yDot;
}

float RigidBody::calcKinecticEnergy() {
	return 0.5f * m * dot(v, v);//+ 0.5f * dot(w, glm::inverse(I_inv) * w);
}

std::vector<float> RigidBody::euler(float t, float h, const std::vector<float>& y0,glm::vec3 addVel,glm::vec3 earthPos,bool moon,bool energy) {
    
	std::vector<float> y1(STATES);
    std::vector<float> dydt0 = dydt(t, y0);
	
	//previous position-velocity
	vec3 vel0 = vec3(y0[7] / m, y0[8] / m, y0[9] / m);
	float  preVel = length(vel0);

	vec3 pos0 = vec3(y0[0], y0[1], y0[2]);

	/*_______main update_________*/
	for (int i = 3; i < STATES; i++) {
        y1[i] = y0[i] + h * dydt0[i];
    }
	y1[0] = y0[0] + h *(addVel.x + y1[7] / m);
	y1[1] = y0[1] + h *(addVel.y + y1[8] / m);
	y1[2] = y0[2] + h *(addVel.z + y1[9] / m);
	
	

	if (energy) {
		//new position-velocity
		vec3 velf = vec3(y1[7] / m, y1[8] / m, y1[9] / m);
		float  finalVel = length(velf);

		vec3 posf = vec3(y1[0], y1[1], y1[2]);

		/*___Energy Conservation____*/
		float Kp = 0.5*m*pow(preVel, 2);
		float Up = -G*(M_earth / distance(earthPos, pos0)) -G*M_sun / distance(vec3(float(2.499573*pow(10, 9)), 0, 0), pos0);

		//float Kf= 0.5*M_moon*pow(length(uLinear_moon), 2);
		float Uf = -G*(M_earth / distance(earthPos, posf)) -G*M_sun / distance(vec3(float(2.499573*pow(10, 9)), 0, 0), posf);

		float Kexpected;
		//valame *pow(10,20) gia na fainetai i metavoli tin taxutita
		Kexpected = Kp + Up/**pow(10,20)*/ - Uf/**pow(10,20)*/;

		float Vexpected = sqrt(Kexpected / (0.5*m));
		//float DV = sqrt(abs(pow(length(velf), 2) - pow(Vexpected, 2)));
		float DV =(pow(length(velf), 2) - pow(Vexpected, 2));
		

		float rad = (80.f / 100.f)*(dot(velf, earthPos - posf) / length(earthPos - posf));

		float perp = sqrt(pow(Vexpected, 2) - pow(rad, 2));
		vec3 newV = rad*normalize(earthPos - posf) - perp*cross(vec3(0, 1, 0), normalize(earthPos - posf));

		y1[7] = m* newV.x;
		y1[8] = m* newV.y;
		y1[9] = m* newV.z;
		
		y1[0] = pos0.x + h *(addVel.x + y1[7] / m);
		y1[1] = pos0.y + h *(addVel.y + y1[8] / m);
		y1[2] = pos0.z + h *(addVel.z + y1[9] / m);
		//posLinear_moon = prePos + (uLinear_earth + newV) *dt;
		
		
	}
    return y1;
}
std::vector<float> RigidBody::euler(float t, float h, const std::vector<float>& y0) {
	std::vector<float> y1(STATES);
	std::vector<float> dydt0 = dydt(t, y0);
	for (int i = 3; i < STATES; i++) {
		y1[i] = y0[i] + h * dydt0[i];
	}
	y1[0] = y0[0] + h *( y1[7] / m);
	y1[1] = y0[1] + h *( y1[8] / m);
	y1[2] = y0[2] + h *( y1[9] / m);

	return y1;
}

std::vector<float> RigidBody::rungeKuta4th(float t, float h, const std::vector<float>& y0) {
    // dydt0 = dydt(x0)
    std::vector<float> dydt0 = dydt(t, y0);

    // y1 = y0 + h * dydt0 / 2
    std::vector<float> y1(STATES);
    for (int i = 0; i < STATES; i++) {
        y1[i] = y0[i] + h * dydt0[i] / 2.0f;
    }
    std::vector<float> dydt1 = dydt(t + h / 2.0f, y1);

    // y2 = y0 + h * dydt1 / 2
    std::vector<float> y2(STATES);
    for (int i = 0; i < STATES; i++) {
        y2[i] = y0[i] + h * dydt1[i] / 2.0f;
    }
    std::vector<float> dydt2 = dydt(t + h / 2.0f, y2);

    // y2 = y0 + h * dydt1 / 2
    std::vector<float> y3(STATES);
    for (int i = 0; i < STATES; i++) {
        y3[i] = y0[i] + h * dydt2[i];
    }
    std::vector<float> dydt3 = dydt(t + h, y3);

    // combine them to estimate the solution.
    std::vector<float> y4(STATES);
    for (int i = 0; i < STATES; i++) {
        y4[i] = y0[i] + h * (dydt0[i] + 2.0f * dydt1[i]
            + 2.0f * dydt2[i] + dydt3[i]) / 6.0f;
    }
    return y4;
}
void RigidBody::advanceState(float t, float h,vec3 addVel,vec3 earthPos,bool m,bool energy) {
    // Task 2e:
    setY(euler(t, h, getY(), addVel, earthPos,m,energy));
    // setY(rungeKuta4th(t, h, getY()));
}
void RigidBody::advanceState(float t, float h) {
	// Task 2e:
	setY(euler(t, h, getY()));
	// setY(rungeKuta4th(t, h, getY()));
}

void RigidBody::advanceStateAdaptive(float t0, float dt, int maxIter, float tol) {
    float t = t0;
    float h = dt;
    int i = 0;
    auto y0 = getY();
    do {
        // euler @ t + h and @ (t + h / 2) + (t + h / 2)
        auto y1 = euler(t, h, y0);
        auto y21 = euler(t, h / 2.0f, y0);
        auto y2 = euler(t, h / 2.0f, y21);

        // find max error
        float error = 0;
        for (int i = 0; i < STATES; i++) {
            float temp = abs(y1[i] - y2[i]);
            if (temp > error) {
                error = temp;
            }
        }

        if (error < tol) { // solution accepted
            setY(y2);
            y0 = y2;
            t += h;
        } else { // step doubling
            h /= 2.0f;
            i++;
            if (i == maxIter) {
                throw std::runtime_error("adaptive time stepping exceeded max_iter\n");
            }
        }
    } while (t < t0 + dt);
}
