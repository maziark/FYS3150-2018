#include "vec3.h"
#include <math.h>
#include <string>
#ifndef CELESTIALOBJ_H
#define CELESTIALOBJ_H

using namespace std;
class CelestialObj
{
public:
    double mass;
    string name = "";
    vec3 position, init_position;
    vec3 velocity, init_velocity;
    vec3 acceleration = vec3(0,0,0);

    vec3 force = vec3(0, 0, 0);

    CelestialObj();
    CelestialObj(double mass, vec3 position, vec3 velocity);
    CelestialObj(string name, double mass, vec3 position, vec3 velocity);


    vec3 calcForce(CelestialObj& other);
    vec3 calcForce(CelestialObj &other, bool Ens);
    vec3 calcForce(CelestialObj& other, double beta);

    void setBack (){
        this->position = this->init_position;
        this->velocity = this->init_velocity;
    }

    friend std::ostream& operator<<(std::ostream& os, const CelestialObj& body); // Allows cout << myVector << endl;

};

#endif // CELESTIALOBJ_H
