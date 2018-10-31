/*
 * This class is used to create Celestial bodies
 * (As its names suggests!)
 * Each Celestial Body has a
 * String name,
 * vec3 initial position, initial velocity
 * The values mentioned about are initiated when the
 * object is created.
 *
 * vec3 force, acceleration
 * are found during the process.
 *
 * This class is whether independantly called with in
 * the main to keep track of certain planets,
 * but mainly used in the SolarSystem class.
 *
*/


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


    /*
     * The function setBack would turn back the time to the point where
     * the object just got defined, so it can used for comparisions
    */
    void setBack (){
        this->position = this->init_position;
        this->velocity = this->init_velocity;
    }

    friend std::ostream& operator<<(std::ostream& os, const CelestialObj& body); // Allows cout << myVector << endl;

};

#endif // CELESTIALOBJ_H
