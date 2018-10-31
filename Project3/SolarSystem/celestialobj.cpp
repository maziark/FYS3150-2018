#include "celestialobj.h"

/* Constructors { */
CelestialObj::CelestialObj()
{
    this->mass = 0;
    this->velocity = vec3(0,0,0);
    this->init_position = vec3(0,0,0);
    this->init_velocity = vec3(0, 0, 0);
    this->position = vec3(0,0,0);
}

CelestialObj::CelestialObj(double mass, vec3 position, vec3 velocity){
    this->mass = mass;
    this->velocity = velocity;
    this->position = position;
    this->init_position = position;
    this->init_velocity = velocity;
}

CelestialObj::CelestialObj(string name, double mass, vec3 position, vec3 velocity){
    this->mass = mass;
    this->velocity = velocity;
    this->position = position;
    this->name = name;
    this->init_position = position;
    this->init_velocity = velocity;
}

/* }Constructors  */

/*
 * This function calculates the force between this body and the other one.
 * Since there is already a function that would do the calculation with
 * beta as a variable, this is just makes a call and returns the value.
 *
*/
vec3 CelestialObj::calcForce(CelestialObj &other){
    return this->calcForce(other, 2.0);
}

/*
 * Calculating gravitational force using the relativity correction.
 *
*/
vec3 CelestialObj::calcForce(CelestialObj &other, bool Ens){

    vec3 f = calcForce(other);

    // 1 + \frac{3l^2}{r^2c^2}
    double c2 = 3993960777.9406304;
    vec3 r = this->position - other.position;
    vec3 v = this->velocity;
    double rel = 3.0 * (r.cross(v).lengthSquared())/(r.lengthSquared()*c2);
    return f * (1 + rel);
}

/*
 * This function calculates the Newtonian gravitational force.
 * and it would return the vector of the force caused by the
 * other planet
 *
*/
vec3 CelestialObj::calcForce(CelestialObj &other, double beta){
    vec3 r = this->position - other.position;
    double d = r.length();
    if (d == 0) return vec3(0, 0, 0);
    vec3 f = this->mass*other.mass*r/pow(d, beta + 1.0);

    return f;
}




/*std::ostream& operator<<(std::ostream& os, const CelestialObj& body)
{
    os << body.position << " " << body.velocity << " ";
    return os;
}*/
