/*
 * The SolarSystem class is the definition for every SolarSystem,
 * it keeps track of every celestial body that belongs to the solar system.
 *
 * This class has few different flags that would enable certain features.
 * If the flag bool En sets to be true, it would uses the general law of
 * relativity and its correction of gravitational force to calculate the
 * forces.
 *
 * double KE, PE (Calculated in every turn, it is a good measure of how circular
 * the orbit is (If KE = PE) or if KE + PE stays constant or not (law of conservation of energy)
 *
 * vec3 momentum (is also another way of testing if the simulation is
 * functioning properly or not)
 *
 * vector<CelestialObj> &CelObj() : This will keep track of every celestial body that
 * belongs to the solar system.
*/

#include "celestialobj.h"
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>

#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

using namespace std;
class SolarSystem
{
public:
    double KE, PE;
    bool En = false;
    vec3 momentum;
    SolarSystem();
    CelestialObj &addCelestialObj (double mass, vec3 position, vec3 velocity);
    CelestialObj &addCelestialObj (string name, double mass, vec3 position, vec3 velocity);

    vector<CelestialObj> &CelObj();

    void calculateForce ();
    void calculateTotalMomentum();
    void writeToFile(string filename);
    void setBack ();
    void zeroTotalMomentum ();
    void zeroTotalMomentum (int n);
    double distanceFromSol (int i);
private:
    vector<CelestialObj> collection;
    ofstream m_file;
};

#endif // SOLARSYSTEM_H
