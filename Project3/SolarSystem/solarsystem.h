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
private:
    vector<CelestialObj> collection;
    ofstream m_file;
};

#endif // SOLARSYSTEM_H
