/*
 * Similar to EulerAlgorithm class, this class is a ODE solver,
 * following the velocity verlet algorithm.
 *
 * double delta_t keeps track of the step size of every calculation.
 *
 * void stepForward, moves the time forward for one unit of time.
 *
 * testEn, runs the system for the defined duration, and keeps track of every
 * perihelion(spelling out this word right, is just IMPOSSIBLE) point during that duration.
 *
*/

#include "solarsystem.h"
#include <string>
#include <fstream>
#include <ctime>
#include <vector>
#include <math.h>
using namespace std;

#ifndef VERLETALGORITHM_H
#define VERLETALGORITHM_H



class VerletAlgorithm
{
public:
    VerletAlgorithm(double dt);
    void stepForward(class SolarSystem& sola);
    double runFor(SolarSystem &sola, int steps, string fileName);
    double runFor(SolarSystem &sola, double duration, string fileName);
    double testEn(SolarSystem &sola, double duration, string fileName, vector<vec3> &perihelion);
private:
    double delta_t;

};

#endif // VERLETALGORITHM_H


