#include "solarsystem.h"
#include <string>
#include <fstream>
#include <ctime>

#ifndef VERLETALGORITHM_H
#define VERLETALGORITHM_H


class VerletAlgorithm
{
public:
    VerletAlgorithm(double dt);
    void stepForward(class SolarSystem& sola);
    double runFor(SolarSystem &sola, int steps, string fileName);
    double runFor(SolarSystem &sola, double duration, string fileName);
private:
    double delta_t;

};

#endif // VERLETALGORITHM_H


