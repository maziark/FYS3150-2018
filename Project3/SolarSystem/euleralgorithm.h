#include "solarsystem.h"
#include <string>
#include <fstream>
#include <ctime>

#ifndef EULERALGORITHM_H
#define EULERALGORITHM_H


class EulerAlgorithm
{
public:
    EulerAlgorithm(double dt);
    void stepForward (class SolarSystem& sola);
    double runFor (class SolarSystem& sola, int steps, string fileName);
    double runFor(SolarSystem &sola, double duration, string fileName);
private:
    double delta_t;
};

#endif // EULERALGORITHM_H
