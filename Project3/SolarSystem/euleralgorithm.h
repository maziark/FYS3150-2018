/*
 * Euler algorithm, one of the solver classes on this assignment,
 * double delta_t defines the step size of every calculation.
 *
 * void stepForward : will only move the system by one unit of delta_t,
 * this can be called by the two runFor functions to run the time for
 * longer durations
 *
*/

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
