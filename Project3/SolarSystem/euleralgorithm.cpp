#include "euleralgorithm.h"
#include <iostream>

using namespace std;
/*
 * constructor {
*/
EulerAlgorithm::EulerAlgorithm(double dt)
{
    this->delta_t = dt;
}
/*
 * }constructor
*/


/*
 * pushes the time forward by one step of size (delta_t)
 * it calculates the force resulted by other bodies in the
 * solar system, and used the Euler algorithm to estimates
 * the position and the velocity.
 *
 * (a_i = f_i/m_i)
 *
 * x_{i+1} = x_i + v_i*dt
 * v_{i+1} = v_i + a_i*dt
 *
 *
*/
void EulerAlgorithm::stepForward(SolarSystem& sola){
    sola.calculateForce();
    for (CelestialObj &celObj : sola.CelObj()){
        celObj.position += celObj.velocity * delta_t;
        celObj.velocity += celObj.force / celObj.mass * delta_t;
        //cout << celObj.velocity.length() << endl;
        //cout << celObj.force.length() << endl;

    }
}

/*
 * calls the stepForward function for a certain number of steps
 * and saves the output in a file (so it can be plotted.
 *
 * I used it as a testing function. It is called by (the other) runFor
 * function.
 *
 * it returns the time it takes to do all the calculations.
 *
*/

double EulerAlgorithm::runFor(SolarSystem &sola, int steps, string fileName){
    ofstream fout (fileName);
    clock_t begin = clock();
    double x, y, z;
    for (int i = 0; i < steps; i++){
        if (i == 1) cout << sola.KE + sola.PE << " "  << sola.momentum.length()<< endl;
        this->stepForward(sola);
        for (CelestialObj &celObj : sola.CelObj()){
            x = celObj.position[0];
            y = celObj.position[1];
            z = celObj.position[2];
            fout << celObj.name << " " << x << " " << y << endl;
        }
    }
    clock_t end = clock();
    cout << sola.KE + sola.PE << " "  << sola.momentum.length()<< endl;
    fout.close();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;

}

/*
 * runs the system for a duration of time, calculates the number of steps it requires and calls
 * the functoin above to do the calculation.
*/
double EulerAlgorithm::runFor(SolarSystem &sola, double duration, string fileName){
    return this->runFor(sola, int(duration/delta_t), fileName);
}
