#include "euleralgorithm.h"
#include <iostream>

using namespace std;
EulerAlgorithm::EulerAlgorithm(double dt)
{
    this->delta_t = dt;
}

void EulerAlgorithm::stepForward(SolarSystem& sola){
    sola.calculateForce();
    for (CelestialObj &celObj : sola.CelObj()){
        celObj.position += celObj.velocity * delta_t;
        celObj.velocity += celObj.force / celObj.mass * delta_t;
        //cout << celObj.velocity.length() << endl;
        //cout << celObj.force.length() << endl;

    }
}

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
            fout << x << " " << y << endl;
        }
    }
    clock_t end = clock();
    cout << sola.KE + sola.PE << " "  << sola.momentum.length()<< endl;
    fout.close();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;

}

double EulerAlgorithm::runFor(SolarSystem &sola, double duration, string fileName){
    return this->runFor(sola, int(duration/delta_t), fileName);
}
