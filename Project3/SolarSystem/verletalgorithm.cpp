#include "verletalgorithm.h"

VerletAlgorithm::VerletAlgorithm(double dt)
{
    this->delta_t = dt;
}

void VerletAlgorithm::stepForward(SolarSystem &sola){
    sola.calculateForce();
    for (CelestialObj &celObj : sola.CelObj()){

        celObj.position += celObj.velocity * delta_t + 0.5 * (celObj.force / celObj.mass) * delta_t * delta_t;
        celObj.velocity += 0.5 * (celObj.force / celObj.mass) * delta_t;
    }

    sola.calculateForce();
    for (CelestialObj &celObj : sola.CelObj()){
        celObj.velocity += 0.5 * (celObj.force / celObj.mass) * delta_t;
    }

}

double VerletAlgorithm::runFor(SolarSystem &sola, int steps, string fileName){
    ofstream fout (fileName);
    clock_t begin = clock();
    double x, y, z;
    for (int i = 0; i < steps; i++){
        this->stepForward(sola);

        for (CelestialObj &celObj : sola.CelObj()){
            x = celObj.position[0];
            y = celObj.position[1];
            z = celObj.position[2];
            fout << celObj.name << " " << x << " " << y << endl;
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    fout.close();
    return elapsed_secs;
}

double VerletAlgorithm::runFor(SolarSystem &sola, double duration, string fileName){
    return this->runFor(sola, int(duration/delta_t), fileName);
}
