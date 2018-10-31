/* I've tried to make sure, every function in EulerAlgorithm has a sibling in this class,
 * so, this way the only thing needs to be changed is the definision of the solver, every-
 * -thing else would remain the same.
*/


#include "verletalgorithm.h"
// Constructor
VerletAlgorithm::VerletAlgorithm(double dt)
{
    this->delta_t = dt;
}


/*
 * Similar to Euler algorithm, this function pushes the system
 * by one unit of time (with step size of delta_t)
 *
 * cosidering the nature of verlet algorithm, every step would take
 * 2ice as long, Since there are two for loops, O(2n*n^2/2) = O(n^3)
 *
 * it calculates the position based on the velocity,
 * then estimates the velocity based on the current acceleration, and then
 * averages the change with the next value of acceleration,
 *
 * this way it reduces the error of calculaiton by a degree
*/
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

/*
 * this function calls the stepForward steps times!
 * and prints the output in a file
 *
 * it returns the time it took to run the algorithm.
 * SolarSystem &sola    : pointer to the solar system
 * int steps            : number of steps with size delta_t
 * string fileName      : where to store the output
 *
*/
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


/*
 * similar to the same function in Euler algorithm this will run for a duration
 * of time instead of number of steps.
 * It calculates the number of steps required and calls the other function with the similar
 * name to do the calculation.
 *
 * it will return the time it takes to run the runFor function.
*/
double VerletAlgorithm::runFor(SolarSystem &sola, double duration, string fileName){
    return this->runFor(sola, int(duration/delta_t), fileName);
}

/*
 * This function is supposed to run the find the perihelion of mercury over a century.
 *
 * Knowing the nature of a Perihelion point, to be able to run the functoin faster, instead
 * of checking for speed and distance, I find the local minima points or to visualize it : \_/
 *
 * This way the results will be pretty accurate even with larger value of delta_t (step sizes)
 *
 * Then it will store the points of perihelion in every Mercury's year, and rights them in a file.
 *
 * SolarSystem &sola        : pointer to the solar system.
 * double duration          : duration (100 earth years!)
 * string fileName          : file to store the data
 * vector<vec3> &perihelion : vector of perihelion points
 *
 * and it returns the time it took to run the program for that duration of time.
 */
double VerletAlgorithm::testEn(SolarSystem &sola, double duration, string fileName, vector<vec3> &perihelion){
    ofstream fout (fileName);
    clock_t begin = clock();
    vec3 temp();
    double steps = int(duration/delta_t);
    double last_distance = sola.distanceFromSol(1);
    double last2_distance = sola.distanceFromSol(1);
    vec3 v = sola.CelObj()[1].position;

    for (int i = 0; i < steps; i++){
        sola.CelObj()[0].position = vec3(0, 0, 0);
        this->stepForward(sola);
        double current = sola.distanceFromSol(1);
        cout << "loading.." << endl;
        // Local minima
        if (last_distance < current && last_distance < last2_distance){
            //cout << last_distance << endl;
            perihelion.push_back(v);}
        last2_distance = last_distance;
        last_distance = current;
        v = sola.CelObj()[1].position;

        //fout << sola.CelObj()[1].name << " " << v[0] << " " << v[1] << endl;
    }
    double x, y;
    cout << perihelion.at(perihelion.size() - 1);
    double angle = atan(perihelion.at(perihelion.size() - 1)[1] / perihelion.at(perihelion.size() - 1)[0]);
    angle = angle * 3600*180/M_PI;
    cout << "angle of change = " << angle << endl;
    for (vec3 &v : perihelion){
        x = v[0];
        y = v[1];
        fout << sola.CelObj()[0].name << " " << 0 << " " << 0 << endl;
        fout << sola.CelObj()[1].name << " " << x << " " << y << endl;

    }



    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    fout.close();
    return elapsed_secs;

}
