#include "solarsystem.h"

SolarSystem::SolarSystem()
{

}

/*
 * Constructors {
*/

// This constructor is just there for testing reasons, I prefer planets with names!
CelestialObj& SolarSystem::addCelestialObj(double mass, vec3 position, vec3 velocity){
    this->collection.push_back(CelestialObj(mass, position, velocity));
    return collection.back();
}


// Will define a celestial object with the values given, and push them in the collection vector.
CelestialObj& SolarSystem::addCelestialObj(string name, double mass, vec3 position, vec3 velocity){
    this->collection.push_back(CelestialObj(name, mass, position, velocity));
    return collection.back();
}

/*
 * }Constructors
*/

/*
 *  return collectoin so the solver functions (The ones defined in Euler,
 *  and Verlet algorithm) get access to planets
*/
vector<CelestialObj> &SolarSystem::CelObj(){
    return collection;
}


/*
 * The heart of the program, and the one that is the most time consuming
 * while returning absoloutly nothing! (void)
 *
 * the idea is that it would set the KE, PE, momentum, and forces
 * all to zero, and then produces pairs of celestiabl bodies, and
 * calculates the forces that are involved.
 *
 * Using the third law of Newton, (I push you, you push me back), it cuts
 * the calculation in half, so for n planets it will run for n(n - 1)/2
 *
 * it then calculates the KE, PE, and the force.
 *
*/
void SolarSystem::calculateForce(){
    double G = 4*M_PI * M_PI;
    this->KE = 0;
    this->PE = 0;
    this->momentum = vec3(0,0,0);
    for (CelestialObj &obj : collection){
        obj.force.zeros();
    }

    for (int i = 0; i < collection.size(); i++){
        CelestialObj &obj1 = collection[i];
        for (int j = i + 1; j < collection.size(); j++){
            CelestialObj &obj2 = collection[j];
            vec3 r = obj1.position - obj2.position;
            vec3 f;
            if (En) f = obj1.calcForce(obj2, true);
            else f = obj1.calcForce(obj2, 2.0);
            obj1.force -= f;
            obj2.force += f;
            this->PE += obj1.mass*obj2.mass*r.length();
        }
        obj1.force *= G;
        this->momentum += obj1.mass * obj1.velocity;
        this->KE += obj1.mass * obj1.velocity.lengthSquared();
    }
    this->KE *= 0.5;
    this->PE *= 2;

}


/*
 * This function had a very bright future, till I realized, I'm clueless
 * when it gets to 3D plotting in python.
 *
 * it basically opens a file, if it is not already open, and then write the
 * values of every celestial body in the solar system in that file
 *
*/
void SolarSystem::writeToFile(string filename){
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_file << collection.size() << endl;
    m_file << "Comment line that needs to be here. Balle." << endl;
    for(CelestialObj &body : collection) {
        m_file << "1 " << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n";
    }
}


/*
 * Calculates the total momentum of the system O(n)
*/
void SolarSystem::calculateTotalMomentum(){
    this->momentum.zeros();
    for(CelestialObj &body : collection) {
        this->momentum += body.velocity * body.mass;
    }
}

/*
 * calls the setBack function of every body in the solar system.
 * O(n)
*/
void SolarSystem::setBack(){
    for (CelestialObj &obj : collection){
        obj.setBack();
    }
}

// zero total momenum of the Solar system with sun
void SolarSystem::zeroTotalMomentum(){
    zeroTotalMomentum(0);
}

/*
 * is used to make sure that the total momentum is zero,
 * if the sun is not the center of orbit.
*/
void SolarSystem::zeroTotalMomentum(int n){
    calculateTotalMomentum();
    int i = 0;
    for (CelestialObj &obj : collection){
        if (i == n) obj.velocity -= this->momentum;
        i++;
    }
}

/*
 * returns the distance of a planet from the sun (considering that the
 * sun is the first body added into the collection.
 *
 * Used for perihelion of Mercury.
*/
double SolarSystem::distanceFromSol(int i){
    return (this->collection.at(i).position - this->collection.at(0).position).length();
}
