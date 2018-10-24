#include "solarsystem.h"

SolarSystem::SolarSystem()
{

}

CelestialObj& SolarSystem::addCelestialObj(double mass, vec3 position, vec3 velocity){
    this->collection.push_back(CelestialObj(mass, position, velocity));
    return collection.back();
}

CelestialObj& SolarSystem::addCelestialObj(string name, double mass, vec3 position, vec3 velocity){
    this->collection.push_back(CelestialObj(name, mass, position, velocity));
    return collection.back();
}

vector<CelestialObj> &SolarSystem::CelObj(){
    return collection;
}



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

void SolarSystem::calculateTotalMomentum(){
    this->momentum.zeros();
    for(CelestialObj &body : collection) {
        this->momentum += body.velocity * body.mass;
    }
}

void SolarSystem::setBack(){
    for (CelestialObj &obj : collection){
        obj.setBack();
    }
}

void SolarSystem::zeroTotalMomentum(){
    // zero total momenum of the Solar system with sun
    zeroTotalMomentum(0);
}

void SolarSystem::zeroTotalMomentum(int n){
    calculateTotalMomentum();
    int i = 0;
    for (CelestialObj &obj : collection){
        if (i == n) obj.velocity -= this->momentum;
        i++;
    }
}

