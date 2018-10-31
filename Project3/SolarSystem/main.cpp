#include <iostream>
#include "solarsystem.h"
#include "celestialobj.h"
#include <cmath>
#include <fstream>
#include "euleralgorithm.h"
#include "verletalgorithm.h"
#include <string>


using namespace std;


/*
 * This function tries to estimate the escape velocity, it keeps increasing the initial
 * velocity of earth, until it gets out of orbit. (Or the distance keeps increasing)
 *
 * It returns a value n so the escape velocity would be n*pi
*/
double escapeVelocity (double epsilon, double step){
    SolarSystem ss;
    CelestialObj &sun = ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    CelestialObj &earth = ss.addCelestialObj(" Earth (399) ", 3e-6, vec3(1, 0, 0), vec3(0, 2*M_PI, 0));
    double count = 0;
    //EulerAlgorithm e(0.01);
    VerletAlgorithm v(.01);
    v.runFor(ss, 1.0, "Escape.txt");
    while (fabs((earth.position - sun.position).length() - 1) < epsilon ){
        count ++;
        ss.setBack();
        earth.velocity = vec3(0, (2 + count * step)*M_PI, 0);
        v.runFor( ss, 1.0, "Escape.txt");
    }
    return 2 + count * step;
}


/*
 * adds all the planet, and simulates the whole solar system.
 *
 * The values are produces with a python code, that gets the initial values with an API.
 * (initialValues.py)
 *
 * And saves the result in build.../Report/SolarSystem.txt
 * using the verlet velocity algorithm
*/
void wholeThing (){
    SolarSystem ss;
    CelestialObj &sun = ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss.addCelestialObj(" Mercury (199) ",  3e-06 , vec3( -0.110692716399 ,  -0.452571796847 ,  -0.0268264553904 ) , vec3( 7.9137485748 ,  -1.92341214333 ,  -0.883162004382 ));
    ss.addCelestialObj(" Venus (299) ",  9.5e-06 , vec3( 0.69659394236 ,  0.198982693119 ,  -0.0374678908737 ) , vec3( -2.05375067135 ,  7.06569310618 ,  0.215451299099 ));
    ss.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));
    ss.addCelestialObj(" Mars (499) ",  2.4500000000000003e-06 , vec3( 1.38156231885 ,  -0.124669709005 ,  -0.0365128741782 ) , vec3( 0.654137997378 ,  5.52344287802 ,  0.0996881387514 ));
    ss.addCelestialObj(" Jupiter (599) ",   0.00027499999999999996 , vec3( -2.64739396282 ,  -4.67305743077 ,  0.0786456882645 ) , vec3( 2.3653772054 ,  -1.22916498934 ,  -0.0478313076139 ));
    ss.addCelestialObj(" Saturn (699) ",  1.65e-07 , vec3( 1.55973004042 ,  -9.94076143511 ,  0.110709244224 ) , vec3( 1.90202860593 ,  0.308411319459 ,  -0.0810363789727 ));
    ss.addCelestialObj(" Uranus (799) ",  4.4e-05 , vec3( 17.1720615577 ,  9.99684099244 ,  -0.185236412458 ) , vec3( -0.730404179716 ,  1.17257139265 ,  0.0137871649976 ));
    ss.addCelestialObj(" Neptune (899) ",  5.15e-05 , vec3( 28.9220401115 ,  -7.72369017482 ,  -0.50755595577 ) , vec3( 0.290653024201 ,  1.11289564521 ,  -0.029775649625 ));
    ss.addCelestialObj(" Pluto (999) ",  6.550000000000001e-09 , vec3( 11.6500836526 ,  -31.58154113 ,  0.00885314582401 ) , vec3( 1.09905629976 ,  0.151051720665 ,  -0.333322090737 ));

    //EulerAlgorithm e(0.00001);
    //e.runFor(ss, 100000, "body2.txt");
    VerletAlgorithm v(1e-3);
    v.runFor(ss, 1, "body2.txt");
    sun.velocity -=  ss.momentum;
    //v.runFor(ss, 1000, "body2.txt");
    v.runFor(ss, 300.0, "Report/SolarSystem.txt");
}

/*
 * Runs the twoBody problem with Euler and Verlet, and compares their result based on
 * conservation of Energy and conservation of momentum.
 *
 * it gets the delta_t and number of years.
 *
 * The output will be stored in build.../Report/twoBodyEnergy_[{1/dt}].txt
 *
*/
void twoBodyEuler (double dt, double years){
    SolarSystem ss_e, ss_v;


    ss_e.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss_e.addCelestialObj(" Earth (399) ",  3e-6, vec3(1, 0, 0), vec3(0, 2*M_PI, 0));

    ss_v.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss_v.addCelestialObj(" Earth (399) ",  3e-6, vec3(1, 0, 0), vec3(0, 2*M_PI, 0));

    ofstream fout ("Report/twoBodyEnergy_" + to_string(int(1/dt)) + ".txt");

    EulerAlgorithm e(dt);
    VerletAlgorithm v(dt);
    for (int i = 1; i < int(years); i++) {
        //e.runFor(ss_e, 1.0, "Report/Eulerbody2_"+ to_string(int(1/dt))+".txt");
        //v.runFor(ss_v, 1.0, "Report/Verletbody2_"+ to_string(int(1/dt))+".txt");
        e.runFor(ss_e, 1.0, "body2.txt");
        v.runFor(ss_v, 1.0, "body2.txt");
        fout << ss_e.PE << " " << ss_e.KE << " " << ss_v.PE << " " << ss_v.KE << endl;//" " << ss.momentum << endl;
    }
    fout.close();
}

/*
 * Compares Euler and Verlet algorithm with different delta_t stamps.
 * Not sure if my groupmate has used this in the report or not, but it is good to have!
*/
void time_report (){
    ofstream fout ("Report/twoBodyTime.txt");
    SolarSystem ss_e, ss_v;
    ss_e.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss_e.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));

    ss_v.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss_v.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));

    for (int i = 1; i < 9; i++){
        double E_T, V_T;

        EulerAlgorithm e(pow(10, double(-1 * i)));
        VerletAlgorithm v(pow(10, double(-1 * i)));

        E_T = e.runFor(ss_e, 1.0, "body2.txt");
        V_T = v.runFor(ss_v, 1.0, "body2.txt");

        fout << pow(10, double(-1 * i)) << " " << E_T << " " << V_T << endl;
        ss_e.setBack();
        ss_v.setBack();

    }
    fout.close();
}


void testEn (){
    //ofstream fout ("Report/testEn.txt");
    SolarSystem ss_e, ss_v;
    //ss_e.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    //ss_e.addCelestialObj(" Mercury (199) ",  3e-06 , vec3( -0.110692716399 ,  -0.452571796847 ,  -0.0268264553904 ) , vec3( 7.9137485748 ,  -1.92341214333 ,  -0.883162004382 ));

    CelestialObj &sun = ss_v.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    CelestialObj &mercury =  ss_v.addCelestialObj(" Mercury (199) ",  3e-06 , vec3( 0.3075 ,  0 ,  0 ) , vec3( 0 ,  12.44 ,  0 ));

    ss_v.En = true;
    ss_v.calculateTotalMomentum();
    sun.velocity -= ss_v.momentum/sun.mass;

    vector<vec3> perihelion;

    VerletAlgorithm v(10e-6);

    cout << "Time it takes : " << v.testEn(ss_v, 100.0, "Report/testEn.txt", perihelion) << endl;
}



void testEnOff (){
    //ofstream fout ("Report/testEnOff.txt");
    SolarSystem  ss_v;
    //ss_e.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    //ss_e.addCelestialObj(" Mercury (199) ",  3e-06 , vec3( -0.110692716399 ,  -0.452571796847 ,  -0.0268264553904 ) , vec3( 7.9137485748 ,  -1.92341214333 ,  -0.883162004382 ));

    CelestialObj &sun = ss_v.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    CelestialObj &mercury =  ss_v.addCelestialObj(" Mercury (199) ",  3e-06 , vec3( 0.3075 ,  0 ,  0 ) , vec3( 0 ,  12.44 ,  0 ));

    ss_v.En = false;
    ss_v.calculateTotalMomentum();
    sun.velocity -= ss_v.momentum/sun.mass;

    vector<vec3> perihelion;

    VerletAlgorithm v(10e-7);

    cout << "Time it takes : " << v.testEn(ss_v, 100.0, "Report/testOff.txt", perihelion) << endl;
}

/* runs the three body simulation, sun, earth and jupiter, with jupiter's actual mass */
void threeBody_1 (){
    SolarSystem ss;
    ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));
    ss.addCelestialObj(" Jupiter (599) ",   0.00027499999999999996 , vec3( -2.64739396282 ,  -4.67305743077 ,  0.0786456882645 ) , vec3( 2.3653772054 ,  -1.22916498934 ,  -0.0478313076139 ));

    VerletAlgorithm v(1e-4);
    v.runFor(ss, 50.0, "threeBody1.txt");
}


/* runs the three body simulation, sun, earth and jupiter, with 10 times jupiter's mass  */
void threeBody_2 (){
    SolarSystem ss;
    ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));
    ss.addCelestialObj(" Jupiter (599) ",   10 * 0.00027499999999999996 , vec3( -2.64739396282 ,  -4.67305743077 ,  0.0786456882645 ) , vec3( 2.3653772054 ,  -1.22916498934 ,  -0.0478313076139 ));

    VerletAlgorithm v(1e-4);
    v.runFor(ss, 50.0, "threeBody10.txt");
}


/* runs the three body simulation, sun, earth and jupiter, with 100 times jupiter's mass  */
void threeBody_3 (){
    SolarSystem ss;
    ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));
    ss.addCelestialObj(" Jupiter (599) ",   100 * 0.00027499999999999996 , vec3( -2.64739396282 ,  -4.67305743077 ,  0.0786456882645 ) , vec3( 2.3653772054 ,  -1.22916498934 ,  -0.0478313076139 ));

    VerletAlgorithm v(1e-4);
    v.runFor(ss, 50.0, "threeBody100.txt");
}


/* M_Jupiter *= 10e3
 * two lovers dancing together and the earth in between cockblocking them!
*/
void threeBody_4 (){
    SolarSystem ss;
    ss.addCelestialObj(" Sun (10) ",  1.0 , vec3( 0.0 ,  0.0 ,  0.0 ) , vec3( 0.0 ,  0.0 ,  0.0 ));
    ss.addCelestialObj(" Earth (399) ",  3.3e-07 , vec3( 0.922303424508 ,  0.378697251855 ,  -2.11210097393e-05 ) , vec3( -2.48970762794 ,  5.78497703731 ,  -3.95468655934e-05 ));
    ss.addCelestialObj(" Jupiter (599) ",   1000 * 0.00027499999999999996 , vec3( -2.64739396282 ,  -4.67305743077 ,  0.0786456882645 ) , vec3( 2.3653772054 ,  -1.22916498934 ,  -0.0478313076139 ));

    VerletAlgorithm v(1e-4);
    v.runFor(ss, 50.0, "threeBody1000.txt");
}



int main(int argc, char *argv[])
{
    twoBodyEuler(10e-1, 10);
    twoBodyEuler(10e-2, 10);
    twoBodyEuler(10e-3, 10);
    twoBodyEuler(10e-4, 10);

    /*SolarSystem ss;
    CelestialObj &sun = ss.addCelestialObj(1.0, vec3(0,0,0), vec3(0,0,0));
    CelestialObj &earth = ss.addCelestialObj(3e-6, vec3(1, 0, 0), vec3(0, 2*M_PI, 0));*/
    //ss.addCelestialObj(9.5e-4, vec3(5.2, 0, 0), vec3(0, 2.75516462, 0));

    //EulerAlgorithm e(0.00001);

    //VerletAlgorithm v(0.001);
    //wholeThing();
    //e.runFor(ss, 100000, "body2.txt");
    //v.runFor(ss, 100000, "body2.txt");


    // To find Escape Velocity
    //cout << escapeVelocity (4.0, 0.01) << endl;
    //wholeThing();
    //time_report();

    //threeBody_1();
    //threeBody_2();
    //threeBody_3();
    //threeBody_4();
    //testEn();
    //testEnOff();
    return 0;
}
