#!/usr/bin/env python

# -*- coding: utf-8 -*-
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

ids = [str(i*100+99) for i in range(1,10)]
ids.insert(0, '10')
#ids.insert(len(ids), '1977 UB')
objs = [Horizons(id=i, id_type= 'majorbody', epochs={'start':'2018-10-16', 'stop':'2018-11-15','step':'1d'}) for i in ids]

vecs = [obj.vectors()[0] for obj in objs]

planets = [[vec['targetname'],
            vec['x'], vec['y'], vec['z'],
            vec['vx'] * 365, vec['vy'] * 365, vec['vz'] * 365] for vec in vecs]

masses = [2e30, 6e24, 1.9e25, 6.6e23, 4.9e24, 5.5e26, 3.3e23, 8.8e25, 1.03e26, 1.31e22]
m_relative = [m/masses[0] for m in masses]
i = -1
for p in planets :
    i += 1;
    print ('ss.addCelestialObj("', p[0], '", ', str(m_relative[i]),
                ', vec3(', str(p[1]), ', ', str(p[2]), ', ', str(p[3]), ')',
                ', vec3(', str(p[4]), ', ', str(p[5]), ', ', str(p[6]), '));')

"""with open('planets.txt', 'w') as file:
    for p in planets"""
