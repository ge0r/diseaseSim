#! /usr/bin/env python2.7

import numpy as np

from my_methods import *


class City:

    def __init__(self, name, population, infected, transport, emergency_policy):
        self.name = name
        self.population = population
        self.infected = infected
        self.susceptible = population - infected
        self.recovered = 0

        # divide each element of the transport vector by the city's population
        self.w = np.matrix(transport)[:]/population

        print self.w

        self.emergency_strategy = emergency_policy[0]
        self.emergency_criterion = emergency_policy[1]
        self.open_airport = True

    def update(self, infected, recovered):
        self.susceptible = self.population - (infected + recovered)
        self.infected = infected
        self.recovered = recovered

    def check_emergency(self, days):

        # if the epidemic lasts for long enough, the city closes its airport
        if self.emergency_strategy == 1:
            if days > self.emergency_criterion:
                if self.open_airport is not False:
                    print "emergency strategy 1"
                    self.open_airport = False

        # if a high enough number of people have been infected in the city, the city closes its airport
        elif self.emergency_strategy == 2:
            if self.infected > self.emergency_criterion:
                if self.open_airport is not False:
                    print "emergency strategy 2"
                    self.open_airport = False

        # in an act of self sacrifice, if a high number of infected travel abroad, the city closes its airport
        elif self.emergency_strategy == 3:
            _sum = 0
            for i in xrange(0, self.w.shape[1]):
                _sum += self.w[0, i]*self.infected

            if _sum > self.emergency_criterion:
                if self.open_airport is not False:
                    print "emergency strategy 3"
                    self.open_airport = False

        return not self.open_airport


# assign parameter values
alpha = 0.017
beta = 0.22
gamma = 0.000526

ratio = 0.05

# h should be a fraction of the smallest duration/largest rate
if beta >= alpha:
    h = ratio*(1/beta)
else:
    h = ratio*(1/alpha)

# duration of the model in days
simulation_duration = 365

steps = int(simulation_duration // h)

# number of airport passengers per day
A_stock_cph = 4011.54
A_cph_oslo = 3936.22
A_oslo_stock = 3672.32

# list containing each city[name, population, initial number of infected]
cities = [City("Stockholm", 2192433, 1,  [0, A_stock_cph, A_oslo_stock], [1, 1]),
          City("Copenhagen", 2016285, 10000, [A_stock_cph, 0, A_cph_oslo], [0, 5000]),
          City("Oslo", 1717900, 1000000,        [A_oslo_stock, A_cph_oslo, 0], [0, 0])]



# SIS.run(alpha, beta, cities[1], h, steps, SIS.euler)
# SIS.run(alpha, beta, cities[1], h, steps, SIS.rk_2)
# SIS.run(alpha, beta, cities[1], h, steps, SIS.rk_4)

# SIR.run(alpha, beta, gamma, cities, h, steps, SIR.euler)
# SIR.run(alpha, beta, gamma, cities, h, steps, SIR.rk_4)
SIR.run(alpha, beta, gamma, cities, h, steps, SIR.sim)

plt.show()
