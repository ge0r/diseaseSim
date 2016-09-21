#! /usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt

class SIS:

    @staticmethod
    def run(alpha, beta, city, h, steps, method):
        # duration of SIS is t = steps*h
        t = np.zeros(steps)
        S = np.zeros(steps)
        I = np.zeros(steps)

        I[0] = city.infected
        S[0] = city.susceptible
        N = city.population

        method(alpha, beta, t, S, I, N, h, steps)

        axes = plt.gca()
        axes.set_ylim([0, N/1000000])
        axes.set_xlim([0, t[steps-1]])

    @staticmethod
    def susceptible(y, args):
        alpha = args[1]
        beta = args[2]
        N = args[3]
        S = y[0]
        I = y[1]
        return -beta*S*I/N + alpha*I

    @staticmethod
    def infected(y, args):
        alpha = args[1]
        beta = args[2]
        N = args[3]
        S = y[0]
        I = y[1]
        return beta*S*I/N - alpha*I

    @staticmethod
    # SIS Euler method
    def euler(alpha, beta, t, S, I, N, h, steps):
        # preparing arguments for functions
        args = [h, alpha, beta, N]

        for i in xrange(0, steps-1):
            y = [S[i], I[i]]

            t[i+1] = t[i] + h
            S[i+1] = S[i] + h * SIS.susceptible(y, args)
            I[i+1] = I[i] + h * SIS.infected(y, args)

        print 'S = '+str(S[27//h]-227882.312249)+'     I = '+ str(I[27/h]-1788402.68775)

        plt.plot(t, S[:] / 1000000, "g", linewidth=2, label='S_sto')
        plt.plot(t, I[:] / 1000000, "r", linewidth=2, label='I_sto')

        plt.title("SIS")
        plt.legend(loc=1)
        plt.xlabel("Time (days)")
        plt.ylabel("Number of persons (millions)")

    @staticmethod
    # SIS Second-Order Runge-Kutta
    def rk_2(alpha, beta, t, S, I, N, h, steps):
        # preparing arguments for functions
        args = [h, alpha, beta, N]
        f = [SIS.susceptible, SIS.infected]

        for i in xrange(0, steps-1):
            y = [S[i], I[i]]
            k_1 = k1(f, y, args, 2)
            k_2 = k2(k_1, f, y, args, 2)

            t[i+1] = t[i] + h
            S[i+1] = S[i] + k_2[0]
            I[i+1] = I[i] + k_2[1]

        plt.plot(t, S, 'r--')
        plt.plot(t, I, 'b--')

    @staticmethod
    # SIS Fourth-Order Runge-Kutta
    def rk_4(alpha, beta, t, S, I, N, h, steps):
        # preparing arguments for functions
        args = [h, alpha, beta, N]
        f = [SIS.susceptible, SIS.infected]

        for i in xrange(0, steps-1):
            y = [S[i], I[i]]

            k_1 = k1(f, y, args, 2)
            k_2 = k2(k_1, f, y, args, 2)
            k_3 = k3(k_2, f, y, args, 2)
            k_4 = k4(k_3, f, y, args, 2)

            t[i+1] = t[i] + h
            S[i+1] = S[i] + k_1[0]/6 + k_2[0]/3 + k_3[0]/3 + k_4[0]/6
            I[i+1] = I[i] + k_1[1]/6 + k_2[1]/3 + k_3[1]/3 + k_4[1]/6

        plt.plot(t, S, 'r')
        plt.plot(t, I, 'b')


class SIR:

    @staticmethod
    def run(alpha, beta, gamma, cities, h, steps, method):
        # duration of SIS is t = steps*h
        t = np.zeros(steps)

        S = np.zeros((steps, len(cities)))
        I = np.zeros((steps, len(cities)))
        R = np.zeros((steps, len(cities)))
        N = np.zeros(len(cities))
        w = np.zeros((len(cities), len(cities)))

        Nmax = 0

        for i in xrange(0, len(cities)):
            I[0][i] = cities[i].infected
            S[0][i] = cities[i].susceptible
            N[i] = cities[i].population

            w[i] = cities[i].w

            # finding maximum population for the y axis of the comparison plot
            if N[i] > Nmax:
                Nmax = N[i]

        method(alpha, beta, gamma, t, S, I, R, N, h, w, steps, cities)

        axes = plt.gca()
        axes.set_ylim([0, Nmax/1000000])
        axes.set_xlim([0, t[steps-1]])


    @staticmethod
    def susceptible(y, args):
        beta = args[2]
        gamma = args[3]
        w = args [4]

        S = y[0]
        I = y[1]
        N = y[3]
        S_vector = y[4]
        city_index = y[7]

        _sum = 0
        for i in xrange(0, w.shape[0]):
            if i != city_index:
                _sum += w[city_index][i] * S_vector[i] - w[i][city_index] * S_vector[city_index]

        return -beta*S*I/N - gamma*S + _sum

    @staticmethod
    def infected(y, args):
        alpha = args[1]
        beta = args[2]
        w = args [4]

        S = y[0]
        I = y[1]
        N = y[3]
        I_vector = y[5]
        city_index = y[7]

        axlad = 0
        _sum = 0
        for i in xrange(0, w.shape[0]):
            if i != city_index:
                _sum += w[city_index][i] * I_vector[i] - w[i][city_index] * I_vector[city_index]

        return beta*S*I/N - alpha*I + _sum

    @staticmethod
    def recovered(y, args):
        alpha = args[1]
        gamma = args[3]
        w = args[4]

        S = y[0]
        I = y[1]
        R_vector = y[6]
        city_index = y[7]

        _sum = 0

        for i in xrange(0, w.shape[0]):
            if i != city_index:
                _sum += w[city_index][i] * R_vector[i] - w[i][city_index] * R_vector[city_index]

        return gamma*S + alpha*I + _sum

    @staticmethod
    # SIR Euler method
    def euler(alpha, beta, gamma, t, S, I, R, N, h, w, steps, cities):
        # preparing arguments for functions
        args = [h, alpha, beta, gamma, w]

        for i in xrange(0, steps-1):
            t[i+1] = t[i] + h

            # for every city
            for j in xrange(0, len(N)):
                y = [S[i][j], I[i][j], R[i][j], N[j], S[i], I[i], R[i]]

                S[i+1][j] = S[i][j] + h * SIR.susceptible(y, args)
                I[i+1][j] = I[i][j] + h * SIR.infected(y, args)
                R[i+1][j] = R[i][j] + h * SIR.recovered(y, args)

        plot(t, S, I, R)


    @staticmethod
    # SIR Second-Order Runge-Kutta
    def rk_2(alpha, beta, gamma, t, S, I, R, N, h, w, steps, cities):
        # preparing arguments for functions
        args = [h, alpha, beta, gamma, w]

        f = [SIR.susceptible, SIR.infected, SIR.recovered]

        for i in xrange(0, steps-1):
            t[i+1] = t[i] + h

            # for every city
            for j in xrange(0, len(N)):
                y = [S[i][j], I[i][j], R[i][j], N[j], S[i], I[i], R[i], j]

                # 3 different derivative values
                k_1 = k1(f, y, args, 3)
                k_2 = k2(k_1, f, y, args, 3)

                S[i+1][j] = S[i][j] + k_2[0]
                I[i+1][j] = I[i][j] + k_2[1]
                R[i+1][j] = R[i][j] + k_2[2]

        plot(t, S, I, R)

    @staticmethod
    # SIR Fourth-Order Runge-Kutta
    def rk_4(alpha, beta, gamma, t, S, I, R, N, h, w, steps, cities):
        # preparing arguments for functions
        args = [h, alpha, beta, gamma, w]

        f = [SIR.susceptible, SIR.infected, SIR.recovered]

        for i in xrange(0, steps-1):
            t[i+1] = t[i] + h

            # for every city
            for j in xrange(0, len(N)):
                y = [S[i][j], I[i][j], R[i][j], N[j], S[i], I[i], R[i], j]

                # 3 different derivative values
                k_1 = k1(f, y, args, 3)
                k_2 = k2(k_1, f, y, args, 3)
                k_3 = k3(k_2, f, y, args, 3)
                k_4 = k4(k_3, f, y, args, 3)

                S[i+1][j] = S[i][j] + k_1[0]/6 + k_2[0]/3 + k_3[0]/3 + k_4[0]/6
                I[i+1][j] = I[i][j] + k_1[1]/6 + k_2[1]/3 + k_3[1]/3 + k_4[1]/6
                R[i+1][j] = R[i][j] + k_1[2]/6 + k_2[2]/3 + k_3[2]/3 + k_4[2]/6

        plot(t, S, I, R)


    @staticmethod
    # SIR emergency policy and deaths
    def sim(alpha, beta, gamma, t, S, I, R, N, h, w, steps, cities):
        # preparing arguments for functions
        args = [h, alpha, beta, gamma, w]

        f = [SIR.susceptible, SIR.infected, SIR.recovered]

        for i in xrange(0, steps-1):
            t[i+1] = t[i] + h


            # for every city
            for j in xrange(0, len(N)):
                y = [S[i][j], I[i][j], R[i][j], N[j], S[i], I[i], R[i], j]

                # 3 different derivative values
                k_1 = k1(f, y, args, 3)
                k_2 = k2(k_1, f, y, args, 3)
                k_3 = k3(k_2, f, y, args, 3)
                k_4 = k4(k_3, f, y, args, 3)

                S[i+1][j] = S[i][j] + k_1[0]/6 + k_2[0]/3 + k_3[0]/3 + k_4[0]/6
                I[i+1][j] = I[i][j] + k_1[1]/6 + k_2[1]/3 + k_3[1]/3 + k_4[1]/6
                R[i+1][j] = R[i][j] + k_1[2]/6 + k_2[2]/3 + k_3[2]/3 + k_4[2]/6

                # update each city with infected and recovered
                cities[j].update(I[i+1][j], R[i+1][j])
                if cities[j].check_emergency(t[i+1]) is True:
                    # close city airport
                    w[j, :] = 0
                    w[:, j] = 0
                    args = [h, alpha, beta, gamma, w]

        plot(t, S, I, R)


def k1(f, y, args, y_num):
    h = args[0]
    k_1 = np.zeros(y_num)

    for i in xrange(0, y_num):
        k_1[i] = h*f[i](y, args)

    return k_1


def k2(k_1, f, y, args, y_num):
    h = args[0]
    k_2 = np.zeros(y_num)

    for i in xrange(0, y_num):
        y[i] += 1/2 * k_1[i]

    for i in xrange(0, y_num):
        k_2[i] = h*f[i](y, args)

    return k_2


def k3(k_2, f, y, args, y_num):
    h = args[0]
    k_3 = np.zeros(y_num)

    for i in xrange(0, y_num):
        y[i] += 1/2 * k_2[i]

    for i in xrange(0, y_num):
        k_3[i] = h*f[i](y, args)

    return k_3


def k4(k_3, f, y, args, y_num):
    h = args[0]
    k_4 = np.zeros(y_num)

    for i in xrange(0, y_num):
        y[i] += k_3[i]

    for i in xrange(0, y_num):
        k_4[i] = h*f[i](y, args)

    return k_4


def plot(t, S, I, R):
    plt.plot(t, S[:, 0] / 1000000, 'r--', linewidth=2, label='S_sto')
    plt.plot(t, I[:, 0] / 1000000, 'r', linewidth=2, label='I_sto')
    plt.plot(t, R[:, 0] / 1000000, 'r', linewidth=1, label='R_sto')
    plt.plot(t, S[:, 1] / 1000000, 'g--', linewidth=2, label='S_cph')
    plt.plot(t, I[:, 1] / 1000000, 'g', linewidth=2, label='I_cph')
    plt.plot(t, R[:, 1] / 1000000, 'g', linewidth=1, label='R_cph')
    plt.plot(t, S[:, 2] / 1000000, 'b--', linewidth=2, label='S_osl')
    plt.plot(t, I[:, 2] / 1000000, 'b', linewidth=2, label='I_osl')
    plt.plot(t, R[:, 2] / 1000000, 'b', linewidth=1, label='R_osl')
    # plt.plot(t, (I[:, 2]+I[:, 1]+I[:, 0]) / 3000000, 'black', linewidth=3, label='I_total')

    Imax = I[:, 2]+I[:, 1]+I[:, 0]

    m = np.amax(Imax/1000000)

    plt.title("SIR, maximum infected = " + str("%.5f" % m) + " million")
    plt.legend(loc = 4)
    plt.xlabel("Time (days)")
    plt.ylabel("Number of persons (millions)")


