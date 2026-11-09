"""Generated Python equations derived from FDA/CiPA; see licenses/."""
import numpy as np
from numpy import exp, log, sqrt, power as pow
from numba import njit

PARAMETER_NAMES = ('A1', 'B1', 'q1', 'A2', 'B2', 'q2', 
                   'A3', 'B3', 'q3', 'A4', 'B4', 'q4', 
                   'A11', 'B11', 'q11', 'A21', 'B21', 
                   'q21', 'A31', 'B31', 'q31', 'A41', 
                   'B41', 'q41', 'A51', 'B51', 'q51', 
                   'A52', 'B52', 'q52', 'A53', 'B53', 
                   'q53', 'A61', 'B61', 'q61', 'A62', 
                   'B62', 'q62', 'A63', 'B63', 'q63', 
                   'Kmax', 'Ku', 'n', 'halfmax', 'Kt', 
                   'Vhalf', 'T', 'timeout', 'starttime')

@njit(cache=True, error_model="numpy")
def evaluate(t, y, p):
    ydot = np.zeros(11)
    yout = np.zeros(0)
    A1 = p[0]
    B1 = p[1]
    q1 = p[2]
    A2 = p[3]
    B2 = p[4]
    q2 = p[5]
    A3 = p[6]
    B3 = p[7]
    q3 = p[8]
    A4 = p[9]
    B4 = p[10]
    q4 = p[11]
    A11 = p[12]
    B11 = p[13]
    q11 = p[14]
    A21 = p[15]
    B21 = p[16]
    q21 = p[17]
    A31 = p[18]
    B31 = p[19]
    q31 = p[20]
    A41 = p[21]
    B41 = p[22]
    q41 = p[23]
    A51 = p[24]
    B51 = p[25]
    q51 = p[26]
    A52 = p[27]
    B52 = p[28]
    q52 = p[29]
    A53 = p[30]
    B53 = p[31]
    q53 = p[32]
    A61 = p[33]
    B61 = p[34]
    q61 = p[35]
    A62 = p[36]
    B62 = p[37]
    q62 = p[38]
    A63 = p[39]
    B63 = p[40]
    q63 = p[41]
    Kmax = p[42]
    Ku = p[43]
    n = p[44]
    halfmax = p[45]
    Kt = p[46]
    Vhalf = p[47]
    T = p[48]
    timeout = p[49]
    starttime = p[50]
    ydot[0] = -(A11*exp(B11*y[10])*y[0]*exp((T-20)*log(q11)/10)-A21*exp(B21*y[10])*y[1]*exp((T-20)*log(q21)/10))+(A51*exp(B51*y[10])*y[2]*exp((T-20)*log(q51)/10)-A61*exp(B61*y[10])*y[0]*exp((T-20)*log(q61)/10))
    ydot[1] = (A11*exp(B11*y[10])*y[0]*exp((T-20)*log(q11)/10)-A21*exp(B21*y[10])*y[1]*exp((T-20)*log(q21)/10))-(A3*exp(B3*y[10])*y[1]*exp((T-20)*log(q3)/10)-A4*exp(B4*y[10])*y[5]*exp((T-20)*log(q4)/10))+(A52*exp(B52*y[10])*y[3]*exp((T-20)*log(q52)/10)-A62*exp(B62*y[10])*y[1]*exp((T-20)*log(q62)/10))
    ydot[2] = -(A1*exp(B1*y[10])*y[2]*exp((T-20)*log(q1)/10)-A2*exp(B2*y[10])*y[3]*exp((T-20)*log(q2)/10))-(A51*exp(B51*y[10])*y[2]*exp((T-20)*log(q51)/10)-A61*exp(B61*y[10])*y[0]*exp((T-20)*log(q61)/10))
    ydot[3] = (A1*exp(B1*y[10])*y[2]*exp((T-20)*log(q1)/10)-A2*exp(B2*y[10])*y[3]*exp((T-20)*log(q2)/10))-(A31*exp(B31*y[10])*y[3]*exp((T-20)*log(q31)/10)-A41*exp(B41*y[10])*y[4]*exp((T-20)*log(q41)/10))-(A52*exp(B52*y[10])*y[3]*exp((T-20)*log(q52)/10)-A62*exp(B62*y[10])*y[1]*exp((T-20)*log(q62)/10))
    ydot[4] = (A31*exp(B31*y[10])*y[3]*exp((T-20)*log(q31)/10)-A41*exp(B41*y[10])*y[4]*exp((T-20)*log(q41)/10))-(A53*exp(B53*y[10])*y[4]*exp((T-20)*log(q53)/10)-A63*exp(B63*y[10])*y[5]*exp((T-20)*log(q63)/10))-(Kmax*Ku*y[9]**n/(y[9]**n+halfmax)*y[4]-Ku*y[7])
    ydot[5] = (A3*exp(B3*y[10])*y[1]*exp((T-20)*log(q3)/10)-A4*exp(B4*y[10])*y[5]*exp((T-20)*log(q4)/10))+(A53*exp(B53*y[10])*y[4]*exp((T-20)*log(q53)/10)-A63*exp(B63*y[10])*y[5]*exp((T-20)*log(q63)/10))-(Kmax*Ku*y[9]**n/(y[9]**n+halfmax)*y[5]-Ku*A53*exp(B53*y[10])*exp((T-20)*log(q53)/10)/(A63*exp(B63*y[10])*exp((T-20)*log(q63)/10))*y[6])
    ydot[6] = (Kmax*Ku*y[9]**n/(y[9]**n+halfmax)*y[5]-Ku*A53*exp(B53*y[10])*exp((T-20)*log(q53)/10)/(A63*exp(B63*y[10])*exp((T-20)*log(q63)/10))*y[6])+(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[6])
    ydot[7] = (Kmax*Ku*y[9]**n/(y[9]**n+halfmax)*y[4]-Ku*y[7])+(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[7])
    ydot[8] = -(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[7])-(Kt/(1+exp(-(y[10]-Vhalf)/6.789))*y[8]-Kt*y[6])
    ydot[9] = 0
    ydot[10] = 0
    return ydot, yout
