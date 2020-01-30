import matplotlib.pyplot as plt
import math

#Constants independent of the rocket
d_t = .01  #the delta time used to perform the numerical simulation
g = 9.8    #gravitational acceleration

#Rocket specific constants
C_d = .4   #coefficient of drag
A = 0.018  #cross-sectional area
m_init =60 #mass on the launch pad
m_min = 50 #dry mass


#Thrust info (used a static values for now, can easly implement them as varying quantity)
m_dot = 2  #nozzle mass flow rate 
u_init = 2000 #exhaust velocity 






#lists to store values for matplotlib 
h_list = list() #altitude
t_list = list() #time
m_list = list() #mass
v_list = list() #velocity

#test_list = list()

def b(h) : #the b coefficient takes is used to find the quadratic drag to which the rocket is subjected,
           #it dependends on atmospheric density, so on altitude
    if h < 11000 :
        T = 15.04 - 0.00649 * h
        p = 101.29  * ((T+273.1)/288.08) ** 5.256
    
    elif h< 25000:
        T = -55.46
        p = 22.65 * math.exp(1.73-0.000157 *h)
    
    else:
        T = -131.31 + .00299 *h
        p = 2.488 * ((T+273.1)/216.6)**(-11.388)

    rho = p/(.2869 * (T + 273.1))
    print (.5 * rho *C_d * A)
    return .5 * rho *C_d * A

def u(t) : #for now, simply returns the constant exhaust velocity as long as there is still fuel in the tanks
    if m(m_dot, t) > m_min :
        return u_init
    else :
        return 0

def m(m_d, t): #decreases the rocket mass according to the (for now) constant mass flow rate, and floors at the rocket dry mass
    m = m_init - m_d * t 
    if m > m_min :
        return m
    else :
        return m_min

def f (v,h,t): #ODE expression
    if v> 0 :
        return m_dot*(u(t)/m(m_dot,t)) - (b(h)*v**2)/m(m_dot,t) -g #as long as ascending, the drag and gravity a pulling the rocket down
    else :
        return -g + (b(h)*v**2)/m(m_dot,t) # when descending, drag is pulling the rocket up,  and no need for thrust term,
                                           #This is mainly to have a nice graph with a nice visible apogee, the descent is very far from reality (no chutes)

def Euler (v,h, t) : #Euler method to numerically solve ODE
    return v + f(v,h,t)*d_t

#def test (t) :
#    return -g*t**2 + 300*t + 600


#initial conditions
t = 0    #time
v = 0.0001 #velocity (slightly positive to get the simulation going, else multiplication by zero stops the velocity to grow)
h = 0      #intial altitude

#Main loop
while t <  60 : #60 seconds simulation
    #test_list.append(test(t))
    

    #Simulate new state
    m_list.append(m(m_dot,t)) 
    h += v*d_t
    v = Euler(v,h,t)

    #record timetep results
    v_list.append(v)
    h_list.append(h)
    t_list.append(t)

    #advance time
    t += d_t

#Matplotlib stuff
#plt.plot(t_list,test_list, label = "lbl")
plt.plot(t_list,m_list, label = "m")
plt.plot(t_list,h_list,label = "h")
plt.plot(t_list,v_list, label = "v")
plt.legend()
plt.show()      