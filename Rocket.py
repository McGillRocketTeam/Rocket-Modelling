from Injector import Injector
from Tank import Tank
from CombustionChamber import CombustionChamber
from Nozzle import Nozzle
import matplotlib.pyplot as plt
import numpy as np

#NOTE REGARDING DEBUG VERBOSITY
# 0 = no debug messages
# 1 = some debug messages
# 2 = as many debug messages as possible, report anything that is odd

class Rocket:

    DEBUG_VERBOSITY = 2

    def __init__(self):
        # The order of initialisation here reflects the hierarchy we are using
        # For now, all properties are hard-coded, rather than inputs
        self.Tank = Tank()
        self.Injector = Injector()
        self.CombustionChamber = CombustionChamber()
        self.Nozzle = Nozzle()
        self.t = 0

        #self.T_tank = 0
        #self.rho_tank_liquid = 0
        self.T_cc = 0
        #self.T_post_comb = 0 #cp = combustion products
        self.P_cc = 0 #cp = combustion products

        self.m_dot_ox = 0
        self.m_dot_fuel = 0
        self.m_dot_choke = 0

    def simulate(self, dt=0.001, max_timesteps=1e6):
        # Simulate until reach zero oxidiser/fuel mass, not based on final time
        loop_ctr = 0
        thrust_curve = []

        self.initiate()

        while self.Tank.m_ox > 0 and self.CombustionChamber.inner_radius < self.CombustionChamber.outer_radius \
                and loop_ctr < max_timesteps:

            if self.DEBUG_VERBOSITY > 0:
                print("t = ", dt*loop_ctr)

            # converge is configured to return -1 in order to end simulation
            self.converge()
            self.update(dt)

            loop_ctr += 1

            thrust_curve.append(self.Nozzle.thrust)

        if self.Tank.m_ox < 1e-5:
            print("[Rocket.Simulate] Tank has emptied of oxidiser, we shall ")
        elif self.CombustionChamber.outer_radius - self.CombustionChamber.inner_radius < 1e-5:
            print("[Rocket.Simulate] Fuel grain has burned away, let us arise soaring from its ashes")
        elif loop_ctr > max_timesteps:
            print("[Rocket.Simulate: ERROR, minor] Simulator has exceeded max timesteps without emptying")
        return thrust_curve

    def initiate(self):
        # Calculate m_dot_ox when ignition occurs
        # Eventually we could change this to simulate the startup transient

        # This assumes combustion is already occurring at steady state and
        # that m_dot_ox is determined by the choke
        self.m_dot_ox = self.Nozzle.get_choked_flow_rate(self.CombustionChamber.pressure, self.CombustionChamber.temperature)


        #self.T_tank = self.Tank.T_tank
        #self.rho_tank_liquid = self.Tank.rho_liquid
        self.T_cc = self.CombustionChamber.temperature
        #self.T_post_comb = 0
        self.P_cc = 0



    def update(self, dt):
        self.Tank.update(dt, self.m_dot_ox) # change m_ox
        # Tank.update is configured to set rho = -1 in certain cases to end simulation
        if self.Tank.rho_liquid == -1:
            return
        #self.Injector.update(dt) # should do NOTHING
        self.CombustionChamber.update(dt, self.m_dot_fuel) # change m_fuel & r_fuel
        #self.Nozzle.update(dt) # should do NOTHING

    def converge(self):
        #The second MAJOR function.  After everything that depends EXPLICITLY on time
        # has changed, call this function to "equilibrate" the various components.  This should
        # be done after every timestep
        epsilon = 1000 # Percent change between steps
        epsilon_min = 10 # margin of error in Pascals
        converge_ctr = 0
        while epsilon > epsilon_min:

            # Does nothing, just for stylistic reasons
            self.Tank.converge()

            # Determine oxidiser mass flow rate based on the injector model and relevant conditions
            self.m_dot_ox = self.Injector.converge(self.Tank.T_tank, self.Tank.rho_liquid, \
                                                   self.CombustionChamber.temperature, self.CombustionChamber.pressure)

            # Determine fuel mass flow rate based on CC model
            self.m_dot_fuel = self.CombustionChamber.converge(self.m_dot_ox, self.m_dot_fuel)

            #With 2 mass flow rates, use nozzle to re-set T_CC
            m_dot_choke = self.m_dot_fuel+self.m_dot_ox
            P_cc_old = self.CombustionChamber.pressure
            P_cc_new = self.Nozzle.converge(self.CombustionChamber.temperature, m_dot_choke)
            self.CombustionChamber.pressure = P_cc_new


            #Judge convergence by change in temperature; should approach zero.
            #Other option is to do by pressure, but Temp is more sensitive
            epsilon = abs(P_cc_new - P_cc_old)

            converge_ctr += 1
            if self.DEBUG_VERBOSITY > 1:
                print("***DEBUG*** [Rocket.converge] Steps to convergence =  ", converge_ctr)
                print("***DEBUG*** [Rocket.converge] Convergence epsilon =  ", epsilon)

        self.Nozzle.set_thrust_from_nozzle(m_dot_choke, self.CombustionChamber.temperature, self.CombustionChamber.pressure)

myRocket = Rocket()
dt = 0.000001
max_steps = 1e6
thrust_curve = myRocket.simulate(dt,)

times = np.linspace(0, endpoint=False, num=len(thrust_curve), step=dt )
plt.plot(times, thrust_curve)
