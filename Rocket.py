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

        self.T_tank = 0
        self.rho_tank = 0
        self.T_cc = 0
        self.rho_cc = 0
        self.T_post_comb = 0 #cp = combustion products
        self.P_cc = 0 #cp = combustion products

        self.m_dot_ox = 0
        self.m_dot_fuel = 0
        self.m_dot_choke = 0

    def simulate(self, dt=0.001, max_timesteps=1e6):
        # Simulate until reach zero oxidiser/fuel mass, not based on final time
        loop_ctr = 0
        thrust_curve = []
        while self.Tank.oxidizer_mass > 0 and self.CombustionChamber.inner_radius < self.CombustionChamber.outer_radius \
                and loop_ctr < max_timesteps:

            if self.DEBUG_VERBOSITY > 0:
                print("t = ", dt*loop_ctr)

            self.converge()
            self.update(dt)

            loop_ctr += 1

            thrust_curve.append(self.Nozzle.thrust)

        if self.Tank.oxidizer_mass < 1e-5:
            print("[Rocket.Simulate] Tank has emptied of oxidiser")
        elif self.CombustionChamber.outer_radius - self.CombustionChamber.inner_radius < 1e-5:
            print("[Rocket.Simulate] Fuel grain has burned away")
        elif loop_ctr > max_timesteps:
            print("[Rocket.Simulate: ERROR, minor] Simulator has exceeded max timesteps without emptying")
        return thrust_curve


    def update(self, dt):
        self.T_tank, self.rho_tank = self.Tank.update(dt, self.m_dot_ox) # change m_ox
        self.Injector.update(dt) # should do NOTHING
        self.CombustionChamber.update(dt, self.m_dot_fuel) # change m_fuel & r_fuel
        self.Nozzle.update(dt) # should do NOTHING

    def converge(self):
        epsilon = 1000 # Percent change between steps
        epsilon_min = 1e-6
        converge_ctr = 1
        while epsilon < epsilon_min:

            self.T_tank, self.rho_tank = self.Tank.converge()

            self.m_dot_ox = self.Injector.converge(self.T_tank, self.rho_tank, self.T_cc, self.rho_cc)

            self.m_dot_fuel, self.T_post_comb = self.CombustionChamber.converge(self.m_dot_ox, self.P_cc)

            m_dot_choke = self.Nozzle.get_mass_flow_rate(self.P_cc, self.T_post_comb)
            m_dot_actual = self.m_dot_ox + self.m_dot_fuel

            self.P_cc = self.Nozzle.converge(m_dot_actual, self.T_post_comb, self.P_cc)

            epsilon = abs(self.m_dot_ox + self.m_dot_fuel - m_dot_choke)

            converge_ctr += 1
        if self.DEBUG_VERBOSITY > 1:
            print("***DEBUG*** [Rocket.converge] Steps to convergence =  ", converge_ctr)

myRocket = Rocket()
dt = 0.001
max_steps = 1e6
thrust_curve = myRocket.simulate()
times = np.linspace(0, endpoint=False, num=len(thrust_curve), step=dt )
plt.plot(times, thrust_curve)
