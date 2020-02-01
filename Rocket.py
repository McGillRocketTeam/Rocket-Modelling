from Injector import Injector
from Tank import Tank
from CombustionChamber import CombustionChamber
from Nozzle import Nozzle


class Rocket:

    def __init__(self, tank_init=[1, 1], injector_init=[1, 1], cc_init=[1, 1], nozzle_init=[1, 1]):
        # The order of initialisation here reflects the hierarchy we are using
        self.Tank = Tank(tank_init)
        self.Injector = Injector(injector_init)
        self.CombustionChamber = CombustionChamber(cc_init)
        self.Nozzle = Nozzle(nozzle_init)
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

    def simulate(self, dt=0.001):
        # Simulate until reach zero oxidiser/fuel mass, not based on final time
        loop_ctr = 0
        max_loops = 1e4
        while self.Tank.oxidizer_mass > 0 and self.CombustionChamber.inner_radius < self.CombustionChamber.outer_radius \
                and loop_ctr<max_loops:
            self.update(dt)
            self.converge()
            loop_ctr += 1

    def update(self, dt):
        self.T_tank, self.rho_tank = self.Tank.update(dt, self.m_dot_ox) # change m_ox
        self.Injector.update(dt) # should do NOTHING
        self.CombustionChamber.update(dt) # change r_fuel
        self.Nozzle.update(dt) # should do NOTHING

    def converge(self):
        epsilon = 1000 # Percent change between steps
        epsilon_min = 1e-3
        while epsilon < epsilon_min:
            self.T_tank, self.rho_tank = self.Tank.converge()
            self.m_dot_ox = self.Injector.converge(self.T_tank, self.rho_tank, self.T_cc, self.rho_cc)
            self.m_dot_fuel, self.T_post_comb = self.CombustionChamber.converge(self.m_dot_ox, self.P_cc)
            m_dot_choke = self.Nozzle.get_mass_flow_rate(self.P_cc, self.T_post_comb)
            m_dot_actual = self.m_dot_ox + self.m_dot_fuel
            self.P_cc = self.Nozzle.converge(m_dot_actual, self.T_post_comb, self.P_cc)
            epsilon = abs(self.m_dot_ox + self.m_dot_fuel - m_dot_choke)