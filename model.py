from dataclasses import dataclass
from typing import List

dt = 1 * (10 ** -3)
dx = 5 * (10 ** -3)
dy = 5 * (10 ** -3)
# Mesenchymal-like cancer cell diffusion coeff
D_M = 1 * (10 ** -4)
# Epithelial-like cancer cell diffusion coeff
D_E = 5 * (10 ** -5)
# Mesenchymal haptotatic sensitivity coeff
phi_M = 5 * (10 ** -4)
# Epithelial haptotatic sensitivity coeff
phi_E = 5 * (10 ** -4)
# MMP-2 diffusion coefficient
D_m = 1 * (10 ** -3)
# MMP-2 production rate
theta = 0.195
# MMP-2 decay rate
MMP2_decay = 0.1
# ECM degradation rate by MT1-MMP
gamma_1 = 1
# ECM degradation rate by MMP-2
gamma_2 = 1
# Time CTC's spend in vasculature
T_v = 0.18
# epithelial doubling time
T_M = 3
# mesenchymal doubling time
T_E = 2
# single CTC survival probability
P_s = 5 * (10 ** -4)
# cluster survival probability
P_C = 2.5 * (10 ** -4)
# extravasion prob to bones
E_1 = 0.5461
# extravasion prob to lungs
E_2 = 0.2553
# extravasion prob to liver
E_3 = 0.1986
# Number of iterations to run
ITERATIONS = 10


@dataclass
class Grid:
    """
    Represents 201 x 201 primary and secondary grids
    -5 dictionaries for mesenchymal, epithelial, MM2, ECM concentration, blood vessels (with keeping track of normal or ruptured)
    -key = location (tuple)
    -value = concentration (int, how many of them there are)
    -Specifc PDES (methods)
    """
    mes: dict = dict()  # mesenchymal
    epi: dict = dict()  # epithelial
    MMP2: dict = dict()  # matrix metalloproteinase-2
    ECM: dict = dict()  # extracellular matrix
    bv: dict = dict()  # blood vessels
    clusters: List[tuple] = []  # clusters leaving

    def initalize_primary():
        """
        adds mesechmmal cancer cells cluster in middle
        add epi
        """

    def initalize_secondary():
        """
        adds normal vessels to grid
        """
        pass

    def update_MMP2(self) -> None:
        pass

    def update_ECM(self) -> None:
        pass

    def update_mesechymal(self) -> None:
        """
        simulate movement
        simulate mitosis
        """
        pass

    def update_epithelial(self) -> None:
        """
        simulate movement
        simulate mitosis
        """
        pass

    def update_all(self) -> None:  # list of cells that move into vasculature
        """
        update MMP2 then ECM then mesechymal then epithelial
        figures out what cells are on blood vessels and removes them from dicts
        adds clusters leaving to self.cluster
        """
        self.update_MMP2()
        self.update_ECM()
        self.update_mesechymal()
        self.update_epithelial()


"""
Set up vascular class 
- List of cell clusters in vascular (#mesenchyml, #epithelia, time)
"""


@dataclass
class Vascular:
    clusters: List[tuple] = []  # (# mes, # epi, time)

    def update_all(self, newClusters) -> None:
        """
        add the new cells to the current clusters
        iterate through the clusters
            increment time
            update the
        """


@dataclass
class Model:
    """
    One primary breast grid object : (grid class)
    One vascular object : (vascular class)
    One secondy bones grid object : (grid class)
    One secondary lungs grid object : (grid class)
    One secondary liver grid object : (grid class)
    move time step method?
    """
    breast: Grid = Grid()
    vascular: Vascular = Vascular()
    bones: Grid = Grid()
    lungs: Grid = Grid()
    liver: Grid = Grid()

    def initialize(self) -> None:
        """
        initalize breast, vascular, bones, lungs, livers
        populations breast with blood vessels, and cancer cells
        """

    def update(self) -> None:
        """
        updates breast, vascular, bones, lungs, liver
        """
        prev = self.breast()
        # prev_vasc = self.vasc()
        self.breast.update_all(prev)
        self.vascular.update_all(prev.clusters)
        self.bones.update_all(prev)
        self.lungs.update_all(prev)
        self.liver.update_all(prev)


def main():
    # create a model
    model = Model()
    # initialize model (populating with blood vessels and cells)
    model.initialize()
    # start with primary grid
    for t in range(ITERATIONS):  # change time later!!! # iterations
        # update the model
        model.update()
        # we are done!!!


if __name__ == "__main__":
    main()
