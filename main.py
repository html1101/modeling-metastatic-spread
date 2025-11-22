from dataclasses import dataclass, field
from typing import List
import random

from sympy.polys.numberfields.subfield import field_isomorphism_factor
length = 201 #board length
width = 201 #board width

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
# Max steps in vasaclature
vasc_time = 6
# Probability cluster disaggregates
P_d = .3 # idk tbh the paper doesn't say...

@dataclass
class PrimaryGrid: 
    """
    Represents 201 x 201 primary and secondary grids
    -5 dictionaries for mesenchymal, epithelial, MM2, ECM concentration, blood vessels (with keeping track of normal or ruptured)
    -key = location (tuple)
    -value = concentration (int, how many of them there are) 
    -Specifc PDES (methods)
    """
    mes: dict = field(default_factory = dict) #mesenchymal
    epi: dict = field(default_factory = dict) #epithelial
    MMP2: dict = field(default_factory = dict) # matrix metalloproteinase-2
    ECM: dict = field(default_factory = dict) # extracellular matrix
    bv: dict = field(default_factory = dict) #blood vessels
    clusters : List[tuple] = field(default_factory = list) #clusters leaving

    #Rishika
    def initalize_primary(self):
        """
        adds mesechmmal cancer cells cluster in middle
        add epi
        """
        #I made my boundry flux thing assuming zero based indicies - Mia 

    #Carmen
    def update_MMP2(self) -> None:
        pass

    #Carmen
    def update_ECM(self) -> None:
        pass

    #Mia
    def update_mesechymal(self) -> None:
        """
        simulate movement 
        simulate mitosis 
        """
        new_mes = (self.mes).copy()
        for (i,j), concentration in self.mes.values(): 
            for _ in range(0,concentration):
                ECM_conc_left = self.MM2[(i-1,j)] if i > 0 else 0
                ECM_conc_right = self.MM2[(i+1,j)] if i < 200 else 0
                ECM_conc_down = self.MM2[(i,j-1)] if j > 0 else 0
                ECM_conc_up = self.MM2[(i,j+1)]  if j < 200 else 0
                z = random.randint() 
                prob_move_left = (dt / (dx)^2) * (D_M - (phi_M/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_right = (dt / (dx)^2) * (D_M + (phi_M/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_down = (dt / (dx)^2) * (D_M + (phi_M/4) *  (ECM_conc_up- ECM_conc_down))
                prob_move_up  = (dt / (dx)^2) * (D_M - (phi_M/4) *  (ECM_conc_up- ECM_conc_down))
                if z < prob_move_left: #cell moves left
                    if i > 0:
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i-1,j)] = new_mes[(i-1,j)] + 1
                    else: 
                        continue #no change if on left boundry 
                elif z < prob_move_right + prob_move_right:  #cell moves right
                    if i < 200: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i+1,j)] = new_mes[(i-1,j)] + 1
                    else: 
                        continue #no change if on right boundry 
                elif z < prob_move_down + prob_move_right + prob_move_left: #cell moves down
                    if j > 0: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i,j-1)] = new_mes[(i,j-1)] + 1
                    else: 
                        continue #no change if on lower boundry 
                elif z < prob_move_down + prob_move_right + prob_move_left + prob_move_up:
                    if j < 200: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i,j+1)] = new_mes[(i,j+1)] + 1   
                    else: 
                        continue #no change if on upper boundry 
        self.mes = new_mes           

    #Mia
    def update_epithelial(self) -> None:
        """
        simulate movement 
        simulate mitosis 
        """
        new_epi = (self.epi).copy()
        for (i,j), concentration in self.epi.values(): 
            for _ in range(0,concentration): 
                ECM_conc_left = self.MM2[(i-1,j)] if i > 0 else 0
                ECM_conc_right = self.MM2[(i+1,j)] if i < 200 else 0
                ECM_conc_down = self.MM2[(i,j-1)] if j > 0 else 0
                ECM_conc_up = self.MM2[(i,j+1)]  if j < 200 else 0
                z = random.randint() 
                prob_move_left = (dt / (dx)^2) * (D_E - (phi_E/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_right = (dt / (dx)^2) * (D_E + (phi_E/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_down = (dt / (dx)^2) * (D_E + (phi_E/4) *  (ECM_conc_up- ECM_conc_down))
                prob_move_up  = (dt / (dx)^2) * (D_E - (phi_E/4) *  (ECM_conc_up- ECM_conc_down))
                if z < prob_move_left: 
                    if i > 0: 
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i-1,j)] = new_epi[(i-1,j)] + 1
                    else: 
                        continue
                elif z < prob_move_right + prob_move_left: 
                    if i < 200:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i+1,j)] = new_epi[(i-1,j)] + 1
                    else: 
                        continue
                elif z < prob_move_down + prob_move_right + prob_move_left: 
                    if j > 0:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i,j-1)] = new_epi[(i,j-1)] + 1
                    else: 
                        continue
                elif z < prob_move_down + prob_move_right + prob_move_left + prob_move_up:
                    if j < 200:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i,j+1)] = new_epi[(i,j+1)] + 1    
                    else: 
                        continue
        self.epi = new_epi  
        
    #Sarah  
    def update_all(self) -> None: # list of cells that move into vasculature
        """
        update MMP2 then ECM then mesechymal then epithelial
        figures out what cells are on blood vessels and removes them from dicts 
        adds clusters leaving to self.cluster 
        """
        self.update_MMP2()
        self.update_ECM()
        self.update_mesechymal()
        self.update_epithelial()
        #need to do: figures out what cells are on blood vessels and removes them from dicts adds clusters leaving to self.cluster 

class SecondaryGrid: 
    """
    Represents 201 x 201 primary and secondary grids
    -5 dictionaries for mesenchymal, epithelial, MM2, ECM concentration, blood vessels)
    -key = location (tuple)
    -value = concentration (int, how many of them there are) 
    -Specifc PDES (methods)
    """
    mes: dict = dict() #mesenchymal
    epi: dict = dict() #epithelial
    MMP2: dict = dict() # matrix metalloproteinase-2
    ECM: dict = dict() # extracellular matrix
    bv: dict = dict() #blood vessels 
    clusters : List[tuple] = [] #clusters leaving 

    #Rishika
    def initalize_secondary(): 
        """
        adds normal vessels to grid
        """
        pass

    #Carmen
    def update_MMP2(self) -> None:
        pass

    #Carmen
    def update_ECM(self) -> None:
        pass

    #Mia
    def update_mesechymal(self) -> None:
        """
        simulate movement 
        simulate mitosis 
        """
        new_mes = (self.mes).copy()
        for (i,j), concentration in self.mes.values(): 
            for _ in range(0,concentration):
                ECM_conc_left = self.MM2[(i-1,j)] if i > 0 else 0
                ECM_conc_right = self.MM2[(i+1,j)] if i < 200 else 0
                ECM_conc_down = self.MM2[(i,j-1)] if j > 0 else 0
                ECM_conc_up = self.MM2[(i,j+1)]  if j < 200 else 0
                z = random.randint() 
                prob_move_left = (dt / (dx)^2) * (D_M - (phi_M/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_right = (dt / (dx)^2) * (D_M + (phi_M/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_down = (dt / (dx)^2) * (D_M + (phi_M/4) *  (ECM_conc_up- ECM_conc_down))
                prob_move_up  = (dt / (dx)^2) * (D_M - (phi_M/4) *  (ECM_conc_up- ECM_conc_down))
                if z < prob_move_left: #cell moves left
                    if i > 0:
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i-1,j)] = new_mes[(i-1,j)] + 1
                    else: 
                        continue #no change if on left boundry 
                elif z < prob_move_right + prob_move_right:  #cell moves right
                    if i < 200: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i+1,j)] = new_mes[(i-1,j)] + 1
                    else: 
                        continue #no change if on right boundry 
                elif z < prob_move_down + prob_move_right + prob_move_left: #cell moves down
                    if j > 0: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i,j-1)] = new_mes[(i,j-1)] + 1
                    else: 
                        continue #no change if on lower boundry 
                elif z < prob_move_down + prob_move_right + prob_move_left + prob_move_up:
                    if j < 200: 
                        new_mes[(i,j)] = new_mes[(i,j)] -1
                        new_mes[(i,j+1)] = new_mes[(i,j+1)] + 1   
                    else: 
                        continue #no change if on upper boundry 
        self.mes = new_mes           

    #Mia
    def update_epithelial(self) -> None:
        """
        simulate movement 
        simulate mitosis 
        """
        new_epi = (self.epi).copy()
        for (i,j), concentration in self.epi.values(): 
            for _ in range(0,concentration): 
                ECM_conc_left = self.MM2[(i-1,j)] if i > 0 else 0
                ECM_conc_right = self.MM2[(i+1,j)] if i < 200 else 0
                ECM_conc_down = self.MM2[(i,j-1)] if j > 0 else 0
                ECM_conc_up = self.MM2[(i,j+1)]  if j < 200 else 0
                z = random.randint() 
                prob_move_left = (dt / (dx)^2) * (D_E - (phi_E/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_right = (dt / (dx)^2) * (D_E + (phi_E/4) *  (ECM_conc_right - ECM_conc_left))
                prob_move_down = (dt / (dx)^2) * (D_E + (phi_E/4) *  (ECM_conc_up- ECM_conc_down))
                prob_move_up  = (dt / (dx)^2) * (D_E - (phi_E/4) *  (ECM_conc_up- ECM_conc_down))
                if z < prob_move_left: 
                    if i > 0: 
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i-1,j)] = new_epi[(i-1,j)] + 1
                    else: 
                        continue
                elif z < prob_move_right + prob_move_left: 
                    if i < 200:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i+1,j)] = new_epi[(i-1,j)] + 1
                    else: 
                        continue
                elif z < prob_move_down + prob_move_right + prob_move_left: 
                    if j > 0:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i,j-1)] = new_epi[(i,j-1)] + 1
                    else: 
                        continue
                elif z < prob_move_down + prob_move_right + prob_move_left + prob_move_up:
                    if j < 200:
                        new_epi[(i,j)] = new_epi[(i,j)] -1
                        new_epi[(i,j+1)] = new_epi[(i,j+1)] + 1    
                    else: 
                        continue
        self.epi = new_epi  
        
    #Sarah    
    def update_all(self, clusters: List[tuple[int, int]]) -> None: # list of cells that move into vasculature
        """
        update MMP2 then ECM then mesechymal then epithelial
        figures out what cells are coming in through the blood vessels and add them to the dicts
        """
        self.update_MMP2()
        self.update_ECM()
        self.update_mesechymal()
        self.update_epithelial()
        #need to do: figures out what cells are coming in through the blood vessels and add them to the dicts


"""
Set up vascular class 
- List of cell clusters in vascular (#mesenchyml, #epithelia, time)
"""
#Funny Joelle 😂😂😂
@dataclass
class Vascular:
    clusters: List[tuple[int, int, int]] = field(default_factory=list) # (# mes, # epi, time)
    bones: List[tuple[int, int]] = field(default_factory=list)
    lungs: List[tuple[int, int]] = field(default_factory=list)
    liver: List[tuple[int, int]] = field(default_factory=list)

    def update_all(self, primary) -> None:
        """
        add the new cells to the current clusters
        iterate through the clusters
        increment time
        assign leaving clusters to new locations, and store in class attribute
        """
        newClusters = primary.clusters
        self.clusters.append(newClusters)
        updatedClusters = []
        leavingVascular = []
        for cluster in self.clusters:
            mes = cluster[0]
            epi = cluster[1]
            time = cluster[2]
            if time == vasc_time:
                leavingVascular.append(cluster)
                continue
            time +=1 #maybe this should be dt? -mia 
            # checking to see if they disaggregate
            if time >= vasc_time/2 and mes+epi >1:
                disaggregate_mes = 0
                for _ in range(mes):
                    r = random.random()
                    if r < P_d:
                        disaggregate_mes += 1
                disaggregate_epi = 0
                for _ in range(epi):
                    r = random.random()
                    if r < P_d:
                        disaggregate_epi += 1
                remaining_mes = mes - disaggregate_mes
                remaining_epi = epi - disaggregate_epi
                if remaining_mes + remaining_epi >=2:
                    updatedClusters.append((remaining_mes, remaining_epi, time))
                for i in range(disaggregate_mes):
                    updatedClusters.append((1, 0, time))
                for i in range(disaggregate_epi):
                    updatedClusters.append((0, 1, time))
            else:
                updatedClusters.append((mes, epi, time))
        bones = []
        lungs = []
        liver = []
        for cluster in leavingVascular:
            if cluster[0] == 0 or cluster[1] == 0:
                prob = P_s # it's a single
            else:
                prob = P_C # it's a cluster
            r = random.random()
            if r < prob:
                continue
            else:
                newLoc = random.randrange(1, 4)
                if newLoc == 1:
                    bones.append((cluster[0], cluster[1]))
                elif newLoc == 2:
                    lungs.append((cluster[0], cluster[1]))
                elif newLoc == 3:
                    liver.append((cluster[0], cluster[1]))

        self.bones = bones
        self.lungs = lungs
        self.liver = liver

@dataclass
class Model: 
    """
    One primary breast grid object : (primary grid class)
    One vascular object : (vascular class)
    One secondy bones grid object : (secondary grid class)
    One secondary lungs grid object : (secondary grid class)
    One secondary liver grid object : (secondary grid class)
    move time step method?
    """
    breast: PrimaryGrid = PrimaryGrid()
    vascular: Vascular = Vascular() 
    bones : SecondaryGrid = SecondaryGrid()
    lungs : SecondaryGrid = SecondaryGrid()
    liver : SecondaryGrid = SecondaryGrid()
    
    #Sarah
    def initialize(self) -> None: 
        """
        initalize breast, vascular, bones, lungs, livers 
        populations breast with blood vessels, and cancer cells         
        """
    
    #Baby
    def update(self) -> None: 
        """
        updates breast, vascular, bones, lungs, liver 
        """
        prev_primary_grid = self.breast
        self.breast.update_all()
        prev_vasc = self.vascular
        self.vascular.update_all(prev_primary_grid)
        # passing in the previous cells that will migrate to bones, lungs, liver
        self.bones.update_all(prev_vasc.bones)
        self.lungs.update_all(prev_vasc.lungs)
        self.liver.update_all(prev_vasc.liver)


    
def main():
    # create a model
    model = Model()
    # initialize model (populating with blood vessels and cells)
    model.initialize()
    # start with primary grid
    for t in range(ITERATIONS): # change time later!!! # I think this shoudl be iterations / dt? # iterations
        #update the model 
        model.update()
        #we are done!!!

if __name__ == "__main__":
    main()