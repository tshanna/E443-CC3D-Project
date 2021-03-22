
from cc3d.core.PySteppables import *
import numpy as np
import random
import sys
import math

small_num = sys.float_info.min
sec_per_mcs = 60
pi = math.pi

class CD8TcellProjectSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        
        model_string = """
        
        E1: -> IR ; lamR1 * ( fAPC ) + ( lamR2 + muDIL2 ) * IRa - ( muAIL2 )*( IL2cm )* IR - ( kR * IR ) ;           // Non-activated IL2R
        E2: -> IRa ; ( muAIL2 )*( IL2cm )* IR - (muDIL2 * IRa) - ( kRa * IRa ) ;                                     // Activated IL2R
        E3: -> Tb ; lamT1 * ( fAPC ) + lamT2 * ( (Tb) / (lamT3 + Tb) ) * Tb - ( kT * Tb ) ;                          // Tbet
        E4: -> Fs ; lamF + ( muDF * Fsa ) - H * ( muAF * Tbcm * Fs ) - ( kF * Fs ) ;                                 // Nn-activated Fas
        E5: -> Fsa ; H * ( muAF * Tbcm * Fs ) - ( muDF * Fsa ) - ( kFa * Fsa ) ;                                     // Activated Fas
        E6: -> C ; lamc1 * ( 1/ (1 + (lamc2*IRa) ) ) * ( 1/ (1 + (lamc3*fAPC) ) ) + ( lamc4 * Fsa ) - ( kC * C ) ;  // Caspase
        
        // Decay rates
        
        kR = 0.0029 ; // 1/min
        kRa = 0.0029 ; // 1/min
        kT = 0.0035 ; // 1/min
        kF = 0.0047 ; // 1/min
        kFa = 0.0047 ; // 1/min
        kC = 0.0038 ; // 1/min
        
        // Feedback strengths
        
        lamR1 = 0.0158 ; // M 1/min
        lamR2 = 0.001 ; // 1/min
        lamT1 = 0.01 ; // M 1/min
        lamT2 = 0.004 ; // 1/min
        lamT3 = 0.01 ; // M
        lamc1 = 0.01 ; // M 1/min
        lamc2 = 100 ; // 1/M
        lamc3 = 0.01 ; // N/A
        lamc4 = 0.004 ; // 1/min
        
        // Association/Dissociation rates
        
        muAIL2 = 6E8 ; // 1/(M*min)
        muDIL2 = 0.006 ; // 1/min
        muAF = 0.86E5 ; // 1/(M*min)
        muDF = 0.006 ; // 1/min
        
        lamF = 3.47E-5 ; // M 1/min
        
        // Initial conditions
        
        IL2cm = 0
        fAPC = 0
        Tbcm = 0
        H = 1
        
        IR0 = 0
        IR = IR0
        IRa0 = 0
        IRa = IRa0
        Tb0 = 0
        Tb = Tb0
        Fs0 = 0
        Fs = Fs0
        Fsa0 = 0
        Fsa = Fsa0
        C0 = 0
        C = C0
        
        """
        
        i = 0
        
        while i < 33:
            x1 = random.sample(range(0,201),1)
            print(x1[0])
            x = x1[0]
            y1 = random.sample(range(0,201),1)
            y = y1[0]
            z = 1
            
            if self.cell_field[x,y,z] is None:
                cell = self.new_cell(self.NAIVE)
                self.cell_field[x, y, z] = cell
                i += 1
        
        
        L = [3,11,20]
        
        for cell in self.cell_list:
            print(cell.id)
            if cell.id in L:
                cell.type = self.APC
                
                if cell.type == self.APC:

                    cell.targetVolume = 250
                    cell.lambdaVolume = 10
                    
                    life = random.sample(range(32,41),1)
                    
                    time = life[0] * sec_per_mcs
                    
                    cell.dict['lifespan'] = time
                
            else:
                
                cell.targetVolume = 25
                cell.lambdaVolume = 10
                
                self.add_antimony_to_cell(model_string=model_string,
                                  model_name='dp',
                                  cell=cell,
                                  step_size=1E-2)
                

    def step(self,mcs):
        
        IL2_secretor = self.get_field_secretor('IL2')
        
        lamR3 = 0.0
        lamR4 = 0.0
        lam1 = 1E-12
        lamT4 = 0.0
        
        
        
        # death of APC
        for cell in self.cell_list_by_type(self.APC):
            
            if mcs % 10 == 0:
                n = random.uniform(-pi,pi)
                
                sign = random.uniform(-3,3)
                abs_sign = abs(sign)
                
                if sign == 0:
                    sign = sign + 1
                    abs_sign = abs(sign)
                
                x_temp = math.cos(n)
                y_temp = math.sin(n)
                
                x = 0.4*(x_temp*x_temp)*(sign/abs_sign)
                y = 0.4*(y_temp*y_temp)*(sign/abs_sign)

                cell.lambdaVecX = x
                cell.lambdaVecY = y
            
            if mcs >= cell.dict['lifespan']:
                
                cell.targetVolume = 0
        
        for cell in self.cell_list_by_type(self.NAIVE, self.EFFECTOR, self.PREACTIVATED, self.ACTIVATED):
            
            if mcs % 10 == 0:
                n = random.uniform(-pi,pi)
                
                sign = random.uniform(-3,3)
                abs_sign = abs(sign)
                
                if sign == 0:
                    sign = sign + 1
                    abs_sign = abs(sign)
                
                x_temp = math.cos(n)
                y_temp = math.sin(n)
                
                x = 3.0*(x_temp*x_temp)*(sign/abs_sign)
                y = 3.0*(y_temp*y_temp)*(sign/abs_sign)
                
                cell.lambdaVecX = x
                cell.lambdaVecY = y
            
            if cell.type == self.PREACTIVATED:
                
                cell.lambdaVecX = 0.0
                cell.lambdaVecY = 0.0
            
            il2_cm = IL2_secretor.amountSeenByCell(cell)
            
            fAPC = 0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                
                if neighbor:
                    if neighbor.type == self.APC:
                        fAPC += 1
            
            #update ode
            cell.sbml.dp['IL2cm'] = il2_cm
            cell.sbml.dp['fAPC'] = fAPC
            
            sbml_simulator = self.get_sbml_simulator(model_name='dp', cell=cell)

            sbml_simulator.timestep()
            
            # Naive -> PA
            if cell.type == self.NAIVE and fAPC > 0:
                
                cell.type = self.cell_type.Preactivated

            # uptake of IL2 and threshold for preactivated -> activated
            if cell.type == self.PREACTIVATED:
                
                if il2_cm > 7:
                    cell.type = self.ACTIVATED

            # second term PDE
            secrete = ( lamR3*(( cell.sbml.dp['IRa'] )/( lamR4 + cell.sbml.dp['IRa'] + small_num)) + lam1 * fAPC ) * ( 1 / (1 + lamT4 * cell.sbml.dp['Tb']) )
            
            # secretion of IL2 by T cells            
            if cell.type is not self.APC or self.NAIVE:
                IL2_secretor.secreteOutsideCellAtBoundary(cell, secrete)
            
    def finish(self):
        
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return



class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list_by_type(self.EFFECTOR, self.ACTIVATED):
            # Effector and activated t cells divide every ~8 hours
            if mcs % 480 == 0:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        
        # one cell inherits k = [0.7,1.0] fraction of parent cell
        # other cell has 2 - k 
        
        # k = random.uniform(0.7,1.0)
        # self.parent_cell.targetVolume = self.cell.volume * k                 

        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        #if self.parent_cell.type==1:
        #    self.child_cell.type=2
        #else:
        #    self.child_cell.type=1
        



