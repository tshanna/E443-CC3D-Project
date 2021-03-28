
from cc3d.core.PySteppables import *
import numpy as np
import random
import sys
import math

ODE_vars = ['IR','IRa','Tb','Fs','Fsa','C']

small_num = sys.float_info.min
sec_per_mcs = 60
# 1E-11 # 1E-15 by calculation # 1.0 to remove
Mvox = 5E-12
#Mvox = 1.0 
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
        
        secrete = 0
        IL2cm = 0
        fAPC = 0
        Tbcm = 0
        H = 0
        
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
        
        #initialize cells throughout the domain
        i = 0
        while i < 33:
            x1 = random.sample(range(0,201),1)
            print(x1[0])
            x = x1[0]
            y1 = random.sample(range(0,201),1)
            y = y1[0]
            z = 1
            
            # make sure cells are not overwriting eachother
            if self.cell_field[x,y,z] is None:
                cell = self.new_cell(self.NAIVE)
                self.cell_field[x, y, z] = cell
                i += 1
        
        
        L = [3,11,20]
        
        for cell in self.cell_list:
            print(cell.id)
            # initialize the 3 starting APC 
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
        
        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
        # concentrationAtCOM = IL2[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
    
        IL2_secretor = self.get_field_secretor('IL2')
        
        IL2 = self.field.IL2
        
        lamR3 = 0.0
        lamR4 = 0.0
        lam1 = 1E-12
        lamT4 = 0.0
          
        # death of APC
        for cell in self.cell_list_by_type(self.APC):
            
            # movement for APC updated every 90 minutes (mcs)
            if mcs % 90 == 0:
                n = random.uniform(-pi,pi)
                # randomize if movement vector is positive or negative
                sign = random.uniform(-3,3)
                abs_sign = abs(sign)
                
                # ensure there is no divide by 0
                if sign == 0:
                    sign = sign + 1
                    abs_sign = abs(sign)
                
                x_temp = math.cos(n)
                y_temp = math.sin(n)
                
                # r = 20 for APC
                x = 20.0*(x_temp*x_temp)*(sign/abs_sign)
                y = 20.0*(y_temp*y_temp)*(sign/abs_sign)

                cell.lambdaVecX = x
                cell.lambdaVecY = y
            
            # APC death after lifespan is reached
            if mcs >= cell.dict['lifespan']:
                
                cell.targetVolume = 0
        
        for cell in self.cell_list_by_type(self.NAIVE, self.EFFECTOR, self.PREACTIVATED, self.ACTIVATED):
            
            # movement for T cells updated every 90 minutes (mcs)
            if mcs % 90 == 0:
                n = random.uniform(-pi,pi)
                
                # randomize if movement vector is positive or negative
                sign = random.uniform(-3,3)
                abs_sign = abs(sign)
                
                # ensure there is no divide by 0
                if sign == 0:
                    sign = sign + 1
                    abs_sign = abs(sign)
                
                x_temp = math.cos(n)
                y_temp = math.sin(n)
                
                x = 150.0*(x_temp*x_temp)*(sign/abs_sign)
                y = 150.0*(y_temp*y_temp)*(sign/abs_sign)
                
                cell.lambdaVecX = x
                cell.lambdaVecY = y            
            
            # count amount of APC a T cell is in contact with
            fAPC = 0
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                
                if neighbor:
                    if neighbor.type == self.APC:
                        fAPC += 1
            
            # keep track of effector-effector or effector-activated contacts (heaviside function H)
            if cell.type == self.EFFECTOR:
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                    if neighbor:
                        if neighbor.type == self.EFFECTOR or self.ACTIVATED:
                            cell.sbml.dp['H'] = 1
                        else:
                            cell.sbml.dp['H'] = 0
            
            # il2_cm = IL2_secretor.amountSeenByCell(cell)/cell.volume
            
            # concentrationAtCOM = IL2[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
            il2_cm = IL2[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]
            # print("IL2CM VALUE IS:", il2_cm)
            
            #update ode
            cell.sbml.dp['IL2cm'] = il2_cm
            cell.sbml.dp['fAPC'] = fAPC
            
            # second term PDE
            secrete = ( lamR3*(( cell.sbml.dp['IRa'] )/( lamR4 + cell.sbml.dp['IRa'] + small_num)) + lam1 * fAPC ) * ( 1 / (1 + lamT4 * cell.sbml.dp['Tb']) )
            cell.sbml.dp['secrete'] = secrete
            
            #if cell.type is not self.cell_type.Naive:
            if cell.type is self.PREACTIVATED or self.ACTIVATED or self.EFFECTOR:
                # secretion of IL2 by T cells 
                IL2_secretor.secreteOutsideCellAtBoundary(cell, cell.sbml.dp['secrete'])
                
                # Caspase threshold for effector, activated, preactivated
                # if cell.sbml.dp['C'] > 2.63*Mvox:
                   # cell.targetVolume = 0.0
            
            # step the simulation
            sbml_simulator = self.get_sbml_simulator(model_name='dp', cell=cell)
            sbml_simulator.timestep()
            
            # Naive -> PA
            if cell.type == self.NAIVE and fAPC > 0:
                
                cell.type = self.cell_type.Preactivated

            # Preactivated cells stop moving until activated
            if cell.type == self.PREACTIVATED:
                
                cell.lambdaVecX = 0.0
                cell.lambdaVecY = 0.0
                
            # if IL2 threshold reached, PA -> A
            if cell.sbml.dp['IL2cm'] > 7*Mvox:
                cell.type = self.ACTIVATED
            
            # if Tbet threshold reached, A -> E
            if cell.type == self.ACTIVATED:
                if cell.sbml.dp['Tb'] > 40 * Mvox:
                    cell.type = self.cell_type.Effector
            

            ###### Moved to line 222
            # second term PDE
            # secrete = ( lamR3*(( cell.sbml.dp['IRa'] )/( lamR4 + cell.sbml.dp['IRa'] + small_num)) + lam1 * fAPC ) * ( 1 / (1 + lamT4 * cell.sbml.dp['Tb']) )
            # cell.sbml.dp['secrete'] = secrete
            
            # secretion of IL2 by T cells            
            # if cell.type is not self.NAIVE:
                # IL2_secretor.secreteOutsideCellAtBoundary(cell, cell.sbml.dp['secrete'])
            
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
        
        # self.parent_cell.targetVolume = self.cell.volume * k                 

        self.clone_parent_2_child()    

        for v in ODE_vars:
            k = random.uniform(0.7,1.0)
            x = self.parent_cell.sbml.dp[v]
            self.child_cell.sbml.dp[v] = x*(2-k)
            self.parent_cell.sbml.dp[v] = x*k

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        #if self.parent_cell.type==1:
        #    self.child_cell.type=2
        #else:
        #    self.child_cell.type=1
        



