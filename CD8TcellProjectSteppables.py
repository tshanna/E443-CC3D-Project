
from cc3d.core.PySteppables import *
import numpy as np
import random

class CD8TcellProjectSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        any code in the start function runs before MCS=0
        """
        
        model_string = """
        
        E1: -> IR ; lamR1 * ( fAPC ) + ( lamR2 + muDIL2 ) * IRa - ( muAIL2 )*( IL2cm )* IR ;          // Non-activated IL2R
        E2: IR -> ;  kR * IR ;                                                                  // Decay of Non-activated IL2R
        E3: -> IRa ; ( muAIL2 )*( IL2cm )* IR - (muDIL2 * IRa) ;                                  // Activated IL2R
        E4: IRa -> ; kRa * IRa ;                                                                // Decay of Activated IL2R
        E5: -> Tb ; lamT1 * ( fAPC ) + lamT2 * ( (Tb) / (lamT3 + Tb) ) * Tb ;                         // Tbet
        E6: Tb -> ; kT * Tb ;                                                                   // Decay of Tbet
        E7: -> Fs ; lamF + ( muDF * Fsa ) - H * ( muAF * Tbcm * Fs ) ;                              // Nn-activated Fas
        E8: Fs -> ; kF * Fs ;                                                                   // Decay of Fas
        E9: -> Fsa ; H * ( muAF * Tbcm * Fs ) - ( muDF * Fsa ) ;                                  // Activated Fas
        E10: Fsa -> ; kFa * Fsa ;                                                               // Decay of Fas
        E11: -> C ; lamc1 * ( 1/ (1 + (lamc2*IRa) ) ) * ( 1/ (1 + (lamc3*fAPC) ) ) + ( lamc4 * Fsa ) ;  // Caspase
        E12: C -> ; kC * C ;                                                                    // Decay of Caspase
        
        // E13: -> IL2 ; D*IL2 + ( lamR3 * ( IRa/ ( lamR4 + IRa)) + (lam1*fAPC)) * ( 1/ (1 + (lamT4*Tb))) - delta*IL2 ;    // concentration of IL2
        
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
        
        // For E13
        
        // D = 0.2 ; // pixels/MCS
        // delta = 0.008 ; // 1/min
        // lamR3 = 0.0 ; // M 1/min
        // lamR4 = 0.0 ; // M
        // lamT4 = 0.0 ; // 1/M
        // lam1 = 10E-12 ; // m 1/min
        
        end"""
        
        #options = {'relative': 1e-10, 'absolute': 1e-12}
        #self.set_sbml_global_options(options)
        
        # generate 3 random numbers between 1,33 inclusive
        # put the 3 numbers in a list
        
        L = [3,11,20]

        #L = random.sample(range(1,34),3)
        
        # self.field.IL2
        
        # IL2_secretor = self.get_field_secretor("IL2")
        # in step()
        
        IL2_secretor = self.field.IL2
        
        for cell in self.cell_list:
            print(cell.id)
            if cell.id in L:
                cell.type = self.APC

            if cell.type == self.APC:

                cell.targetVolume = 250
                cell.lambdaVolume = 10
                
            else:
                
                cell.targetVolume = 25
                cell.lambdaVolume = 10
                self.add_antimony_to_cell(model_string=model_string,
                                  model_name='dp',
                                  cell=cell,
                                  step_size=1E-2)
                

    def step(self,mcs):
        
        self.timestep_sbml()

        IL2_secretor = self.field.IL2
        
        for cell in self.cell_list:
            # loop over cells to calculate secretion (2nd term PDE) 
            # ...
            pass
            
        for cell in self.cell_list_by_type(self.NAIVE):
            # naive -> PA on contact with APC
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.APC:
                        cell.type = 2

            
            
                
            
    def finish(self):
        
        """
        Finish Function is called after the last MCS
        """

    def on_stop(self):
        # this gets called each time user stops simulation
        return



