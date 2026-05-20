from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
import random
import os
import winsound

# changed: IFNhm, NS1_del, initial_infected_fraction, MOI type, MOI_PFU, IFNn, secretion_rate_IFN_baseline, secretion_rate_IFN
# changed: endocytosis_prob_, GDP
# change between SIP and without SIP: virus secretion, IFN secretion, MOI_mixed=MOI_PFU, Gene defect prob, and file name
# 1 pixel = 1.26 microns^2
# cell diameter = 10um, or 7.94 pixels
# 6 MCS = 1 h
pixel_size = 1.26
min_mcs = 10.0  # min/mcs
hours_mcs = min_mcs / 60.0  # hours/mcs
days_mcs = min_mcs / 1440.0  # day/mcs
secretion_rate_V = 71.6/6. # 1900/6 (fit with SIPs). 1850/6 (fit without SIPs). we mean "virus release rate" #estimated value = 71.6/6.0 = 71.6 viruses released per hour (converted to 10 minutes). 
secretion_rate_IFN = 2.5 # 4.55 (fit with SIPs). 6.8 (fit without SIPs). estimated = 2.5 IFN released [micrograms/mL] every MCS (every 10 minutes). 

secretion_rate_IFN_baseline = 0.0031*10.0 # (estimates 0.0031*0.01 from data) IFN molecules released every MCS (every 10 minutes), decay*IFNic

IRF7_secretion_multiplier = 10 # 10x more IFN molecules released every MCS (every 10 minutes)
IFITM_effectiveness = 0.96 # I estimated 96% effectiveness against virus endocytosis
IFN_stimulation_threshold = 0.00001 # this is to avoid activating network with negligible IFN concentration values, saving time
endocytosis_prob_ = 0.1 # chance of endocytosis every MCS (10 steps after infection)
death_prob = (1./6)*(1./5) # 3.3% chance of death every MCS (30 steps to die, 5h of virus production before death on avg)
infection_initialization_type = "mixed virus sample"    # "test" = 1 infected cell in the middle, 
                                                        # "test_baseline" = no virus, to test baseline IFN secretion, 
                                                        # "mixed virus sample" = increments the field 
IFNn = 4
IFNhm = 1.0 # default 1 [micrograms/mL]
NS1_del = 0 # 0 or 1
Wait_50_steps = 1
#initial_infected_fraction = 0.01
MOI_PFU = 0.0012 # fit value = 0.0012 (0.00006 without defects)
MOI_mixed = MOI_PFU*(100+1.5)/1.5 # infectios units per cell
file_name = f'longNODEFECT2_noPKRapop_IS{secretion_rate_IFN:.3f}_VS{secretion_rate_V:.3f}_MOI{MOI_PFU}_ISB{secretion_rate_IFN_baseline:.6f}'

# GDP = genome defect probabilities
GDP = {"PB1g":0.567211899,
       "PB2g":0.568645964,
       "PAg":0.547128424,
       "NS1g":0.224665344,
       "NPg":0.424879137,
       "HAg":0.464785577,
       "NAg":0.402172369,
       "M1g":0.243308949,
       "M2g":0.101758182,
       "NEPg":0.125294801}
'''GDP = {"PB1g":0.0,
       "PB2g":0.0,
       "PAg":0.0,
       "NS1g":0.0,
       "NPg":0.0,
       "HAg":0.0,
       "NAg":0.0,
       "M1g":0.0,
       "M2g":0.0,
       "NEPg":0.0}'''

virus_IFN_mbs = """

// viral replication phases
// TODO : make innate immunity defenses imperfect to allow disease spread

node vRNPn
{
    logic = $internalized_virus > 0;
    rate_up = @logic ? $vRNP_import : 0.0;
    rate_down = @logic ? 0.0 : 0.0;
}

node Producing
{
    logic = vRNPn && RC;
    rate_up = @logic ? $virus_replication : 0.0;
    rate_down = @logic ? 0.0 : $virus_decay;
}

node Releasing
{
    logic = Producing && HA && NA && M1 && M2 && NEP;
    rate_up = @logic ? 1000.0 : 0.0; // instantaneous transition
    rate_down = @logic ? 0.0 : 1.0;
}

// viral proteins, replication capacity

node PB1
{
    logic = $PB1g > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node PB2
{
    logic = $PB2g > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node PA
{
    logic = $PAg > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node NP
{
    logic = $NPg > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node RC
{
    logic = PB1 && PB2 && PA && NP;
    rate_up = @logic ? $RC_formation : 0.0;
    rate_down = @logic ? 0.0 : $RC_degradation;
}

node NS1
{
    logic = $NS1g > 0 && vRNPn && !(OAS || PKR); //ns1del
    rate_up = @logic*!($NS1_del) ? $protein_production : 0.0;
    rate_down = @logic*!($NS1_del) ? 0.0 : $protein_degradation;
}

node HA
{
    logic = $HAg > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node NA
{
    logic = $NAg > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node M1
{
    logic = $M1g > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node M2
{
    logic = $M2g > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

node NEP
{
    logic = $NEPg > 0 && vRNPn && !(OAS || PKR);
    rate_up = @logic ? $protein_production : 0.0;
    rate_down = @logic ? 0.0 : $protein_degradation;
}

// Sensing viruses, producing IFN

node RIGI
{
    logic = 1; // TODO : NS1 downregulates RIGI
    rate_up = @logic ? 0.0 : 0.0;
    rate_down = @logic ? 0.0 : 0.0;
}

node IRF3
{
    logic = 1; // TODO : RIGI eventually goes low when both vRNPn and RIGI are activated
    rate_up = @logic ? 0.0 : 0.0;
    rate_down = @logic ? 0.0 : 0.0;
}

node TLR7
{
    logic = 1; // TODO : TLR7 goes to 1 immediately after internalized virus and decays with vRNP_import rate, block TLR7 pathway when endocytosis is blocked
    rate_up = @logic ? 0.0 : 0.0;
    rate_down = @logic ? 0.0 : 0.0;
}

node IFNmRNA 
{
    logic = (((IRF3 || IRF7) && (RIGI && vRNPn)) || (TLR7 && IRF7 && vRNPn)) and ! NS1; // NS1 does not allow this to go up, but also not pull it down
    // TODO : TLR7 pathway being active all the time could make NS1 inhibition of RIGI unimportant, fix?
    rate_up = @logic ? $mRNA_upregulation : 0.0;
    rate_down = @logic ? 0.0 : $mRNA_degradation; // TODO : maybe NS1 does bring this down by shutting down host transcription
}

// IFN stimulation and ISGs
// *(1.0-$IFN_hill_fuction)

node IFITM
{
    logic = $IFNc > 0 && !(NS1);
    rate_up = @logic ? $IFITM_upregulation*$IFN_hill_fuction : 0.0; 
    //rate_up = @logic ? 0.0 : 0.0; //change here
    rate_down = @logic ? $IFITM_downregulation : $IFITM_downregulation;
}

node OAS
{
    logic = $IFNc > 0 && !(NS1);
    rate_up = @logic ? $mRNA_upregulation*$IFN_hill_fuction : 0.0; 
    //rate_up = @logic ? 0.0 : 0.0; //change here
    rate_down = @logic ? $mRNA_degradation : $mRNA_degradation;
}

node PKR
{
    logic = $IFNc > 0 && !(NS1);
    rate_up = @logic ? $mRNA_upregulation*$IFN_hill_fuction : 0.0; 
    //rate_up = @logic ? 0.0 : 0.0; //change here
    rate_down = @logic ? $mRNA_degradation : $mRNA_degradation;
}

node IRF7
{
    logic = $IFNc > 0 && !(NS1); 
    rate_up = @logic ? $mRNA_upregulation*$IFN_hill_fuction : 0.0; 
    rate_down = @logic ? $mRNA_degradation : $mRNA_degradation;
}


"""

virus_IFN_mbs_config_unformated = """

$PB1g = 0;
$PB2g = 0;
$PAg = 0;
$NPg = 0;
$NS1g = 0;
$HAg = 0;
$NAg = 0;
$M1g = 0;
$M2g = 0;
$NEPg = 0;

$vRNP_import = 1.0/(<<vRNP_import_time>>*6.0); //1h
$virus_decay = 1.0/(<<intracellular_virus_decay_time>>*6.0); //10h
$protein_degradation = 1.0/(<<protein_degradation_time>>*6.0); //1h
$RC_formation = 1.0/(<<RC_formation_time>>*6.0); //instantaneous
$RC_degradation = 1.0/(<<RC_degradation_time>>*6.0); //1h
$mRNA_degradation = 1.0/(<<mRNA_degradation_time>>*6.0); //12h
$IFITM_downregulation = 1.0/(<<IFITM_downregulation_time>>*6.0); //24h

$virus_replication = 1.0/(<<virus_replication_time>>*6.0); //5h
$protein_production = 1.0/(<<protein_production_time>>*6.0); //6min
$mRNA_upregulation = 1.0/(<<mRNA_upregulation_time>>*6.0); //6h
$IFITM_upregulation = 1.0/(<<IFITM_upregulation_time>>*6.0); //12h

$IFN_hill_fuction = 0.0;
$internalized_virus = 0;
$IFNc = 0.0;
$NS1_del = 0.0;
//$IFNh = 10.0;
//$IFNn = 4.0;

RIGI.istate = 1;
IRF3.istate = 1;
TLR7.istate = 1;
IRF7.istate = 0;
IFNmRNA.istate = 0;
PB1.istate = 0;
PB2.istate = 0;
PA.istate = 0;
NP.istate = 0;
NS1.istate = 0;
HA.istate = 0;
NA.istate = 0;
M1.istate = 0;
M2.istate = 0;
NEP.istate = 0;
RC.istate = 0;
IFITM.istate = 0;
OAS.istate = 0;
PKR.istate = 0;
vRNPn.istate = 0;
Producing.istate = 0;
Releasing.istate = 0;

"""

virus_IFN_mbs_config = """

$PB1g = 0;
$PB2g = 0;
$PAg = 0;
$NPg = 0;
$NS1g = 0;
$HAg = 0;
$NAg = 0;
$M1g = 0;
$M2g = 0;
$NEPg = 0;

$vRNP_import = 1.0/(3*6.0); //3h
$virus_decay = 1.0/(10*6.0); //10h
$RC_formation = 1.0/(0.0001*6.0); //instantaneous
$RC_degradation = 1.0/(10*6.0); //10h
$mRNA_degradation = 1.0/(12*6.0); //12h
$IFITM_downregulation = 1.0/(24*6.0); //24h
$protein_degradation = 1.0/(24*6.0); //24h

$protein_production = 1.0/(1*6.0); //1h
$virus_replication = 1.0/(1*6.0); //1h
$mRNA_upregulation = 1.0/(4*6.0); //4h
$IFITM_upregulation = 1.0/(12*6.0); //12h

$IFN_hill_fuction = 0.0;
$internalized_virus = 0;
$IFNc = 0.0;
$NS1_del = 0.0;
//$IFNh = 10.0;
//$IFNn = 4.0;

RIGI.istate = 1;
IRF3.istate = 1;
TLR7.istate = 1;
IRF7.istate = 0;
IFNmRNA.istate = 0;
PB1.istate = 0;
PB2.istate = 0;
PA.istate = 0;
NP.istate = 0;
NS1.istate = 0;
HA.istate = 0;
NA.istate = 0;
M1.istate = 0;
M2.istate = 0;
NEP.istate = 0;
RC.istate = 0;
IFITM.istate = 0;
OAS.istate = 0;
PKR.istate = 0;
vRNPn.istate = 0;
Producing.istate = 0;
Releasing.istate = 0;

"""

class Viral_Replication_ModelSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        
        #fields
        SteppableBasePy.__init__(self,frequency)
        self.create_scalar_field_cell_level_py("vRNPn")
        self.create_scalar_field_cell_level_py("Producing")
        self.create_scalar_field_cell_level_py("Releasing")
        self.create_scalar_field_cell_level_py("NS1")
        self.create_scalar_field_cell_level_py("IRF3")
        self.create_scalar_field_cell_level_py("RIGI")
        self.create_scalar_field_cell_level_py("IFNmRNA")
        self.create_scalar_field_cell_level_py("IFITM")
        self.create_scalar_field_cell_level_py("IRF7")

    def start(self):

        self.endocytosis_prob = 1

        #fields
        self.vRNPn = self.field.vRNPn
        self.Producing = self.field.Producing
        self.Releasing = self.field.Releasing
        self.NS1 = self.field.NS1
        self.IRF3 = self.field.IRF3
        self.RIGI = self.field.RIGI
        self.IFNmRNA = self.field.IFNmRNA
        self.IFITM = self.field.IFITM
        self.IRF7 = self.field.IRF7
        
        self.secretor_virus = self.get_field_secretor("virus")
        self.secretor_IFN = self.get_field_secretor("IFN")
        self.vfield = self.field.virus
        
        for cell in self.cell_list:
            
            cell.dict["total_virus_released"] = 0
            cell.targetVolume = 50
            cell.lambdaVolume = 10
            
            
        rs = self.simulator.getRandomSeed()
        base_dir = os.path.dirname(os.path.abspath(__file__))
        #output_folder = r'C:\Trabalhos\APL_Project\VRPDM_Continuous_Activation_ISG\output'
        output_folder = os.path.join(base_dir, "output")
        os.makedirs(output_folder, exist_ok=True)
        vknockout = 'WT' if NS1_del==0 else 'NS1del'
        #MOI_ = 'high' if MOI_ else 'low'
        #file_name = file_name 
        global file_name
        count = count_files_with_string(output_folder, file_name)
        file_name = file_name+"_"+str(count)+".txt"
        self.output_path= os.path.join(output_folder, file_name)
        with open(self.output_path, 'w') as file:
            file.write(f"time \t totIFNc \t dead \t infected \t totVirus \t IFNmRNA \t NS1 \t IRF7 \n")
            
    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """
        
        if mcs == 60:
            self.endocytosis_prob = endocytosis_prob_
            
        if (secretion_rate_IFN_baseline or Wait_50_steps) and mcs==50:        
            
            # one infected cell in the middle
            if infection_initialization_type=="test":
                cell = self.cell_field[self.dim.x/2, self.dim.y/2, 0]
                self.infect_cell_PFU(cell)
            if infection_initialization_type=="test_baseline":
                pass
            if infection_initialization_type=="mixed virus sample":
                virus_field = self.field.virus
                #for x, y, z in self.every_pixel():
                for i in range(int(MOI_mixed*len(self.cell_list))):
                    x = np.random.randint(0, self.dim.x)
                    y = np.random.randint(0, self.dim.y)
                    if virus_field[x, y, 0] == 0:
                        virus_field[x, y, 0] += 1.99
                    elif virus_field[x, y, 0] > 0:
                        virus_field[x, y, 0] += 1.00
                    
        #TODO : model medium wash
        
        #dead_cells = len(self.cell_list_by_type(self.D))
        #infected_cells = len(self.cell_list_by_type(self.I))
        #total_cells = len(self.cell_list)
        #fieldIFN = self.field.IFN
        #total_IFN = 0
        #secretor_i = self.get_field_secretor("IFN")
        #total_IFN = secretor_i.totalFieldIntegral()
        #secretor_v = self.get_field_secretor("virus")
        #total_virus = secretor_v.totalFieldIntegral()
        #total_IFN_concentration = total_IFN/(self.dim.x*self.dim.y*pixel_size**3)
        
        #time in MCS, total IFN in micromols/microns^3, total virus per cell
        #print(f"{mcs/6.}\t {total_IFN_concentration} \t {len(self.cell_list_by_type(self.D))/len(self.cell_list)}\t {total_virus/len(self.cell_list)}")
        
        #with open(self.output_path, 'a') as file:
            #file.write(f"{mcs/6.}\t {total_IFN_concentration} \t {len(self.cell_list_by_type(self.D))/len(self.cell_list)}\t {total_virus/len(self.cell_list)}\n")
        
        self.timestep_maboss()
        
        #fields
        self.vRNPn.clear() 
        self.Producing.clear() 
        self.Releasing.clear() 
        self.NS1.clear() 
        self.IRF3.clear()  
        self.RIGI.clear() 
        self.IFNmRNA.clear() 
        self.IFITM.clear() 
        self.IRF7.clear() 
        
        IFNmRNA_cell = 0
        IRF7_cell = 0
        NS1_cell = 0
        
        for cell in self.cell_list_by_type(self.I):
            
            #fields
            vRNPn = int(cell.maboss.VModel['vRNPn'].state) 
            Producing = int(cell.maboss.VModel['Producing'].state) 
            Releasing = int(cell.maboss.VModel['Releasing'].state)
            NS1 = int(cell.maboss.VModel['NS1'].state)
            IRF3 = int(cell.maboss.VModel['IRF3'].state)
            RIGI = int(cell.maboss.VModel['RIGI'].state)
            IFNmRNA = int(cell.maboss.VModel['IFNmRNA'].state)
            IFITM = int(cell.maboss.VModel['IFITM'].state)
            PKR = int(cell.maboss.VModel['PKR'].state)
            IRF7 = int(cell.maboss.VModel['IRF7'].state)
            
            #fields
            self.vRNPn[cell] = vRNPn 
            self.Producing[cell] = Producing 
            self.Releasing[cell] = Releasing 
            self.NS1[cell] = NS1 
            self.IRF3[cell] = IRF3 
            self.RIGI[cell] = RIGI 
            self.IFNmRNA[cell] = IFNmRNA 
            self.IFITM[cell] = IFITM 
            self.IRF7[cell] = IRF7 
            
            IFNc = self.secretor_IFN.amountSeenByCell(cell)/cell.volume
            cell.maboss.VModel.network.symbol_table["IFNc"] = IFNc
            cell.maboss.VModel.network.symbol_table["IFN_hill_fuction"] = IFNc**IFNn / (IFNc**IFNn+IFNhm**IFNn)
            #cell.maboss.VModel.network.symbol_table["IFN_hill_fuction"] = np.heaviside(IFNc-1, 1)
            
            viral_exposure = self.secretor_virus.amountSeenByCell(cell)
            # virus endocytosis
            # TODO increase chance of endocytosis depending on viral concentration or endocytose more than 1 virus
            if viral_exposure > 1 and np.random.uniform() < self.endocytosis_prob and not Releasing:
                self.infect_cell(cell)
                # no extra internalization because network is boolean
                self.secretor_virus.secreteInsideCell(cell, -1./cell.volume)
            # virus secretion
            if Releasing:
                secreted_virus = Producing * secretion_rate_V
                res = self.secretor_virus.secreteInsideCellTotalCount(cell, secreted_virus/cell.volume)
                cell.dict["total_virus_released"] += res.tot_amount
                #print(res.tot_amount, cell.dict["total_virus_released"])
                # no reduction of internal virus on secretion because network is boolean
            # IFN secretion
            if IFNmRNA:
                # secreted_IFN = secretion_rate_IFN # OLD
                secreted_IFN = secretion_rate_IFN * (1 + IRF7 * (IRF7_secretion_multiplier-1))
                self.secretor_IFN.secreteInsideCell(cell, secreted_IFN/cell.volume)
            # cell death
            #if (Producing or (PKR and not NS1)) and np.random.uniform() < death_prob:
            if Producing and np.random.uniform() < death_prob:
                self.delete_maboss_from_cell(cell=cell, model_name="VModel")
                cell.type = self.D
                
            if IFNmRNA:
                IFNmRNA_cell += 1
            if IRF7:
                IRF7_cell += 1
            if NS1:
                NS1_cell += 1
                
        for cell in self.cell_list_by_type(self.U):
            # virus endocytosis
            # TODO increase chance of endocytosis depending on viral concentration or endocytose more than 1 virus
            viral_exposure = self.secretor_virus.amountSeenByCell(cell)
            if viral_exposure > 1 and np.random.uniform() < self.endocytosis_prob:
                self.infect_cell(cell)
                self.secretor_virus.secreteInsideCell(cell, -1./cell.volume)
            
            # Baseline IFN secretion
            self.secretor_IFN.secreteInsideCell(cell, secretion_rate_IFN_baseline/cell.volume)
            
            # IFN stimulation
            IFN_exposure = self.secretor_IFN.amountSeenByCell(cell)
            if IFN_exposure > IFN_stimulation_threshold:
                self.stimulate_cell(cell)
                
        for cell in self.cell_list_by_type(self.E):
            
            #fields
            vRNPn = int(cell.maboss.VModel['vRNPn'].state)
            Producing = int(cell.maboss.VModel['Producing'].state) 
            Releasing = int(cell.maboss.VModel['Releasing'].state)
            NS1 = int(cell.maboss.VModel['NS1'].state)
            IRF3 = int(cell.maboss.VModel['IRF3'].state)
            RIGI = int(cell.maboss.VModel['RIGI'].state)
            IFNmRNA = int(cell.maboss.VModel['IFNmRNA'].state)
            IFITM = int(cell.maboss.VModel['IFITM'].state)
            IRF7 = int(cell.maboss.VModel['IRF7'].state)
            PKR = int(cell.maboss.VModel['PKR'].state)
            #fields
            self.vRNPn[cell] = vRNPn 
            self.Producing[cell] = Producing 
            self.Releasing[cell] = Releasing 
            self.NS1[cell] = NS1 
            self.IRF3[cell] = IRF3 
            self.RIGI[cell] = RIGI 
            self.IFNmRNA[cell] = IFNmRNA 
            self.IFITM[cell] = IFITM 
            self.IRF7[cell] = IRF7 
            
            IFNc = self.secretor_IFN.amountSeenByCell(cell)/cell.volume
            cell.maboss.VModel.network.symbol_table["IFNc"] = IFNc
            cell.maboss.VModel.network.symbol_table["IFN_hill_fuction"] = IFNc**IFNn / (IFNc**IFNn+IFNhm**IFNn)
            #print(IFNc, IFNhm, IFNc**IFNn / (IFNc**IFNn+IFNhm**IFNn))
            #cell.maboss.VModel.network.symbol_table["IFN_hill_fuction"] = np.heaviside(IFNc-1,1)
            
            viral_exposure = self.secretor_virus.amountSeenByCell(cell)
            # TODO : Cells with IFITM still endocytose?
            if viral_exposure > 1 and np.random.uniform() < self.endocytosis_prob * (1 - IFITM_effectiveness*IFITM):
                self.infect_cell(cell)
                self.secretor_virus.secreteInsideCell(cell, -1./cell.volume)
                
            # cytotoxicity of IFN
            # if PKR and np.random.uniform() < death_prob:
            #     self.delete_maboss_from_cell(cell=cell, model_name="VModel")
            #     cell.type = self.D
                
            # Baseline IFN secretion
            self.secretor_IFN.secreteInsideCell(cell, secretion_rate_IFN_baseline/cell.volume)
                
            if IFNmRNA:
                IFNmRNA_cell += 1
            if IRF7:
                IRF7_cell += 1
                
        # TODO : model medium wash
        dead_cells = len(self.cell_list_by_type(self.D))
        infected_cells = len(self.cell_list_by_type(self.I))
        total_cells = len(self.cell_list)
        fieldIFN = self.field.IFN
        total_IFN = 0
        secretor_i = self.get_field_secretor("IFN")
        total_IFN = secretor_i.totalFieldIntegral()
        secretor_v = self.get_field_secretor("virus")
        total_virus = secretor_v.totalFieldIntegral()
        total_IFN_concentration = total_IFN/(self.dim.x*self.dim.y)
        # time in MCS, total IFN in micromols/microns^3, total virus per cell
        #print(f"{mcs/6.}\t {total_IFN_concentration} \t {len(self.cell_list_by_type(self.D))/len(self.cell_list)}\t {total_virus/len(self.cell_list)}")
        
        with open(self.output_path, 'a') as file:
            file.write(f"{mcs/6.-50/6.} \t {total_IFN_concentration} \t {len(self.cell_list_by_type(self.D))/len(self.cell_list)} \t {len(self.cell_list_by_type(self.I))/len(self.cell_list)} \t {total_virus/len(self.cell_list)} \t {IFNmRNA_cell/len(self.cell_list)} \t {NS1_cell/len(self.cell_list)} \t {IRF7_cell/len(self.cell_list)}\n")
        
        
    def finish(self):
        
        for cell in self.cell_list:
            print(cell.dict["total_virus_released"])
        
        winsound.Beep(440, 500)
        pass

    def on_stop(self):
        
        pass
        
    def stimulate_cell(self, cell):
        if cell.type == self.U:
            cell.type = self.E
            self.add_maboss_to_cell(cell=cell,model_name="VModel",
                                    bnd_str=virus_IFN_mbs,
                                    cfg_str=virus_IFN_mbs_config,
                                    time_step=1.0,time_tick=1.0,
                                    seed=random.randint(0, int(1E9)))
    def infect_cell(self, cell):
        if cell.type == self.U:
            self.add_maboss_to_cell(cell=cell,model_name="VModel",
                                    bnd_str=virus_IFN_mbs,
                                    cfg_str=virus_IFN_mbs_config,
                                    time_step=1.0,time_tick=1.0,
                                    seed=random.randint(0, int(1E9)))
        cell.type = self.I
        cell.maboss.VModel.network.symbol_table["internalized_virus"] = 1
        cell.maboss.VModel.network.symbol_table["PB1g"] += np.random.uniform()>GDP["PB1g"]
        cell.maboss.VModel.network.symbol_table["PB2g"] += np.random.uniform()>GDP["PB2g"]
        cell.maboss.VModel.network.symbol_table["PAg"] += np.random.uniform()>GDP["PAg"]
        cell.maboss.VModel.network.symbol_table["NPg"] += np.random.uniform()>GDP["NPg"]
        cell.maboss.VModel.network.symbol_table["NS1g"] += np.random.uniform()>GDP["NS1g"]
        cell.maboss.VModel.network.symbol_table["HAg"] += np.random.uniform()>GDP["HAg"]
        cell.maboss.VModel.network.symbol_table["NAg"] += np.random.uniform()>GDP["NAg"]
        cell.maboss.VModel.network.symbol_table["M1g"] += np.random.uniform()>GDP["M1g"]
        cell.maboss.VModel.network.symbol_table["M2g"] += np.random.uniform()>GDP["M2g"]
        cell.maboss.VModel.network.symbol_table["NEPg"] += np.random.uniform()>GDP["NEPg"]
        cell.maboss.VModel.network.symbol_table["NS1_del"] = NS1_del
    
        # if not np.prod(list(cell.dict["GP"].values())):
            # # estimate viral protein defects and update viral protein integrity inside the cell
            # cell.dict["GP"] = {key : cell.dict["GP"].get(key,0)+1*(np.random.uniform()>GDP.get(key,0))
                                # for key in set(cell.dict["GP"])}
        # if np.prod(list(cell.dict["GP"].values())):
            # cell.maboss.VModel.network.symbol_table['internalized_virus'] = 1
            
    def infect_cell_PFU(self, cell):
        if cell.type == self.U:
            self.add_maboss_to_cell(cell=cell,model_name="VModel",
                                        bnd_str=virus_IFN_mbs,
                                        cfg_str=virus_IFN_mbs_config,
                                        time_step=1.0,time_tick=1.0,
                                        seed=random.randint(0, int(1E9)))
        cell.type = self.I
        cell.maboss.VModel.network.symbol_table["internalized_virus"] = 1
        cell.maboss.VModel.network.symbol_table["PB1g"] += 1
        cell.maboss.VModel.network.symbol_table["PB2g"] += 1
        cell.maboss.VModel.network.symbol_table["PAg"] += 1
        cell.maboss.VModel.network.symbol_table["NPg"] += 1
        cell.maboss.VModel.network.symbol_table["NS1g"] += 1
        cell.maboss.VModel.network.symbol_table["HAg"] += 1
        cell.maboss.VModel.network.symbol_table["NAg"] += 1
        cell.maboss.VModel.network.symbol_table["M1g"] += 1
        cell.maboss.VModel.network.symbol_table["M2g"] += 1
        cell.maboss.VModel.network.symbol_table["NEPg"] += 1
        cell.maboss.VModel.network.symbol_table["NS1_del"] = NS1_del
        
        if NS1_del: cell.maboss.VModel.network.symbol_table["NS1g"] = 0
        
def count_files_with_string(folder_path, search_string):
    count = 0

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path):
            if search_string in filename:
                count += 1

    return count
