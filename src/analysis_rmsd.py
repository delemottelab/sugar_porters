import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms

### Here we define a function to calculate rmsd for each sugar porter, the residue numbers for each group need to be adjusted accordingly

###GLUT3
def get_RMSDglut3(ref, u, what_group):
    what_group_dict = {'all':'all',
                       'backbone':'backbone',
                      'skip_ICH5':'backbone and not resid 454-470',
                      'only_ICH5':'resid 454-470',
                       'ICH no ICH5':'backbone and resid 206-263',
                      'skip_all_ICH':'backbone and not resid 206-263+454-470',
                      'sugar':'resname BGLC'}
    
    selection=what_group_dict[what_group]

    R = rms.RMSD(reference = ref, 
                    atomgroup = u, 
                    select = selection, 
                    center = True
                   )
        
    R.run()
    
    return R.results.rmsd



###PfHT1
def get_RMSDpfht1(ref, u, what_group):
    what_group_dict = {'all':'all',
                       'backbone':'backbone',
                      'skip_ICH5':'backbone and not resid 480-498',
                      'only_ICH5':'resid 480-498',
                       'ICH no ICH5':'backbone and resid 230-288',
                      'skip_all_ICH':'backbone and not resid 230-288+480-498',
                'TM':'backbone and resid 25-56+73-101+104-122+127-154+159-185+200-229+289-324+328-353+357-380+386-417+422-452+455-475',
                     'TM6':'backbone and resid 200-229',
                     'TM4':'backbone and resid 127-154',
                     'TM10':'backbone and resid 386-417',
                     'noTM6':'backbone and not resid 200-229',
                     'noTM6noICH5':'backbone and not resid 200-229+480-489',
                     'sugar':'resname BGLC'}
    
    selection=what_group_dict[what_group]

    R = rms.RMSD(reference = ref, 
                    atomgroup = u, 
                    select = selection, 
                    center = True
                   )
        
    R.run()
    
    return R.results.rmsd