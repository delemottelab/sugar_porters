import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances 
from numpy.linalg import norm

##Get gate distances as an array. Gate_EC and gate_IC shape should be : \
##gate_EC = [(TM1b start,TM1b end), (TM7b start, TM7b end)] \
##gate_IC = [(TM4b start, TM4b end), (TM10b start, TM10b end)]\
##Refs for some sugar porters (make sure to check numbering in ref structure): \
##GLUT1: gate_EC = [(29,37), (288,295)], gate_IC = [(137,146), (385,394)]\
##GLUT3:gate_EC = [(29,35), (288,294)], gate_IC = [(135,144), (383,391)]\
##GLUT5:gate_EC = [(30,37), (289,295)], gate_IC = [(136,145), (386,394)]\
##PfHT1: gate_EC = [(45, 51), (311, 317)],gate_IC = [(148, 154), (412, 418)]\
##XylE: 

##Get distances to calculate TR angle as an array. For TR angle one should take i and i+4 residues at TM2 start, TM2 end, TM8 start, TM8 end. So the shape should be:
##TR_angle=[(TM2start,TM2start+4),(TM2end-4,TM2end),(TM8start,TM8start+4),(TM8end-4,TM8end)]
##Refs for some sugar porters (make sure to check numbering in ref structure): \
##GLUT1:
##GLUT3:TR_angle=[(58,62),(82,86),(305,309),(325,329)]
##GLUT5:
##PfHT1:TR_angle=[(72,76),(96,100),(327,331),(351,355)]
##XylE:


def make_gate_arr (md_uni, gate_EC, gate_IC):


    gate_EC_dists = []
    gate_IC_dists = []
    
    for timestep in md_uni.trajectory:
        tm1 = md_uni.select_atoms('resid %i-%i' %(gate_EC[0][0], gate_EC[0][1])).center_of_mass()
        tm7 = md_uni.select_atoms('resid %i-%i' %(gate_EC[1][0], gate_EC[1][1])).center_of_mass()
        tm4 = md_uni.select_atoms('resid %i-%i' %(gate_IC[0][0], gate_IC[0][1])).center_of_mass()
        tm10 = md_uni.select_atoms('resid %i-%i' %(gate_IC[1][0], gate_IC[1][1])).center_of_mass()    


        gate_EC_dists.append(float(distances.distance_array(tm1, tm7)))
        gate_IC_dists.append(float(distances.distance_array(tm4, tm10)))
   # print("returning EC gate, IC gate dists")
    gate_EC_dists = np.array(gate_EC_dists)
    gate_IC_dists = np.array(gate_IC_dists)
    return gate_EC_dists, gate_IC_dists


def TRangle(md_uni,TR_angle):
    Theta=[]
    for timestep in md_uni.trajectory:
    
        """Calculate TR angle"""
        A = md_uni.select_atoms('resid %i-%i and (name CA)'%(TR_angle[0][0],TR_angle[0][1])).center_of_geometry()
        B = md_uni.select_atoms('resid %i-%i and (name CA)'%(TR_angle[1][0],TR_angle[1][1])).center_of_geometry()
        C = md_uni.select_atoms('resid %i-%i and (name CA)'%(TR_angle[2][0],TR_angle[2][1])).center_of_geometry()
        D = md_uni.select_atoms('resid %i-%i and (name CA)'%(TR_angle[3][0],TR_angle[3][1])).center_of_geometry()


        AB = A - B
        CD = C - D

        theta = np.arccos(np.dot(AB, CD)/(norm(AB)*norm(CD)))
        Theta.append(theta)
    return np.rad2deg(Theta)

def make_dist_arr(md_uni,res1,res2):
    distCC = []

    residue1 = md_uni.select_atoms("resid "+str(res1))
    residue2 = md_uni.select_atoms("resid "+str(res2))
    for timestep in md_uni.trajectory:
        distCC.append(np.min(distances.distance_array(residue1.positions, residue2.positions)))
        
        
    return distCC

def classify(md_uni,gate_EC,gate_IC,TR_angle):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tr_angle=TRangle(md_uni,TR_angle)
    ec_com,ic_com=make_gate_arr(md_uni, gate_EC, gate_IC)

    tr_angle = np.asarray(tr_angle)
    ec_com = np.asarray(ec_com)
    ic_com = np.asarray(ic_com)
    ax1.plot(tr_angle, '-',color='whitesmoke', label='TR angle')
    ax2.plot(ec_com,'--',color='lightgrey',alpha=0)
    ax2.plot(ic_com,'.',color='gray',alpha=0)
    ax1.plot(ec_com,'--',color='lightgrey',label='EC gate')
    ax1.plot(ic_com,'.',color='gray',label='IC gate')

    ax1.plot([0,np.shape(tr_angle)[0]],[30]*2,'k--',alpha=0.2)
    ax1.plot([0,np.shape(tr_angle)[0]],[40]*2,'k-.',alpha=0.2)

    outopen=np.where((ec_com>14) & (tr_angle<30) & (ic_com<=16))
    inwopen=np.where((ic_com>16) & (tr_angle>30) & (ec_com<=12)& (tr_angle<40))
    outocc=np.where((ec_com<=14) & (tr_angle<30) & (ic_com<=16) & (ec_com>11))
    inwocc=np.where((ic_com<=16) & (tr_angle>30) & (ec_com<=14) & (ic_com>13))    
    occ=np.where((ic_com<=13) & (ec_com<=11))
    
    unclass=np.where((tr_angle>40))

    ax1.plot(outopen[0],tr_angle[outopen[0]],'b*',label='out.op.')
    ax1.plot(outocc[0],tr_angle[outocc[0]],'c.',label='out.occ.')
    ax1.plot(occ[0],tr_angle[occ[0]],'kd',label='occ.')

    ax1.plot(inwocc[0],tr_angle[inwocc[0]],'+',color='pink',label='inw.occ.')
    ax1.plot(inwopen[0],tr_angle[inwopen[0]],'r^',label='inw.op.')
    ax1.legend(loc='center left', bbox_to_anchor=(1.5, 0.5))  # Legend outside plot box
    
    ax1.set_xlabel('Frames')
    ax2.set_ylabel('Gates distances ($\AA$)')
    ax1.set_ylabel('TR angle (°)')
    
    ax1.set_ylim(8,48)
    ax2.set_ylim(8,48)

    colorp=tr_angle.copy()
    colorp[outopen[0]]=0
    colorp[outocc[0]]=1
    colorp[occ[0]]=2
    colorp[inwocc[0]]=3
    colorp[inwopen[0]]=4
    colorp[unclass[0]]=5
    #colorp[channellike[0]]=5
    colorp[np.where(colorp>5)[0]]=5
    return colorp