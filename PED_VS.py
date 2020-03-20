# ********************************************************************************************************************************
# *************************** A VIRTUAL SCREENING SOFTWARE BASED ON PED MOLECULAR DESCRIPTORS ************************************
# ********************************************************************************************************************************
# Antonio Oliver Gelabert ( ORCID : http://orcid.org/0000-0001-8571-2733 )
# March, 2020
# For more details and citation : https://www.nature.com/articles/srep43738
# *******************************************************************************************************************************
# PARAMETERS
#
# EXAMPLE USAGE

import os
import numpy as np
import pandas as pd
import argparse

directory_in_str="."
directory = os.fsencode(directory_in_str)

#t0= tm.clock()

# Parsing optional arguments of the program  
ap = argparse.ArgumentParser()
ap.add_argument("-iq", "--iq", required=False, default='query.mol2', help="3D query filename file")
ap.add_argument("-fmin", "--filtmin", required=False, default='1', help="Mininum distance filter (in Amstrongs). Default value: 1 Amstrong")
ap.add_argument("-fmax", "--filtmax", required=False, default='99999999', help="maximum area filter (in Amstrongs). Default value:no limit")
args = vars(ap.parse_args())

# Assign args to program variables

fquery=args["iq"]
fmax=float(args["filtmax"])
fmin=float(args["filtmin"])

mdesc=[]
sim=[]
desf = open('PED_comp_'+fquery+'_vs_mol2files.shd', 'w')
scores= open('PED_comp_'+fquery+'_vs_mol2files_scores.shd', 'w')
desf.write('u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,label,S\n')
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    file=fquery
    
    fat = open("atoms.txt", 'w')
    fat.write("Id Atype x y z Otype Str1 Str2 Charge\n")
    with open(file) as f:
        for line in f:
            if line.rstrip() == "@<TRIPOS>ATOM":
                print("")
                for line in f:
                    if line.rstrip() == "@<TRIPOS>BOND":
                        break
                    fat.write(line.rstrip()+'\n')
    fat.close()
    
    df=pd.read_csv('atoms.txt', delimiter=r"\s+")
    Eij=[]
    for index, row in df.iterrows():
    #    print(row['x'], row['y'],row['z'],row['Charge'])
        for index2, row2 in df.iterrows():
            if(index2>index):
                dij=((row['x']-row2['x'])**2+(row['y']-row2['y'])**2+(row['z']-row2['z'])**2)**0.5
                if(dij > 1.0):
                    Eij.append(row['Charge']*row2['Charge']/dij)

    PED=[]
    Eij.sort(reverse = True)
    PED.extend(Eij[0:6])
    PED.extend(Eij[len(Eij)-6:len(Eij)])

    if filename.endswith(".mol2"):
        file=filename
        fat = open("atoms.txt", 'w')
        fat.write("Id Atype x y z Otype Str1 Str2 Charge\n")
        with open(file) as f:
            for line in f:
                if line.rstrip() == "@<TRIPOS>ATOM":
                    for line in f:
                        if line.rstrip() == "@<TRIPOS>BOND":
                            break
                        fat.write(line.rstrip()+'\n')
        fat.close()

        df=pd.read_csv('atoms.txt', delimiter=r"\s+")

        Eij=[]
        for index, row in df.iterrows():
            for index2, row2 in df.iterrows():
                if(index2>index):
                    dij=((row['x']-row2['x'])**2+(row['y']-row2['y'])**2+(row['z']-row2['z'])**2)**0.5
                    if(dij > fmin) & (dij < fmax):
                        Eij.append(row['Charge']*row2['Charge']/dij)
                
        PED2=[]
        Eij.sort(reverse = True)
        PED2.extend(Eij[0:6])
        PED2.extend(Eij[len(Eij)-6:len(Eij)])
        
        MhD=0.0
        
        for i in range(0,len(PED2)-1):
            MhD=MhD+1/12*np.absolute(PED[i]-PED2[i])
        S=1.0/(1.0+MhD)
        scores.write(filename+','+str(np.round(S,4))+'\n')
        fat.close()
        sim.append(S)
        mdesc.extend(PED)
        desf.write(str(round(PED2[0],4))+','+str(round(PED2[1],4))+','+str(round(PED2[2],4))+','+str(round(PED2[3],4))+','+str(round(PED2[4],4))+','+str(round(PED2[5],4))+','+str(round(PED2[6],4))+','+str(round(PED2[7],4))+','+str(round(PED2[8],4))+','+str(round(PED2[9],4))+','+str(round(PED2[10],4))+','+str(round(PED2[11],4))+', '+filename+','+str(round(S,4))+'\n')

desf.close()
df=pd.read_csv('PED_comp_'+fquery+'_vs_mol2files.shd')
dfs=df.sort_values(by='S', ascending=False) 

desfs = open('PED_comp_'+fquery+'_vs_ranked.shd', 'w')
#for index, row in desfs.iterrows():
desfs.write(str(dfs))
desfs.close()

scores.close()
