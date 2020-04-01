# ********************************************************************************************************************************
# *************************** COMPUTING OF PED MOLECULAR DESCRIPTORS *************************************************************
# ********************************************************************************************************************************
# Antonio Oliver Gelabert ( ORCID : http://orcid.org/0000-0001-8571-2733 )
# March, 2020
# For more details and citation : https://www.nature.com/articles/srep43738
# *******************************************************************************************************************************
# INPUT : a set of mol2 structure files in the same directory
# OUTPUT : a set of 12 descriptors provided by PED methodology that uses spatial coordinates and molecular partial charges
# PARAMETERS : fmin : filter of minimum distance between atoms (by default, 1 Amstrong)
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
desf = open('PED_comp_'+fquery+'.shd', 'w')
desf.write('u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,label,S\n')
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".mol2"):
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
        #    print(row['x'], row['y'],row['z'],row['Charge'])
            for index2, row2 in df.iterrows():
                if(index2>index):
                    dij=((row['x']-row2['x'])**2+(row['y']-row2['y'])**2+(row['z']-row2['z'])**2)**0.5
                    if(dij > fmin)&(dij < fmax):
                        Eij.append(row['Charge']*row2['Charge']/dij*14.4)

        Eij.sort(reverse = False)
        PED=[]
        PED.extend(Eij[0:6])
        PED.extend(Eij[len(Eij)-6:len(Eij)])
        desf.write(str(round(PED[0],4))+','+str(round(PED[1],4))+','+str(round(PED[2],4))+','+str(round(PED[3],4))+','+str(round(PED[4],4))+','+str(round(PED[5],4))+','+str(round(PED[6],4))+','+str(round(PED[7],4))+','+str(round(PED[8],4))+','+str(round(PED[9],4))+','+str(round(PED[10],4))+','+str(round(PED[11],4))+', '+filename+'\n')

desf.close()

