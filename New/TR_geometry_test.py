import argparse
from pathlib import Path

import numpy as np
from numpy.linalg import norm
import pandas as pd
import pymol
from pymol import cmd

import Bio
from Bio import PDB
from Bio.PDB import Polypeptide
from Bio.PDB.PDBParser import PDBParser

import scipy
from scipy.spatial.transform import Rotation
import sklearn
from sklearn.decomposition import PCA
from skimage.measure import CircleModel,EllipseModel

def circle_fit(C):
    def to_minimize(center,points=C):
        Ri=np.mean([norm(center-point) for point in points])
        return np.array([np.abs(norm(center-point)-Ri) for point in points])
    x0=(C[0]+C[-1])/2
    c=scipy.optimize.leastsq(to_minimize,x0)[0]
    del to_minimize
    return c

def orthogonalize(v,n):
    n_orth=n/norm(n)
    v_orth = v-np.dot(v,n_orth)*n_orth
    v_orth /= norm(v_orth)
    return v_orth

def dihedral_angle(v1,v2,n):
    n /= norm(n)
    v1_to_orth=v1/norm(v1)
    v2_to_orth=v2/norm(v2)
    v1_orth=orthogonalize(v1_to_orth,n)
    v2_orth=orthogonalize(v2_to_orth,n)
    dirn=np.cross(v1_orth,v2_orth)
    dirn /= norm(dirn)
    return get_angle(v1_orth,v2_orth),np.sign(np.dot(n,dirn))

def get_angle(v1,v2):
    return np.arccos(np.dot(v1,v2)/(norm(v1)*norm(v2)))

def create_list(list_indexes):   #Used to process the unit_def argument
    list_indexes=list_indexes.split(',')
    L=[]
    for item in list_indexes:
        item=item.strip().split('_')
        item=(int(item[0]),int(item[1]))
        L.append(item)
    return L
        


#---------------------------------------------------------------------------------------------------------------------------------

arg_parser=argparse.ArgumentParser()   #Parse command line
arg_parser.add_argument('filepath',action='store',help='Path to input file')
arg_parser.add_argument('chain',action='store',help='Chain')
arg_parser.add_argument('unit_def',action='store',help='Unit limits, written as s1_e1,s2_e2,...')
arg_parser.add_argument('-ins',action='store',help='Starts and ends of insertions, formatted like the units')
arg_parser.add_argument('-o',action='store',help='Output file if desired, ex. Outputs\my_pdb.csv')
arg_parser.add_argument('--draw',action='store_true',help='Use if Pymol drawing is desired')
args=arg_parser.parse_args()
filepath=args.filepath
list_indexes=args.unit_def

#---------------------------------------------------------------------------------------------------------------------------------

parser=PDBParser(QUIET=True)   #Parse .pdb file, define units, calculate center of each unit
structure=parser.get_structure('structure',Path(filepath))

chain=structure[0][args.chain]
units_ids=create_list(list_indexes)
try:   #If we have insertions, we make sure to remove them from the structure
    ins_ids=create_list(args.ins)
    units=[]
    for limits in units_ids:
        a,b=limits
        unit_ins=[ins for ins in ins_ids if a<=ins[0]<=b or a<=ins[1]<=b]
        unit=[]
        for residue in chain:
            res_id=residue.get_id()[1]
            res_in_ins=np.any([ins[0]<=res_id<=ins[1] for ins in unit_ins])
            if a<=res_id<=b and not res_in_ins and Polypeptide.is_aa(residue):
                unit.append(residue)
        units.append(unit)
                   
except Exception:
    units=[]    
    for limits in units_ids:
        a,b=limits
        unit=[]
        for residue in chain:
            if a<=residue.get_id()[1]<=b:
                if Polypeptide.is_aa(residue):
                    unit.append(residue)
        units.append(unit)

units_coords=[]   #For each unit, store coordinate of each CA atom
region=[]
for unit in units:
    ca_coords=[residue['CA'].get_coord() for residue in unit]
    units_coords.append(ca_coords)
    region.extend(ca_coords)

geometric_centers=[sum(coords)/len(coords) for coords in units_coords]  #Define geometric center for each unit
N=len(geometric_centers)
#----------------------------------------------------------------------------------------------------------------------------------

#CURVATURE
def rotation_fit(C=geometric_centers,window=6):
    data=np.asarray(C)
    N=len(data)
    index_list=[]
    rotation_list=[]
    centers_list=[]
    res_list=[]
    for i in range(max(N-window+1,1)):  #Fit circle to each window
        indexes=[*range(i,min(i+window,N))]
        data_to_fit=data[indexes]
        pca=PCA()
        pca.fit(data_to_fit)
        data_transformed=pca.transform(data_to_fit)


        #Delete third axis, then fit ellipse
        circle_data=np.delete(data_transformed,2,1)
        circle=CircleModel() 
        circle.estimate(circle_data)
        res=circle.residuals(circle_data)
        c=list(circle.params[0:2])
        #c=circle_fit(circle_data).tolist()
        c.append(0)
        centers_list.append(pca.inverse_transform(c))
        index_list.append(indexes)
        res_list.append(np.sum(np.square(res)))
        
    def_centers=np.empty((N-1,3))
    best_res=np.full(N-1,np.inf)
    for center,indexes,res in zip(centers_list,index_list,res_list):   #For each unit pair, select circle with best fit
        act_indexes=indexes[:-1]
        mask=best_res[act_indexes] > res
        def_centers[act_indexes]=center
        best_res[act_indexes]=res
    print(np.mean(best_res))
    return def_centers

if N>=3:
    rot_centers=rotation_fit()
    rot_angles=[get_angle(geometric_centers[i]-rot_centers[i],geometric_centers[i+1]-rot_centers[i]) for i in range(N-1)]
    if np.mean(rot_angles)<0.1:  # When there is no curvature, we fix the pitch axis
        rot_centers=[rot_centers[0] for i in range(len(rot_centers))]    
    
    

#---------------------------------------------------------------------------------------------------------------------------------

pitch_axis=[]
for i in range(N-1):
    vec_1=rot_centers[i]-geometric_centers[i]
    vec_1/=norm(vec_1)
    vec_2=rot_centers[i]-geometric_centers[i+1]
    vec_2/=norm(vec_2)
    pitch_axis.append((vec_1,vec_2))
    
twist_axis=[]
for i in range(N-1):
    twist_vect=geometric_centers[i+1]-geometric_centers[i]
    vec_1=orthogonalize(twist_vect,pitch_axis[i][0])
    vec_2=orthogonalize(twist_vect,pitch_axis[i][1])
    twist_axis.append((vec_1,vec_2))
    
rots=[Rotation.align_vectors([twist_axis[i][0],pitch_axis[i][0]],[twist_axis[i][1],pitch_axis[i][1]])[0] for i in range(N-1)]

units_rots=[]
pca=PCA()
units_comps=[]

for i in range(N):
    unit=units_coords[i]-geometric_centers[i]
    pca.fit(unit)
    comps=pca.components_
    if i > 0:   # We need to do 2 things
        def_comps=[]
        rot_comps=rots[i-1].apply(comps)
        ref_comps=units_comps[i-1]
        
        for j in [0,1]:
            coeffs=[abs(np.dot(vec,ref_comps[j])) for vec in rot_comps]
            index=np.argmax(coeffs)
            def_comps.append(comps[index])
            comps=np.delete(comps,index,0)
            rot_comps=np.delete(rot_comps,index,0)
        comps=def_comps
        rot_comps=rots[i-1].apply(comps)
        
        for j in [0,1]:
            if np.dot(rot_comps[j],ref_comps[j])<0:
                comps[j]=-comps[j]
        
                     
        comps=np.append(comps,[np.cross(comps[0],comps[1])],axis=0)   #Third component is fully determined by the first 2 
    else:
        comps[2]=np.cross(comps[0],comps[1])
    units_comps.append(comps) 
    
for i in range(N-1):
    dir_unit1=units_comps[i]
    dir_unit2=units_comps[i+1]
    dir_unit2=rots[i].apply(dir_unit2)
    units_rots.append(Rotation.align_vectors(dir_unit2,dir_unit1)[0])

pitchlist=[]
twistlist=[]
handednesslist=[]
for i in range(N-1):
    ref_pitch=units_rots[i].apply(twist_axis[i][0])
    pitchlist.append(dihedral_angle(twist_axis[i][0],ref_pitch,pitch_axis[i][0])[0])
    
    ref_twist=units_rots[i].apply(pitch_axis[i][0])
    res=dihedral_angle(pitch_axis[i][0],ref_twist,twist_axis[i][0])
    twistlist.append(res[0])
    handednesslist.append(res[1])

#---------------------------------------------------------------------------------------------------------------------------------
rot_angles.append(np.mean(rot_angles))
twistlist.append(np.mean(twistlist))
pitchlist.append(np.mean(pitchlist))
handednesslist.append(np.mean(handednesslist))
rot_angles.append(np.std(rot_angles))
twistlist.append(np.std(twistlist))
pitchlist.append(np.std(pitchlist))
handednesslist.append(np.std(handednesslist))

rot_angles.insert(0,0)
twistlist.insert(0,0)
handednesslist.insert(0,0)
pitchlist.insert(0,0)

starts=[unit[0] for unit in units_ids]
starts.append('mean')
starts.append('std deviation')
ends=[unit[1] for unit in units_ids]
ends.append('-')
ends.append('-')
pdbs=[filepath.split('\\')[-1][:-11] for i in range(N+2)]
chains=[args.chain for i in range(N+2)]

d={'pdb_id':pdbs,'chain':chains,'unit start':starts,'unit end':ends,'curvature':rot_angles,'twist':twistlist,'handedness':handednesslist,'pitch':pitchlist}
df=pd.DataFrame(data=d)



with pd.option_context('display.max_rows', None,'display.max_columns', None):
    print(df)


try:
    out_path=args.o
    
    df.to_csv(Path(out_path+'\out_'+pdbs[0]+'.csv'))
except Exception:
    pass
    
#-----------------------------------------------------------------------------------------------------------------------------------
#Pymol drawing
#Pymol drawing
if args.draw:
    pymol.finish_launching()
    cmd.load(filepath,format='pdb')
    cmd.hide('all')
    for i in range(N):   #Place pseudoatoms to draw distances and angles
        cmd.pseudoatom('geo_centers',pos=tuple(geometric_centers[i]))
        cmd.pseudoatom('first_component',pos=tuple(geometric_centers[i]+12*units_comps[i][0]))
        cmd.pseudoatom('second_component',pos=tuple(geometric_centers[i]+12*units_comps[i][1]))
        cmd.pseudoatom('third_component',pos=tuple(geometric_centers[i]+12*units_comps[i][2]))

    for i in range(len(rot_centers)):
        cmd.pseudoatom('rot_centers',pos=tuple(rot_centers[i]))


    for i in range(N-1):  #Draw rotation angles and protein geometry
        cmd.select('point1',selection='model geo_centers and name PS{}'.format(str(i+1)))
        cmd.select('point2',selection='model geo_centers and name PS{}'.format(str(i+2)))
        cmd.select('rot_center',selection='model rot_centers and name PS{}'.format(str(i+1)))
        cmd.angle('rot_angle',selection1='point1',selection2='rot_center',selection3='point2')
        cmd.distance('superaxis',selection1='point1',selection2='point2')


    for i in range(N):   #Draw reference system for each unit and principal components
        cmd.select('geo_center',selection='model geo_centers and name PS{}'.format(str(i+1)))
        cmd.select('ref_first',selection='model first_component and name PS{}'.format(str(i+1)))
        cmd.select('ref_second',selection='model second_component and name PS{}'.format(str(i+1)))
        cmd.select('ref_third',selection='model third_component and name PS{}'.format(str(i+1)))
        cmd.distance('pca_first',selection1='geo_center',selection2='ref_first')
        cmd.distance('pca_second',selection1='geo_center',selection2='ref_second')
        cmd.distance('pca_third',selection1='geo_center',selection2='ref_third')
        
    cmd.color('red','pca_first')
    cmd.color('orange','pca_second')
    cmd.color('white','pca_third')
    cmd.hide('labels')
    cmd.deselect()