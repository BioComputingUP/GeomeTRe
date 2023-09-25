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
from Bio.PDB.MMCIFParser import MMCIFParser

import scipy
from scipy.spatial.transform import Rotation
import sklearn
from sklearn.decomposition import PCA

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
        item=item.strip().split('-')
        item=[int(item[0]),int(item[1])]
        L.append(item)
    list_indexes=L
    return L

def find_rotation(A,B):
    n=len(A)
    m=len(B)
    j=abs(n-m)
    best_fit=np.inf
    flag=n >= m
    rot=None
    for i in range(j+1):
        if flag:
            clipped1=A[i:i+m]
            clipped2=B
        else:
            clipped1=A
            clipped2=B[i:i+n]
        new_fit=Rotation.align_vectors(clipped2,clipped1)
        if new_fit[1] < best_fit:
            rot=new_fit[0]
    return rot
        


#---------------------------------------------------------------------------------------------------------------------------------

arg_parser=argparse.ArgumentParser()   #Parse command line
arg_parser.add_argument('filepath',action='store',help='Path to input file')
arg_parser.add_argument('chain',action='store',help='Chain')
arg_parser.add_argument('unit_def',action='store',help='Unit limits, written as s1-e1,s2-e2,...')
arg_parser.add_argument('-o',action='store',help='Output file if desired, ex. Outputs\my_pdb.csv')
args=arg_parser.parse_args()
filepath=args.filepath
list_indexes=args.unit_def

#---------------------------------------------------------------------------------------------------------------------------------

parser=MMCIFParser(QUIET=True)   #Parse .cif file, define units, calculate center of each unit
structure=parser.get_structure('structure',Path(filepath))

chain=structure[0][args.chain]
units_ids=create_list(list_indexes)

units=[]     #Add residues to each unit
for limits in units_ids:
    a,b=limits
    unit=[]
    for residue in chain:
        if a<=residue.get_id()[1]<=b:
            if Polypeptide.is_aa(residue):
                unit.append(residue)
    units.append(unit)

units_coords=[]   #For each unit, store coordinate of each CA atom
for unit in units:
    units_coords.append([residue['CA'].get_coord() for residue in unit])

geometric_centers=[sum(coords)/len(coords) for coords in units_coords]  #Define geometric center for each unit
N=len(geometric_centers)
#----------------------------------------------------------------------------------------------------------------------------------

#CURVATURE
def rotation_fit(C,window=6):
    data=np.asarray(C)
    N=len(data)
    index_list=[]
    rotation_list=[]
    centers_list=[]
    for i in range(N-window+1):
        indexes=[*range(i,min(i+window,N))]
        data_to_fit=data[min(i,N-window):min(i+window,N)]
        pca=PCA()
        pca.fit(data_to_fit)
        data_transformed=pca.transform(data_to_fit)

        #Delete third axis, then fit ellipse
        circle_data=np.delete(data_transformed,2,1)
        c=circle_fit(circle_data).tolist()
        c.append(0)
        centers_list.append(pca.inverse_transform(c))
        index_list.append(indexes)

    av_centers=[]
    m=len(centers_list)
    for i in range(N):
        indexes=[j for j in range(m) if i in index_list[j]]
        to_av=[centers_list[j] for j in indexes]
        av_centers.append(np.mean(to_av,axis=0))
    return av_centers

rot_centers=rotation_fit(geometric_centers,6)
curv_centers=[(rot_centers[i]+rot_centers[i+1])/2 for i in range(N-1)]
rot_angles=[get_angle(geometric_centers[i]-curv_centers[i],geometric_centers[i+1]-curv_centers[i]) for i in range(N-1)]
if np.mean(rot_angles)<0.1:  # When there is no curvature, we fix the pitch axis
    rot_centers=[rot_centers[0] for i in range(len(rot_centers))]

#---------------------------------------------------------------------------------------------------------------------------------

new_pitch_axis=[]
for i in range(N):
    new_pitch_axis.append(rot_centers[i]-geometric_centers[i])
    new_pitch_axis[-1] /= norm(new_pitch_axis[-1])

new_twist_axis=[geometric_centers[i+2]-geometric_centers[i] for i in range(N-2)]
for i in range(N-2):
    new_twist_axis[i] /= norm(new_twist_axis[i])

rots=[Rotation.align_vectors([new_twist_axis[i],new_pitch_axis[1:-1][i]],[new_twist_axis[i+1],new_pitch_axis[1:-1][i+1]])[0] for i in range(N-3)]

new_twist_axis.insert(0,rots[0].apply(new_twist_axis[0]))
new_twist_axis.append(rots[-1].apply(new_twist_axis[-1],inverse=True))
for i in range(N):
    new_twist_axis[i]=orthogonalize(new_twist_axis[i],new_pitch_axis[i])
    
rots.insert(0,Rotation.align_vectors([new_twist_axis[0],new_pitch_axis[0]],[new_twist_axis[1],new_pitch_axis[1]])[0])
rots.append(Rotation.align_vectors([new_twist_axis[-2],new_pitch_axis[-2]],[new_twist_axis[-1],new_pitch_axis[-1]])[0])

units_rots=[]
for i in range(N-1):
    unit1=units_coords[i]-geometric_centers[i]
    unit2=units_coords[i+1]-geometric_centers[i+1]
    unit2=rots[i].apply(unit2)
    units_rots.append(find_rotation(unit1,unit2))

pitchlist=[]
twistlist=[]
handednesslist=[]
for i in range(N-1):
    ref_pitch=units_rots[i].apply(new_twist_axis[i])
    pitchlist.append(dihedral_angle(new_twist_axis[i],ref_pitch,new_pitch_axis[i])[0])
    
    ref_twist=units_rots[i].apply(new_pitch_axis[i])
    res=dihedral_angle(new_pitch_axis[i],ref_twist,new_twist_axis[i])
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
pdbs=[filepath[-8:-4] for i in range(N+2)]
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
pymol.finish_launching()
cmd.load(filepath,format='cif')
cmd.hide('all')
start_vector=units_coords[0][-1]-units_coords[0][0]
for i in range(N):   #Place pseudoatoms to draw distances and angles
    cmd.pseudoatom('geo_centers',pos=tuple(geometric_centers[i]))
    cmd.pseudoatom('twist_ref',pos=tuple(geometric_centers[i]+6*new_twist_axis[i]))
    cmd.pseudoatom('pitch_ref',pos=tuple(geometric_centers[i]+6*new_pitch_axis[i]))
    cmd.pseudoatom('start_ref',pos=tuple(geometric_centers[i]-0.5*start_vector))
    cmd.pseudoatom('end_ref',pos=tuple(geometric_centers[i]+0.5*start_vector))
    if i < N-1:
        start_vector=units_rots[i].apply(rots[i].apply(start_vector,inverse=True))
    
for i in range(len(curv_centers)):
    cmd.pseudoatom('rot_centers',pos=tuple(curv_centers[i]))
    

for i in range(N-1):  #Draw rotation angles and protein geometry
    cmd.select('point1',selection='model geo_centers and name PS{}'.format(str(i+1)))
    cmd.select('point2',selection='model geo_centers and name PS{}'.format(str(i+2)))
    cmd.select('rot_center',selection='model rot_centers and name PS{}'.format(str(i+1)))
    cmd.angle('rot_angle',selection1='point1',selection2='rot_center',selection3='point2')
    cmd.distance('superaxis',selection1='point1',selection2='point2')
        
        
for i in range(N):   #Draw reference system for each unit
    cmd.select('geo_center',selection='model geo_centers and name PS{}'.format(str(i+1)))
    cmd.select('twist',selection='model twist_ref and name PS{}'.format(str(i+1)))
    cmd.select('pitch',selection='model pitch_ref and name PS{}'.format(str(i+1)))
    cmd.distance('twist_axis',selection1='geo_center',selection2='twist')
    cmd.distance('pitch_axis',selection1='geo_center',selection2='pitch')
    
for i in range(N-1): #Draw vectors representing unit orientations
    cmd.select('start_point_1',selection='model start_ref and name PS{}'.format(str(i+1)))
    cmd.select('end_point_1',selection='model end_ref and name PS{}'.format(str(i+1)))
    cmd.select('start_point_2',selection='model start_ref and name PS{}'.format(str(i+2)))
    cmd.select('end_point_2',selection='model end_ref and name PS{}'.format(str(i+2)))
    cmd.distance('unit_vectors',selection1='start_point_1',selection2='end_point_1')
    cmd.distance('unit_vectors',selection1='start_point_2',selection2='end_point_2')
    cmd.distance('units_axis',selection1='start_point_1',selection2='start_point_2')
    cmd.distance('units_axis',selection1='end_point_1',selection2='end_point_2')
        
cmd.color('green','twist_axis')
cmd.color('blue','pitch_axis')
cmd.color('orange','unit_vectors')
cmd.color('red','units_axis')
cmd.hide('labels')
cmd.deselect()