import argparse
from pathlib import Path
import copy

import numpy as np
from numpy.linalg import norm
import pandas as pd
import pymol
from pymol import cmd

import Bio
from Bio import PDB
from Bio.PDB import Polypeptide
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.cealign import CEAligner

import scipy
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
import sklearn
from sklearn.decomposition import PCA
from skimage.measure import CircleModel,EllipseModel

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

def widest_circle(c,data):   # Widest circular crown that's within units
    nearest=np.NINF
    farthest=np.inf
    for unit in data:
        distances=[norm(ca-c) for ca in unit]
        is_farthest=max(distances)
        if is_farthest < farthest:
            farthest=is_farthest
        is_nearest=min(distances)
        if is_nearest > nearest:
            nearest=is_nearest
    return nearest-farthest  # We use scipy.minimize, so we return the opposite of the width
        
def get_unit_rotation(unit1,unit2,rot):  # Align 2 units using CEalign, and return rotation
    chain1=Chain('A1')
    chain2=Chain('A2')
    unit1_toalign=copy.deepcopy(unit1)
    unit1_copy=copy.deepcopy(unit1)
    unit2_copy=copy.deepcopy(unit2)
    
    for residue in unit1_toalign:
        chain1.add(residue)
    for residue in unit2_copy:
        chain2.add(residue)
        
    chain2.transform(rot.inv().as_matrix(),np.zeros(3))  #Align reference axes of the units
    
    model1=Model('U1')
    model2=Model('U2')
    model1.add(chain1)
    model2.add(chain2)
    win_size=min([8,len(chain1)//2,len(chain2)//2])
    try:
        aligner=CEAligner(window_size=win_size,max_gap=20)
        aligner.set_reference(model2)
        aligner.align(model1,transform=True)
        coord1=np.asarray([residue['CA'].get_coord() for residue in unit1_copy])
        coord2=np.asarray([residue['CA'].get_coord() for residue in model1['A1']])

        coord1 -= np.mean(coord1,axis=0)
        coord2 -= np.mean(coord2,axis=0)
        return Rotation.align_vectors(coord2,coord1)[0]
    except:
        return None
    
def widest_circle_fit(units,centers,window=6):   # Alternative method for curvature
    N=len(units)
    index_list=[]
    centers_list=[]
    score_list=[]
    
    for i in range(max(N-window+1,1)):
        min_index=i
        max_index=min(i+window,N)
        data_to_fit=units[min_index:max_index]
        pca_centers=centers[min_index:max_index]
        
        pca=PCA(n_components=2)  # Find plane of rotation of units, and project them onto it
        pca.fit(pca_centers)
        pca_centers=pca.transform(pca_centers)
        data_transformed=[]
        for unit in data_to_fit:
            data_transformed.append(pca.transform(unit))
        circle=CircleModel()
        circle.estimate(pca_centers)
        res=minimize(widest_circle,circle.params[0:2],args=(data_transformed))  # Find widest crown in the 2d plane
        
        centers_list.append(pca.inverse_transform(res.x))
        index_list.append([*range(min_index,max_index)])
        score_list.append(res.fun)
    
    def_centers=np.empty((N-1,3))
    best_score=np.full(N-1,np.NINF)
    for center,indexes,score in zip(centers_list,index_list,score_list):   #For each unit pair, select center corresponding to the widest crown
        act_indexes=indexes[:-1]
        mask=best_score[act_indexes] < score
        def_centers[act_indexes]=center
        best_score[act_indexes]=score
    return def_centers

def build_ref_axes(geometric_centers,rot_centers):
    N=len(geometric_centers)
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
    return pitch_axis,twist_axis,rots
    
def Pymol_drawing(filepath,geometric_centers,rot_centers,twist_axis,rots,units_rots,unit_vector):
    N=len(geometric_centers)
    pymol.finish_launching()
    cmd.load(filepath,format='pdb')
    cmd.hide('all')
    for i in range(N):   #Place pseudoatoms to draw distances and angles
        cmd.pseudoatom('geo_centers',pos=tuple(geometric_centers[i]))
        cmd.pseudoatom('ref_1',pos=tuple(geometric_centers[i]+6*unit_vector))
        cmd.pseudoatom('ref_2',pos=tuple(geometric_centers[i]-6*unit_vector))
        cmd.select('unit_1',selection='model ref_1 and name PS{}'.format(str(i+1)))
        cmd.select('unit_2',selection='model ref_2 and name PS{}'.format(str(i+1)))
        cmd.distance('unit_vector',selection1='unit_1',selection2='unit_2')
        if i < N-1:
            unit_vector=units_rots[i].apply(rots[i].apply(unit_vector,inverse=True))


    for i in range(len(rot_centers)):
        cmd.pseudoatom('rot_centers',pos=tuple(rot_centers[i]))


    for i in range(N-1):  #Draw rotation angles and protein geometry
        cmd.pseudoatom('twist_ref',pos=tuple(geometric_centers[i]+6*twist_axis[i][0]))
        cmd.select('point1',selection='model geo_centers and name PS{}'.format(str(i+1)))
        cmd.select('point2',selection='model geo_centers and name PS{}'.format(str(i+2)))
        cmd.select('rot_center',selection='model rot_centers and name PS{}'.format(str(i+1)))
        cmd.select('twist_point',selection='model twist_ref and name PS{}'.format(str(i+1)))
        cmd.angle('rot_angle',selection1='point1',selection2='rot_center',selection3='point2')
        cmd.distance('superaxis',selection1='point1',selection2='point2')
        cmd.distance('twist_axis',selection1='point1',selection2='twist_point')
        
        cmd.select('point1_1',selection='model ref_1 and name PS{}'.format(str(i+1)))
        cmd.select('point1_2',selection='model ref_2 and name PS{}'.format(str(i+1)))
        cmd.select('point2_1',selection='model ref_1 and name PS{}'.format(str(i+2)))
        cmd.select('point2_2',selection='model ref_2 and name PS{}'.format(str(i+2)))
        cmd.distance('dist_1',selection1='point1_1',selection2='point2_1')
        cmd.distance('dist_2',selection1='point1_2',selection2='point2_2')
        
        
        
        
    cmd.color('orange','unit_vector')
    cmd.color('red','dist_1')
    cmd.color('red','dist_2')
    cmd.color('green','twist_axis')
    cmd.hide('labels')
    cmd.deselect()
    


def compute_geometry(filepath,chain,units_ids,ins_ids=[],draw=False):
    parser=PDBParser(QUIET=True)   #Parse .pdb file, define units, calculate center of each unit
    structure=parser.get_structure('structure',Path(filepath))
    chain_s=structure[0][chain]
    if len(ins_ids)>0:   #If we have insertions, we make sure to remove them from the structure
        units=[]
        for limits in units_ids:
            a,b=limits
            unit_ins=[ins for ins in ins_ids if a<=ins[0]<=b or a<=ins[1]<=b]
            unit=[]
            for residue in chain_s:
                res_id=residue.get_id()[1]
                for atom in residue:
                    if atom.is_disordered() and atom.get_altloc() != "A":   # Remove disordered / duplicated atoms from residue
                        residue.detach_child(atom.get_id())   
                res_in_ins=np.any([ins[0]<=res_id<=ins[1] for ins in unit_ins])
                if a<=res_id<=b and not res_in_ins and Polypeptide.is_aa(residue):
                    unit.append(residue)
            units.append(unit)

    else:
        units=[]    
        for limits in units_ids:
            a,b=limits
            unit=[]
            for residue in chain_s:
                for atom in residue:
                    if atom.is_disordered() and atom.get_altloc() != "A":   # Remove disordered / duplicated atoms from residue
                        residue.detach_child(atom.get_id())
                if a<=residue.get_id()[1]<=b:
                    if Polypeptide.is_aa(residue):
                        unit.append(residue)
            units.append(unit)
    units_coords=[]   #For each unit, store coordinate of each CA atom
    region=[]
    for unit in units:
        ca_coords=[residue['CA'].get_coord() for residue in unit]
        if len(ca_coords) > 0:
            units_coords.append(ca_coords)
            region.extend(ca_coords)
        else:
            del(units_ids[i])
            del(units[i])

    geometric_centers=[sum(coords)/len(coords) for coords in units_coords]  #Define geometric center for each unit
    N=len(geometric_centers)
    assert N>=3,'At least 3 units needed'
    rot_centers=widest_circle_fit(units_coords,geometric_centers)
    rot_angles=[get_angle(geometric_centers[i]-rot_centers[i],geometric_centers[i+1]-rot_centers[i]) for i in range(N-1)]

    pitch_axis,twist_axis,rots=build_ref_axes(geometric_centers,rot_centers)


    units_rots=[]
    for i in range(N-1):
        res=get_unit_rotation(units[i],units[i+1],rots[i])
        if res is not None:
            units_rots.append(res)
        else:
            units_rots.append('FAIL')
            print('CEAlign failure for units {}-{}'.format(i,i+1))
    
    pitchlist=[]
    twistlist=[]
    handednesslist=[]
    for i in range(N-1):   # Decompose rotation into pitch and twist
        rotation=units_rots[i]
        if rotation != 'FAIL':
            ref_pitch=units_rots[i].apply(twist_axis[i][0])
            pitchlist.append(dihedral_angle(twist_axis[i][0],ref_pitch,pitch_axis[i][0])[0])

            ref_twist=units_rots[i].apply(pitch_axis[i][0])
            res=dihedral_angle(pitch_axis[i][0],ref_twist,twist_axis[i][0])
            twistlist.append(res[0])
            handednesslist.append(res[1])
        else:
            pitchlist.append(np.NAN)
            twistlist.append(np.NAN)
            handednesslist.append(np.NAN)
                  
    if draw:
        draw_pca=PCA()
        draw_pca.fit(units_coords[0])
        unit_vector=draw_pca.components_[0]
        Pymol_drawing(filepath,geometric_centers,rot_centers,twist_axis,rots,units_rots,unit_vector)
        
    return rot_angles,twistlist,pitchlist,handednesslist

def Repeats_geometry(filepath,chain,units_ids,ins_ids='',o_path=None,draw=False,):
    units_ids=create_list(units_ids)
    rot_angles,twistlist,pitchlist,handednesslist=compute_geometry(filepath,args.chain,units_ids,ins_ids,draw)
    
    # DataFrame output
    rot_angles.extend([np.nanmean(rot_angles),np.nanstd(rot_angles)])   
    twistlist.extend([np.nanmean(twistlist),np.nanstd(twistlist)])
    pitchlist.extend([np.nanmean(pitchlist),np.nanstd(pitchlist)])
    handednesslist.extend([np.nanmean(handednesslist),np.nanstd(handednesslist)])
    #rmsds.extend([np.nanmean(rmsds),np.nanstd(rmsds)])

    rot_angles.insert(0,0)
    twistlist.insert(0,0)
    handednesslist.insert(0,0)
    pitchlist.insert(0,0)
    #rmsds.insert(0,0)

    N=len(rot_angles)-2
    starts=[unit[0] for unit in units_ids]
    starts.append('mean')
    starts.append('std deviation')
    ends=[unit[1] for unit in units_ids]
    ends.append('-')
    ends.append('-')
    pdb=filepath.split('\\')[-1]
    pdbs=[pdb for i in range(N+2)]
    chains=[chain for i in range(N+2)]
    d={'pdb_id':pdbs,'chain':chains,'unit start':starts,'unit end':ends,'curvature':rot_angles,'twist':twistlist,'handedness':handednesslist,'pitch':pitchlist}
    df=pd.DataFrame(data=d)
    if o_path:
        df.to_csv(Path(o_path+'\out_'+pdb+'.csv'))
    else:
        with pd.option_context('display.max_rows', None,'display.max_columns', None):
            print(df)
            
#---------------------------------------------------------------------------------------------------------------------------------

if __name__=='__main__':
    arg_parser=argparse.ArgumentParser()   #Parse command line
    arg_parser.add_argument('filepath',action='store',help='Path to input file')
    arg_parser.add_argument('chain',action='store',help='Chain')
    arg_parser.add_argument('unit_def',action='store',help='Unit limits, written as s1_e1,s2_e2,...')
    arg_parser.add_argument('-ins',action='store',help='Starts and ends of insertions, formatted like the units')
    arg_parser.add_argument('-o',action='store',help='Output file if desired, ex. Outputs\my_pdb.csv')
    arg_parser.add_argument('--draw',action='store_true',help='Use if Pymol drawing is desired')
    args=arg_parser.parse_args()
    ins=args.ins
    if ins is None:
        ins=''
    o_path=args.o
    Repeats_geometry(args.filepath,args.chain,args.unit_def,ins,o_path,args.draw)