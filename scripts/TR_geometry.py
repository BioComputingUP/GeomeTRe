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
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.cealign import CEAligner
import Bio.SeqUtils
from Bio.SeqUtils import seq1
import tmtools

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
    n = n / norm(n)
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
    list_indexes=list(dict.fromkeys(list_indexes))
    L=[]
    for item in list_indexes:
        item=item.strip().split('_')
        item=(int(item[0]),int(item[1]))
        if len(L)>0:
            if item[0] > L[-1][1]:
                L.append(item)
        else:
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
        
def get_unit_rotation(coords,seqs,rotations):  # Align 2 units using CEalign, and return rotation
    coords_1=coords[0]
    coords_2=coords[1]
    coords_1=rotations[0].apply(coords_1)
    coords_2=rotations[1].apply(coords_2)
    alignment=tmtools.tm_align(coords_1,coords_2,seqs[0],seqs[1])
    return alignment
    
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
        
        
        center=res.x
        centers_list.append(pca.inverse_transform(center))
        index_list.append([*range(min_index,max_index)])
        score_list.append(np.std([norm(center-geo_center) for geo_center in pca_centers]))
    
    def_centers=np.empty((N-1,3))
    best_score=np.full(N-1,np.inf)
    for center,indexes,score in zip(centers_list,index_list,score_list):   #For each unit pair, select center corresponding to the widest crown
        act_indexes=indexes[:-1]
        score_to_confront=best_score[act_indexes]
        centers_to_confront=def_centers[act_indexes]
        mask=score_to_confront > score
        centers_to_confront[mask]=np.vstack([center for i in range(len(centers_to_confront[mask]))])
        score_to_confront[mask]=score
        best_score[act_indexes]=score_to_confront
        def_centers[act_indexes]=centers_to_confront
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
    
    rots=[]
    for i in range(N-1):
        try:
            rots.append((Rotation.from_matrix(np.array([twist_axis[i][0],pitch_axis[i][0],np.cross(twist_axis[i][0],pitch_axis[i][0])])),Rotation.from_matrix(np.array([twist_axis[i][1],pitch_axis[i][1],np.cross(twist_axis[i][1],pitch_axis[i][1])]))))
                        
        except:
            rots.append(Rotation.from_matrix(np.eye(3)))
    return pitch_axis,twist_axis,rots
    
def Pymol_drawing(filepath,geometric_centers,rot_centers,twist_axis,pitch_axis,rots,units_rots,unit_vector):
    N=len(geometric_centers)
    pymol.finish_launching()
    cmd.load(filepath)
    cmd.hide('all')
    for i in range(N):   #Place pseudoatoms to draw distances and angles
        cmd.pseudoatom('geo_centers',pos=tuple(geometric_centers[i]))
        #cmd.pseudoatom('ref_1',pos=tuple(geometric_centers[i]+6*unit_vector))
        #cmd.pseudoatom('ref_2',pos=tuple(geometric_centers[i]-6*unit_vector))
        #cmd.select('unit_1',selection='model ref_1 and name PS{}'.format(str(i+1)))
        #cmd.select('unit_2',selection='model ref_2 and name PS{}'.format(str(i+1)))
        #cmd.distance('unit_vector',selection1='unit_1',selection2='unit_2')
        if i < N-1:
            unit_vector=units_rots[i] @ (rots[i][0].apply(unit_vector))
            unit_vector=rots[i][1].apply(unit_vector,inverse=True)


    #for i in range(len(rot_centers)):
        #cmd.pseudoatom('rot_centers',pos=tuple(rot_centers[i]))


    for i in range(N-1):  #Draw rotation angles and protein geometry
        cmd.pseudoatom('twist_ref',pos=tuple(geometric_centers[i]+9*twist_axis[i][0]))
        cmd.pseudoatom('pitch_ref',pos=tuple(geometric_centers[i]+9*pitch_axis[i][0]))
        cmd.pseudoatom('yaw_ref',pos=tuple(geometric_centers[i]+9*np.cross(pitch_axis[i][0],twist_axis[i][0])))
        cmd.select('point1',selection='model geo_centers and name PS{}'.format(str(i+1)))
        #cmd.select('point2',selection='model geo_centers and name PS{}'.format(str(i+2)))
        #cmd.select('rot_center',selection='model rot_centers and name PS{}'.format(str(i+1)))
        cmd.select('twist_point',selection='model twist_ref and name PS{}'.format(str(i+1)))
        cmd.select('pitch_point',selection='model pitch_ref and name PS{}'.format(str(i+1)))
        cmd.select('yaw_point',selection='model yaw_ref and name PS{}'.format(str(i+1)))
        #cmd.angle('rot_angle',selection1='point1',selection2='rot_center',selection3='point2')
        #cmd.distance('superaxis',selection1='point1',selection2='point2')
        cmd.distance('twist_axis',selection1='point1',selection2='twist_point')
        cmd.distance('pitch_axis',selection1='point1',selection2='pitch_point')
        cmd.distance('yaw_axis',selection1='point1',selection2='yaw_point')
        
    cmd.pseudoatom('twist_ref',pos=tuple(geometric_centers[-1]+6*twist_axis[-1][1]))
    cmd.pseudoatom('pitch_ref',pos=tuple(geometric_centers[-1]+6*pitch_axis[-1][1]))
    cmd.pseudoatom('yaw_ref',pos=tuple(geometric_centers[-1]+6*np.cross(pitch_axis[-1][1],twist_axis[-1][1])))
    cmd.select('point1',selection='model geo_centers and name PS{}'.format(str(N)))
    #cmd.select('point2',selection='model geo_centers and name PS{}'.format(str(i+2)))
    #cmd.select('rot_center',selection='model rot_centers and name PS{}'.format(str(i+1)))
    cmd.select('twist_point',selection='model twist_ref and name PS{}'.format(str(N)))
    cmd.select('pitch_point',selection='model pitch_ref and name PS{}'.format(str(N)))
    cmd.select('yaw_point',selection='model yaw_ref and name PS{}'.format(str(N)))
    #cmd.angle('rot_angle',selection1='point1',selection2='rot_center',selection3='point2')
    #cmd.distance('superaxis',selection1='point1',selection2='point2')
    cmd.distance('twist_axis',selection1='point1',selection2='twist_point')
    cmd.distance('pitch_axis',selection1='point1',selection2='pitch_point')
    cmd.distance('yaw_axis',selection1='point1',selection2='yaw_point')
        
        
        #cmd.select('point1_1',selection='model ref_1 and name PS{}'.format(str(i+1)))
        #cmd.select('point1_2',selection='model ref_2 and name PS{}'.format(str(i+1)))
        #cmd.select('point2_1',selection='model ref_1 and name PS{}'.format(str(i+2)))
        #cmd.select('point2_2',selection='model ref_2 and name PS{}'.format(str(i+2)))
        #cmd.distance('dist_1',selection1='point1_1',selection2='point2_1')
        #cmd.distance('dist_2',selection1='point1_2',selection2='point2_2')
        
        
        
        
    #cmd.color('orange','unit_vector')
    #cmd.color('red','dist_1')
    #cmd.color('red','dist_2')
    cmd.color('red','twist_axis')
    cmd.color('orange','pitch_axis')
    cmd.color('white','yaw_axis')
    cmd.hide('labels')
    cmd.deselect()
    
    
#--------------------------------------------------------------------------------------------------------------------------------

def Repeats_geometry(filepath,chain,units_ids,ins_ids='',o_path='',draw=False,batch=False):
    units_ids=create_list(units_ids)
    mmcif_parser=MMCIFParser(QUIET=True)   #Parse .pdb file, define units, calculate center of each unit
    pdb_parser=PDBParser(QUIET=True) #call pdb parser
    file_type = filepath[-3:]
    if file_type == 'cif':
        structure=mmcif_parser.get_structure('structure',Path(filepath)) 
    elif file_type == 'pdb':
        structure=pdb_parser.get_structure('structure',Path(filepath))

    
    chain_s=structure[0][chain]
    if ins_ids:
        ins_ids=create_list(ins_ids)
    else:
        ins_ids=[]   
    if len(ins_ids)>0:   #If we have insertions, we make sure to remove them from the structure
        units=[]
        to_remove=[]
        for limits in units_ids:
            a,b=limits
            unit_ins=[ins for ins in ins_ids if a<=ins[0]<=b or a<=ins[1]<=b]
            unit=[]
            for residue in chain_s:
                res_id=residue.get_id()[1]
                res_in_ins=np.any([ins[0]<=res_id<=ins[1] for ins in unit_ins])
                if a<=res_id<=b and not res_in_ins and Polypeptide.is_aa(residue):
                    unit.append(residue)
            if len(unit)>0:
                units.append(unit)
            else:
                to_remove.append(limits)
        for ids in to_remove:
            units_ids.remove(ids)

    else:
        units=[]
        to_remove=[]
        for limits in units_ids:
            a,b=limits
            unit=[]
            for residue in chain_s:
                if a<=residue.get_id()[1]<=b:
                    if Polypeptide.is_aa(residue):
                        unit.append(residue)
            if len(unit)>0:
                units.append(unit)
            else:
                to_remove.append(limits)
        for ids in to_remove:
            units_ids.remove(ids)
    units_seqs=[]
    for unit in units:
        seq=seq1(''.join([res.get_resname() for res in unit]))
        units_seqs.append(seq)
                
    units_coords=[]   #For each unit, store coordinate of each CA atom
    for unit in units:
        ca_coords=[residue['CA'].get_coord() for residue in unit]
        units_coords.append(ca_coords)
        
    geometric_centers=[sum(coords)/len(coords) for coords in units_coords]  #Define geometric center for each unit
    N=len(geometric_centers)
    rot_centers=widest_circle_fit(units_coords,geometric_centers)
    rot_angles=[get_angle(geometric_centers[i]-rot_centers[i],geometric_centers[i+1]-rot_centers[i]) for i in range(N-1)]

    pitch_axis,twist_axis,rots=build_ref_axes(geometric_centers,rot_centers)


    units_rots=[]
    tmscores=[]
    for i in range(N-1):
        alignment=get_unit_rotation(units_coords[i:i+2],units_seqs[i:i+2],rots[i])
        units_rots.append(alignment.u)
        tmscores.append(alignment.tm_norm_chain1)
    
    pitchlist=[]
    twistlist=[]
    twist_handednesslist=[]
    pitch_handednesslist=[]
    yawlist=[]
    for i in range(N-1):   # Decompose rotation into pitch, twist and yaw
        rotation=units_rots[i]
        
        ref_twist=rotation @ pitch_axis[i][0]
        res=dihedral_angle(pitch_axis[i][0],ref_twist,np.array([1,0,0]))
        
        
        twist,pitch,yaw=Rotation.from_matrix(rotation).as_euler('xyz')
        twistlist.append(abs(twist))
        twist_handednesslist.append(np.sign(twist))
        pitchlist.append(abs(pitch))
        pitch_handednesslist.append(np.sign(pitch))
        yawlist.append(abs(yaw))  
        
    stats=[np.nanmean(rot_angles),np.nanstd(rot_angles),np.nanmean(twistlist),np.nanstd(twistlist),
           np.nanmean(twist_handednesslist),np.nanstd(twist_handednesslist),np.nanmean(pitchlist),
           np.nanstd(pitchlist),np.nanmean(pitch_handednesslist),np.nanstd(pitch_handednesslist),
           np.nanmean(tmscores),np.nanstd(tmscores),np.nanmean(yawlist),np.nanstd(yawlist)]
    
    
    # DataFrame output
    rot_angles.extend(stats[0:2])   
    twistlist.extend(stats[2:4])
    twist_handednesslist.extend(stats[4:6])
    pitchlist.extend(stats[6:8])
    pitch_handednesslist.extend(stats[8:10])
    tmscores.extend(stats[10:12])
    yawlist.extend(stats[12:14])

    rot_angles.insert(0,0)
    twistlist.insert(0,0)
    twist_handednesslist.insert(0,0)
    pitch_handednesslist.insert(0,0)
    pitchlist.insert(0,0)
    tmscores.insert(0,0)
    yawlist.insert(0,0)

    N=len(rot_angles)-2
    starts=[unit[0] for unit in units_ids]
    starts.append('mean')
    starts.append('std deviation')
    ends=[unit[1] for unit in units_ids]
    ends.append('-')
    ends.append('-')
    pdb=filepath.split('/')[-1][0:4]
    pdbs=[pdb for i in range(N+2)]
    chains=[chain for i in range(N+2)]

    d={'pdb_id':pdbs,'chain':chains,'unit start':starts,'unit end':ends,'curvature':rot_angles,'twist':twistlist,'twist_hand':twist_handednesslist,'pitch':pitchlist,'pitch_hand':pitch_handednesslist,'TM-score':tmscores,'yaw':yawlist}
    df=pd.DataFrame(data=d)
    if not batch:
        if o_path:
            df.to_csv(Path(o_path+'\out_'+pdb+'.csv'))
        else:
            with pd.option_context('display.max_rows', None,'display.max_columns', None):
                print(df)
    else:
        out_list=[pdb,chain,units_ids[0][0],units_ids[-1][1]]+stats
        print('\t'.join(str(x) for x in out_list))
            
    if draw:
        draw_pca=PCA()
        draw_pca.fit(units_coords[0])
        unit_vector=draw_pca.components_[0]
        Pymol_drawing(filepath,geometric_centers,rot_centers,twist_axis,pitch_axis,rots,units_rots,unit_vector)
        
    return df,stats
            
#---------------------------------------------------------------------------------------------------------------------------------

if __name__=='__main__':
    arg_parser=argparse.ArgumentParser()   #Parse command line
    arg_parser.add_argument('filepath',action='store',help='Path to input file')
    arg_parser.add_argument('chain',action='store',help='Chain')
    arg_parser.add_argument('unit_def',action='store',help='Unit limits, written as s1_e1,s2_e2,...')
    arg_parser.add_argument('-ins',action='store',default='',help='Starts and ends of insertions, formatted like the units')
    arg_parser.add_argument('-o',action='store',default='',help='Output file if desired, ex. Outputs\my_pdb.csv')
    arg_parser.add_argument('--draw',action='store_true',help='Use if Pymol drawing is desired')
    arg_parser.add_argument('--batch',action='store_true',help='Batch mode, used in run_TR_geometry.sh')
    args=arg_parser.parse_args()
    Repeats_geometry(args.filepath,args.chain,args.unit_def,args.ins,args.o,args.draw,args.batch)
