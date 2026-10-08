import numpy as np
from itertools import combinations

def cart1d_to_distances1d(vec):
    vec = vec.reshape(-1,3)
    n = len(vec)
    distance_matrix = np.zeros((n,n))
    for i,j in combinations(range(len(vec)),2):
        R = np.linalg.norm(vec[i]-vec[j])
        distance_matrix[j,i] = R
    distance_vector = distance_matrix[np.tril_indices(len(distance_matrix),-1)]
    return distance_vector


def read_xyz_blocks(fname):
    """
    Returns
    -------
    geoms : list of ndarray (6,3)  Cartesian coordinates
    elems : list of list[str]      element symbols per geometry
    """
    geoms = []
    with open(fname) as fh:
        lines = [ln.rstrip() for ln in fh]
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "6":                  # start of a geometry block
            i += 2                       # skip “6” line + blank/comment line
            coords = []
            for _ in range(6):
                parts = lines[i].split()
                if len(parts) < 4:
                    raise ValueError(f"Malformed atom line @ line {i+1}")
                coords.append([float(p) for p in parts[1:4]])
                i += 1
            geoms.append(np.asarray(coords, dtype=float))
        else:
            i += 1
                                  # skip other lines
    return geoms

def calculate_com_velocity(data, masses):
    total_mass = sum(masses)
    velocity_sum = np.zeros(3)
    for velocity, mass in zip(data, masses):
        velocity_sum += velocity
    return velocity_sum / total_mass


def calculate_cosine_theta(v_rel_initial, v_rel_final):
    cos_theta = np.dot(v_rel_initial, v_rel_final) / (np.linalg.norm(v_rel_initial) * np.linalg.norm(v_rel_final))
    return cos_theta

def calculate_theta(v_rel_initial, v_rel_final):
    cos_theta = calculate_cosine_theta(v_rel_initial, v_rel_final)
    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # Ensure the value is within valid range
    return np.degrees(theta)  # Convert to degrees
#===========================================
n = 6
masses=[14.007,1.008,1.008,35.453,1.008,15.999]
order=[4,1,2,3,0,5]
inver_list=[]



idp=np.loadtxt('traj_13_boltzmann/NZPA_analysis/direct_sn2.txt',skiprows=1)
#for set in range(1,16):
    #for fold in range(1,10001):
for info in range(len(idp)):
        set=int(idp[info,0])
        fold=int(idp[info,1])
        n=6
        if set>10:
            set1=set-10
            base=f'/home/raj/tapish/traj_13_boltzmann_2/set{set1}/{fold}'
        else:
            base=f'/home/raj/tapish/traj_13_boltzmann/set{set}/{fold}'
        #base='/home/raj/tapish/traj_13_boltzmann/set1/3418'
        order=[4,1,2,3,0,5]
        cart=read_xyz_blocks(f"{base}/movie.xyz")
        reordered_geometries = [geom[order] for geom in cart]
        inter_nuclie= [cart1d_to_distances1d(g) for g in reordered_geometries]
        inter=np.copy(inter_nuclie)
        cart=np.array(cart)        
        #print(cart.shape) 
        inver=np.loadtxt(f"{base}/sign")#idx,(10 or 15),value,(0 or 1 or -1)
        initial_value=inver[0,2]
        bond_breaking_parm=False #for checking of NHCl is formed or not as in that region sign is -1
        current_sign=0
        info_list=[0,0,0,0,0,0,0]#invertion due to bond break,invertion due to OH interaction,invertion outside before sn2,inversion by cl after sn2,invertion outside after sn2,retention(near TS),overall retention (1 if written thing present else 0) 
        for i in range(len(inter)):
            if inver[i,3]==-1:
                bond_breaking_parm=True
            else:
                if inver[i,1]==10:
                    if inver[i,3]==current_sign:    
                        if bond_breaking_parm==True:
                            bond_breaking_parm=False
                    else:
                        if ((inter[i, 6]<1.5 and inter[i,10]<1.5) or (inter[i, 7]<1.5 and inter[i,11]<1.5) or (inter[i, 8]<1.5 and inter[i,12]<1.5)) and bond_breaking_parm==True:
                         #inversion after proton transfer and return in inverted position 
                            info_list[0]=1
                        elif ((inter[i, 6]<1.5 and inter[i,10]<1.5) or (inter[i, 7]<1.5 and inter[i,11]<1.5) or (inter[i, 8]<1.5 and inter[i,12]<1.5)) and bond_breaking_parm==False:
                         #inversion have because of OH in interaction region
                            info_list[1]=1
                        else:
                            #inversion happen outside interaction region
                            info_list[2]=1
                        current_sign=inver[i,3]
                if inver[i,1]==15:
                    if inver[i,3]!=current_sign:
                        if inter[i,9]<2.5 and inter[i,14]<2.5:
                                #inversion after sn2
                            if (initial_value/inver[i,2])>0:
                                #retention
                                info_list[5]=1
                            #else: inversion
                        elif ((inter[i, 6]<1.5 and inter[i,3]<2.3) or (inter[i, 7]<1.5 and inter[i,4]<2.3) or (inter[i, 8]<1.5 and inter[i,5]<2.3)):
                            info_list[3]=1
                        else:
                                #inversion happen outside interaction region
                            info_list[4]=1
                        current_sign=inver[i,3]
            if i==(len(inter)-2):
                if (initial_value*inver[i,2])>0:
                    info_list[6]=1
        inver_list.append([set,fold,info_list[0],info_list[1],info_list[2],info_list[3],info_list[4],info_list[5],info_list[6]])

np.savetxt('double_inver_dirsn2.txt',np.array(inver_list),fmt='%d',header='bond_break inversion_OH inversion_before_SN2 inversion_by_Cl_after_SN2 inversion_after_SN2 retention_near_TS overall_retention')                    
            

















 
