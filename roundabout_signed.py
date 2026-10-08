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


#import numpy as np

def signed_angle(a, b, n):
    a = np.array(a)
    b = np.array(b)
    n = np.array(n)
    
    return np.arctan2(
        np.dot(np.cross(a, b),n),
        np.dot(a, b)
    )
#===========================================
n = 6
masses=[14.007,1.008,1.008,35.453,1.008,15.999]
order=[4,1,2,3,0,5]
mech_list=[]


idp=np.loadtxt('tapish/traj_13_boltzmann/NZPA_analysis/indirect_proton.txt',skiprows=1)
for info in range(len(idp)):
#for info in range(1,2):
    set=int(idp[info,0])
    fold=int(idp[info,1])
    #set=1
    #fold=8234
    n=6
    if idp[info,0]>10:
        set1=set-10
        base=f'traj_13_boltzmann_2/set{set1}/{fold}'
    else:
        base=f'traj_13_boltzmann/set{set}/{fold}'

    #base='/home/raj/tapish/traj_13_boltzmann/set2/4174'

    order=[4,1,2,3,0,5]
    cart=read_xyz_blocks(f"{base}/movie.xyz")
    reordered_geometries = [geom[order] for geom in cart]
    inter_nuclie= [cart1d_to_distances1d(g) for g in reordered_geometries]
    inter=np.copy(inter_nuclie)
    
    indices = [i for i in range(len(inter)) if inter[i, 9] < 3.4 and inter[i, 14] < 4.2 and inter[i,8]<1.3 and inter[i,7]<1.3]
    start = indices[0]
    last = indices[-1]
    angle=[]
    cart=np.array(cart)
            
    count=False
    mech=[0,8]
    first=0
    for i in range(start,last):
                                                                                                              
        if inter[i,7]<1.2 and inter[i,8]<1.2 and inter[i,9]<2.2:
            if i==first+1 and first>0:
                L_vec_in=np.cross(cart_first,((14.007*cart[i,0])+(1.008*cart[i,1])+(1.008*cart[i,2]))/(14.007+1.008+1.008) - cart[i,3])
            if i>first and first>0:
                tan_angle=signed_angle(cart_first, ((14.007*cart[i,0])+(1.008*cart[i,1])+(1.008*cart[i,2]))/(14.007+1.008+1.008) - cart[i,3], L_vec_in)
                theta_deg = np.degrees(tan_angle)
                theta = (theta_deg + 360) % 360 # theta_0_360
                                                                             
            if count:
                if (inter[i,14]<2.5 and np.max(inter[first:i,14])>4.2) or (inter[i,14]<3 and np.max(inter[first:i, 14])>4.2 and inter[i-1,14]-inter[i,14]>0 and inter[i-20,14]-inter[i,14]>0 and \
                    inter[i,14]-inter[i+1,14]<0 and inter[i,14]-inter[i+20,14]<0):#END OF ONE CYCLE- CHECK IF ITS GOES OUT INTERACTION REGION BEFORE RETURNING AND IF ITS A TURNING POINT OR NEAR TURING POINT                                                  
                    angle_list=[bin_0_45, bin_45_90, bin_90_135, bin_135_180, bin_180_225, bin_225_270, bin_270_315, bin_315_360]
                    zero_bins=sum(v == 0 for v in angle_list)
                    #-----ANALYSIS------------
                    if zero_bins<3:    
                        mech=[1,zero_bins]
                    elif np.max(inter[first:i, 14])>4.2:
                        if mech[0]==0:
                            mech=[2,zero_bins]
                    else:
                        if mech[0]==0:
                            mech=[0,zero_bins] 
                    
                    first=i
                    theta=0
                    bin_0_45=0
                    bin_45_90=0
                    bin_90_135=0
                    bin_135_180=0
                    bin_180_225=0
                    bin_225_270=0
                    bin_270_315=0
                    bin_315_360=0
                    if inter[i,14]<3 and inter[i-1,14]-inter[i,14]>0 and inter[i-20,14]-inter[i,14]>0 and \
                        inter[i,14]-inter[i+1,14]<0 and inter[i,14]-inter[i+20,14]<0:#IF ITS A TURNING POINT AND INSIDE RANGE   
                        count=True
                        first=i
                        theta=0
                        cart_first=((14.007*cart[first,0])+(1.008*cart[first,1])+(1.008*cart[first,2]))/(14.007+1.008+1.008) - cart[first,3]
                        bin_0_45=0
                        bin_45_90=0
                        bin_90_135=0
                        bin_135_180=0
                        bin_180_225=0
                        bin_225_270=0
                        bin_270_315=0
                        bin_315_360=0
                            #print('--',i)
                    else:
                        count=False
            if inter[i,14]<3 and inter[i-1,14]-inter[i,14]>0 and inter[i-20,14]-inter[i,14]>0 and\
                inter[i,14]-inter[i+1,14]<0 and inter[i,14]-inter[i+20,14]<0 and count==False:#TO FIND TURNING POINT
                count=True
                first=i
                theta=0
                cart_first=((14.007*cart[first,0])+(1.008*cart[first,1])+(1.008*cart[first,2]))/(14.007+1.008+1.008) - cart[first,3]
                bin_0_45=0
                bin_45_90=0
                bin_90_135=0
                bin_135_180=0
                bin_180_225=0
                bin_225_270=0
                bin_270_315=0
                bin_315_360=0
            if count:
                if 2<theta<45:
                    bin_0_45+=1
                elif 45<=theta<90:
                    bin_45_90+=1
                elif 90<=theta<135:
                    bin_90_135+=1
                elif 135<=theta<180:
                    bin_135_180+=1
                elif 180<=theta<225:
                    bin_180_225+=1
                elif 225<=theta<270:
                    bin_225_270+=1
                elif 270<=theta<315:
                    bin_270_315+=1
                elif 315<=theta<360:
                    bin_315_360+=1
        if first>0:         
            angle.append([i,theta])
    np.savetxt(f'{base}/angle',np.array(angle),fmt='%f')

    mech_list.append([set,fold,mech[0],mech[1],np.max(inter[start:last, 14]),last-start])
    
#h-RA:3
#RA:1
#roam:2
#nan:0
np.savetxt('roundabout_proton_new.txt',np.array(mech_list),fmt="%d %d %d %d %f %d")    
            
