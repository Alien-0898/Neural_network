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



idp=np.loadtxt('/home/raj/tapish/traj_13_boltzmann/NZPA_analysis/sn2.txt',skiprows=1)
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
        fs=0#1 if its front side 
        for i in range(len(inter)):
            if inter[i,9]<2.5 and inter[i,14]<2.5:
                theta=calculate_theta((cart[i,3]-cart[i,0]),(cart[i,5]-cart[i,0]))
                if theta<100:
                    fs=1
                else:
                    fs=0
        inver_list.append([set,fold,fs,theta])

np.savetxt('front_side_sn2.txt',np.array(inver_list),fmt='%d')                    
            

















 
