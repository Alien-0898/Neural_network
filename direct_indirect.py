import os
import numpy as np
from itertools import combinations
from scipy.ndimage import gaussian_filter1d


path=np.loadtxt('all_product_analysis_new.txt')
# ------------ READ XYZ ----------------
def read_xyz_blocks(fname):
    geoms = []
    if not os.path.isfile(fname):
        return None   # File not found (Nan later)

    with open(fname) as fh:
        lines = [ln.rstrip() for ln in fh]

    i = 0
    while i < len(lines):
        if lines[i].strip() == "6":
            i += 2
            coords = []
            for _ in range(6):
                parts = lines[i].split()
                coords.append([float(p) for p in parts[1:4]])
                i += 1
            geoms.append(np.asarray(coords))
        else:
            i += 1

    return geoms if geoms else None


# -------- DISTANCE VECTOR -------------
def cart_to_distances(vec):
    vec = vec.reshape(-1,3)
    n = len(vec)
    distance_matrix = np.zeros((n,n))
    for i,j in combinations(range(len(vec)),2):
        R = np.linalg.norm(vec[i]-vec[j])
        distance_matrix[j,i] = R
    distance_vector = distance_matrix[np.tril_indices(len(distance_matrix),-1)]
    return distance_vector

# ===============================
# MAIN ANALYSIS FOR ONE TRAJECTORY
# ===============================
def analyze_traj(xyzfilei,path):
    # Load trajectory geometries
    cart_geoms = read_xyz_blocks(xyzfile)
    if cart_geoms is None:
        return np.nan, np.nan, np.nan

    # Reorder atoms
    new_order = [4, 1, 2, 3, 0, 5]
    reordered = [g[new_order] for g in cart_geoms]

    # Convert to internuclear distances
    inter = np.array([cart_to_distances(g) for g in reordered])

    # If geometry count < 3 → cannot detect turning points
    if len(inter) < 3:
        return np.nan, np.nan, np.nan

    # ----------------------------------------
    #TIME + TURNING
    # ----------------------------------------
    indices = [i for i in range(len(inter))
                   if inter[i, 9] < 3.4 and inter[i, 14] < 4.2]#int1 N-O:4.16 and int5 N-Cl:3.31
    inter_15=gaussian_filter1d(inter[:,14], sigma=14,truncate=1)#N-O dist
    inter_10=gaussian_filter1d(inter[:,9], sigma=14,truncate=1)#N-Cl dist
    if len(indices) == 0:
        time = np.nan
        turning = np.nan
        turning_1=np.nan
    else:
        first = indices[0]
        last = indices[-1]
        time = 0.24 * (last - first + 1)

        turning = 0
        #print('path',path)
        for i in indices:
            if inter_15[i]>2.26 or path==2 or path==3 or  path==5:#if path is proton transfer or system without NH2OH formation
                if 1 <= i < len(inter) - 1:
                    if (inter_15[i - 1] - inter_15[i]) >0 and (inter_15[i] - inter_15[i + 1]) < 0:
                        turning += 1
            else:#if NH2OH formed and path is sn2
                if 1 <= i < len(inter) - 1 and inter_15[i]<=2.26:
                    if (inter_10[i - 1] - inter_10[i]) > 0 \
                     and (inter_10[i] - inter_10[i + 1]) < 0:
                        turning += 1
    return turning, time


# ======================================================
# LOOP FOR 10,000 TRAJECTORIES AND SAVE RESULTS
# ======================================================
base = "traj_13_boltzmann"
N = 10000

results = []
for j in range(1,11):
    print("Processed:",j)
    for i in range(1, 10001):
        #print(i,j)
        # Print progress every 1000 iterations
        if i % 1000 == 0:
            print("Processed:", i)

        xyzfile = f"{base}/set{j}/{i}/movie.xyz"

        if not os.path.isfile(xyzfile):
            # store NaN for missing data
            results.append([i, np.nan, np.nan,  np.nan])
            continue
        p=path[(j-1)*10000+(i-1),1]
        #print(p)
        sn2_tn,sn2_tm = analyze_traj(xyzfile,p)
        results.append([i, sn2_tn, sn2_tm])

# convert to array
results = np.array(results, dtype=float)

# save result
np.savetxt(
    f"{base}/all_turning_time_new1.txt",
    results,
    fmt="%7.0f %10.3f" "%10.3f"
)

print("Saved:", results.shape)
