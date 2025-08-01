import numpy as np
from numba import njit, prange
import math
from sklearn_extra.cluster import KMedoids
from sklearn.preprocessing import MinMaxScaler

geom=np.loadtxt('/home/Tapish/active_learning_2.1/active_learning_2.1.0/nn_peslearn/pes_data_active_learn.txt')
new_geom=np.loadtxt('/home/Tapish/active_learning_2.1/active_learning_2.1.0/selected_geom/try_r_0.1_std_0.01/final_geom.txt')
geom1=geom[:,:-1]
scaler = MinMaxScaler()  # default scales to [0, 1]
geom1 = scaler.fit_transform(geom1)

# Step 2: Compute your custom distance matrix (example: squared Euclidean)
@njit

def custom_distance_matrix(X):
    n = X.shape[0]
    D = np.zeros((n, n))
    bond_list = np.array([
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14],
        [1, 0, 2, 3, 5, 4, 6, 8, 7, 9,10,12,11,13,14],
        [0, 2, 1, 4, 3, 5, 7, 6, 8, 9,11,10,12,13,14],
        [2, 0, 1, 4, 5, 3, 7, 8, 6, 9,11,12,10,13,14],
        [1, 2, 0, 5, 3, 4, 8, 6, 7, 9,12,10,11,13,14],
        [2, 1, 0, 5, 4, 3, 8, 7, 6, 9,12,11,10,13,14]
    ])

    for i in range(n):
        for j in range(i+1, n):  # only compute upper triangle
            f_list = np.empty(len(bond_list))
            for k in range(len(bond_list)):
                diff = X[i, bond_list[k]] - X[j]
                f_list[k] = math.sqrt(np.sum(diff**2))
            d = np.min(f_list)
            D[i, j] = d
            D[j, i] = d  # symmetric

    return D

D_custom = custom_distance_matrix(geom1)


# kmedoids = KMedoids(n_clusters=3, metric='euclidean', init='random', random_state=345,max_iter=1000)
# kmedoids.fit(geom1)

kmedoids = KMedoids(n_clusters=3, metric='precomputed', init='random', random_state=0)
kmedoids.fit(D_custom)
exemplars = geom1[kmedoids.medoid_indices_]
print(exemplars)
list=np.array(exemplars)
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/info.txt",list,fmt=("%5.8f"))
@njit
def cluster_func_numba(row1, row2):
    bond_list = [
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14],
        [1, 0, 2, 3, 5, 4, 6, 8, 7, 9,10,12,11,13,14],
        [0, 2, 1, 4, 3, 5, 7, 6, 8, 9,11,10,12,13,14],
        [2, 0, 1, 4, 5, 3, 7, 8, 6, 9,11,12,10,13,14],
        [1, 2, 0, 5, 3, 4, 8, 6, 7, 9,12,10,11,13,14],
        [2, 1, 0, 5, 4, 3, 8, 7, 6, 9,12,11,10,13,14]
    ]
    f_list=np.empty(len(bond_list))

    for j in prange(len(bond_list)):
        index_list=bond_list[j]
        f = 0.0
        for i in range(len(row1)):
            f += (row1[index_list[i]] - row2[i]) ** 2
        f_list[j]=math.sqrt(f)
    
    return min(f_list)


sub_geom1=np.empty((0,16))
sub_geom2=np.empty((0,16))
sub_geom3=np.empty((0,16))
sub_geom4=np.empty((0,16))
sub_geom5=np.empty((0,16))
sub_geom6=np.empty((0,16))
sub_geom7=np.empty((0,16))
sub_geom8=np.empty((0,16))
sub_geom9=np.empty((0,16))
for i in range(len(geom)):
    # l1=[0,0,0,0,0,0,0,0,0]
    l1=[0,0,0]
    for j in range(len(list)):
        l1[j]=cluster_func_numba(geom1[i],list[j])
        # print(l1)
    min_index = l1.index(min(l1))
    if min_index==0:
        sub_geom1=np.vstack((sub_geom1,geom[i]))
    elif min_index==1:
        sub_geom2=np.vstack((sub_geom2,geom[i]))
    elif min_index==2:
        sub_geom3=np.vstack((sub_geom3,geom[i]))
    # elif min_index==3:
    #     sub_geom4=np.vstack((sub_geom4,geom[i]))
    # elif min_index==4:
    #     sub_geom5=np.vstack((sub_geom5,geom[i]))
    # elif min_index==5:
    #     sub_geom6=np.vstack((sub_geom6,geom[i]))
    # elif min_index==6:
    #     sub_geom7=np.vstack((sub_geom7,geom[i]))
    # elif min_index==7:
    #     sub_geom8=np.vstack((sub_geom8,geom[i]))
    # elif min_index==8:
    #     sub_geom9=np.vstack((sub_geom9,geom[i]))
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom1",sub_geom1,fmt=("%5.8f"))
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom2",sub_geom2,fmt=("%5.8f"))
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom3",sub_geom3,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom4",sub_geom4,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom5",sub_geom5,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom6",sub_geom6,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom7",sub_geom7,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom8",sub_geom8,fmt=("%5.8f"))
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom9",sub_geom9,fmt=("%5.8f"))

np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom1_informat",sub_geom1,fmt=("%5.8f"),delimiter=',')
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom2_informat",sub_geom2,fmt=("%5.8f"),delimiter=',')
np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.3/geom3_informat",sub_geom3,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom4_informat",sub_geom4,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom5_informat",sub_geom5,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom6_informat",sub_geom6,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom7_informat",sub_geom7,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom8_informat",sub_geom8,fmt=("%5.8f"),delimiter=',')
# np.savetxt(f"/home/Tapish/active_learning_2.1/cluster_fitting/cluster_1.1/geom9_informat",sub_geom9,fmt=("%5.8f"),delimiter=',')
