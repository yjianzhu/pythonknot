import numpy as np

def bond_length(xyz:np.ndarray,type = "open"):
    """计算轨迹的键长分布,输入应该为N_frames, N_atoms, 3的numpy array. type="open"表示开链，type="ring"表示闭链."""
    #print("original shape:", xyz.shape)
    temp = xyz.transpose(1,2,0)
    #print("temp shape:", temp.shape)
    if(type == "ring"):
        bond_length = np.linalg.norm(temp[1:] - temp[:-1], axis=1)
        # 闭链需要计算首尾原子的键长
        end_to_start_bond = np.linalg.norm(temp[-1] - temp[0], axis=0)
        #print("end_to_start_bond shape:", end_to_start_bond.shape)
        # 从(N_frames) 变为 (1, N_frames)
        end_to_start_bond = end_to_start_bond.reshape(1,-1)
        bond_length = np.concatenate((bond_length, end_to_start_bond), axis=0)
    elif(type == "open"):
        bond_length = np.linalg.norm(temp[1:] - temp[:-1], axis=1)
    print("bond_length shape:", bond_length.shape)
    bond_length = bond_length.flatten()
    return bond_length

def radius_of_gyration(trajectory:np.ndarray):
    """计算轨迹的质心半径分布,输入应该为N_frames, N_atoms, 3的numpy array."""
    # 计算每帧的质心
    centroids = np.mean(trajectory, axis=1)
    
    # 广播质心回到trajectory的形状，以便进行逐元素操作
    # centroids[:, np.newaxis, :] 使得 centroids 从 (N_frames, 3) 扩展为 (N_frames, 1, 3)
    # 这样可以与 trajectory (N_frames, N_atoms, 3) 对齐
    displacement = trajectory - centroids[:, np.newaxis, :]
    
    # 计算每个原子到质心的距离的平方
    distances_squared = np.sum(displacement**2, axis=2)
    
    # 计算每帧的回转半径
    gyration_radius = np.sqrt(np.mean(distances_squared, axis=1))
    
    return gyration_radius

def end_to_end_distance(trajectory:np.ndarray):
    """计算轨迹的端到端距离分布,输入应该为N_frames, N_atoms, 3的numpy array."""
    # 计算每帧的端到端距离
    end_to_end = trajectory[:,-1] - trajectory[:,0]
    end_to_end_distance = np.linalg.norm(end_to_end, axis=1)
    
    return end_to_end_distance