#Import standard packages
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
#from scipy import io
#from scipy import stats
import pickle
import scipy
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
# If you would prefer to load the '.h5' example file rather than the '.pickle' example file. You need the deepdish package
# import deepdish as dd 

#Import function to get the covariate matrix that includes spike history from previous bins
#from Neural_Decoding.preprocessing_funcs import get_spikes_with_history


import os
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, LSTM, SimpleRNN, GRU, Activation, Dropout
from tensorflow.keras import optimizers

from sklearn import svm, datasets
#from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
import sklearn.metrics
from sklearn.utils.multiclass import unique_labels
import scipy.io
import scipy.signal
import pandas as pd


###
def process_mat_file(beh_fp,Fc_fp,Pc_fp,v_thresh,reward_thresh,start_f,
                     force_remove=None, 
                     recalculate_velocity=False,
                     velocity_scale=0.3,# Might waht to twik it.
                     all_neuron=True,filterdata=True):
    
    beh_data = scipy.io.loadmat(beh_fp)
    Fc_data = scipy.io.loadmat(Fc_fp)
    Pc_data = scipy.io.loadmat(Pc_fp)

    pos=beh_data['behavior']['ybinned'][0,0]
    velocity=beh_data['behavior']['velocity'][0,0]
    
    reward=beh_data['behavior']['reward'][0,0]
    switch_point=start_f#beh_data['behavior']['startframe'][0]
    startframe=switch_point#[0,-1]
    plt.figure
    plt.plot(velocity[0])

    # force calculate velocity
    if recalculate_velocity:
        velocity = np.diff(pos, axis=1)/velocity_scale
        velocity[velocity<0] = 0.0    

    activity = Fc_data['data']['Fc3'][0,0]
    if filterdata:
        activity=Fc3_filter(activity)

    
    PF_num=Pc_data['number_of_PFs'][0]
    PC_id = (PF_num==1)
    if all_neuron:
        PC_activity = activity
    else:
        PC_activity = activity[:,PC_id==True]

    frame_num = pos.shape
    
    reward_id = np.where(reward[0]>reward_thresh)
    #print(reward_id)
    
    remove_both, reward_id = compute_remove_id_v1(reward_id, velocity, v_thresh, reward_thresh)
    
    if force_remove:
        force_remove_mask = np.zeros(reward.shape[1])
        force_remove_mask[np.arange(force_remove[0], force_remove[1])] = 1
        remove_both=np.logical_or(remove_both, force_remove_mask)
    #print(startframe)
    wrong_lap=np.where(reward_id<startframe[0][0])
    #print('here',wrong_lap)
    reward_id=np.delete(reward_id,wrong_lap)
    #print(wrong_lap)
    #print(reward_id)
    return {
        'remove_id': remove_both,
        'PC_activity': PC_activity,
        'PC_id':PC_id,
        'pos': pos,
        'start_frame': startframe,
        'reward_id': reward_id,
        'velocity': velocity
    }

    ###

def Fc3_filter(activity):

    filter_activity=np.zeros(activity.shape)
    for i in range(activity.shape[1]):
        cur_activity=activity[:,i]
        peaks, _ = scipy.signal.find_peaks(cur_activity,distance=150)

        new_range=find_range(cur_activity)
        range_dict=filter_peak(peaks, new_range,cur_activity)
        new_peaks = list(range_dict.keys())
        #print(new_peaks)
        for p in new_peaks:
            #print(p)
            filter_activity[range_dict[p][0]:p,i]=1
            
    return filter_activity

def find_range(activity):
    
    iszero=np.equal(activity, 0)
    absdiff = np.abs(np.diff(iszero))
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    
    new_ranges = []
    i = 0
    while i < len(ranges)-1:
        r = ranges[i]
        #print(i, r)
        j = 1
        r_curr = r
        conn = True
        while conn and i+j+1 < len(ranges):
            r_next = ranges[i+j]
            if r_next[0] - r_curr[1] <= 2:
                #print(r_curr, r_next, 'connect')
                r_curr = r_next
                r_next = ranges[i+j+1]
                j += 1
            else:
                conn = False
        i += j
        new_ranges.append((r[0], r_curr[1]))
        
   # pprint(new_ranges)
    return new_ranges

def filter_peak(peaks,ranges,activity):
    range_dict={}
    for p in peaks:
        for r in ranges:
            if p>r[0] and p<r[1]:
                range_dict[p]=r
    keys = list(range_dict.keys())
    
    for k1 in keys:
        for k2 in keys:
            if k1 not in range_dict.keys() or k2 not in range_dict.keys():
                continue
            elif k1 == k2:
                continue
            #elif (range_dict[k1] == range_dict[k2]).all():
            elif range_dict[k1] == range_dict[k2]:
                #print('dup')
                del range_dict[k2]
        
    filter_r = [tuple(r) for r in range_dict.values()]
    original_r = [tuple(r) for r in ranges]
    
    miss_r = [r for r in original_r if r not in filter_r]
    #print('miss r', miss_r)
    
    for m in miss_r:
        #add_peak=max(activity[m[0]:m[1]])
        add_peak, _ = scipy.signal.find_peaks(activity[m[0]:m[1]],distance=150)
        #print('add',add_peak)
        if add_peak:
            range_dict[add_peak[0]+m[0]]=tuple(m)
    #for r in ranges:
    #miss_range=[r for r in ranges if r not in filter_range]
        
    #print(range_dict)
    return range_dict



def compute_remove_id_v1(reward_id, velocity, v_thresh, reward_thresh,force_remove=None,):
#     print('01_reward_id:', reward_id)
    diff_reward_id = np.diff(reward_id[0])
#     print('diff_reward_id:', diff_reward_id)
    one_diffs = np.where(diff_reward_id == 1)[0]+ 1
#     print('one_diffs:', one_diffs)
    #print(reward_id)
    reward_id = np.delete(reward_id, one_diffs)
    #print(reward_id)    
#     print('reward_id in func:', reward_id)
    freeze_id=np.where(velocity[0]<v_thresh)
#     print(freeze_id)
    #print(velocity.shape)

    remove_end=np.zeros(velocity.shape[1])
    for i in reward_id:
        for j in range(17):
            remove_end[i+j]=1


    remove_freeze=np.zeros(velocity.shape[1])
    # for i in freeze_id
    remove_freeze[freeze_id[0]] = 1


    remove_both=np.logical_or(remove_freeze,remove_end)
    
    if force_remove:
        force_remove_mask = np.zeros(velocity.shape[1])
        force_remove_mask[np.arange(force_remove[0], force_remove[1])] = 1
        remove_both=np.logical_or(remove_both, force_remove_mask)
    
    
    return remove_both, reward_id

###$$ GET_SPIKES_WITH_HISTORY #####
def get_spikes_with_history(neural_data,bins_before,bins_after,bins_current=1):
    """
    Function that creates the covariate matrix of neural activity

    Parameters
    ----------
    neural_data: a matrix of size "number of time bins" x "number of neurons"
        the number of spikes in each time bin for each neuron
    bins_before: integer
        How many bins of neural data prior to the output are used for decoding
    bins_after: integer
        How many bins of neural data after the output are used for decoding
    bins_current: 0 or 1, optional, default=1
        Whether to use the concurrent time bin of neural data for decoding

    Returns
    -------
    X: a matrix of size "number of total time bins" x "number of surrounding time bins used for prediction" x "number of neurons"
        For every time bin, there are the firing rates of all neurons from the specified number of time bins before (and after)
    """

    num_examples=neural_data.shape[0] #Number of total time bins we have neural data for
    num_neurons=neural_data.shape[1] #Number of neurons
    surrounding_bins=bins_before+bins_after+bins_current #Number of surrounding time bins used for prediction
    X=np.empty([num_examples,surrounding_bins,num_neurons]) #Initialize covariate matrix with NaNs
    X[:] = np.NaN
    #Loop through each time bin, and collect the spikes occurring in surrounding time bins
    #Note that the first "bins_before" and last "bins_after" rows of X will remain filled with NaNs, since they don't get filled in below.
    #This is because, for example, we cannot collect 10 time bins of spikes before time bin 8
    start_idx=0
    for i in range(num_examples-bins_before-bins_after): #The first bins_before and last bins_after bins don't get filled in
        end_idx=start_idx+surrounding_bins; #The bins of neural data we will be including are between start_idx and end_idx (which will have length "surrounding_bins")
        X[i+bins_before,:,:]=neural_data[start_idx:end_idx,:] #Put neural data from surrounding bins in X, starting at row "bins_before"
        start_idx=start_idx+1;
    return X


###get X_train
def get_act_pos(experiment,bins_before,bins_after,session_type=1 , bin_num=50 ,lap_range=(0, 1)):
    familiar_laps = experiment['reward_id'][experiment['reward_id'] < experiment['start_frame'][0][1]]
    novel_laps = experiment['reward_id'][experiment['reward_id'] >= experiment['start_frame'][0][1]] 
    # 1. find lap reward ids
    if session_type==1:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(count)
            

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            print(lap_reward_ids)
            if experiment['reward_id'][0]<200:#experiment['start_frame']:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                pass
            else:
                lap_reward_ids.insert(0,bins_before+1)
                print('test',lap_reward_ids)
    elif session_type==2 :
        novel_start=familiar_laps.shape[0]+1
        print(novel_start)
        if lap_range[0]>novel_start:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(lap_reward_ids)

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(lap_reward_ids)
            lap_reward_ids[0]=experiment['start_frame'][0][1]
            
    elif session_type==3:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]          

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            print(lap_reward_ids)
            if experiment['reward_id'][0]>experiment['start_frame'][0][0]:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                lap_reward_ids.insert(0,experiment['start_frame'][0][0])
            else:
                lap_reward_ids[0]=experiment['start_frame'][0][0]
                #lap_reward_ids.insert(0,bins_before+1)
                print('test',lap_reward_ids)
        
#             print(lap_reward_ids)

    # 2. find laps for those laps
    #print(lap_reward_ids)
    lap_mask = np.logical_not(experiment['remove_id'])
    
    
    lap_mask[0:lap_reward_ids[0]-1-bins_before] = False
    lap_mask[lap_reward_ids[0]-1-bins_before:lap_reward_ids[0]-1] = True#
#    plt.figure
#    plt.plot(lap_mask)    
    lap_mask[lap_reward_ids[-1]+1:] = False
    lap_mask[lap_reward_ids[-1]+1:lap_reward_ids[-1]+1+bins_after] = True
    plt.figure(figsize=(20,10))
    plt.plot( experiment['pos'][0])
    plt.plot(lap_mask)
    
    print(experiment['remove_id'].shape)
    print(lap_mask.shape)
    
    masked_pos = experiment['pos'][0, lap_mask]
    masked_PC_activity = experiment['PC_activity'][lap_mask] 
#     print('before',experiment['pos'])
#     print('after',masked_pos)
    start_pos=np.min(masked_pos)
    end_pos=np.max(masked_pos)+0.01
    bin_edges=np.linspace(start_pos,end_pos,bin_num+1)
    binned_pos=np.digitize(masked_pos,bin_edges)

    # post process
#     binary_masked_PC_activity = (masked_PC_activity > 1).astype(np.float32)
    print('masked shape:', masked_PC_activity.shape)   
    print('binned pos shape:', binned_pos.shape)
#     convolved_PC_activity = np.stack([
#         np.convolve(masked_PC_activity[:,i], np.ones(14), mode='same')                      
#         for i in range(X_train.shape[1])], axis=1)
#     print('convolved', convolved_PC_activity.shape)
    return masked_PC_activity, binned_pos
###
def get_act_pos_group2(experiment,bins_before,
                      bins_after,group_bin=90,session_type=1 , bin_num=50 ,lap_range=(0, 1)):
    familiar_laps = experiment['reward_id'][experiment['reward_id'] < experiment['start_frame'][0][1]]
    novel_laps = experiment['reward_id'][experiment['reward_id'] >= experiment['start_frame'][0][1]] 
    # 1. find lap reward ids
    if session_type==1:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(count)
            

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            print(lap_reward_ids)
            if experiment['reward_id'][0]<200:#experiment['start_frame']:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                pass
            else:
                lap_reward_ids.insert(0,bins_before+1)
                print('test',lap_reward_ids)
    elif session_type==2 :
        novel_start=familiar_laps.shape[0]+1
        print('start',novel_start)
        #print(lap_range[0])
        if lap_range[0]>novel_start:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print((lap_reward_ids))

        elif lap_range[0]==novel_start:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1]-1)]
            
            lap_reward_ids=np.insert(lap_reward_ids,0,experiment['start_frame'][0][1])
            #print((lap_reward_ids))
        #print('is_num',lap_reward_ids) 

    elif session_type==3:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]          

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(lap_reward_ids)
            if experiment['reward_id'][0]>experiment['start_frame'][0][0]:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                lap_reward_ids.insert(0,experiment['start_frame'][0][0])
            else:
                lap_reward_ids[0]=experiment['start_frame'][0][0]
                #lap_reward_ids.insert(0,bins_before+1)
                #print('test',lap_reward_ids)
        
             #print(lap_reward_ids)

    # 2. find laps for those laps
    #print('ids',lap_reward_ids)
    lap_mask = np.logical_not(experiment['remove_id'])
    
    
    #lap_mask[0:lap_reward_ids[0]-1-bins_before] = False
    #lap_mask[lap_reward_ids[0]-1-bins_before:lap_reward_ids[0]-1] = True#  
    #lap_mask[lap_reward_ids[-1]+1:] = False
    #lap_mask[lap_reward_ids[-1]+1:lap_reward_ids[-1]+1+bins_after] = True
    #print(lap_reward_ids)
    lap_mask[0:lap_reward_ids[0]-bins_before] = False
    lap_mask[lap_reward_ids[0]-bins_before:lap_reward_ids[0]] = True#  
    lap_mask[lap_reward_ids[-1]:] = False
    lap_mask[lap_reward_ids[-1]:lap_reward_ids[-1]+bins_after] = True
    
    
    
    plt.figure(figsize=(20,10))
    plt.plot( experiment['pos'][0])#[lap_reward_ids[0]-1000:lap_reward_ids[1]])
    plt.plot(lap_mask)#[lap_reward_ids[0]-1000:lap_reward_ids[1]])
    
    #print(experiment['remove_id'].shape)
    #print(lap_mask.shape)
    #print('where',np.where(lap_mask==True)[0:10])
    
    masked_pos = experiment['pos'][0, lap_mask]
    masked_PC_activity = experiment['PC_activity'][lap_mask] 
#     print('before',experiment['pos'])
#     print('after',masked_pos)
    start_pos=np.min(masked_pos)-1e-5
    end_pos=np.max(masked_pos)+1e-5
    print('poses', start_pos, end_pos)
    bin_edges=np.linspace(start_pos,end_pos,bin_num+1)
    binned_pos=np.digitize(masked_pos,bin_edges)
    
    group_lap_num=lap_range[1]-lap_range[0]
    total_group_bin=group_lap_num*group_bin+bins_before+bins_after+1
    group_activity2=np.zeros([total_group_bin,masked_PC_activity.shape[1]])
    
    
    group_edges=np.linspace(start_pos,end_pos,group_bin+1)
    group_pos=np.digitize(masked_pos,group_edges)
    group_activity=np.zeros([group_bin,masked_PC_activity.shape[1]])
    id_seq=groupSequence(group_pos)

    sequence_dict = {}  
    curr_lap = 0
    prev_p = -1

    for i, p in enumerate(group_pos[bins_before:-bins_after]):
        if prev_p-p>80:
        #p == 1 and prev_p == 90:
            curr_lap += 1
        sequence_dict[i] = {
            'id': i,
            'pos': p,
            'lap': curr_lap,
            'activity': masked_PC_activity[i+bins_before, :],   
            'real_pos': masked_pos[i+bins_before],
        }
        prev_p = p    
    laps = np.unique([v['lap'] for v in sequence_dict.values()])
    #print(sequence_dict[1]['lap'])
    remapped_dict = {}
    remapped_pos_dict={}
    neuron_num = masked_PC_activity.shape[1]
    df = pd.DataFrame(sequence_dict.values())
    #print(df)
    
    gdf = df.groupby(['lap', 'pos'], as_index=True).agg({'activity': list})
    
    gdf_pos= df.groupby(['lap', 'pos'], as_index=True).agg({'real_pos': list})
    #print(gdf_pos.loc)
    for l in laps:
        remapped_activity = -1 * np.ones((group_bin, neuron_num))
        remapped_pos = -1 * np.ones((group_bin, 1))
        #if l == 0:
         #   continue
        for p in range(1, group_bin+1):
            if (l, p) in gdf.index:
                remapped_activity[p-1, :] = np.stack(gdf.loc[(l, p)]['activity'], axis=1).mean(axis=1)
                remapped_pos[p-1] = np.mean(gdf_pos.loc[(l, p)]['real_pos'])
                #print(remapped_activity[p-1, :].shape)
            else:
                pass

        
        remapped_dict[l] = remove_minus_one(remapped_activity, axis=0)
        #remapped_pos[l]= remove_minus_one(remapped_pos, axis=0)
    
        #print(len(remapped_pos))
        minus_one_id=np.where(remapped_pos==-1)[0]
        #print(minus_one_id)
        for x in minus_one_id:
            #print(x)
            if x!=len(remapped_pos)-1:
                j=x-1
                k=x+1
                while remapped_pos[j]==-1:
                    j-=1

                while remapped_pos[k]==-1:
                    k+=1
                remapped_pos[x]=(remapped_pos[j]+remapped_pos[k])/2
            else:
                remapped_pos[x]=remapped_pos[x-1]

        remapped_pos_dict[l]=remapped_pos
       
    full_activity = np.concatenate(list(remapped_dict.values()), axis=0)
    full_group_pos=np.concatenate(list(remapped_pos_dict.values()), axis=0).ravel()
    binned_group_pos=np.digitize(full_group_pos,bin_edges)
    
    #print(binned_group_pos)
    
    before_act=np.zeros((bins_before,neuron_num))
    after_act=np.zeros((bins_after,neuron_num))
    
    before_pos=np.ones((bins_before,))*binned_group_pos[0]
    after_pos=np.ones((bins_after,))*binned_group_pos[-1]
    
    full_activity=np.append(before_act,full_activity,axis=0)
    full_activity=np.append(full_activity,after_act,axis=0)
    #print(binned_group_pos.shape)
    #print(before_pos.shape)
    binned_group_pos=np.append(before_pos,binned_group_pos,axis=0)
    binned_group_pos=np.append(binned_group_pos,after_pos,axis=0)
    
    
    #print(binned_group_pos.shape)
    return masked_PC_activity, binned_pos, full_activity, binned_group_pos

def create_group_bin(lap_num,bins_before,bins_after, group_bin=90):
    lap_bin=list(range(group_bin))
    pos=lap_bin[-bins_before:]
    print(lap_bin)
    for i in range(lap_num):
        pos=np.append(pos,lap_bin)

    pos=np.append(pos,lap_bin[0:bins_after])
    pos=pos+1
    
    return pos

def remove_minus_one(arr, axis=0):
    minus_one_idx = np.where(arr == -1)
    #print('minus one rows', minus_one_idx)
    #print(len(arr))
    for x, y in zip(*minus_one_idx):
        if x!=len(arr)-1:
            j = x - 1
            k = x + 1
            while arr[j, y] == -1:
                j -= 1
            while arr[k, y] == -1:
            

                k += 1
            

            arr[x, y] = (arr[j, y] + arr[k, y]) / 2
        else:
            arr[x,y]=0
                
    return arr

def produce_cross_validationset(G_CA3,validation_unit,train_group_X,train_group_pos,test_lap_num,session_type,
                            bins_before,bins_after,bins_current,group_bin,bin_num,trail):

    valid_group_X=np.array([])
    for key in G_CA3.keys():

        
        valid_range=(G_CA3[key]['familiar_lap']+test_lap_num+1+trail*validation_unit,G_CA3[key]['familiar_lap']+test_lap_num+(trail+1)*validation_unit)
        #print(valid_range)


        exp=G_CA3[key]['exp']
        #print(exp['reward_id'][45])

        X_valid,y_valid,validgroup_X,valid_group_pos = get_act_pos_group2(exp, bins_before,bins_after,group_bin,session_type,bin_num,valid_range)

        valid_group_X = np.hstack((valid_group_X, validgroup_X))if valid_group_X.size else validgroup_X

    train_group_X=np.delete(train_group_X, np.s_[bins_before+(trail)*group_bin+1:bins_before+(validation_unit+trail)*group_bin],0)
    train_group_pos=np.delete(train_group_pos, np.s_[bins_before+(trail)*group_bin+1:bins_before+(validation_unit+trail)*group_bin],0)                        

    return train_group_X, train_group_pos,valid_group_X,valid_group_pos


def get_act_pos_group(experiment,bins_before,
                      bins_after,group_bin=50,session_type=1 , bin_num=50 ,lap_range=(0, 1)):
    familiar_laps = experiment['reward_id'][experiment['reward_id'] < experiment['start_frame'][0][1]]
    novel_laps = experiment['reward_id'][experiment['reward_id'] >= experiment['start_frame'][0][1]] 
    # 1. find lap reward ids
    if session_type==1:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(count)
            

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            print(lap_reward_ids)
            if experiment['reward_id'][0]<200:#experiment['start_frame']:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                pass
            else:
                lap_reward_ids.insert(0,bins_before+1)
                print('test',lap_reward_ids)
    elif session_type==2 :
        novel_start=familiar_laps.shape[0]
        print(novel_start)
        if lap_range[0]>novel_start:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(lap_reward_ids)

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            #print(lap_reward_ids)
            lap_reward_ids[0]=experiment['start_frame'][0][1]
            
    elif session_type==3:
        if lap_range[0]>1:
            
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]          

        else:
            lap_reward_ids = [experiment['reward_id'][i] for i in range(lap_range[0]-1, lap_range[1])]
            print(lap_reward_ids)
            if experiment['reward_id'][0]>experiment['start_frame'][0][0]:
                #lap_reward_ids.insert(0,0)
                #lap_reward_ids.pop(0) ### NEED refine
                #print(lap_reward_ids)
                lap_reward_ids.insert(0,experiment['start_frame'][0][0])
            else:
                lap_reward_ids[0]=experiment['start_frame'][0][0]
                #lap_reward_ids.insert(0,bins_before+1)
                print('test',lap_reward_ids)
        
#             print(lap_reward_ids)

    # 2. find laps for those laps
    #print(lap_reward_ids)
    lap_mask = np.logical_not(experiment['remove_id'])
    
    
    lap_mask[0:lap_reward_ids[0]-1-bins_before] = False
    lap_mask[lap_reward_ids[0]-1-bins_before:lap_reward_ids[0]-1] = True#
#    plt.figure
#    plt.plot(lap_mask)    
    lap_mask[lap_reward_ids[-1]+1:] = False
    lap_mask[lap_reward_ids[-1]+1:lap_reward_ids[-1]+1+bins_after] = True
#    plt.figure(figsize=(20,10))
#    plt.plot( experiment['pos'][0])
#    plt.plot(lap_mask)
    
    print(experiment['remove_id'].shape)
    print(lap_mask.shape)
    
    masked_pos = experiment['pos'][0, lap_mask]
    masked_PC_activity = experiment['PC_activity'][lap_mask] 
#     print('before',experiment['pos'])
#     print('after',masked_pos)
    start_pos=np.min(masked_pos)
    end_pos=np.max(masked_pos)+0.01
    bin_edges=np.linspace(start_pos,end_pos,bin_num+1)
    binned_pos=np.digitize(masked_pos,bin_edges)
    
    group_edges=np.linspace(start_pos,end_pos,group_bin+1)
    group_pos=np.digitize(masked_pos,group_edges)
    group_activity=np.zeros([group_bin,masked_PC_activity.shape[1]])
    id_seq=groupSequence(group_pos)
    print((id_seq))
    for i in range(max(id_seq)-1):
        cur_bin=i+1
        #grouppos_id=np.where(group_pos==cur_bin)
        cur_id=np.where(id_seq==cur_bin)
        if cur_bin==1:
            cur_activity=np.mean(masked_PC_activity[cur_id,:],axis=1)
            cur_pos=np.mean(masked_pos[cur_id])
            group_activity=cur_activity
            g_pos=cur_pos
            #print(group_activity.shape)
        else:
            cur_activity=np.mean(masked_PC_activity[cur_id,:],axis=1)
            cur_pos=np.mean(masked_pos[cur_id])
            #print(cur_activity.shape)
            group_activity=np.append(group_activity,cur_activity,axis=0)
            g_pos=np.append(g_pos,cur_pos)
    #group_activity=group_activity.T
    new_binned_pos=np.digitize(g_pos,bin_edges)            

    return masked_PC_activity, binned_pos,new_binned_pos,group_activity

def groupSequence(lst): 
    res2=[[]]
    count=0
    for i in range(len(lst)): 
        if lst[i-1]== lst[i]: 
            res2[-1].append(count)
  
        else: 
            count=count+1
            res2[-1].append(count)
            
    res=np.array(res2[0])
    return res

def get_random_sample(X,X_v,X_t,N_neuron):

    total_N = X.shape[2]
    #print(total_N)
    random_neuron_set = np.random.choice(np.arange(total_N), N_neuron, replace=False)
    X_random=X[:,:,random_neuron_set]
    X_v_random=X_v[:,:,random_neuron_set]
    X_t_random=X_t[:,:,random_neuron_set]
    
    return X_random, X_v_random,X_t_random
###
def get_lap_id(y,bins_before):
    pos_diffs=np.diff(y)
    # plt.plot(pos_diffs)
    lap_id=np.where(pos_diffs< -7)[0]
    if lap_id[0]!=5:
        lap_id=np.insert(lap_id,0,bins_before-1)

    #print(lap_id)
    diff_id=np.diff(lap_id)
    #print(diff_id)
    one_diffs = np.where(diff_id < 5)[0]+ 1
    lap_id = np.delete(lap_id, one_diffs)
    lap_id = np.append(lap_id, lap_id[-1]+90)
    print(lap_id)
    return lap_id
###
def compute_lap_score(X_b_test, y_test,LSTM_model,bins_before,bins_after):
#     pos_diffs=np.diff(y_test)
#     # plt.plot(pos_diffs)
#     lap_id=np.where(pos_diffs< -5)[0]
#     lap_id = np.append(lap_id,-1, y_test.shape[0])
#     print(y_test.shape)
    le = LabelEncoder()
    ohe = OneHotEncoder()
    lap_id = get_lap_id(y_test,bins_before)
    #lap_id.insert(0,bins_before-1)
    
    print(lap_id)
    lap_score=np.zeros(lap_id.shape[0]-1)
    lap_diff = np.zeros(lap_id.shape[0]-1)
    lap_abs_diff = np.zeros(lap_id.shape[0]-1)
    
        
#   #print('lap_id', lap_id.shape)
    for i in range(lap_id.shape[0]-1):
#        if i==0:
#           lap_y_test=y_test[0:lap_id[i]]
#            lap_X_b_test=X_b_test[0:lap_id[i],:]
#        else:
#            lap_y_test=y_test[lap_id[i-1]+1:lap_id[i]]
#            lap_X_b_test=X_b_test[lap_id[i-1]+1:lap_id[i],:]
        
        lap_y_test=y_test[lap_id[i]+1-bins_before:lap_id[i+1]+bins_after]
        lap_X_b_test=X_b_test[lap_id[i]+1-bins_before:lap_id[i+1]+bins_after,:,:]
        print(lap_X_b_test.shape)
        y_test_lst=[[i] for i in y_test]
        Y_test=ohe.fit_transform(y_test_lst).toarray()

        y_test_predicted=LSTM_model.predict(lap_X_b_test[bins_before:-bins_after-1,:,:])
        
        y_test_pre_inverse=ohe.inverse_transform(y_test_predicted)
        #if i==lap_id.shape[0]-2:
        #plt.figure
        #print(len(lap_y_test[bins_before:-bins_after-1]))
        #print(len(y_test_pre_inverse))
        lap_score[i] = sklearn.metrics.r2_score(lap_y_test[bins_before:-bins_after-1], y_test_pre_inverse)
        lap_diff[i] = sklearn.metrics.mean_squared_error(lap_y_test[bins_before:-bins_after-1], y_test_pre_inverse)
        lap_abs_diff[i] = sklearn.metrics.mean_absolute_error(lap_y_test[bins_before:-bins_after-1], y_test_pre_inverse)
    return lap_score, lap_diff,lap_abs_diff
###
def fit_LSTM(X_train, y_train, X_v,Y_valid,epochs, filepath, model_params={'unit': 512, 'layer': 2},force=True):
    assert 'unit' in model_params and 'layer' in model_params
    if os.path.exists(filepath) and not force:
        model = tf.keras.models.load_model(filepath, compile=True)
    else:
        model = Sequential()  # Declare model
        for _ in range(model_params['layer']-1):
            model.add(LSTM(model_params['unit'],
                input_shape=(X_train.shape[1], X_train.shape[2]), return_sequences=True))  # Within recurrent layer, include dropout    
        model.add(LSTM(model_params['unit'],
            input_shape=(X_train.shape[1], X_train.shape[2]), return_sequences=False))  # Within recurrent layer, include dropout
#         model.add(Dropout(0.2))
        
        #model.add(LSTM(X_train.shape[2],
                       #input_shape=(X_train.shape[1], X_train.shape[2])))
        #model.add(LSTM(X_train.shape[2],
        #               input_shape=(X_train.shape[1], X_train.shape[2])))  # Within recurrent layer, include dropout
        model.add(Dense(y_train.shape[1]))
        #sgd = optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
        #sgd=tf.keras.optimizers.SGD(learning_rate=0.01, momentum=0.0, nesterov=False)
        #model.compile(loss='mse', optimizer=sgd, metrics=['accuracy'])  # Set loss function and optimizer

        model.compile(loss='mse', optimizer='rmsprop', metrics=['accuracy'])  # Set loss function and optimizer
    model.fit(X_train, 
              y_train,
              batch_size=None,
              epochs=epochs, 
              verbose=1,
              validation_data=(X_v, Y_valid),
              validation_freq=10)
    #model.evaluate(X_v,Y_valid)
    tf.keras.models.save_model(model, filepath)
    # Fit the model
    # evaluate the model
#     scores = model.evaluate(X_train, y_train)
    scores = model.evaluate(X_v, Y_valid, verbose=0)

    print("\nAccuracy of valid set: %.2f%%" % (scores[1] * 100))

    return model
    ###

def plot_confusion_matrix(y_true, y_pred, classes,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    #print(cm)

    fig, ax = plt.subplots()
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax)
    # We want to show all ticks...
    #ax.set(xticks=np.arange(cm.shape[1]),
           #yticks=np.arange(cm.shape[0]),
           # ... and label them with the respective list entries
           #xticklabels=classes, yticklabels=classes,
           #xticklabels=np.arange(1,51,5), yticklabels=classes,
           #title=title,
           #ylabel='True label',
           #xlabel='Predicted label')
    ax.set(ylabel='True label',xlabel='Predicted label')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    #for i in range(cm.shape[0]):
        #for j in range(cm.shape[1]):
            #ax.text(j, i, format(cm[i, j], fmt),
                    #ha="center", va="center",
                    #color="white" if cm[i, j] > thresh else "black")
    fig.tight_layout()
    return ax

