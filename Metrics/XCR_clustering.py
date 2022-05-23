import numpy as np
import pandas as pd
from Levenshtein import hamming
import markov_clustering as mc
from pathlib import Path
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelmin

    
def dist_to_nearest(df):
    '''
    Calculate distances to the nearest neighbours normed on cdr3nt length
    
    :param df: dataframe with clonotypes
    :return nn_dist: list with normed distances 
    '''
    
    nn_dist = []
    for vj_len in df['vj_len'].unique():
        df_subset_vj_len = df[df.vj_len == vj_len]
        cdr3seq_length = int(df_subset_vj_len.cdr3seq_length.values[0])
        for cdr3seq_1 in df_subset_vj_len.cdr3seq:
            h_dist = df_subset_vj_len.cdr3seq.apply(lambda cdr3seq_2: hamming(cdr3seq_1, cdr3seq_2)/cdr3seq_length).to_list()
            h_dist.remove(0)
            if len(h_dist) != 0:
                nn_dist.append(min(h_dist))
                
    return nn_dist


def find_threshold(nn_dist):
    '''
    Calculate threshold normed on cdr3nt length for bimodal distribution of distances
    
    :param nn_dist: list with normed distances
    :return threshold: float, normed threshold
    '''
    
    bug_iterator = 0
    bandwidth = 0.03
    number = 10
    nn_dist += list(-1*np.array(nn_dist))
    Matrix = np.array(nn_dist)[:, np.newaxis]
    Matrix_plot = np.linspace(0, max(nn_dist), 150)[:, np.newaxis]
    while number != 1:
        if bug_iterator != 200:
            kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(Matrix)
            aprx = np.exp(kde.score_samples(Matrix_plot))
            number = len(argrelmin(aprx)[0])
            if number == 0:
                threshold = 0.1
                break
            threshold = Matrix_plot[int(argrelmin(aprx)[0][0])][0]
            bandwidth += 0.005
            bug_iterator += 1
        else:
            print('Bad data! Problems with finding threshold.')
            break
            
    return threshold
  
    
def create_adjacency_matrix(df, h_threshold):
    '''
    Calculate adjacency matrix (binary matrix: 0 - edge between nodes is absent, else - 1)
    
    :param df: dataframe with clonotypes
    :param h_threshold: int, distance threshold showing whether we should 
                        draw an edge between clonotypes or not
    :return distances: np.array, 2-dimensional adjacency matrix
    '''
    
    distances = []
    for cdr3seq_1 in df.cdr3seq:
        distances_tmp = df.cdr3seq.apply(lambda cdr3seq_2: 1 if 0 < hamming(cdr3seq_1, cdr3seq_2) <= h_threshold else 0)
        distances.append(distances_tmp.to_list())
    distances = np.array(distances)
    
    return distances
 
    
def markov_clustering_subset(df, inflation, expansion, pruning_threshold, h_threshold, n_slippage, xcr):
    '''
    Markov clustering for subset of clonotypes in the same V/J/length group
    
    :param df: dataframe with clonotypes
    :param inflation: float, MCL inflation factor
    :param expansion: int, MCL expansion factor
    :param pruning_threshold: float, threshold below which matrix elements will be set to 0
    :param h_threshold: int, distance threshold showing whether we should 
                        draw an edge between clonotypes or not
    :param n_slippage: int, number which we should add to each cluster number
    :param xcr: str, type of receptor (BCR or TCR)
    :return df: dataframe with clonotypes and assigned clusters
    '''
    
    if xcr == 'BCR':
        # find threshold for bimodal distribution
        cdr3seq_length = int(df.cdr3seq_length.values[0])
        nn_dist = dist_to_nearest(df=df)
        h_threshold_normed = find_threshold(nn_dist=nn_dist.copy())
        h_threshold = h_threshold_normed*cdr3seq_length   
    adjacency_matrix = create_adjacency_matrix(df=df, h_threshold=h_threshold)
    resulting_matrix = mc.run_mcl(matrix=adjacency_matrix, 
                                  inflation=inflation,
                                  expansion=expansion, 
                                  pruning_threshold=pruning_threshold)
    clusters = mc.get_clusters(resulting_matrix)
    
    # retrieval of cluster labels 
    clusters_obtained = {i+n_slippage+1:tup for i, tup in enumerate(clusters)}
    labels_assigned = {}
    for numb, cluster in clusters_obtained.items():
        for member in cluster:
            labels_assigned[member] = numb
    labels_assigned = dict(sorted(labels_assigned.items(), key=lambda x: x[0]))   
    df['n_cluster'] = labels_assigned.values()
    
    return df


def MCL(df_sample, chain, inflation, expansion, pruning_threshold, h_threshold=1, xcr='TCR', type_seq='aa'):
    '''
    Markov clustering algorithm for TCR/BCR
    
    :param df_sample: dataframe with clonotypes
    :param chain: str, name of the chain
    :param inflation: float, MCL inflation factor
    :param expansion: int, MCL expansion factor
    :param pruning_threshold: float, threshold below which matrix elements will be set to 0
    :param h_threshold: int, distance threshold showing whether we should 
                        draw an edge between clonotypes or not
    :param xcr: str, type of receptor (BCR or TCR)
    :param type_seq: type of cdr3 sequence (aa or nt)
    :return df_clust: dataframe with clonotypes and assigned clusters
    '''
    
    df = df_sample.copy()
    df.rename(columns={f'cdr3{type_seq}': 'cdr3seq'}, inplace=True)
    df = df[df.chain == chain]
    df['cdr3seq_length'] = df.cdr3seq.apply(lambda x: str(len(x)))
    df['vj_len'] = df[['v_region', 'j_region', 'cdr3seq_length']].apply(lambda x: ' '.join(x), axis=1)
    df_clust = pd.DataFrame()
    n_slippage = 0 
    for vj_len in df.vj_len.unique():
        df_subset_vj_len = df[df.vj_len == vj_len]
        if df_subset_vj_len.shape[0] > 1:
            df_clust_subset = markov_clustering_subset(df=df_subset_vj_len,
                                                       inflation=inflation,
                                                       expansion=expansion,
                                                       pruning_threshold=pruning_threshold,
                                                       h_threshold=h_threshold,
                                                       n_slippage=n_slippage,
                                                       xcr=xcr)
            df_clust = pd.concat([df_clust, df_clust_subset])
        else:
            df_subset_vj_len['n_cluster'] = [n_slippage+1]
            df_clust = pd.concat([df_clust, df_subset_vj_len])
        n_slippage = df_clust['n_cluster'].max()
    
    if not df_clust.empty:
        df_clust.drop('cdr3seq_length', axis=1, inplace=True)
        df_clust.rename(columns={'cdr3seq': f'cdr3{type_seq}'}, inplace=True)
        df_clust.sort_values(by=['n_cluster'], inplace=True)
        df_clust.reset_index(drop=True, inplace=True)
        
    return df_clust
