import pandas as pd
import numpy as np
from pathlib import Path
from XCR_clustering import *
# from Levenshtein import distance

# DB_ASSOCIATIONS = pd.read_pickle('/uftp/Blood/prototype/local_databases/tcr-ab-database.pkl')
DB_ASSOCIATIONS = pd.read_csv('/uftp/Blood/prototype/local_databases/TCR_all.csv')

def is_valid_hla_file(path_to_hla_file):
    '''
    Check if file with patient's HLA exists and it is not broken
    
    :param path_to_hla_file: str, path to file with patient's HLA
    :return: bool value (True if file can be read and it is not empty) 
    '''
    
    if not Path(path_to_hla_file).exists():
        print(f'{path_to_hla_file} does not exist!')
        return False
    try:
        df_hla = pd.read_csv(path_to_hla_file, sep='\t')
    except:
        return False
    return not df_hla.empty


def is_valid_sample_file(path_to_sample_file, chain):
    '''
    Check if file with patient's sample exists and it is not broken and not empty
    
    :param path_to_sample_file: str, path to file with patient's repertoire
    :param chain: str, name of the chain
    :return: bool value (True if file can be read and it is not empty) 
    '''
    
    if not Path(path_to_sample_file).exists():
        print(f'{path_to_sample_file} does not exist!')
        return False
    try:
        df_sample = pd.read_csv(path_to_sample_file, sep='\t')
    except:
        return False
    return not df_sample[df_sample.chain == chain].empty


def filter_cohorts_df(df_cohorts, chain, use_hla):
    '''
    Check dataframe with cohorts and leave only valid samples 
    
    :param df_cohorts: dataframe with 2 columns: first - path to sample folder, second - category
    :param chain: str, name of the chain
    :param use_hla: bool value (True - if we want to intersect database with patient's HLA)
    :return df_cohorts_filtered: filtered dataframe with cohorts
    '''
    
    indexes = [False]*df_cohorts.shape[0]
    for i, path_to_sample in enumerate(df_cohorts.Sample):
        path_to_sample_file = Path(path_to_sample) / 'results/TCR_qc_pass.txt'
        path_to_hla_file = Path(path_to_sample) / 'hla/hla-I.tsv'
        if use_hla:
            indexes[i] = (is_valid_sample_file(path_to_sample_file=path_to_sample_file, chain=chain) and
                          is_valid_hla_file(path_to_hla_file=path_to_hla_file) and use_hla)
        else:
            indexes[i] = is_valid_sample_file(path_to_sample_file=path_to_sample_file, chain=chain)
            
    df_cohorts_filtered = df_cohorts[indexes]   
    df_cohorts_filtered.reset_index(drop=True, inplace=True)
    
    return df_cohorts_filtered

    
def make_db_slice(df_hla, df_sample, chain, category, use_hla):
    '''
    Clean associations database and then intersect it with patient's repertoire
    by V/J/length groups and HLA (optional). 
    
    :param df_hla: dataframe with patient's HLA information
    :param df_sample: dataframe with patient's clonotypes
    :param chain: str, name of the chain
    :param category: str, name of the category to subset database with associations 
                     (e.g consider only Autoimmune)
    :param use_hla: bool value (True - if we want to intersect database with patient's HLA)
    :return db_assoc_slice: subset of the initial database with associations obtained after intersection
    '''
 
    #clean associations database
    db_assoc = DB_ASSOCIATIONS.copy()
    if not category is None:
        db_assoc = db_assoc[db_assoc.Category == category]
    columns_to_subset = ['Disease.name', 'Epitope.sequence', 'HLA', f'CDR3.{chain}.aa', 
                         f'V{chain}', f'J{chain}', 'Pubmed.id', 'Grade']
    db_assoc = db_assoc[columns_to_subset]
    new_cols = {f'CDR3.{chain}.aa': 'cdr3aa', f'V{chain}': 'v_region', f'J{chain}': 'j_region',
                'Epitope.sequence': 'Epitope', 'Disease.name': 'Disease', 'Pubmed.id': 'Pubmed_id'}
    db_assoc = db_assoc.rename(columns=new_cols)
    db_assoc = db_assoc.replace('-', np.nan)
    db_assoc = db_assoc.replace('No', np.nan)
    db_assoc.dropna(inplace=True)
    
    db_assoc = db_assoc.apply(lambda x: x.replace('∗', '*'))
    db_assoc = db_assoc.apply(lambda x: x.replace('.', '*'))
    db_assoc = db_assoc.apply(lambda x: x.replace('•', '*'))
    db_assoc = db_assoc.apply(lambda x: x.replace("'", '*'))
        
    db_assoc['v_region'] = db_assoc['v_region'].apply(lambda x: x.split('*')[0])
    db_assoc['j_region'] = db_assoc['j_region'].apply(lambda x: x.split('*')[0])
    db_assoc['HLA'] = db_assoc['HLA'].apply(lambda x: x.replace('HLA-', ''))
#     db_assoc['HLA'] = db_assoc['HLA'].apply(lambda x: x.split(':')[0])

    db_assoc['cdr3aa'] = db_assoc['cdr3aa'].apply(lambda x: x.upper())
    db_assoc.drop_duplicates(inplace=True, ignore_index=True)
    
    # add new columns
    db_assoc['cdr3aa_length'] = db_assoc.cdr3aa.apply(lambda x: str(len(x)))
    db_assoc['vj_len'] = db_assoc[['v_region', 'j_region', 'cdr3aa_length']].apply(lambda x: ' '.join(x), axis=1)
    db_assoc['chain'] = [chain]*db_assoc.shape[0]
    db_assoc['Source'] = ['base']*db_assoc.shape[0]
    
    if use_hla:
        # get list of patient's hla
        if len(df_hla) == 1:
            df_hla = df_hla.loc[:, ~df_hla.columns.isin(['Unnamed: 0', 'Reads', 'Objective'])]
            df_hla = df_hla.dropna(axis = 1)
            hla_list = list(df_hla.iloc[0, :].unique())
        else:
            hla_list = [df_hla.iloc[:, 0].name]
            hla_list += list(df_hla.iloc[:, 0])
            hla_list = [x.replace('HLA-', '') for x in hla_list]
        # intersect associations databse with patient
        db_assoc = db_assoc.loc[db_assoc['HLA'].isin(hla_list), :]
    db_assoc_slice = db_assoc.loc[db_assoc['vj_len'].isin(df_sample.vj_len.unique()), :]
    
#     diseases = ['Cytomegalovirus', 'Epstein_Barr_Virus', 'Influenza_A']
#     epitopes = ['VSDGGPNLY', 'CTELKLSDY ', 'GILGFVFTL', 'FMYSDFHFI', 'CLGGLLTMV',
#                 'GLCTLVAML', 'NLVPMVATV', 'KTGGPIYKR', 'RVLSFIKGTK', 'ILRGSVAHK',
#                 'RVRAYTYSK', 'RLRAEAQVK', 'SIIPSGPLK', 'AVFDRKSDAK', 'IVTDFSVIK',
#                 'ATIGTAMYK', 'DYCNVLNKEF', 'LPFDKTTVM', 'RPPIFIRRL', 'ELRSRYWAI',
#                 'RAKFKQLL', 'FLRGRAYGL', 'QAKWRLQTL', 'SDEEEAIVAYTL', 'SRYWAIRTR',
#                 'ASCMGLIY', 'RRIYDLIEL', 'YPLHEQHGM', 'IPSINVHHY', 'EENLLDFVRF',
#                 'EFFWDANDIY', 'TPRVTGGGAM']
    
#     db_assoc_slice = db_assoc_slice[db_assoc_slice['Disease'].isin(diseases)].reset_index(drop=True)
#     db_assoc_slice = db_assoc_slice[db_assoc_slice['Epitope'].apply(
#         lambda x: any([distance(x, y) < 3 for y in epitopes]))].reset_index(drop=True)
    
    return db_assoc_slice
    

def make_sample_slice(db_assoc_slice, df_sample, chain, sample_name):
    '''
    Prettify patient's sample and intersect it with associations database by V/J/length groups.
    
    :param db_assoc_slice: cleaned and intersected association database
    :param df_sample: dataframe with patient's repertoire
    :param chain: str, name of the chain
    :param sample_name: str, name of the sample
    :return df_sample_slice: subset of patient's sample obtained after intersection
    '''
    
    df_sample = df_sample[df_sample.chain == chain]
    df_sample = df_sample[['fraction', 'cdr3aa', 'v_region', 'j_region', 'vj_len', 'chain']]
    df_sample['Source'] = [sample_name]*df_sample.shape[0]
    df_sample_slice = df_sample.loc[df_sample['vj_len'].isin(db_assoc_slice.vj_len.unique()), :]

    return df_sample_slice


def combine_db_assoc_with_sample(df_sample, df_hla, chain, sample_name, category, use_hla):
    '''
    Clean association database and patient's sample, intersect them 
    and then combine into one dataframe
    
    :param df_sample: dataframe with patient's repertoire
    :param df_hla: dataframe with patient's HLA information
    :param chain: str, name of the chain
    :param sample_name: str, name of the sample
    :param category: str, name of the category to subset database with associations 
                 (e.g consider only Autoimmune)
    :param use_hla: bool value (True - if we want to intersect database with patient's HLA)
    :return df_sample_assoc: dataframe with clonotypes from patient's sample and associations database
    '''
    
    df_sample['cdr3aa_length'] = df_sample.cdr3aa.apply(lambda x: str(len(x)))
    df_sample['vj_len'] = df_sample[['v_region', 'j_region', 'cdr3aa_length']].apply(lambda x: ' '.join(x), axis=1)
    
    db_assoc_slice = make_db_slice(df_hla=df_hla, 
                                   df_sample=df_sample, 
                                   chain=chain, 
                                   category=category, 
                                   use_hla=use_hla)
    
    df_sample_slice = make_sample_slice(db_assoc_slice=db_assoc_slice, 
                                        df_sample=df_sample,
                                        chain=chain,
                                        sample_name=sample_name)
    
    df_sample_assoc = pd.concat([db_assoc_slice, df_sample_slice])
    df_sample_assoc.reset_index(drop=True, inplace=True)
    
    return df_sample_assoc


def assoc_search(df_sample, df_hla, chain, sample_name, inflation, expansion, pruning_threshold, use_hla=True, category=None):
    '''
    Cluster patient's clonotypes (TCR) with clonotypes from association database
    
    :param df_sample: dataframe with patient's repertoire
    :param df_hla: dataframe with patient's HLA information
    :param chain: str, name of the chain
    :param sample_name: str, name of the sample
    :param inflation: float, MCL inflation factor
    :param expansion: int, MCL expansion factor
    :param pruning_threshold: float, threshold below which matrix elements will be set to 0
    :param use_hla: bool value (True - if we want to intersect database with patient's HLA)
    :param category: str, name of the category to subset database with associations 
                    (e.g consider only Autoimmune)
    :return df_clust: dataframe with mixed clonotypes and assigned clusters
    '''
    
    df_sample_assoc = combine_db_assoc_with_sample(df_sample=df_sample, 
                                                   df_hla=df_hla, 
                                                   chain=chain, 
                                                   sample_name=sample_name, 
                                                   category=category, 
                                                   use_hla=use_hla)
    if df_sample_assoc.empty:
        print(f'There is no intersection between database (category = {category}) and {sample_name}. Abort caclulations!')
        return df_sample_assoc
    
    df_clust = MCL(df_sample=df_sample_assoc, 
                   chain=chain,
                   inflation=inflation,
                   expansion=expansion,
                   pruning_threshold=pruning_threshold, 
                   xcr='TCR',
                   type_seq='aa')
    
    df_clusters_intersected = df_clust.groupby(['n_cluster', 'Source'], 
                                            as_index=False).count()[['n_cluster', 'Source', 'chain']]
    
    df_clusters_intersected = pd.DataFrame(df_clusters_intersected ['n_cluster'].value_counts())
    df_clusters_intersected.reset_index(inplace=True)
    df_clusters_intersected.columns = ['n_cluster', 'counts']
    clusters_intersected = df_clusters_intersected[df_clusters_intersected.counts == 2].n_cluster.values
    if len(clusters_intersected) == 0:
        print(f'{sample_name} clonotypes does not cluster with database (category = {category}). Abort caculations!')
    df_clust = df_clust[df_clust.n_cluster.isin(clusters_intersected)]
    df_clust.reset_index(inplace=True, drop=True)
    
    # find associations
    df_clust['associations'] = [' ']*df_clust.shape[0]
    for cluster in df_clust.n_cluster:
        df_clust_subset = df_clust[df_clust.n_cluster == cluster]
        df_clust_subset_sample = df_clust_subset[df_clust_subset.Source != 'base']
        df_clust_subset_database = df_clust_subset[df_clust_subset.Source == 'base']
        for ind, sample_clonotype in zip(df_clust_subset_sample.index, df_clust_subset_sample.cdr3aa):
            h_dist = df_clust_subset_database.cdr3aa.apply(lambda db_clonotype: hamming(sample_clonotype, db_clonotype))
            diseases = df_clust_subset_database[h_dist == min(h_dist)].Disease.unique().tolist()
            df_clust.at[ind, 'associations'] = diseases
    df_clust['associations'] = df_clust.associations.replace(' ', np.nan)
    
    return df_clust
    

def calc_metrics(df_sample, df_clust, chain, sample_name):
    '''
    Calculate such metrics as clusters_proportion, clonotypes_proportion, cum_fraction, n_epitopes
    
    :param df_sample: dataframe with patient's repertoire
    :param df_clust: dataframe with clustering results
    :param chain: str, name of the chain
    :param sample_name: str, name of the sample
    
    :return metrics_dict: dict with calculated metrics
    '''
    
    n_clonotypes_total = df_sample[df_sample.chain == chain].shape[0]
    clusters_prop = df_clust.n_cluster.nunique() / n_clonotypes_total
    clonotypes_prop = df_clust[df_clust.Source == sample_name].shape[0] / n_clonotypes_total
    cum_fraction = df_clust.fraction.sum()
    n_epitopes = df_clust.Epitope.nunique()
    metrics_dict = {'clusters_proportion': [clusters_prop], 
                    'clonotypes_proportion': [clonotypes_prop],
                    'cum_fraction': [cum_fraction],
                    'n_epitopes': [n_epitopes]}
    
    return metrics_dict
    
    
def cohorts_anaysis(df_cohorts, chain, use_hla, inflation, expansion, pruning_threshold):
    '''
    Clusterize samples from cohort and calculate metrics for clustering results
    
    :param df_cohorts: dataframe with 2 columns: first - path to sample folder, second - category
    :param chain: str, name of the chain
    :param use_hla: bool value (True - if we want to intersect database with patient's HLA)
    :param inflation: float, MCL inflation factor
    :param expansion: int, MCL expansion factor
    :param pruning_threshold: float, threshold below which matrix elements will be set to 0
    :return df_res: dataframe with clusterization metrics
    '''
    
    df_cohorts = filter_cohorts_df(df_cohorts=df_cohorts, chain=chain, use_hla=use_hla)
    df_metrics = pd.DataFrame(columns=['clusters_proportion', 'clonotypes_proportion', 'cum_fraction', 'n_epitopes'])
    df_hla = pd.DataFrame()
    for path_to_sample, category in zip(df_cohorts.Sample, df_cohorts.Category):
        df_sample = pd.read_csv(Path(path_to_sample) / 'results/TCR_qc_pass.txt', sep='\t')
        sample_name = Path(path_to_sample).name
        if use_hla:
            df_hla = pd.read_csv(Path(path_to_sample) / 'hla/hla-I.tsv', sep='\t')

        df_clust = assoc_search(df_sample=df_sample, 
                                df_hla=df_hla, 
                                chain=chain, 
                                sample_name=sample_name,
                                inflation=inflation, 
                                expansion=expansion,
                                pruning_threshold=pruning_threshold,
                                use_hla=use_hla,
                                category=category)
        
        metrics_dict = {'clusters_proportion': [0], 
                        'clonotypes_proportion': [0],
                        'cum_fraction': [0],
                        'n_epitopes': [0]}
            
        if not df_clust.empty:
            metrics_dict = calc_metrics(df_sample=df_sample, 
                                        df_clust=df_clust, 
                                        chain=chain, 
                                        sample_name=sample_name)

        df_metrics = pd.concat([df_metrics, pd.DataFrame(metrics_dict)])
    df_metrics.reset_index(drop=True, inplace=True)
    df_res = pd.concat([df_cohorts, df_metrics], axis=1)
    
    return df_res
        
    