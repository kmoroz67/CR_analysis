def calc_Morisita(df_A, df_B, chain, cdr_type):
    '''
    Calculate Morisita index
    :param a_i(b_i): the T cell counts of clone i in samples A and B 
    :param A(B): the total number of T cells in samples A and B, 
    :param S: list of unique clones in the union of samples A and B.
    :param df: dataframe with patient's clonotype repertoire
    :param chain: name of the chain
    
    :return Morisita index
    '''
    if type(df_A) == str:
        df_A = pd.read_csv(df_A, sep='\t')
        df_B = pd.read_csv(df_B, sep='\t')
    if chain == 'light':
        df_A = df_A[df_A.chain.isin(['kappa', 'lambda'])]
        df_B = df_B[df_B.chain.isin(['kappa', 'lambda'])]
    elif 'Ig' in chain:
        df_A = df_A[df_A.isotype == chain]
        df_B = df_B[df_B.isotype == chain]
    else:
        df_A = df_A[df_A.chain == chain]
        df_B = df_B[df_B.chain == chain]

    A = len(df_A)
    B = len(df_B)

    overlap_cdr = list(set(df_A[cdr_type]).intersection(set(df_B[cdr_type])))

    df_A_overlap = df_A[df_A[cdr_type].isin(overlap_cdr)]
    df_B_overlap = df_B[df_B[cdr_type].isin(overlap_cdr)]
    if cdr_type == 'cdr3aa':
        df_A_overlap = df_A_overlap.groupby(cdr_type)['clonotype_count'].apply(
            sum).reset_index(name='clonotype_count')
        df_B_overlap = df_B_overlap.groupby(cdr_type)['clonotype_count'].apply(
            sum).reset_index(name='clonotype_count')
    df_A_overlap.sort_values(by=[cdr_type], inplace=True)
    df_B_overlap.sort_values(by=[cdr_type], inplace=True)
    df_A_overlap.reset_index(drop=True, inplace=True)
    df_B_overlap.reset_index(drop=True, inplace=True)
    up = sum(2 * df_A_overlap['clonotype_count'] * df_B_overlap['clonotype_count'])

    ai_2 = df_A['clonotype_count'].apply(lambda x: x**2).sum()
    bi_2 = df_B['clonotype_count'].apply(lambda x: x**2).sum()
    down = (A * B * (ai_2 / (A**2) + bi_2 / (B**2)))

    Moris = up / down
    return Moris

def morisita_heatmap(patient_chain_data,
                     chain,
                     cbar_chains,
                     xticklabel_list,
                     yticklabel_list,
                     DPI=150,
                     figsize=(7, 5),
                     fontsize=20):
    if chain in cbar_chains:
        cbar = True
    else:
        cbar = False
    fig, ax = plt.subplots(figsize=figsize, facecolor='white', dpi=DPI)
    ax = sns.heatmap(patient_chain_data,
                     vmin=0,
                     vmax=1,
                     square=True,
                     cbar=cbar)
    if chain not in xticklabel_list:
        ax.set_xticklabels([])
        ax.get_xaxis().set_ticks([])
    else:
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=fontsize)
    if chain not in yticklabel_list:
        ax.set_yticklabels([])
        ax.get_yaxis().set_ticks([])
    else:
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize)
    if 'Ig' not in chain:
        ax.set_title(chain.capitalize(), fontsize=fontsize)
    else:
        ax.set_title(chain, fontsize=fontsize)
    if cbar:
        cbar_obj = ax.collections[0].colorbar
        cbar_obj.ax.tick_params(labelsize=fontsize)
    plt.xlabel('')
    plt.ylabel('')
    plt.gca().invert_yaxis()
    return fig, ax

def morisita_multyplot(sample_names,
                       sample_pathes,
                       chains,
                       cdr_type,
                       cbar_chains,
                       xticklabel_list,
                       yticklabel_list,
                       DPI=150,
                       figsize=(10, 7),
                       fontsize=20):
    samples_dict = dict(zip(sample_names, sample_pathes))
    pair_list = [[x, y] for x in sample_names for y in sample_names]
    for chain in chains:
        morisita_patient_data = {
            'check_list': [],
            'DX': [],
            'DY': [],
            'morisita': []
        }
        for pair in pair_list:
            if pair[::-1] in morisita_patient_data['check_list']:
                morisita_patient_data['morisita'].append(
                    morisita_patient_data['morisita'][
                        morisita_patient_data['check_list'].index(pair[::-1])])
            else:
                df_A = pd.read_csv(samples_dict[pair[0]], sep='\t')
                df_B = pd.read_csv(samples_dict[pair[1]], sep='\t')
                morisita_patient_data['morisita'].append(
                    calc_Morisita(df_A, df_B, chain, cdr_type))
            morisita_patient_data['DX'].append(pair[0])
            morisita_patient_data['DY'].append(pair[1])
            morisita_patient_data['check_list'].append(pair)
        patient_chain_data = pd.DataFrame(morisita_patient_data)
        patient_chain_data = patient_chain_data.pivot('DX', 'DY', 'morisita')
        morisita_heatmap(patient_chain_data,
                         chain,
                         cbar_chains,
                         xtiklabel_list,
                         ytiklabel_list,
                         DPI=DPI,
                         figsize=figsize,
                         fontsize=fontsize)
        
sample_names = [
#     'PBMC spring GM', 'PBMC spring LK', 'Sorted spring GM', 'Sorted spring LK',
#     'Tumor autumn GM', 'PBMC autumn GM', 'Sorted autumn GM', 'Tumor autumn LK',
#     'PBMC autumn LK', 'Sorted autumn LK', '1663 PBMC AB', '1665 CD8 AB',
#     '1666 CD4 AB', 'CP000046 11', '(CP000046) BG001082',
#     '(CP000046) BG001082.2', 'CP000046 14', 'CP000046 PBMC GM',
#     'CP000046 B cells GM', 'BG001082, Cryo PBMCS 4/5/21 GM',
#     'BG001082.2, Cryo PBMCs 8/24/21 GM', 'CP000046, Cryo PBMCs 11/29/21 GM',
    'CP000046 PBMC LK', 'CP000046 B cells LK',
    'BG001082, Cryo PBMCS 4/5/21 LK', 'BG001082.2, Cryo PBMCs 8/24/21 LK',
    'CP000046, Cryo PBMCs 11/29/21 LK'
]

sample_pathes = [
#     '/uftp/Blood/db_calc_pipeline/BG001082_sum/PBMC_Target_GM_spring.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082_sum/PBMC_Target_LK_spring.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082_sum/Sorted_Target_GM_spring.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082_sum/Sorted_Target_LK_spring.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_Tumor_GM/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_PBMC_GM/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_B_Cells_GM/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_Tumor_LK/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_PBMC_LK/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001422_BCR_B_Cells_LK/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001082_WB_TCR_R0001663_PBMC_AB/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001082_WB_TCR_R0001663_PBMC_AB/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/BG001082_WB_TCR_R0001663_PBMC_AB/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211204_NovaC_Sample_11/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211204_NovaC_Sample_12/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211204_NovaC_Sample_13/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211204_NovaC_Sample_14/TCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_12/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_13/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_14/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_15/BCR_processed.txt',
#     '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_16/BCR_processed.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_17/BCR_processed.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_18/BCR_processed.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_19/BCR_processed.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_20/BCR_processed.txt',
    '/uftp/Blood/db_calc_pipeline/BG001082/211202_NovaB_Sample_21/BCR_processed.txt'
]
chains = ['light']
cbar_chains = ['light']
xticklabel_list = ['light']
yticklabel_list = ['light']
cdr_type = 'cdr3nt'
morisita_multyplot(sample_names,
                   sample_pathes,
                   chains,
                   cdr_type,
                   cbar_chains,
                   xticklabel_list,
                   yticklabel_list,
                   DPI=150,
                   figsize=(10, 7),
                   fontsize=20)