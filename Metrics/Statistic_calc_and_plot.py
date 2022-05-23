import numpy as np
import pandas as pd
import random
import math
import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
import plotly.io as pio
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as patches
import itertools
import os

def calc_clonality_aa(df, chain, clonality_type):
    '''
    Calculate clonality (with loge) of unique aa sequences
    
    :param df: dataframe with patient's clonotype repertoire
    :param chain: name of the chain
    :param clonality_type: type of clonality. Available values:
    1) conality_type = loge 
    2) clonality_type = log2
    3) clonality_type = simpson
    :return clonality
    '''
    
    if chain == 'light':
        df = df[df.chain.isin(['kappa', 'lambda'])]
        sum_i = df['clonotype_count'].sum()
        df['fraction'] = df['clonotype_count'].apply(lambda x: x/sum_i)
    else:
        df = df[df.chain == chain]
    df_aa = df[['cdr3aa', 'fraction']].groupby('cdr3aa', as_index=False).sum()
    if not df_aa.empty:
        if df_aa.shape[0] == 1:
            clonality = 1
        else:
            if clonality_type == 'loge':
                shannon_entropy = -np.sum(list(df_aa.fraction.apply(lambda p: p*np.log(p))))
                clonality = round(1-shannon_entropy/np.log(df_aa.shape[0]), 4)
            elif clonality_type == 'log2':
                shannon_entropy = -np.sum(list(df_aa.fraction.apply(lambda p: p*np.log2(p))))
                clonality = round(1-shannon_entropy/np.log2(df_aa.shape[0]), 4)
            elif clonality_type == 'log10':
                shannon_entropy = -np.sum(list(df_aa.fraction.apply(lambda p: p*np.log10(p))))
                clonality = round(1-shannon_entropy/np.log10(df_aa.shape[0]), 4)
            elif clonality_type == 'simpson':
                clonality = round(np.sqrt(np.sum([p**2 for p in df_aa.fraction])), 4)
        return clonality
    return np.nan

def calc_richness(df, chain):
    '''
    Calculate richness(%) - the ratio (expressed as a percentage) between 
    the number of observed rearrangements in a sample and the number of possible 
    theoretical rearrangements between V families and J genes.
    The number of possible theoretical rearrangements for chains (from imgt):
    alpha: v = 61, j = 61, sum = 3721
    beta: sum = 938
    lambda: v = 79, j = 7, sum = 553
    kappa: v = 77, j = 5, sum = 385
    heavy: v = 110, j = 9, sum = 990
    
    :param df: dataframe with combined clonotypes
    :param chain: name of the chain
    
    :return richness
    '''
    
    if chain == 'light':
        df = df[df.chain.isin(['kappa', 'lambda'])]
        sum_i = df['clonotype_count'].sum()
        df['fraction'] = [i/sum_i for i in df['clonotype_count'].to_list()]
    else:
        df = df[df.chain == chain]
    df['vj'] = df[['v_region', 'j_region']].apply(lambda x: x[0] + ' ' + x[1], axis=1)
    df_vj = df[['vj', 'fraction']].groupby('vj', as_index=False).sum()
    if not df_vj.empty:
        if chain == 'alpha':
            richness = round(100*df_vj.shape[0] / 3721, 4)
        elif chain == 'beta':
            richness = round(100*df_vj.shape[0] / 938, 4)
        elif chain == 'lambda':
            richness = round(100*df_vj.shape[0] / 553, 4)
        elif chain == 'kappa':
            richness = round(100*df_vj.shape[0] / 385, 4)
        elif chain == 'light':
            richness = round(100*df_vj.shape[0] / 938, 4)
        elif chain == 'heavy':
            richness = round(100*df_vj.shape[0] / 990, 4)
        return richness
    return np.nan

def calc_evenness(df, chain):
    '''
    Calculate evenness(%)
    An evenness value for each patient in this study was calculated 
    as the ratio of how many rearrangements among the most frequent 
    were necessary to account for 50 % of the global map intensity 
    (cumulative sum of each rearrangement’s frequency) divided by 
    the total number of rearrangements present. 
    
    :param df: dataframe with patient's clonotype repertoire
    :param chain: name of the chain
    
    :return evenness
    '''
    
    if chain == 'light':
        df = df[df.chain.isin(['kappa', 'lambda'])]
        sum_i = df['clonotype_count'].sum()
        df['fraction'] = [i/sum_i for i in df['clonotype_count'].to_list()]
    else:
        df = df[df.chain == chain]
    df['vj'] = df[['v_region', 'j_region']].apply(lambda x: x[0] + ' ' + x[1], axis=1)
    df_vj = df[['vj', 'fraction']].groupby('vj', as_index=False).sum()
    df_vj = df_vj.sort_values(by='fraction',ascending=False)
    if not df_vj.empty:
        observed_rearrang_numb = 0
        cumulative_fraction = 0
        for fraction in df_vj.fraction:
            if cumulative_fraction >= 0.5:
                break
            cumulative_fraction += fraction
            observed_rearrang_numb += 1
        evenness = round(100*observed_rearrang_numb/df_vj.shape[0], 4)
        return evenness
    else:
        return np.nan

def clonality_plot(patient_data, healthy_data, save_path, nbins=5):
    """
    функция для отрисовки графика клональности;
    :patient_data: dataframe с данными пациента
    :healthy_data: dataframe с данными здоровой когорты
    :save_path: путь сохранения графика
    :samples_name: список с названиями образцов 
    :chain: строка с названием цепи, по которой будет строиться график ('alpha', 'beta', 'heavy', etc)"""

    chains = patient_data['chain'].unique()
    for chain in chains:
        # отфильртровываем данные по нужной цепи
        patient_data_chain = patient_data[patient_data['chain']==chain]
        patient_clonality = list(patient_data_chain['clonality'])
        samples_name = list(patient_data_chain['sample_name'])
        patient_clonality_list = dict(zip(samples_name, patient_clonality))
        healthy_clonality = healthy_data['clonality'].loc[healthy_data['chain'] ==
                                                          chain]

        # healthy reference plot
        fig = px.histogram(healthy_clonality,
                           x='clonality',
                           marginal='violin',
                           histnorm='probability density',
                           color_discrete_sequence=['#cbcbcb'],
                           nbins=nbins,
                           width=485,
                           height=285)

        # patient line
        patient_colors = ['darkred', 'darkgreen', 'darkblue']
        patient_colors_list = dict(zip(samples_name, patient_colors))
        for sample in samples_name:
            fig.add_vline(x=patient_clonality_list[sample],
                          line_width=2.7,
                          line_dash="dash",
                          line_color=patient_colors_list[sample])

            # дублирует команды выше, нужно для корреткного отображения в легенде (vline не отображаются)
            fig.add_trace(
                go.Scatter(x=[patient_clonality_list[sample], patient_clonality_list[sample]],
                           y=[0, 0],
                           mode='lines',
                           name=sample + '    ',
                           line=dict(color=patient_colors_list[sample],
                                     width=2.2,
                                     dash='dash')))

        # выбираем значения оси х
        tick_values = list(
            np.arange(
                round(min(healthy_clonality), 1) - 0.1,
                round(max(max(healthy_clonality), max(patient_clonality)), 2) +
                0.1, 0.1))
        tick_values = [round(i, 2) for i in tick_values] + [
            round(patient_clonality[i], 2) for i in range(len(patient_clonality))
        ]

        fig.update_layout(plot_bgcolor=('#F4F5F7'),
                          margin=dict(l=5, r=5, t=5, b=5),
                          xaxis=dict(tickmode='array',
                                     tick0=0,
                                     tickvals=tick_values,
                                     ticktext=tick_values,
                                     tickangle=0,
#                                      ticklabeloverflow='hide past div',
                                     tickfont={'size': 9}),
                          xaxis2=dict(tickmode='array',
                                      tick0=0,
                                      tickvals=tick_values,
                                      ticktext=tick_values,
                                      tickangle=0,
#                                       ticklabeloverflow='hide past div',
                                      tickfont={'size': 9.5}))

        # legend
        fig.update_layout(font_family="Source Sans Pro",
                          legend=dict(yanchor="top",
                                      y=0.83,
                                      xanchor="left",
                                      x=1,
                                      bgcolor=('white')))

        plot_borders = max(max(patient_clonality), max(healthy_clonality)) / 20

        # выбираем масштаб оси
        min_value = round(min(min(patient_clonality), min(healthy_clonality)),
                          2) - plot_borders
        max_value = round(max(max(patient_clonality), max(healthy_clonality)),
                          2) + plot_borders

        fig.update_xaxes(range=[min_value, max_value], gridwidth=2.3)

        # plot improvements
        fig['layout']['xaxis']['title'] = 'Clonality, ' + chain + ' chain'
        fig['data'][0]['showlegend'] = True
        fig['data'][0]['name'] = 'Healthy Donor  '
        fig['data'][1]['marker']['size'] = 1
        fig.update_annotations(font_size=12)
        fig.write_image(save_path+'_C_{}.svg'.format(chain), width=460, height=290)

def richness_evenness_boxplot(data_patient, data_healthy, save_path):
    """
    data_patient: dataframe (columns: sample_name, chain, evenness, richness)
    data_healthy: reference data
    save_path: path to save the plot 
    patient_code: for legend
    timepoints_num: def=1, number of the timepoints (e.g pre, post = 2 timepoints)
    """

    chains = list(data_patient['chain'].unique())
    chains.sort()
    sample_names = list(set(data_patient['sample_name']))
    sample_names.sort()

    data_healthy = data_healthy.sort_values(by=['chain']).reset_index(
        drop=True)

    palette = {'richness': '#ffe26e', 'evenness': '#ffb38c'}
    counter = 1
    color_list = ['darkred', 'darkgreen', 'darkblue']
    timepoint_palette = dict(zip(sample_names, color_list[:len(sample_names)]))

    fig = make_subplots(rows=1,
                        cols=2,
                        subplot_titles=("Richness", "Evenness"),
                        horizontal_spacing=0.05)

    # plot healthy reference box
    for metric in ["richness", "evenness"]:
        fig.add_trace(go.Box(name='Reference {}   '.format(metric),
                             x=data_healthy['chain'],
                             y=data_healthy[metric],
                             marker_color=palette[metric],
                             boxpoints=False),
                      row=1,
                      col=counter)
        if len(chains) <= 2:
            x0 = 0.1
            width = 0.3
            step = 0.5
        elif len(chains) == 3:
            x0 = 0.06666666667
            width = 0.2
            step = 0.33333333333
        elif len(chains) == 4:
            x0 = 0.05
            width = 0.15
            step = 0.25
#         x0 = 0.1
        for chain in chains:
            for sample_name in sample_names:
                y0 = data_patient.loc[
                    (data_patient['sample_name'] == sample_name)
                    & (data_patient['chain'] == chain)].reset_index(
                        drop=True).at[0, metric]
                fig.add_hline(y=y0,
                              x0=x0,
                              x1=x0 + width,
                              row=1,
                              col=counter,
                              line_dash='dashdot',
                              name=sample_name + '       ',
                              line_color=timepoint_palette[sample_name])

            x0 += step

        counter += 1

    # for legend
    for sample_name in sample_names:
        fig.add_trace(
            go.Scatter(x=[chain, chain],
                       y=[y0, y0],
                       mode='lines',
                       name=sample_name + '       ',
                       line=dict(color=timepoint_palette[sample_name],
                                 width=2.2,
                                 dash='dashdot')))

    # plot updates
    fig.update_layout(title='',
                      plot_bgcolor=('#F4F5F7'),
                      font_family="Source Sans Pro",
                      margin=dict(l=0, r=0, t=30, b=0),
                      legend=dict(yanchor="top", y=1, xanchor="left", x=1))
    fig.update_annotations(font_size=12)

    fig.write_image(save_path+'_RE.svg', width=780, height=340)

def statistic_reference_calc(seq_type, material_type, diagnosis,
                                     samples_base, xcr, save_path):
    '''
    Create statistic df by required data.
    seq_type: BULK or Target
    material_type: WB, PBMC, T_cells, etc.
    diagnosis: Healthy, etc.
    samples_base: Lab, Blood, Sorted, Tissue, TCGA
    xcr: TCR, BCR
    save_path: /*/*/*/
    '''
    data_reference_dict = {'chain': [], 'clonality': [], 'evenness': [], 'richness': []}
    if xcr == 'tcr_ab':
        chains = ['alpha', 'beta']
    else:
        chains = ['heavy', 'kappa', 'lambda', 'light']
    XCR_exist_chain = pd.read_csv(
        '/uftp/Blood/prototype/local_databases/XCR_exist_chain_22_02_22.txt',
        sep='\t')
    XCR_exist_chain = XCR_exist_chain[
        (XCR_exist_chain['Type'] == seq_type)
        & (XCR_exist_chain['Cell_type'] == material_type) &
        (XCR_exist_chain['Diagnosis'] == diagnosis) &
        (XCR_exist_chain['Base'] == samples_base)]
    for chain in chains:
        if chain != 'light':
            XCR_exist_chain = XCR_exist_chain[XCR_exist_chain[chain]]
    if 'tcr' in xcr:
        file = 'TCR_qc_pass.txt'
    else:
        file = 'BCR_qc_pass.txt'
    for sample in XCR_exist_chain['Path']:
        df = pd.read_csv(sample + file, sep='\t')
        for chain in chains:
            data_reference_dict['chain'].append(chain)
            data_reference_dict['clonality'].append(calc_clonality_aa(df, chain, clonality_type='loge'))
            data_reference_dict['evenness'].append(calc_richness(df, chain))
            data_reference_dict['richness'].append(calc_evenness(df, chain))
    data_reference = pd.DataFrame(data_reference_dict)
    data_reference.to_csv(save_path, sep='\t', index=False)

def statistic_plot(samples, names, data_healthy, save_path, chains):
    '''
    Plot clonality/richness/evenness.
    samples: list
    names: list
    data_healthy: df
    save_path: str
    chains: list
    '''
#     xcr_dict = {'TCR_qc_pass.txt': [],
#                 'BCR_qc_pass.txt': []}
#     xcr_dict['tcr'].append([x for x in chains if x in ['alpha', 'beta', 'gamma', 'delta']])
#     xcr_list.append([x for x in chains if x in ['heavy', 'kappa', 'lambda', 'light']])
#     if xcr == 'tcr_ab':
#         chains = ['alpha', 'beta']
#     elif xcr == 'tcr_gd':
#         chains = ['gamma', 'delta']
#     else:
#         chains = ['heavy', 'kappa', 'lambda']
#     for file in ['TCR_qc_pass.txt', 'BCR_qc_pass.txt']:
#         if len(xcr_dict[file])==0:
#             continue
#         else:
    samples_dict = dict(zip(names, samples))
    data_patient_dict = {'sample_name': [], 'chain': [], 'clonality': [], 'evenness': [], 'richness': []}
    for name in names:
        df = pd.read_csv(samples_dict[name], sep='\t')
        for chain in chains:
            data_patient_dict['sample_name'].append(name)
            data_patient_dict['chain'].append(chain)
            data_patient_dict['clonality'].append(calc_clonality_aa(df, chain, clonality_type='loge'))
            data_patient_dict['evenness'].append(calc_richness(df, chain))
            data_patient_dict['richness'].append(calc_evenness(df, chain))
    data_patient = pd.DataFrame(data_patient_dict)
    data_healthy = pd.read_csv(data_healthy, sep='\t')
    data_healthy = data_healthy[data_healthy['chain'].isin(chains)].reset_index(drop=True)
    data_patient['chain'] = data_patient['chain'].apply(
        lambda x: x.capitalize() if 'Ig' not in x else x)
    data_healthy['chain'] = data_healthy['chain'].apply(
        lambda x: x.capitalize() if 'Ig' not in x else x)
    clonality_plot(data_patient, data_healthy, save_path)
    richness_evenness_boxplot(data_patient, data_healthy, save_path)