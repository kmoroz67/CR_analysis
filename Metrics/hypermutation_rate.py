import pandas as pd
import sys
import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
sys.path.append('/home/jovyan/Metrics/')
from XCR_clustering import *

def calc_reference_hypermutation_metrics(save_path, diagnosis, material_type, seq_type,
                                         samples_base, age, chains=['heavy', 'lambda', 'kappa']):
    XCR_exist_chain = pd.read_csv(
        '/uftp/Blood/prototype/local_databases/XCR_exist_chain_04_03_22.txt',
        sep='\t')
    if age != None:
        XCR_exist_chain = XCR_exist_chain[XCR_exist_chain['Age'].apply(
            lambda x: True in [
                y in x
                for y in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
            ])].reset_index(drop=True)
        XCR_exist_chain['Age'] = XCR_exist_chain['Age'].astype(float)
        XCR_exist_chain = XCR_exist_chain.query('Age>({}-5)&Age<({}+5)'.format(
            age, age))
    XCR_exist_chain = XCR_exist_chain[
        (XCR_exist_chain['Type'] == seq_type)
        & (XCR_exist_chain['Cell_type'] == material_type) &
        (XCR_exist_chain['Diagnosis'] == diagnosis) &
        (XCR_exist_chain['Base'] == samples_base) &
        (XCR_exist_chain['Stimulation'] == 'Without_stimulation') &
        (XCR_exist_chain['Therapy'] == 'Without_therapy')]
    df_hypermutation_metrics = pd.DataFrame()
    for chain in chains:
        pathes = list(XCR_exist_chain[XCR_exist_chain[chain]]['Path'].apply(lambda x: x+'BCR_qc_pass.txt'))
        names = list(XCR_exist_chain[XCR_exist_chain[chain]]['Path'])
        df_hypermutation_metrics = pd.concat([
            df_hypermutation_metrics,
            calc_hypermutation_metrics(pathes, names, [chain])
        ],
                                             ignore_index=True)
    df_hypermutation_metrics.to_csv('{}hypermutatin_refernce.tsv'.format(save_path), sep='\t', index=False)

def calc_hypermutation_metrics(pathes, names, chains=['heavy', 'lambda', 'kappa']):
    hypermutation_metrics = []
    samples_dict = dict(zip(names, pathes))
    if type(chains[0]) == str:
        chains_dict = {x: chains for x in names}
    else:
        chains_dict = dict(zip(names, chains))
    for name in names:
        df_sample = pd.read_csv(samples_dict[name], sep='\t')
        sample_name = name
        for chain in chains_dict[name]:
            df_clust = MCL(df_sample=df_sample,
                           chain=chain,
                           inflation=2,
                           expansion=2,
                           pruning_threshold=0.001,
                           xcr='BCR',
                           type_seq='nt')
            df_stats = pd.DataFrame(
                df_clust.n_cluster.value_counts()).reset_index()
            df_stats.columns = ['n_cluster', 'cluster_size']
            df_stats = df_stats[df_stats.cluster_size > 1]

            proportion = df_clust.loc[df_clust.n_cluster.isin(
                df_stats.n_cluster)].shape[0] / df_clust.shape[0] * 100
            number = round(10000 / df_clust.shape[0] *
                           df_stats.n_cluster.nunique())
            hypermutation_metrics.append([name, chain, proportion, number])
    df_hypermutation_metrics = pd.DataFrame(hypermutation_metrics)
    df_hypermutation_metrics.columns = [
        'Sample', 'chain', 'fraction',
        'clusters'
    ]
    return df_hypermutation_metrics

def hypermutation_rate_plot(pathes,
                            names,
                            save_path,
                            material_type,
                            seq_type,
                            data_healthy = None,
                            chains=['heavy', 'lambda', 'kappa'],
                            cut_interval_fraction=0,
                            cut_interval_clusters=0):

    if data_healthy == None:
        data_healthy = pd.read_csv('{}hypermutatin_refernce.tsv'.format(save_path), sep='\t')
    else:
        data_healthy = pd.read_csv(data_healthy, sep='\t')
    data_patient = calc_hypermutation_metrics(pathes, names, chains)
    data_healthy = data_healthy[data_healthy['chain'].isin(chains)]
    colors = ['darkred', 'darkgreen', 'darkblue', '#F28500']
    cut_interval_fr = cut_interval_fraction
    cut_interval_cl = cut_interval_clusters

    data_patient.sort_values(by='chain', inplace=True)
    data_healthy.sort_values(by='chain', inplace=True)

    if cut_interval_fr == 0 and cut_interval_cl == 0:
        rows = 1
    else:
        rows = 2

    fig = make_subplots(
        rows=rows,
        cols=2,
        vertical_spacing=0.05,
        subplot_titles=(
            "Fraction of all hypermutated BCRs",
            "Hypermutated BCR clusters per 10 000 unique BCR-clonotypes", "",
            ""))

    # healthy fraction
    for i in range(rows):
        fig.add_trace(
            go.Box(
                name='Healthy donors',
                x=data_healthy["chain"],
                y=data_healthy["fraction"],
                #                 boxpoints="all",
                #                              showlegend=False,
                marker_color='#47B6BC',
                width=.5),
            row=i + 1,
            col=1)

        # healthy clusters
        fig.add_trace(
            go.Box(
                name='Healthy donors',
                x=data_healthy["chain"],
                y=data_healthy["clusters"],
                # boxpoints="all",
                # showlegend=False,
                marker_color='#F4AE63',
                width=.5),
            row=i + 1,
            col=2)

    names = list(data_patient['Sample'].unique())
    pallete = dict(zip(names, colors))
    # patient fraction
    for name in names:
        sample_chains = list(
            data_patient[data_patient['Sample'] == name].chain)
        if len(sample_chains) == 1:
            x0 = 0.025
            width = 0.95
            step = 0.95
        elif len(sample_chains) == 2:
            x0 = 0.02
            width = 0.46
            step = 0.48
        elif len(sample_chains) == 3:
            x0 = 0.015
            width = 0.25
            step = 0.36
        x = 0
        for chain in sample_chains:
            fig.add_hline(y=data_patient[(data_patient['Sample'] == name) & (
                data_patient['chain'] == chain)]['fraction'].iloc[0],
                          x0=x0 + step * x,
                          x1=x0 + width + step * x,
                          row=1,
                          col=1,
                          line_dash='dash',
                          line_color=pallete[name])
            x += 1
        fig.add_trace(
            go.Scatter(
                x=[sample_chains[-1], sample_chains[-1]],
                y=[
                    data_patient[(data_patient['Sample'] == name)
                                 & (data_patient['chain'] == sample_chains[-1]
                                    )]['fraction'].iloc[0],
                    data_patient[(data_patient['Sample'] == name)
                                 & (data_patient['chain'] == sample_chains[-1]
                                    )]['fraction'].iloc[0]
                ],
                mode='lines',
                name=name,
                line=dict(color=pallete[name], width=2.2, dash='dash')))

        # patient clusters
        x = 0
        for chain in sample_chains:
            fig.add_hline(y=data_patient[(data_patient['Sample'] == name) & (
                data_patient['chain'] == chain)]['fraction'].iloc[0],
                          x0=x0 + step * x,
                          x1=x0 + width + step * x,
                          row=1,
                          col=2,
                          line_dash='dash',
                          line_color=pallete[name])
            x += 1

    # fraction updates
    if cut_interval_fr != 0:
        fig.update_yaxes(range=[
            cut_interval_fr[1],
            int(data_healthy['fraction'].max() * 1.1)
        ],
                         row=1,
                         col=1)
        fig.update_xaxes(visible=False, row=1, col=1)
        fig.update_yaxes(range=[0, cut_interval_fr[0]], row=2, col=1)

    # clusters updates
    if cut_interval_cl != 0:
        fig.update_yaxes(range=[
            cut_interval_cl[1],
            int(data_healthy['clusters'].max() * 1.1)
        ],
                         row=1,
                         col=2)
        fig.update_xaxes(visible=False, row=1, col=2)
        fig.update_yaxes(range=[0, cut_interval_cl[0]], row=2, col=2)

    fig.update_layout(title='',
                      plot_bgcolor=('#F4F5F7'),
                      font_family="Source Sans Pro",
                      margin=dict(l=0, r=0, t=30, b=0),
                      legend=dict(yanchor="top", y=1, xanchor="left", x=1))
    fig.update_annotations(font_size=12)
    fig.show()
    # fig.write_image('{}{}_{}_hypermutation.svg'.format(save_path,
    #                                              material_type,
    #                                              seq_type),
    #                  width=780,
    #                  height=340)