import plotly.express as px
from dash import html
from dash import dcc
import dash
import pandas as pd
import numpy as np
from dash.dependencies import Input, Output, State
import datetime as dt


DATA_DIR = 'GSE185381_TCR/clonotypes'

# fig1: umap, Author's cell types
clinical_meta = pd.read_csv('GSE185381_TCR/clinical_metadata.tsv', sep='\t', index_col=0)

# fig2: single, ambigous clonotypes proportion
single_ambiguos_clonotypes = pd.read_csv(DATA_DIR + '/fraction_clo_df.tsv', sep='\t', index_col=0)
reformated_single_ambiguos_clonotypes_dct = {
    'subject id': [],
    'clonotypes types': [],
    'ratio': []
}
for i in single_ambiguos_clonotypes.index:
    for j in single_ambiguos_clonotypes.columns:
        reformated_single_ambiguos_clonotypes_dct['subject id'].append(i)
        reformated_single_ambiguos_clonotypes_dct['clonotypes types'].append(j)
        reformated_single_ambiguos_clonotypes_dct['ratio'].append(
            single_ambiguos_clonotypes.loc[i, j]
        )
reformated_single_ambiguos_clonotypes_df = pd.DataFrame(reformated_single_ambiguos_clonotypes_dct)

# fig3: clonotypes diversity estimation: Shannon entropy
shannon_df = pd.read_csv(DATA_DIR + '/shannon_df.tsv', sep='\t', index_col=0)

# fig4: clonotypes enrichment:
clonotypes_enrichment = pd.read_csv(DATA_DIR + '/clono_proportion_df.tsv', sep='\t', index_col=0)
reformated_clonotypes_enrichment_dct = {
    'subject id': [],
    'Clonotypes size': [],
    'ratio': []
}
for i in clonotypes_enrichment.index:
    for j in clonotypes_enrichment.columns:
        reformated_clonotypes_enrichment_dct['subject id'].append(i)
        reformated_clonotypes_enrichment_dct['Clonotypes size'].append(j)
        reformated_clonotypes_enrichment_dct['ratio'].append(
            clonotypes_enrichment.loc[i, j]
        )
reformated_clonotypes_enrichment_df = pd.DataFrame(reformated_clonotypes_enrichment_dct)

# fig5: cdr3 alpha, beta length:
cdr3_alpha_length = pd.read_csv(DATA_DIR + '/cdr3alpha_length_df.tsv', sep='\t', index_col=0)
reformated_cdr3_alpha_length_dct = {
    'Condition': [],
    'CDR3 length': [],
    'Ratio': []
}
for i in cdr3_alpha_length.index:
    for j in cdr3_alpha_length.columns:
        reformated_cdr3_alpha_length_dct['Condition'].append(i)
        reformated_cdr3_alpha_length_dct['CDR3 length'].append(j)
        reformated_cdr3_alpha_length_dct['Ratio'].append(
            cdr3_alpha_length.loc[i, j]
        )
reformated_cdr3_alpha_length_df = pd.DataFrame(reformated_cdr3_alpha_length_dct)

cdr3_beta_length = pd.read_csv(DATA_DIR + '/cdr3beta_length_df.tsv', sep='\t', index_col=0)
reformated_cdr3_beta_length_dct = {
    'Condition': [],
    'CDR3 length': [],
    'Ratio': []
}
for i in cdr3_beta_length.index:
    for j in cdr3_beta_length.columns:
        reformated_cdr3_beta_length_dct['Condition'].append(i)
        reformated_cdr3_beta_length_dct['CDR3 length'].append(j)
        reformated_cdr3_beta_length_dct['Ratio'].append(
            cdr3_beta_length.loc[i, j]
        )
reformated_cdr3_beta_length_df = pd.DataFrame(reformated_cdr3_beta_length_dct)

# # fig6: UMAP simplified cell type + UMAP condition proportions
# UMAP_scatter_pie = pd.read_csv(DATA_DIR + '/UMAP_porportion.tsv', sep='\t')


app = dash.Dash(__name__)
app.config.suppress_callback_exceptions = True

fig1 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 1',
    color="Author's cell type",
    labels={'color': "Author's cell type"},
    height=500,
    width=900,
)
fig1.update_traces(marker={'size': 5})
fig1.update_layout(margin_t=50, margin_r=300)

fig11 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 2',
    color="simplified_celltype",
    labels={'color': "simplified_celltype"},
    height=500,
    width=900,
)
fig11.update_traces(marker={'size': 5})
fig11.update_layout(margin_t=50, margin_r=300)

fig12 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 3',
    color="Condition",
    labels={'color': "Condition"},
    height=500,
    width=900,
)
fig12.update_traces(marker={'size': 5})
fig12.update_layout(margin_t=50, margin_r=300)

fig13 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 4',
    color="Subject ID",
    labels={'color': "Subject ID"},
    height=500,
    width=900,
    category_orders={"Subject ID": sorted(clinical_meta['Subject ID'].unique())[::-1]}
)
fig13.update_traces(marker={'size': 5})
fig13.update_layout(margin_t=50, margin_r=300)


fig2 = px.bar(
    reformated_single_ambiguos_clonotypes_df,
    x='subject id',
    y="ratio",
    title='''Ratio of Single / Ambiguous Clonotypes<br><br>
    <sup>A proper clonotype is defined as having only 1 TRA and 1 TRB. However, other types of ambiguous clonotypes might also show some clinical significance.<br>
    From the picture, the number of T cells with only a single chain TRB is abundant across patients.
    </sup>
    ''',
    color="clonotypes types",
    barmode='stack',
    height=700,
    width=1000,
)
fig2.update_layout(margin_t=150)

fig3 = px.box(
    shannon_df,
    x='Condition',
    y='entropy per patient',
    title='''Clonotypes Diversity Estimation<br><br>
    <sup>A higher Shannon entropy value indicates a more diverse population. <br>
    A lower Shannon entropy value in AML disease suggests that the body favors specific types of clonotypes to combat the disease.
    </sup>
    ''',
    color='Condition',
    height=600,
    width=500,
)
fig3.update_layout(margin_t=150)
fig3.update_traces(
    boxpoints='all'
)

fig4 = px.bar(
    reformated_clonotypes_enrichment_df,
    x="subject id",
    y="ratio",
    title='''Clononal Expansion<br><br>
    <sup> TCR clonal expansion across normal and Acute Myeloid Leukocytes patients.<br>
    Although most cells contains unique TCR, AML patients show higher degree of expanded clonotype - clonal expansion (n>5)
    </sup>
    ''',
    color="Clonotypes size",
    barmode = 'stack',
    height=700,
    width=1000,
)
fig4.update_layout(margin_t=150)

fig5 = px.line(
    reformated_cdr3_alpha_length_df,
    x="CDR3 length",
    y="Ratio",
    color="Condition",
    height=500,
    width=800,
    title='CDR3 - Alpha Length across Conditions',
    markers=True
)
fig5.update_layout(yaxis_range=[0, 0.35])
fig5.update_layout(margin_t=50)

fig6 = px.line(
    reformated_cdr3_beta_length_df,
    x="CDR3 length",
    y="Ratio",
    color="Condition",
    height=500,
    width=800,
    title='CDR3 - Beta Length across Conditions',
    markers=True
)
fig6.update_layout(yaxis_range=[0, 0.35])
fig6.update_layout(margin_t=50)

app.layout = html.Div(
    children=[
        html.H1(
            'GSE185381 - TCR ANALYSIS',
            style={'textAlign': 'center', 'color': '#503D36', 'font-size': 40}
        ),
        html.Div([
            dcc.Graph(figure=fig1),
            dcc.Graph(figure=fig11),
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig12),
            dcc.Graph(figure=fig13),
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig2)
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig3)
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig4)
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig5),
            dcc.Graph(figure=fig6),
        ], style={'display': 'flex', 'justify-content': 'center'}),
    ]
)


if __name__ == '__main__':
    app.run_server()