import plotly.express as px
from dash import html
from dash import dcc
import dash
import pandas as pd
import numpy as np
from dash.dependencies import Input, Output, State
import datetime as dt


DATA_DIR = 'GSE185381_TCR/clonotypes'
PLOTLY_DATA_DIR = DATA_DIR + '/plotly_data'

# fig1: umap, Author's cell types
clinical_meta = pd.read_csv('GSE185381_TCR/clinical_metadata.tsv', sep='\t')

# fig2: single, ambigous clonotypes proportion
single_ambiguos_clonotypes_df = pd.read_csv(PLOTLY_DATA_DIR + '/fraction_clo_df.tsv', sep='\t')

# fig3: clonotypes diversity estimation: Shannon entropy
shannon_df = pd.read_csv(PLOTLY_DATA_DIR + '/shannon_df.tsv', sep='\t')

# fig4: clonotypes expansion:
clonotypes_expansion_df = pd.read_csv(PLOTLY_DATA_DIR + '/clono_expansion_df.tsv', sep='\t')

# fig5: cdr3 alpha, beta length:
cdr3_alpha_length_df = pd.read_csv(PLOTLY_DATA_DIR + '/cdr3alpha_length_df.tsv', sep='\t')
cdr3_beta_length_df = pd.read_csv(PLOTLY_DATA_DIR + '/cdr3beta_length_df.tsv', sep='\t')

# # fig6: clonotypes tracing
clonotype_tracing_bar_df = pd.read_csv(PLOTLY_DATA_DIR + '/clono_tracing_bar_df.tsv', sep='\t')
clonotype_tracing_scatter_df = pd.read_csv(PLOTLY_DATA_DIR + '/clono_tracing_scatter_df.tsv', sep='\t')


app = dash.Dash(__name__)
app.config.suppress_callback_exceptions = True

fig1_1 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 1',
    color="Author's cell type",
    labels={'color': "Author's cell type"},
    height=500,
    width=900,
)
fig1_1.update_traces(marker={'size': 3})
fig1_1.update_layout(margin_t=50, margin_r=300)

fig1_2 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 2',
    color="simplified_celltype",
    labels={'color': "simplified_celltype"},
    category_orders={"simplified_celltype": ['naive / central memory T cells', 'transitional T cells', 'terminal effector T cells']},
    height=500,
    width=900,
)
fig1_2.update_traces(marker={'size': 3})
fig1_2.update_layout(margin_t=50, margin_r=300)

fig1_3 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP 3',
    color="Condition",
    labels={'color': "Condition"},
    height=500,
    width=900,
)
fig1_3.update_traces(marker={'size': 3})
fig1_3.update_layout(margin_t=50, margin_r=300)

fig1_4 = px.scatter(
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
fig1_4.update_traces(marker={'size': 3})
fig1_4.update_layout(margin_t=50, margin_r=300)

###

fig2 = px.bar(
    single_ambiguos_clonotypes_df,
    x='Subjects',
    y="Ratio",
    title='''Ratio of Single / Ambiguous Clonotypes<br><br>
    <sup>A proper clonotype is defined as having only 1 TRA and 1 TRB. However, other types of ambiguous clonotypes might also show some clinical significance.<br>
    From the picture, the number of T cells with only a single chain TRB is abundant across patients.
    </sup>
    ''',
    color="Clonotypes types",
    barmode='stack',
    height=700,
    width=1000,
)
fig2.update_layout(margin_t=150)

###

fig3 = px.box(
    shannon_df,
    x='Condition',
    y='Entropy per Subject ID',
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

###

fig4 = px.bar(
    clonotypes_expansion_df,
    x="Groups",
    y="Ratio",
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

####

max_cdr3_ratio = max(
    max(cdr3_alpha_length_df['Ratio'].values),
    max(cdr3_beta_length_df['Ratio'].values),
)
fig5_1 = px.line(
    cdr3_alpha_length_df,
    x="CDR3 length",
    y="Ratio",
    color="Condition",
    height=500,
    width=800,
    title='CDR3 - Alpha Length across Conditions',
    markers=True
)
fig5_1.update_layout(yaxis_range=[0, max_cdr3_ratio + 0.05])
fig5_1.update_layout(margin_t=50)

fig5_2 = px.line(
    cdr3_beta_length_df,
    x="CDR3 length",
    y="Ratio",
    color="Condition",
    height=500,
    width=800,
    title='CDR3 - Beta Length across Conditions',
    markers=True
)
fig5_2.update_layout(yaxis_range=[0, max_cdr3_ratio + 0.05])
fig5_2.update_layout(margin_t=50)

###

fig6_1 = px.bar(
    clonotype_tracing_bar_df,
    x="Groups",
    y="Ratio",
    title='''Clononal Tracing across Groups<br><br>
    <sup> Ratio of TCR clonotypes across groups<br>
    Some diseases / celltypes favors specific clonotype sequences. These 5 clonotypes are the most common across them.
    </sup>
    ''',
    color="Clonotypes",
    category_orders={"Groups": ['naive / central memory T cells', 'transitional T cells', 'terminal effector T cells']},
    barmode='stack',
    height=700,
    width=1000,
)
fig6_1.update_layout(margin_t=150)

fig6_3 = px.scatter(
    clonotype_tracing_scatter_df,
    x='X_UMAP',
    y='Y_UMAP',
    color="cdr3",
    labels={'color': "cdr3"},
    height=500,
    width=900,
    size='umis',
)
fig6_4 = px.scatter(
    clinical_meta,
    x='X_UMAP',
    y='Y_UMAP',
    title='UMAP',
    color_discrete_sequence=['gray'],
    height=500,
    width=900,
    opacity=0.01
)
fig6_4.update_traces(hoverinfo='skip', hovertemplate=None)
fig6_3.add_trace(fig6_4.data[0])
fig6_3.update_layout(margin_t=50, margin_r=300)

###

app.layout = html.Div(
    children=[
        html.H1(
            'GSE185381 - TCR ANALYSIS',
            style={'textAlign': 'center', 'color': '#503D36', 'font-size': 40}
        ),
        html.Div([
            dcc.Graph(figure=fig1_1),
            dcc.Graph(figure=fig1_2),
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig1_3),
            dcc.Graph(figure=fig1_4),
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
            dcc.Graph(figure=fig5_1),
            dcc.Graph(figure=fig5_2),
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig6_1)
        ], style={'display': 'flex', 'justify-content': 'center'}),
        html.Div([
            dcc.Graph(figure=fig1_2),
            dcc.Graph(figure=fig6_3),
        ], style={'display': 'flex', 'justify-content': 'center'}),
    ]
)


if __name__ == '__main__':
    app.run_server()