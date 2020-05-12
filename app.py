import dash
import dash_cytoscape as cyto
import dash_html_components as html
from dash.dependencies import Input, Output,State
import dash_core_components as dcc
import pickle
import dash_table
import datetime
import json
import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.figure_factory as ff
import numpy as np
###   Build default table

from special_func import SpecialFun

color_scheme={'up':'#CD6155', 'down':'#48C9B0', 'no change':'#A6ACAF'}

SF=SpecialFun()

genes_def=['PBANKA_1342200.1'] # these are genes we are going to display by default


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__,external_stylesheets=external_stylesheets)
server = app.server


# load steady state and time series data and hierachial clustering
data_df=pickle.load(open("final_cluster_time_SS_data_gene_name.pkl", "rb"))

#load sex related data
sex_df=pickle.load(open("final_sex_df.pkl", "rb"))

# print(data_df.head())

print(sex_df.head())


##  by default to see

# print(data_df.loc['PBANKA_103780.1',:])

##  get median expression level each cluster




med_df=SF.getMedianCluster(data_df)
sub_df=SF.getMediansubCluster(data_df)

pre_style = {"backgroundColor": "#ddd", "fontSize": 20, "padding": "10px", "margin": "10px"}
hidden_style = {"display": "none"}
hidden_inputs = html.Div(id="hidden-inputs", style=hidden_style, children=[])

app.layout = html.Div(children=[
                html.Div(
                [
                        html.Img(
                        src=app.get_asset_url("umu-logo.svg"),
                        className="logo2",
                    ),
                ],
                className="row",
            ),

            ## we put list of publications
            html.H6('Clustering analysis using AP2 KO and AP2-G over-expression data set',className='p'),
            html.Div('Publications:',style={'color': 'blue', 'fontSize': 12},className='paper'),
            html.A('(1) A Knockout Screen of ApiAP2 Genes Reveals Networks of Interacting Transcriptional Regulators Controlling the Plasmodium Life Cycle Cell Host Microbe. 2017 Jan 11; 21(1): 11–22.' ,href='https://www.sciencedirect.com/science/article/pii/S1931312816305145?via%3Dihub',style={'fontSize': 9},className='paper'),
            html.Br([]),
            html.A('(2) Inducible developmental reprogramming redefines commitment to sexual development in the malaria parasite Plasmodium berghei Nat Microbiol. 2018 Nov;3(11):1206-1213.',href='https://www.nature.com/articles/s41564-018-0223-6',style={'fontSize': 9},className='paper'),
            html.Br([]),

            ## we are going to search genes

            html.P('Search genes (e.g., PBANKA_1342200)',className='p'),
            html.Div([
                        html.Div(
                                [dcc.Textarea(
                                id='input_text_gene',
                                placeholder='Enter Genes',
                                rows=20,
                                style={'width': '60%'}
                                 )],className="six columns",style={'marginLeft': 15}),

                                html.Div(html.Button('Submit', id='gene_button'),className="six columns",style={'marginLeft': 15}),

                                ], className="row"),
            # paste table here
            html.Br([]),
            # load_default_table(genes_def),
            html.Div(dash_table.DataTable(id='data-table-cluster'),className='table'),
            # html.Div(id="gene-table_id"),
            html.Br([]),
            # plot dendogram
            html.Div('Level-wise clustering',style={'color': 'blue', 'fontSize': 16},className='paper'),





            # dcc.Graph(id='datatable-dendogram'),
            html.Div([

                # plot Level-1 clustering
            html.Div([

                # dcc.Graph(id='datatable-dendogram',config={'displayModeBar': False})
                dcc.Graph(id='datatable-dendogram')

                    ], style={'width': '25%', 'display': 'inline-block'}, id='Level-1-div'),

            # plot Level-2 clustering
            html.Div([

                # dcc.Graph( id='dendogram-level2',config={'displayModeBar': False})
                dcc.Graph( id='dendogram-level2')

             ], style={'width': '25%', 'display': 'inline-block'}, id='Level-2-div'),

             # plot Level-3 clustering
            html.Div([

                # dcc.Graph( id='dendogram-level3',config={'displayModeBar': False})
                dcc.Graph( id='dendogram-level3')

             ], style={'width': '50%', 'display': 'inline-block'}, id='Level-3-div')

            ], style={ 'padding': '5px 5px' }),

            # ## hidden input box
            # dcc.Input(id="input-time1"),
            # ## hidden input box
            # dcc.Input(id="input-time2"),
            # ## hidden input box
            # dcc.Input(id="input-time3"),
            html.Pre(id='update-on-click-data', style=pre_style),
            hidden_inputs

            ],




            className="page")


@app.callback(
    [Output('data-table-cluster', 'data'),
     Output('data-table-cluster', 'columns'),
     Output('data-table-cluster', 'style_cell_conditional'),
     Output('data-table-cluster', 'style_data_conditional'),
     Output('data-table-cluster', 'style_header'),
     Output('data-table-cluster', 'row_selectable')],
    [Input('gene_button', 'n_clicks')],
    [State("input_text_gene", 'value')]
)

def plot_cluster_dendogram(nclick,val):
    # get  data frame by PBANKA id

    if not val==None:

        genes = [x.strip() for x in val.split(',')]

        found=[]
        not_found=[]

        for gene in genes:
            g=gene+'.1'
            if g in data_df.index:
                found.append(g)
            else:
                not_found.append(g)

        # identify which genes are found and which genes are not found

        if len(found)==0:
            res_val=html.Div([

                html.P('Not a single gene is found in our data',className="p"),
            ])
        else:
            found_genes=found
            data,columns,style_cell_conditional,style_data_conditional,style_header,row_selectable=SF.load_default_table(data_df,found)
    else:
        found_genes=genes_def
        data,columns,style_cell_conditional,style_data_conditional,style_header,row_selectable=SF.load_default_table(data_df,genes_def)


    # now we are going to plot dendogram based on found genes for clustering
    # We are going to gnerate deafult figure
    # first we generate  high level clusters which was computed based on coexpression and TOM similarity
    # then we used median expression values to see association based on clusters
    # color_labels=[data_df.loc[gene,'cluster'] for gene in found_genes]
    # # print('color', color_labels)
    # height=1000
    # orig_labels=med_df.index
    #
    # fig=SF.load_dendogram_cluster(med_df,orig_labels,height,color_labels)
    # fig['layout']['yaxis']['showticklabels']=False
    return data,columns,style_cell_conditional,style_data_conditional,style_header,row_selectable

# # slect rows and perform another dendogram
# @app.callback(
#     [Output('datatable-dendogram', "figure"),Output('dendogram-level2', "figure"),Output('dendogram-level3', "figure")],
#     [Input('data-table-cluster', "data"),
#      Input('data-table-cluster', "selected_rows")])
# def plot_dendogram_subClusterLevel(data, rows):
#
#     if rows is None:
#
#         # find cluster and sub_cluster
#
#         cluster=data_df.loc[genes_def,'cluster'].to_list()[0]
#         sub_cluster=data_df.loc[genes_def,'sub_cluster'].to_list()[0]
#         target_gene=genes_def[0]  # target genes
#     else:
#         sub_cluster=data[rows[0]]['sub_cluster']
#         cluster=data[rows[0]]['cluster']
#         target_gene=data[rows[0]]['Gene ID']
#
#
#
#      # now we are going to plot dendogram based on found genes for clustering
#     # We are going to gnerate deafult figure
#     # first we generate  high level clusters which was computed based on coexpression and TOM similarity
#     # then we used median expression values to see association based on clusters
#
#     found_genes=[target_gene]
#     color_labels=[data_df.loc[gene,'cluster'] for gene in found_genes]
#     # print('color', color_labels)
#     height=1000
#     orig_labels=med_df.index
#
#     fig=SF.load_dendogram_cluster(med_df,orig_labels,height,color_labels)
#
#
#     # print('here', data,cluster,sub_cluster)
#     # get dendogram for sub cluster levels
#     sub_clusters=data_df[data_df['cluster']==cluster]['sub_cluster'].unique()
#     sub_inds=[ cluster+item for item in sub_clusters]
#     tmp=sub_df.loc[sub_inds,:].copy()
#
#
#     # fig1= ff.create_dendrogram(tmp.to_numpy(), orientation='right')
#     color_labels=[sub_cluster]
#     # print(rows,cluster,sub_cluster,sub_clusters,tmp,subcluster_clickData)
#
#     # print('color', color_labels)
#     height=1000
#     orig_labels=sub_clusters
#     fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
#
#     # find all sub cluster genes
#     tmp=data_df[data_df['sub_cluster']==sub_cluster]
#
#     # get rapamycin yes columns
#     rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
#     tmp=tmp.loc[:,rapa_yes_cols].copy()
#
#     # get original labels
#
#     tmp_sex=sex_df.loc[tmp.index,:].copy()
#     orig_labels=sex_df.loc[tmp.index,'label']
#     height=1000
#     target_gene_label=sex_df.loc[target_gene,'label']
#     # print(target_gene_label,target_gene)
#     fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)
#
#
#     # ################ we are going to add naviagation through the clicking figures
#     #
#     # # print(cluster_clickData,subcluster_clickData,gene_clickData)
#     #
#     # # get cluster level
#     #
#     # if (not cluster_clickData==None) and 'text' in cluster_clickData['points'][0].keys():
#     #     # we will update figure at sub-cluster level
#     #     cluster_name=cluster_clickData['points'][0]['text']
#     #     # update figure
#     #     sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
#     #     # import pdb;pdb.set_trace()
#     #     sub_inds=[ cluster_name+item for item in sub_clusters]
#     #     tmp=sub_df.loc[sub_inds,:].copy()
#     #
#     #     color_labels=[sub_clusters[0]]
#     #     # print('color', color_labels)
#     #     height=1000
#     #     orig_labels=sub_clusters
#     #     fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
#     #     cluster_clickData=None
#     #
#     # if (not subcluster_clickData==None) and 'text' in subcluster_clickData['points'][0].keys():
#     #      # find all sub cluster genes
#     #     sub_cluster=subcluster_clickData['points'][0]['text']
#     #     tmp=data_df[data_df['sub_cluster']==sub_cluster]
#     #
#     #     # get rapamycin yes columns
#     #     rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
#     #     tmp=tmp.loc[:,rapa_yes_cols].copy()
#     #
#     #     # get original labels
#     #
#     #     tmp_sex=sex_df.loc[tmp.index,:].copy()
#     #     orig_labels=sex_df.loc[tmp.index,'label']
#     #     # print('orig',orig_labels,sub_cluster)
#     #     height=1000
#     #     target_gene_label=sex_df.loc[tmp.index[0],'label']
#     #     # print(target_gene_label,target_gene)
#     #     fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)
#     #
#     #     subcluster_clickData=None
#
#     #### End of the navigatipns
#
#
#
#      # # re adjust height
#     if len(orig_labels)>70:
#         height=int(len(orig_labels)*(height/70))
#
#     fig2.update_layout(height=height)
#     fig1.update_layout(height=height)
#     fig.update_layout(height=height)
#     # import pdb;pdb.set_trace()
#     # l1_fig.update_layout(height=height)
#     # fig.show()
#
#
#     return fig,fig1,fig2


#

# # slect rows and perform another dendogram
# @app.callback(
#     [Output('datatable-dendogram', "figure"),Output('dendogram-level2', "figure"),Output('dendogram-level3', "figure")],
#     [Input('data-table-cluster', "data"),
#      Input('data-table-cluster', "selected_rows"),Input('datatable-dendogram', "clickData"),Input('dendogram-level2', "clickData"),Input('dendogram-level3', 'clickData')])
# def plot_dendogram_subClusterLevel(data, rows,cluster_clickData,subcluster_clickData,gene_clickData):
#


# # get time stamp
# @app.callback(Output('input-time1', "value"),
#               [Input('datatable-dendogram', "clickData")])
#
# def time1(clickData1):
#     now = datetime.datetime.now()
#     print(now,clickData1)
#     return now

# # slect rows and perform another dendogram
# @app.callback(
#     [Output('datatable-dendogram', "figure")],
#     [Input('data-table-cluster', "data"),
#      Input('data-table-cluster', "selected_rows")])
# def plot_dendogram_subClusterLevel(data, rows):
#
#     if rows is None:
#
#         # find cluster and sub_cluster
#
#         cluster=data_df.loc[genes_def,'cluster'].to_list()[0]
#         sub_cluster=data_df.loc[genes_def,'sub_cluster'].to_list()[0]
#         target_gene=genes_def[0]  # target genes
#     else:
#         sub_cluster=data[rows[0]]['sub_cluster']
#         cluster=data[rows[0]]['cluster']
#         target_gene=data[rows[0]]['Gene ID']
#
#
#
#      # now we are going to plot dendogram based on found genes for clustering
#     # We are going to gnerate deafult figure
#     # first we generate  high level clusters which was computed based on coexpression and TOM similarity
#     # then we used median expression values to see association based on clusters
#
#     found_genes=[target_gene]
#     l1_color=[data_df.loc[gene,'cluster'] for gene in found_genes]
#     # print('color', color_labels)
#     height=1000
#     l1_df=med_df.copy()    ######## first level of clustering
#     l1_label=med_df.index
#
#     # fig=SF.load_dendogram_cluster(med_df,orig_labels,height,color_labels)
#
#
#     # print('here', data,cluster,sub_cluster)
#     # get dendogram for sub cluster levels
#     sub_clusters=data_df[data_df['cluster']==cluster]['sub_cluster'].unique()
#     sub_inds=[ cluster+item for item in sub_clusters]
#     l2_df=sub_df.loc[sub_inds,:].copy()  ######## second level of clustering
#
#
#     # fig1= ff.create_dendrogram(tmp.to_numpy(), orientation='right')
#     color_labels=[sub_cluster]
#     # print(rows,cluster,sub_cluster,sub_clusters,tmp,subcluster_clickData)
#
#     # print('color', color_labels)
#
#     orig_labels=sub_clusters
#     # fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
#
#     # find all sub cluster genes
#     tmp=data_df[data_df['sub_cluster']==sub_cluster]
#
#     # get rapamycin yes columns
#     rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
#     l3_df=tmp.loc[:,rapa_yes_cols].copy()
#
#     # get original labels
#
#     tmp_sex=sex_df.loc[tmp.index,:].copy()
#     orig_labels=sex_df.loc[tmp.index,'label']
#     height=1000
#     target_gene_label=sex_df.loc[target_gene,'label']
#
#     fig=SF.subplot_dendogram(l1_df,l2_df,l3_df,tmp_sex,l1_label,l1_color)
#     # print(target_gene_label,target_gene)
#     # fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)
#
#
#     # ################ we are going to add naviagation through the clicking figures
#     #
#     # # print(cluster_clickData,subcluster_clickData,gene_clickData)
#     #
#     # # get cluster level
#     #
#     # if (not cluster_clickData==None) and 'text' in cluster_clickData['points'][0].keys():
#     #     # we will update figure at sub-cluster level
#     #     cluster_name=cluster_clickData['points'][0]['text']
#     #     # update figure
#     #     sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
#     #     # import pdb;pdb.set_trace()
#     #     sub_inds=[ cluster_name+item for item in sub_clusters]
#     #     tmp=sub_df.loc[sub_inds,:].copy()
#     #
#     #     color_labels=[sub_clusters[0]]
#     #     # print('color', color_labels)
#     #     height=1000
#     #     orig_labels=sub_clusters
#     #     fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
#     #     cluster_clickData=None
#     #
#     # if (not subcluster_clickData==None) and 'text' in subcluster_clickData['points'][0].keys():
#     #      # find all sub cluster genes
#     #     sub_cluster=subcluster_clickData['points'][0]['text']
#     #     tmp=data_df[data_df['sub_cluster']==sub_cluster]
#     #
#     #     # get rapamycin yes columns
#     #     rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
#     #     tmp=tmp.loc[:,rapa_yes_cols].copy()
#     #
#     #     # get original labels
#     #
#     #     tmp_sex=sex_df.loc[tmp.index,:].copy()
#     #     orig_labels=sex_df.loc[tmp.index,'label']
#     #     # print('orig',orig_labels,sub_cluster)
#     #     height=1000
#     #     target_gene_label=sex_df.loc[tmp.index[0],'label']
#     #     # print(target_gene_label,target_gene)
#     #     fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)
#     #
#     #     subcluster_clickData=None
#
#     #### End of the navigatipns
#
#
#
#      # # re adjust height
#     # if len(orig_labels)>70:
#     #     height=int(len(orig_labels)*(height/70))
#     #
#     # fig2.update_layout(height=height)
#     # fig1.update_layout(height=height)
#     # fig.update_layout(height=height)
#     # import pdb;pdb.set_trace()
#     # l1_fig.update_layout(height=height)
#     # fig.show()
#
#
#     return fig


# # slect rows and perform another dendogram
# @app.callback(
#     [Output('datatable-dendogram', "figure"),Output('dendogram-level2', "figure"),Output('dendogram-level3', "figure")],
#     [Input('datatable-dendogram', "clickData"),Input('dendogram-level2', "clickData"),Input('dendogram-level3', 'clickData')])
# def navigate_clusters(cluster_clickData,subcluster_clickData,gene_clickData):
#
#      ################ we are going to add naviagation through the clicking figures
#
#     # print(cluster_clickData,subcluster_clickData,gene_clickData)
#
#     # get cluster level
#
#     if (not cluster_clickData==None) and 'text' in cluster_clickData['points'][0].keys():
#         # we will update figure at sub-cluster level
#         cluster_name=cluster_clickData['points'][0]['text']
#         # update figure
#         sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
#         # import pdb;pdb.set_trace()
#         sub_inds=[ cluster_name+item for item in sub_clusters]
#         tmp=sub_df.loc[sub_inds,:].copy()
#
#         color_labels=[sub_clusters[0]]
#         # print('color', color_labels)
#         height=1000
#         orig_labels=sub_clusters
#         fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
#
#
#     if (not subcluster_clickData==None) and 'text' in subcluster_clickData['points'][0].keys():
#          # find all sub cluster genes
#         sub_cluster=subcluster_clickData['points'][0]['text']
#         tmp=data_df[data_df['sub_cluster']==sub_cluster]
#
#         # get rapamycin yes columns
#         rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
#         tmp=tmp.loc[:,rapa_yes_cols].copy()
#
#         # get original labels
#
#         tmp_sex=sex_df.loc[tmp.index,:].copy()
#         orig_labels=sex_df.loc[tmp.index,'label']
#         # print('orig',orig_labels,sub_cluster)
#         height=1000
#         target_gene_label=sex_df.loc[tmp.index[0],'label']
#         # print(target_gene_label,target_gene)
#         fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)
#
#
#
#     #### End of the navigatipns
#
#
#
#      # # re adjust height
#     if len(orig_labels)>70:
#         height=int(len(orig_labels)*(height/70))
#
#     fig2.update_layout(height=height)
#     fig1.update_layout(height=height)
#     fig.update_layout(height=height)
#     # import pdb;pdb.set_trace()
#     # l1_fig.update_layout(height=height)
#     # fig.show()
#
#
#     return fig,fig1,fig2


def last_clicked(*dash_input_keys):
    """ Get the clickData of the most recently clicked graph in a list of graphs.
    The `value` you will receive as a parameter in your callback will be a dict. The keys you will want to
    pay attention to are:
        - "last_clicked": the id of the graph that was last clicked
        - "last_clicked_data": what clickData would usually return
    This function working depends on a `hidden_inputs` variable existing in the global / file scope. It should be an
    html.Div() input with styles applied to be hidden ({"display": "none"}).
    but why, I hear you ask?
    clickData does not get set back to None after you've used it. That means that if a callback needs the latest
    clickData from two different potential clickData sources, if it uses them both, it will get two sets of clickData
    and no indication which was the most recent.
    :type dash_input_keys: list of strings representing dash components
    :return: dash.dependencies.Input() to watch value of
    """
    # import pdb;pdb.set_trace()
    dash_input_keys = sorted(list(dash_input_keys))
    str_repr = str(dash_input_keys)
    last_clicked_id = str_repr + "_last-clicked"
    existing_child = None
    for child in hidden_inputs.children:
        if child.id == str_repr:
            existing_child = child
            break

    if existing_child:
        return Input(last_clicked_id, 'value')

    # If we get to here, this is the first time calling this function with these inputs, so we need to do some setup
    # make feeder input/outputs that will store the last time a graph was clicked in addition to it's clickdata
    if existing_child is None:
        existing_child = html.Div(id=str_repr, children=[])
        hidden_inputs.children.append(existing_child)

    input_clicktime_trackers = [str_repr + key + "_clicktime" for key in dash_input_keys]
    existing_child.children.append(dcc.Input(id=last_clicked_id, style=hidden_style, value=None))
    for hidden_input_key in input_clicktime_trackers:
        existing_child.children.append(dcc.Input(id=hidden_input_key, style=hidden_style, value=None))

    # set up simple callbacks that just append the time of click to clickData
    for graph_key, clicktime_out_key in zip(dash_input_keys, input_clicktime_trackers):

        @app.callback(Output(clicktime_out_key, 'value'),
                      [Input(graph_key, 'clickData')],
                      [State(graph_key, 'id')])
        def update_clicktime(clickdata, graph_id):

            result = {
                "click_time": datetime.datetime.now().timestamp(),
                "click_data": clickdata,
                "id": graph_id
            }
            return json.dumps(result)


    cb_output = Output(last_clicked_id, 'value')
    cb_inputs = [Input(clicktime_out_key, 'value') for clicktime_out_key in input_clicktime_trackers]
    cb_current_state = State(last_clicked_id, 'value')

    # use the outputs generated in the callbacks above _instead_ of clickData
    @app.callback(cb_output, cb_inputs, [cb_current_state])
    def last_clicked_callback(*inputs_and_state):
        clicktime_inputs = inputs_and_state[:-1]
        last_state = inputs_and_state[-1]
        if last_state is None:
            last_state = {
                "last_clicked": None,
                "last_clicked_data": None,
            }
        else:
            largest_clicktime = -1
            largest_clicktime_input = None
            for clicktime_input in clicktime_inputs:
                clicktime_input = json.loads(clicktime_input)
                click_time = int(clicktime_input['click_time'])
                if clicktime_input['click_data'] and click_time > largest_clicktime:
                    largest_clicktime_input = clicktime_input
                    largest_clicktime = click_time
            if largest_clicktime:
                # largest_clicktime_input=json.loads(largest_clicktime_input)

                last_state=json.loads(last_state)
                last_state['last_clicked'] = largest_clicktime_input["id"]
                last_state['last_clicked_data'] = largest_clicktime_input["click_data"]
        return json.dumps(last_state)

    return Input(last_clicked_id, 'value')





# slect rows and perform another dendogram
@app.callback(
    [Output('datatable-dendogram', "figure"),Output('dendogram-level2', "figure"),Output('dendogram-level3', "figure")],
    [Input('data-table-cluster', "data"),
     Input('data-table-cluster', "selected_rows"),last_clicked('datatable-dendogram', 'dendogram-level2','dendogram-level3')])
def plot_dendogram_subClusterLevel(data, rows,last_clickdata):

    if rows is None:

        # find cluster and sub_cluster

        cluster=data_df.loc[genes_def,'cluster'].to_list()[0]
        sub_cluster=data_df.loc[genes_def,'sub_cluster'].to_list()[0]
        target_gene=genes_def[0]  # target genes
    else:
        sub_cluster=data[rows[0]]['sub_cluster']
        cluster=data[rows[0]]['cluster']
        target_gene=data[rows[0]]['Gene ID']



     # now we are going to plot dendogram based on found genes for clustering
    # We are going to gnerate deafult figure
    # first we generate  high level clusters which was computed based on coexpression and TOM similarity
    # then we used median expression values to see association based on clusters

    found_genes=[target_gene]
    color_labels=[data_df.loc[gene,'cluster'] for gene in found_genes]
    # print('color', color_labels)
    height=1000
    orig_labels=med_df.index
    n=len(orig_labels)
    # print('n1',n)
    fig=SF.load_dendogram_cluster(med_df,orig_labels,height,color_labels)


    # print('here', data,cluster,sub_cluster)
    # get dendogram for sub cluster levels
    sub_clusters=data_df[data_df['cluster']==cluster]['sub_cluster'].unique()
    sub_inds=[ cluster+item for item in sub_clusters]
    tmp=sub_df.loc[sub_inds,:].copy()


    # fig1= ff.create_dendrogram(tmp.to_numpy(), orientation='right')
    color_labels=[sub_cluster]
    # print(rows,cluster,sub_cluster,sub_clusters,tmp,subcluster_clickData)

    # print('color', color_labels)
    height=1000
    orig_labels=sub_clusters
    fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)

    # find all sub cluster genes
    tmp=data_df[data_df['sub_cluster']==sub_cluster]

    # get rapamycin yes columns
    rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
    tmp=tmp.loc[:,rapa_yes_cols].copy()

    # get original labels

    tmp_sex=sex_df.loc[tmp.index,:].copy()
    orig_labels=sex_df.loc[tmp.index,'label']
    if n<len(orig_labels):
        n=len(orig_labels)
    # print('n2',n)
    height=1000
    target_gene_label=sex_df.loc[target_gene,'label']
    # print(target_gene_label,target_gene)
    fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)



    last_clickdata=json.loads(last_clickdata)

    click_data = last_clickdata["last_clicked_data"]
    clicked_id = last_clickdata["last_clicked"]

    # print(clicked_id,click_data)
    # ################ we are going to add naviagation through the clicking figures
    # if cluster 1 is clicked then cluster 2 and three will change
    if clicked_id =='datatable-dendogram' and 'text' in click_data['points'][0]:
        cluster_name=click_data['points'][0]['text']
        sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
        # import pdb;pdb.set_trace()
        sub_inds=[ cluster_name+item for item in sub_clusters]
        tmp=sub_df.loc[sub_inds,:].copy()
        sub_cluster=cluster_name+'_sub_0'
        color_labels=[sub_cluster]
        # print('color', color_labels)
        height=1000

        orig_labels=sub_clusters
        if n<len(orig_labels):
            n=len(orig_labels)
        # print('n3',n)
        fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)

        ### now we want to upload level 3
        # for this we need to use take first sub cluster and then plot genes

        tmp=data_df[data_df['sub_cluster']==sub_cluster]

        # get rapamycin yes columns
        rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
        tmp=tmp.loc[:,rapa_yes_cols].copy()

        # get original labels

        tmp_sex=sex_df.loc[tmp.index,:].copy()
        orig_labels=sex_df.loc[tmp.index,'label']
        if n<len(orig_labels):
            n=len(orig_labels)
        # print('nlast',n)
        # print('orig',orig_labels,sub_cluster)
        height=1000
        target_gene_label=sex_df.loc[tmp.index[0],'label']
        # print(target_gene_label,target_gene)
        fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)


    if clicked_id =='dendogram-level2' and 'text' in click_data['points'][0]:
    # if cluster 2 is clicked then cluster 3 will change
     # find all sub cluster genes

        sub_cluster=click_data['points'][0]['text']
        cluster_name=sub_cluster.split('_')[0]+'_'+sub_cluster.split('_')[1]
        tmp=data_df[data_df['sub_cluster']==sub_cluster]

        # get rapamycin yes columns
        rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
        tmp=tmp.loc[:,rapa_yes_cols].copy()

        # get original labels

        tmp_sex=sex_df.loc[tmp.index,:].copy()
        orig_labels=sex_df.loc[tmp.index,'label']
        if n<len(orig_labels):
            n=len(orig_labels)
        # print('n4',n)
        # print('orig',orig_labels,sub_cluster)
        height=1000
        target_gene_label=sex_df.loc[tmp.index[0],'label']
        # print(target_gene_label,target_gene)
        fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)



        sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
        # import pdb;pdb.set_trace()
        sub_inds=[ cluster_name+item for item in sub_clusters]
        tmp=sub_df.loc[sub_inds,:].copy()

        color_labels=[sub_cluster]
        # print('color', color_labels)
        height=1000

        orig_labels=sub_clusters
        if n<len(orig_labels):
            n=len(orig_labels)
        # print('n5',n)
        fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
    #### End of the navigatipns
    # print(len(orig_labels),n)
     # # re adjust height
    if n>70:
        height=int(n*(height/70))*1.3

    fig2.update_layout(height=height)
    fig1.update_layout(height=height)
    fig.update_layout(height=height)
    # import pdb;pdb.set_trace()
    # l1_fig.update_layout(height=height)
    # fig.show()


    return fig,fig1,fig2


@app.callback(Output('update-on-click-data', 'children'),
              [last_clicked('datatable-dendogram', 'dendogram-level2','dendogram-level3')])
def update_onclick_callback(last_clickdata):
    last_clickdata=json.loads(last_clickdata)

    click_data = last_clickdata["last_clicked_data"]
    clicked_id = last_clickdata["last_clicked"]

    # # if cluster 1 is clicked then cluster 2 and three will change
    # if clicked_id =='datatable-dendogram':
    #     cluster_name=click_data['points'][0]['text']
    #     sub_clusters=data_df[data_df['cluster']==cluster_name]['sub_cluster'].unique()
    #     # import pdb;pdb.set_trace()
    #     sub_inds=[ cluster_name+item for item in sub_clusters]
    #     tmp=sub_df.loc[sub_inds,:].copy()
    #
    #     color_labels=[sub_clusters[0]]
    #     # print('color', color_labels)
    #     height=1000
    #     orig_labels=sub_clusters
    #     fig1=SF.load_dendogram_cluster(tmp,orig_labels,height,color_labels)
    #
    # elif clicked_id =='dendogram-level2':
    # # if cluster 2 is clicked then cluster 3 will change
    #  # find all sub cluster genes
    #     sub_cluster=click_data['points'][0]['text']
    #     tmp=data_df[data_df['sub_cluster']==sub_cluster]
    #
    #     # get rapamycin yes columns
    #     rapa_yes_cols=data_df.columns[data_df.columns.str.contains('Rapa_yes')]
    #     tmp=tmp.loc[:,rapa_yes_cols].copy()
    #
    #     # get original labels
    #
    #     tmp_sex=sex_df.loc[tmp.index,:].copy()
    #     orig_labels=sex_df.loc[tmp.index,'label']
    #     # print('orig',orig_labels,sub_cluster)
    #     height=1000
    #     target_gene_label=sex_df.loc[tmp.index[0],'label']
    #     # print(target_gene_label,target_gene)
    #     fig2=SF.load_dendogram_genes(tmp,tmp_sex,orig_labels,color_scheme,target_gene_label,height)



    return "{} was last clicked and contains clickdata:\n{}".format(clicked_id, click_data)

if __name__ == '__main__':
    app.run_server(debug=True)

