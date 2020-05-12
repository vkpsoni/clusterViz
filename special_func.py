import pandas as pd
import plotly.graph_objs as go
import plotly.figure_factory as ff
import numpy as np

from plotly.subplots import make_subplots
import plotly.graph_objects as go


class SpecialFun:

  def __init__(self):
      pass

  def load_default_table(self,data_df,found):

    # We are going to load by deafult some genes for clustering and visualization of the data
    # get  data frame by PBANKA id  genes_def=['PBANKA_1342200.1','PBANKA_1242600.1']


    if len(found)>0:

        # construct table
        tmpdf=data_df.loc[found,['cluster','sub_cluster','Product Description','Gene Name or Symbol']]
        tmpdf=tmpdf.reset_index()
        tmpdf=tmpdf.rename(columns={'index':'Gene ID'})



        columns=[{"name": i, "id": i} for i in tmpdf.columns]
        data=tmpdf.to_dict('records')
        row_selectable="single"

            # do style

        style_cell_conditional=[
            {
            'if': {'column_id': c},
            'textAlign': 'left'
            } for c in ['Gene ID', 'cluster','sub_cluster','Product Description','Gene Name or Symbol']
            ]

        style_data_conditional=[
            {
            "if": {"row_index": "odd"},
                    "backgroundColor": "rgb(248, 248, 248)",
                    # 'color': 'white'
            }
            ]

        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
                'fontWeight': 'bold'
                }


    else:
        data=[{}]
        columns=[]
        style_cell_conditional=[]
        style_data_conditional=[]
        style_header={}
    return data,columns,style_cell_conditional,style_data_conditional,style_header,row_selectable

  def subplot_dendogram(self,l1_df,l2_df,l3_df,l1_label,l1_color):
    # we will pass three dataframe 1) high level clustering
    # 2) second level clustering
    # 3) third level clustering
    # l1_df, data_frame for first level clustering
    # l2_df, data frame for second level clustering
    # l3_df, dataframe for 3rd level clustering
    # this function we will use to plot dendogram for cluster and sub-clusterLevel


    fig = make_subplots(rows=1, cols=3,subplot_titles=("Level-1", "Level-2", "Level-3" ))

    if len(l1_label)>1:
        fig= ff.create_dendrogram(l1_df.to_numpy(), orientation='right', labels= l1_label)
        y=list(fig['layout']['yaxis']['tickvals'])
        new_labels=list(fig['layout']['yaxis']['ticktext'])
    else:
        fig=go.Figure()
        new_labels=list(l1_label)
        y=[1]
    x=np.full((1, len(y)), .1)[0]
    size=np.full((1, len(y)), 10)[0]
    color=np.full((1, len(y)), '#AED6F1')[0]

    ## we want to change color of sub cluster which contains that genes
    inds=[]

    if not l1_color==None:

        try:
            for label in l1_color:

                inds.append(new_labels.index(label))

        except Exception:
         pass

    color[inds]='#F1C40F'
    # print('inds',inds,new_labels,color_labels)
    fig.add_trace(go.Scatter(mode='markers',x=x,y=y,text=new_labels,hoverinfo='text',marker=dict(size=size,color=color)))


    ####

    return fig

  def load_dendogram_cluster(self,med_df,orig_labels,height=900,color_labels=None):
    # this function we will use to plot dendogram for cluster and sub-clusterLevel
    if len(orig_labels)>1:
        fig= ff.create_dendrogram(med_df.to_numpy(), orientation='right', labels= orig_labels)
        y=list(fig['layout']['yaxis']['tickvals'])
        new_labels=list(fig['layout']['yaxis']['ticktext'])
    else:
        fig=go.Figure()
        new_labels=list(orig_labels)
        y=[1]
    x=np.full((1, len(y)), .1)[0]
    size=np.full((1, len(y)), 10)[0]
    color=np.full((1, len(y)), '#AED6F1')[0]

    ## we want to change color of sub cluster which contains that genes
    inds=[]

    if not color_labels==None:

        try:
            for label in color_labels:

                inds.append(new_labels.index(label))

        except Exception:
         pass

    color[inds]='#F1C40F'
    # print('inds',inds,new_labels,color_labels)
    fig.add_trace(go.Scatter(mode='markers',x=x,y=y,text=new_labels,hoverinfo='text',marker=dict(size=size,color=color)))
    # fig.update_layout(width=600, height=900)
    fig.update_layout(height=height,clickmode='event+select')
    fig['layout']['yaxis']['side']='right'
    fig['layout']['margin']['r']=5
    fig['layout']['margin']['l']=5

    return fig

  def load_dendogram_genes(self,tmp,tmp_sex,orig_labels,color_scheme,target_gene,height=900,xshift=1.5):

    # we are going to plot sub cluster of genes
    # then we will map based on sex ratios
    ## now apply dedogram function of plotly

    if len(orig_labels)>1:

        fig2= ff.create_dendrogram(tmp.to_numpy(), orientation='right', labels=orig_labels)

        new_labels=list(fig2['layout']['yaxis']['ticktext'])

        y=list(fig2['layout']['yaxis']['tickvals'])

        # fig2['layout']['xaxis']['type']='log'

    else:
        ## when only one element found in the tree
        fig2=go.Figure()
        new_labels=list(orig_labels)
        y=[1]

    x=np.full((1, len(y)), .1)[0]
    size=np.full((1, len(y)), 8)[0]

    # ## change colors of genes
    # print(new_labels.index(target_gene))
    # color[new_labels.index(target_gene)]=1

    # set index as a label
    tmp_sex_c=tmp_sex.copy()
    tmp_sex_c.set_index('label',inplace=True)
    tmp_sex_c=tmp_sex_c.replace(color_scheme)

    # get index of the target gene
    tid=new_labels.index(target_gene)
    # we are going to plot for male
    symbol=np.full((1, len(y)),  'triangle-up')[0]
    # symbol[tid]='circle'
    # size[tid]=10
    color=tmp_sex_c.loc[new_labels,'reg_m'].to_list()
    # color=np.full((1, len(y)), 0)[0]
    fig2.add_trace(go.Scatter(mode='markers',x=x,y=y,opacity=1,name='male',text=new_labels,hoverinfo='text',marker=dict(size=size,color=color,symbol=symbol)))



    fig2.add_trace(go.Scatter(mode='markers',x=[xshift],y=[y[tid]],opacity=1,marker=dict(color=['#F1C40F'],symbol=['star-dot'],size=[12])))
    # fig = ff.create_dendrogram(X, color_threshold=1.5)
    fig2['layout']['yaxis']['side']='right'

    # fig2.update_layout(height=height,clickmode='event+select')


    # we are going to plot for female
    # x=np.full((1, len(y)), -xshift)[0]
    symbol=np.full((1, len(y)),  'circle-open')[0]
    size=np.full((1, len(y)), 13)[0]
    # symbol[tid]='circle'
    ##
    # opacity=np.full((1, len(y)),  0.5)[0]

    ##
    color=tmp_sex_c.loc[new_labels,'reg_f'].to_list()
    # color=np.full((1, len(y)), 0)[0]
    fig2.add_trace(go.Scatter(mode='markers',x=x,y=y,opacity=1,name='female',text=new_labels, hoverinfo='text',marker=dict(size=size,color=color,symbol=symbol,line=dict(width=3,color=color))))
    # fig = ff.create_dendrogram(X, color_threshold=1.5)
    # color=np.full((1, len(y)), 0)[0]
    #fig2.add_trace(go.Scatter(mode='markers',x=[x[tid]],y=[y[tid]],opacity=1,name='female',text=new_labels[tid], hoverinfo='text',marker=dict(size=size[tid],color=color[tid],symbol=symbol[tid])))

    fig2['layout']['yaxis']['side']='right'
    fig2['layout']['margin']['r']=5
    fig2['layout']['margin']['l']=5
    # # re adjust height
    # if len(y)>100:
    #     height=int(len(y)*(height/100))

    # legend=go.layout.Legend(
    #     x=0,
    #     y=1,
    #     traceorder="normal",
    #     font=dict(
    #         family="sans-serif",
    #         size=12,
    #         color="black"
    #     ),
    #     bgcolor="LightSteelBlue",
    #     bordercolor="Black",
    #     borderwidth=2
    # )
    ## update legend
    fig2.update_layout(height=height,clickmode='event+select')

    # ###  addd male female legends
    #
    # fig2.add_trace(go.Scatter(mode='markers',x=[-xshift,0.1],y=[y[-1]+15,y[-1]+15],text=['Female','Male'],hoverinfo='text',opacity=1,marker=dict(color=['#7B7D7D','#7B7D7D'],symbol=['triangle-down-open','triangle-up-open'],size=[10,10])))


    ###
    return fig2


  def load_dendogram_genes_old(self,tmp,tmp_sex,orig_labels,color_scheme,target_gene,height=900,xshift=1):

    # we are going to plot sub cluster of genes
    # then we will map based on sex ratios
    ## now apply dedogram function of plotly

    if len(orig_labels)>1:

        fig2= ff.create_dendrogram(tmp.to_numpy(), orientation='right', labels=orig_labels)

        new_labels=list(fig2['layout']['yaxis']['ticktext'])

        y=list(fig2['layout']['yaxis']['tickvals'])

        # fig2['layout']['xaxis']['type']='log'

    else:
        ## when only one element found in the tree
        fig2=go.Figure()
        new_labels=list(orig_labels)
        y=[1]

    x=np.full((1, len(y)), .1)[0]
    size=np.full((1, len(y)), 13)[0]

    # ## change colors of genes
    # print(new_labels.index(target_gene))
    # color[new_labels.index(target_gene)]=1

    # set index as a label
    tmp_sex_c=tmp_sex.copy()
    tmp_sex_c.set_index('label',inplace=True)
    tmp_sex_c=tmp_sex_c.replace(color_scheme)

    # get index of the target gene
    tid=new_labels.index(target_gene)
    # we are going to plot for male
    symbol=np.full((1, len(y)),  'triangle-up')[0]
    # symbol[tid]='circle'
    # size[tid]=10
    color=tmp_sex_c.loc[new_labels,'reg_m'].to_list()
    # color=np.full((1, len(y)), 0)[0]
    fig2.add_trace(go.Scatter(mode='markers',x=x,y=y,opacity=1,name='male',text=new_labels,hoverinfo='text',marker=dict(size=size,color=color,symbol=symbol)))



    fig2.add_trace(go.Scatter(mode='markers',x=[0.5],y=[y[tid]],opacity=1,marker=dict(color=['#F1C40F'],symbol=['star-dot'])))
    # fig = ff.create_dendrogram(X, color_threshold=1.5)
    fig2['layout']['yaxis']['side']='right'

    # fig2.update_layout(height=height,clickmode='event+select')


    # we are going to plot for female
    x=np.full((1, len(y)), -xshift)[0]
    symbol=np.full((1, len(y)),  'triangle-down')[0]
    # symbol[tid]='circle'
    ##
    # opacity=np.full((1, len(y)),  0.5)[0]

    ##
    color=tmp_sex_c.loc[new_labels,'reg_f'].to_list()
    # color=np.full((1, len(y)), 0)[0]
    fig2.add_trace(go.Scatter(mode='markers',x=x,y=y,opacity=1,name='female',text=new_labels, hoverinfo='text',marker=dict(size=size,color=color,symbol=symbol)))
    # fig = ff.create_dendrogram(X, color_threshold=1.5)
    # color=np.full((1, len(y)), 0)[0]
    #fig2.add_trace(go.Scatter(mode='markers',x=[x[tid]],y=[y[tid]],opacity=1,name='female',text=new_labels[tid], hoverinfo='text',marker=dict(size=size[tid],color=color[tid],symbol=symbol[tid])))

    fig2['layout']['yaxis']['side']='right'

    # # re adjust height
    # if len(y)>100:
    #     height=int(len(y)*(height/100))

    # legend=go.layout.Legend(
    #     x=0,
    #     y=1,
    #     traceorder="normal",
    #     font=dict(
    #         family="sans-serif",
    #         size=12,
    #         color="black"
    #     ),
    #     bgcolor="LightSteelBlue",
    #     bordercolor="Black",
    #     borderwidth=2
    # )
    ## update legend
    fig2.update_layout(height=height,clickmode='event+select')

    ###  addd male female legends

    fig2.add_trace(go.Scatter(mode='markers',x=[-xshift,0.1],y=[y[-1]+15,y[-1]+15],text=['Female','Male'],hoverinfo='text',opacity=1,marker=dict(color=['#7B7D7D','#7B7D7D'],symbol=['triangle-down-open','triangle-up-open'],size=[10,10])))


    ###
    return fig2




  def getMediansubCluster(self,data_df):

    res=data_df.copy()
    res=res.drop(columns=['cluster','sub_cluster','Product Description','Gene Name or Symbol'])
    med_df=pd.DataFrame(columns=res.columns)

    for k,vals in data_df.groupby(['cluster','sub_cluster']).indices.items():
        med_df.loc[k[0]+k[1],res.columns]=res.iloc[vals,:].median(axis=0)
    return med_df

  def getMedianCluster(self,data_df):

    res=data_df.copy()
    res=res.drop(columns=['cluster','sub_cluster','Product Description','Gene Name or Symbol'])
    med_df=pd.DataFrame(columns=res.columns)
    for k,vals in data_df.groupby(['cluster']).indices.items():
        med_df.loc[k,res.columns]=res.iloc[vals,:].median(axis=0)
    return med_df







