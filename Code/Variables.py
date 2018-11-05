import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def get_idList(df):
    idList = []
    for i in range(max(df.event)):
        new_df= df.loc[df['event']==i+1]
        x = np.array(new_df.x).reshape(-1,1)
        y = np.array(new_df.z).reshape(-1,1)
        if x.shape==(0,1) or y.shape==(0,1):
            continue
        idList.append(new_df['id'[0]])
        return idList
                            

def get_nHits(df):
    "returns nHits for each event"
    nHits=[]
    for i in range(max(df.event)):
        new_df= df.loc[df['event']==i+1]
        x = np.array(new_df.x).reshape(-1,1)
        y = np.array(new_df.z).reshape(-1,1)
        if x.shape==(0,1) or y.shape==(0,1):
            continue
        nHits.append(len(new_df))
    return nHits

def RSS(df):
    '''Residual sum of squares for a fitted linear regression should be higher for showers than tracks'''
    rssArray = []
    for i in range(max(df.event)):
        new_df = df.loc[df['event']==i+1]
        x = np.array(new_df.x).reshape(-1,1)
        y = np.array(new_df.z).reshape(-1,1)
        if x.shape==(0,1) or y.shape==(0,1):
            continue
        reg = LinearRegression().fit(x, y)
        y_pred = reg.predict(x)
        rssArray.append(mean_squared_error(y, y_pred))
    return rssArray

def pc2VarExplained(df):
    '''2nd pc should explain more variance for showers then tracks'''
    pcaVar = []
    scaler = StandardScaler()
    count = 1
    #nHits = get_nHits(df)
    for i in range(max(df.event)):
        print(count)
        count= count+1
        new_df = df.loc[df['event']==i+1]
        x = np.array(new_df.x).reshape(-1,1)
        y = np.array(new_df.z).reshape(-1,1)
        if x.shape==(0,1) or y.shape==(0,1):
            continue
        scaler.fit(new_df[['x','z']])
        new_df[['x', 'z']] = scaler.transform(new_df[['x','z']])
        data = [[np.array(new_df.x)[i], np.array(new_df.z)[i]] for i in range(len(new_df.x))]
        pca = PCA(n_components=1)
        pca.fit_transform(data)
        pcaVar.append(1- pca.explained_variance_ratio_)
    return pcaVar
