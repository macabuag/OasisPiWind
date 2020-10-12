# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 10:12:04 2020

@author: silvia Bertelli
"""

############################ LOAD LIBRARIES ###################################
import os
import sys
from pathlib import Path
import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt 
from shapely.geometry import Polygon, MultiPoint, Point, box
import pprint

from matplotlib.colorbar import Colorbar
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits.axes_grid1 import AxesGrid

#call the python code for the windstorm estimation
import windfield_TEST_SPYDER as wnf
import GuadelopeTest as gt

########################### SET REFERENCES ####################################
#print(os.getcwd()) #this is my current working directory

#PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname('__file__'), '..'))
#print("Current working directory is:" + PACKAGE_DIR)
#sys.path.append(PACKAGE_DIR)
#FILE_DIR = os.path.dirname('__file__')
#FILE_DATA_DIR = os.path.join(TEST_DIR, 'inputs/')
########################  FILE DIRECTORIES ####################################

# input files
fp = Path.cwd() / 'inputs' / 'STORM' / 'storm370805.csv' #STORM database single event 
admin = Path.cwd() / 'inputs'/ 'gadm36_GLP_shp' / 'gadm36_GLP_0.shp'

fp_storm = Path.cwd() / 'inputs' / 'STORM' / 'STORM_DATA_IBTRACS_NA_1000_YEARS_0.txt' #STORM database

# output files
fp_storm_csv = Path.cwd() / 'inputs' / 'STORM' / 'STORM_DATA_IBTRACS_NA_1000_YEARS_0.csv'
fp_storm_selection = Path.cwd() / 'inputs' / 'STORM' / 'selection.csv'
fp_storm_shp = Path.cwd() / 'outputs' / 'shapefiles' / 'storm.shp'

############################## FUNCTIONS ######################################

# load a textfile the from storm database

def load_STORM(path):
    """ import txt files, add the header, add an ID column and export as csv"""
    df = pd.read_csv(path, sep=",", header = None, skiprows=[0],
                     names=("Year", "Month", "TCnumber", "TimeStep", "BasinID", 
                            "Lat", "Lon","Press_hPa", "Vmax_ms", "RMW_km",
                            "Category", "Landfall", "Dist_km")) 
    
    #add wind category
    df['V'] = df['Vmax_ms']
    df.loc[:,'category'] = df.apply(wnf.get_wind_category,axis=1)
    
    #set ID column
    df['ID_event'] = df.Year.astype(str) + '_' +  df.Month.astype(str) + '_' +  df.TCnumber.astype(str)
    
    # clean the dataframe   
    usecols = ['ID_event','TimeStep', 'Lat', 'Lon','Vmax_ms', 'RMW_km', 'category']  
    df_clean = df[usecols].copy() 
    
    # save as csv
    df_clean.to_csv(fp_storm_csv, index=False)
    return df_clean

def centroid_admin(gdf_admin):
    gdf_admin_centroid = gdf_admin.copy()
    gdf_admin_centroid['geometry'] = gdf_admin_centroid['geometry'].centroid

    #reproject in decimal degree
    gdf_admin_centroid = wnf.reproject(gdf_admin_centroid, EPSG = 4326)  #centroid of the country in decimal degree 


    # save lat and lon in two different columns
    gdf_admin_centroid['Lat_admin'] = gdf_admin_centroid.geometry.y
    gdf_admin_centroid['Lon_admin'] = gdf_admin_centroid.geometry.x
    return gdf_admin_centroid

def get_distance(gdf_storm, gdf_centroid_admin):
    df_wind = (gdf_storm.assign(key=1)                  #
          .merge(gdf_centroid_admin.assign(key=1), on="key")
          .drop("key", axis=1))

    # estimate the distance with the harvesine formula
    df_wind['r_km'] = wnf.haversine(df_wind['Lon'], df_wind['Lat'], df_wind['Lon_admin'], df_wind['Lat_admin'])
    return df_wind

def select_events(df_wind, distance, category):
    #select rows with r_km minor than the distance
    df = df_wind.copy()
    df = df.loc[df['r_km'] < distance]
    
    # get the corresponding list of unique values in ID_events
    events_ID = df.ID_event.unique()
    #print(events_ID)
    
    #select the events
    df1 = df_wind.copy()
    df1 = df1.loc[df1['ID_event'].isin(events_ID)]
    
    #select rows with category greater than X
    df2 = df1.loc[df1['category'] >= category]
    # get the corresponding list of unique values in ID_events
    events_ID2 = df2.ID_event.unique()
       
    #clean dataframe
    columns_selection = ['ID_event', 'TimeStep', 'Lat', 'Lon', 'Vmax_ms',
                         'RMW_km', 'category', 'r_km']
    
    df3 = df1[columns_selection].copy()
    df3 = df3.loc[df3['ID_event'].isin(events_ID2)]
    
    print(df3.info())
    #save as csv
    df3.to_csv(fp_storm_selection, index=False)
    #save shapefile
    gdf = wnf.geolocalization(df3)
    gdf.to_file(fp_storm_shp)
    return gdf
  
def load_STORM_single(fp):
    df = pd.read_csv(fp, sep=",")

    #estimate B_shape factor
    df['B'] = wnf.B_P05(df['Vmax_ms'], df['Lat'])
    
    #delete r_km column
    df.drop('r_km', axis=1, inplace=True)
    print(df.info())
    return df
    
################################ MAIN #########################################

def main():
    #load storm events
    #df_storm = load_STORM(fp_storm)
    #print(df_storm.info())
    #gdf_storm = wnf.geolocalization(df_storm)
    
    #load and estimate the admn centroid
    #gdf_admin = gt.get_admin(admin, projected_crs = "EPSG:2970" )
    #gdf_admin_centroid = centroid_admin(gdf_admin)
    #print(gdf_admin_centroid)
    
    #merge the two dataframe
    #storm_dist = get_distance(gdf_storm, gdf_admin_centroid)
    #print(storm_dist.info())
    #print(storm_dist.head())
    
    #selection = select_events(storm_dist, distance = 250, category = 4) #add max distance in km
    
    # load single event after cleaning
    df = load_STORM_single(fp)
    
##############################################################################
if __name__ == "__main__":
    main()