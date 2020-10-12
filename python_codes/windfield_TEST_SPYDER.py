# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 16:55:43 2020

@author: silvia Bertelli
"""
############################ LOAD LIBRARIES ###################################
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt 
from geopandas import GeoDataFrame
from shapely.geometry import Point, Polygon 
import re
import geoplot
import folium
from glob import glob #loading multiple files



########################### SET REFERENCES ####################################
#print(os.getcwd()) #this is my current working directory

PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname('__file__'), '..'))
print("Current working directory is:" + PACKAGE_DIR)
sys.path.append(PACKAGE_DIR)
FILE_DIR = os.path.dirname('__file__')
#FILE_DATA_DIR = os.path.join(TEST_DIR, 'inputs/')
########################  FILE DIRECTORIES ####################################

#

## Input files

# wind track
fp = Path.cwd() / 'inputs' / 'IBTrACS' / 'George1998.csv' #IBTrACS Database
#fp = Path.cwd() / 'inputs' / 'STORM' / 'test.csv' #STORM database
path = Path.cwd() / 'inputs' / 'STORM' 
#fp_B =  Path.cwd() / 'inputs' / 'STORM' / 'STORM_DATA_IBTRACS_NA_1000_YEARS_1.txt'

# exposure - locations
fp1 = Path.cwd() / 'inputs' /'exposure' /'Barbados_cities.csv'

## Output files

#6-hourly storm track points saved as excel file (includes radius and B values) and shapefile:
wind_dataset = Path.cwd() / 'outputs' /'excel' / 'wind' / 'cod_wind_storm_test.xlsx' 
wind_dataset_shp = Path.cwd() / 'outputs' / 'shapefiles' / 'wind' /'cod_wind_storm_test.shp'

#Exposure estimation save as shapefile:
exposure_points = Path.cwd() / 'outputs' /'shapefiles' / 'exposure' /'cod_exposure_Barbados.shp'

#Exposure & distance to windtrack and estimated velocity save as excel spredsheet:
combination_dataset = Path.cwd() / 'outputs' /'excel' / 'combined'/ 'cod_Barbados&storm_test.xlsx'

#Groupby of the exposure & distance to windtrack and estimated velocity save as excel spredsheet and shapefile:
groupby_dataset = Path.cwd() / 'outputs' /'excel' / 'groupby' / 'cod_Barbados&storm_test_groupby.xlsx'
groupby_dataset_shp = Path.cwd() / 'outputs' /'shapefiles' / 'groupby'/ 'cod_Barbados&storm_test_groupby.shp'

#Windtrack map (grid) save as excel file with all possible combinations: 
windtrack_grid = Path.cwd() / 'outputs' /'excel' / 'grid' / 'cod_grid_storm_test.xlsx'

#Groupby of the windtrack map (grid) save as shapefile with the AreaPeril_ID and max Vc value:
groupby_windtrack_grid_shp = Path.cwd() / 'outputs' /'shapefiles' / 'grid'/ 'cod_grid_storm_test_groupby.shp'

############################## FUNCTIONS ######################################

## 1 Data pre-processing functions

### 1.1 Unit conversion

# Function to convert latitude and longitude to decimal degree:
def dms2dd(s):
    """convert lat and long to decimal degrees"""
    direction = s[-1]
    degrees = s[0:4]
    dd = float(degrees) 
    if direction in ('S','W'):
        dd*= -1
    return dd

# Function to convert from nautical miles to kilometers:
def nmiles_to_km(N):
    """convert nautical miles to km"""
    N = N * 1.852
    return N

#Function to convert velocity from knot to m/s (NOTE: double check this formula)  
def knot_to_msec(Velocity):
    """convert Knots (= natucal miles per hour) to m/s:
        1 Knot = 1852 meters per hour; 1 h = 3600s"""
    Velocity = Velocity * 1852 / 3600
    return Velocity

### 1.2 Distance estimation

#Function to estimate the distance between two points in decimal degree based 
#on the __Harvesine law__.
#Check Wikipedia for more information on its formulation: 
    #https://en.wikipedia.org/wiki/Haversine_formula

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1.values)
    lat1 = np.radians(lat1.values)
    lon2 = np.radians(lon2.values)
    lat2 = np.radians(lat2.values)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    dlat = np.subtract(lat2, lat1)

    a = np.add(np.power(np.sin(np.divide(dlat, 2)), 2), 
               np.multiply(np.cos(lat1), 
                           np.multiply(np.cos(lat2),
                                       np.power(np.sin(np.divide(dlon, 2)), 
                                                2))))
    c = np.multiply(2, np.arcsin(np.sqrt(a)))
    r = 6372.795477598 # approximate radius of earth in km

    return c*r #this is in km

### 1.3 Database normalization

#The following function load the data from the __Hurdat database__ and apply 
#the unit conversions:
    #-  Convert latitude and longitude in decimal degree;
    #-  convert the wind velocity from [knot] to [m/s];
    #-  Estimate the Radius of maximum wind according to Willoughby and Rahn 
        #(2006) (you can select a different function from the ones available 
        #on section 4.1)
    #-  Estimate the Holland parameter B according to Powell et al (2005);
    
def load_data_hurdat(fp):
    """ inport the excel file of the windstorm database and apply unit 
    convertions
    
    NOTE: I NEED TO DIVIDE THE TIME HEADER TO NORMALISED WITH OTHER DATASETS"""
    
    df = pd.read_excel(fp, header = None, skiprows=[0], 
                       names=("YYMMDD", "Hours", "Status", "Lat", "Lon",
                              "Vmax_kts", "Pmax_mb", 
                              "34_r_NE", "34_r_SE", "34_r_SW", "34_r_NW",
                              "50_r_NE", "50_r_SE", "50_r_SW", "50_r_NW", 
                              "64_r_NE", "64_r_SE", "64_r_SW", "64_r_NW"))
    
    df_event_clean = df_event[['YYMMDD','Hours','Lat','Lon', 'Vmax_ms', 
                               'Rmax_km', 'B', 'category']].copy()
    
    df['Lat'] =  df['Lat'].apply(dms2dd)
    df['Lon'] =  df['Lon'].apply(dms2dd)
    df['Vmax_ms'] = df['Vmax_kts'].apply(knot_to_msec)
    df['RMW_km'] = Rmax_W06(df['Vmax_ms'], df['Lat'])
    df['B'] = B_P05(df['Vmax_ms'], df['Lat'])
    
    df_clean = df[['YYMMDD','Hours','Lat','Lon', 'Vmax_ms', 'Rmax_km', 
                   'B']].copy()
    return df_clean

#The following function load the data from the IBTrACS database and 
#apply the unit conversions.

def load_data_IBTrACS(fp):
    """ inport the csv file of the windstorm database and apply unit 
    convertions"""  
    
    df = pd.read_csv(fp, sep = ',', header = 0)
    
    usecols = ['Year','Month', 'Day','Hour', 'Lat', 'Lon', 'Vmax_ms', 
               'RMW_km', 'B']
        
    #convert the column 'ISO_TIME' in data time format
    df['ISO_TIME'] = pd.to_datetime(df['ISO_TIME'], errors='coerce')
    df['Year'] = df['ISO_TIME'].dt.year
    df['Month'] = df['ISO_TIME'].dt.month
    df['Day'] = df['ISO_TIME'].dt.day
    df['Hour'] = df['ISO_TIME'].dt.hour
    
    #replace missing values with a zero and convert to float
    df['USA_RMW'] = df['USA_RMW'].fillna(0)
    df['USA_RMW'] = pd.to_numeric(df['USA_RMW'], errors='coerce')
    
    #apply convertions
    df['Lat'] = df['LAT']
    df['Lon'] = df['LON']
    df['Vmax_ms'] = df['USA_WIND'].apply(knot_to_msec)
    df['RMW_km'] = df['USA_RMW'].apply(nmiles_to_km)
    
    #estimate the B Holland parameter
    df['B'] = B_P05(df['Vmax_ms'], df['Lat'])
    
    #clean the dataframe
    df_clean = df[usecols].copy()   
    return df_clean

#The following function load the data from the __STORM database__ in csv 
#format and apply the unit conversions.

def load_data_STORM(fp):
    """ import the txt file of the STORM database, add header,
    and apply unit convertions"""  
    
    usecols = ['Year','Month',"TCnumber",'TimeStep', 'Lat', 'Lon', 
               'Vmax_ms', 'RMW_km', 'B']
    
    df = pd.read_csv(fp, sep=",", header = None, skiprows=[0],
                     names=("Year", "Month", "TCnumber", "TimeStep", "BasinID", 
                            "Lat", "Lon","Press_hPa", "Vmax_ms", "RMW_km",
                            "Category", "Landfall", "Dist_km" ))
    
    df['B'] = B_P05(df['Vmax_ms'], df['Lat'])
        
    df_clean = df[usecols].copy()   
    return df_clean

# function to load multiple files and join in a unique database

def load_multiple_STORM(path):
    """ import multiple txt files into a unique dataframe"""
    all_files = glob(os.path.join(path, "*.txt"))
    df_from_each_file = (pd.read_csv(f, sep=",", header = None, skiprows=[0],
                                     names=("Year", "Month", "TCnumber", 
                                            "TimeStep", "BasinID", 
                                            "Lat", "Lon","Press_hPa", 
                                            "Vmax_ms", "RMW_km",
                                            "Category", "Landfall", "Dist_km")
                                     ) for f in all_files)
    df   = pd.concat(df_from_each_file, ignore_index=True)
    
    usecols = ['Year','Month',"TCnumber",'TimeStep', 'Lat', 'Lon', 
               'Vmax_ms', 'RMW_km', 'B']
    df['B'] = B_P05(df['Vmax_ms'], df['Lat'])
        
    df_clean = df[usecols].copy() 
    
    return df_clean


#Load the exposure dataset. The following excel file has been created from 
#scratch and contains the following fields: ID, Latitude, Longitude, 
#Location name.

def load_data_exposure(fp1):
    """import the csv file of the exposure data points"""
    df1 = pd.read_csv(fp1, 
                      names=("id", "Lon", "Lat", "City"),
                      index_col = 'id')
    
    # convert lat and lon type to float
    df1['Lon'] = pd.to_numeric(df1['Lon'], errors='coerce')
    df1['Lat'] = pd.to_numeric(df1['Lat'], errors='coerce')
    
    return df1

### 1.4 Geolocalization

#Convert the dataframe in a geodataframe, and apply the coordinate system WGS84
def geolocalization(df):
    """function to convert the database in a geodatabase and add the coordiate
    system WGS84 (("EPSG:4326")"""
    gdf = gpd.GeoDataFrame( df, geometry=gpd.points_from_xy(
        x=df.Lon, y=df.Lat),
        crs = "EPSG:4326")
    return gdf

#Re-projcet the geodataframe from WGS 84 coordinate system to EPSG XXXX 
#coordinate system. Please specify the new coordinate system with the following
#format when calling this function: "EPSG:XXXX"

def reproject(gdf, EPSG):
    """ convert a geodataframe from WGS84 ("EPSG:4326") coordinate system to 
    another EPSG coordinate system of your choice"""
    gdf = gdf.to_crs(EPSG)
    return gdf

### 1.5 Plotting

#Function for plotting points-like features on a folium map:
def plotDot(gdf, color):
    '''input: series that contains a numeric named latitude and a numeric named longitude
    this function creates a CircleMarker and adds it to the map'''
    return folium.CircleMarker(location=[gdf.Lat, gdf.Lon], radius=2, weight=3,
                               fill=True, fill_opacity=0.1, color = color
                               ).add_to(gdf_map)

## 2 Windfield estimation functions

### 2.1 Radius of maximum wind (Rmax)

def Rmax_W04(Vmax, Lat):
    """ Estimation of the radius of maximum wind according to the formula 
    proposed by Willoughby and Rahn (2004), eq. (12.1)
    Note: possibly under-estimating Rmax when compared to IBTrACS"""
    Rmax = 46.29 * (np.exp(-0.0153*Vmax + 0.0166*Lat))
    return Rmax #this is ok if the formula is in km

def Rmax_W06(Vmax, Lat):
    """ Estimation of the radius of maximum wind according to the formula 
    proposed by Willoughby and Rahn (2006), equation 7a"""
    Rmax = 46.4 * (np.exp(-0.0155*Vmax + 0.0169*Lat))
    return Rmax #this is ok if the formula is in km

def Rmax_Q11(Vmax):
    """ Estimation of the radius of maximum wind according to the formula proposed
    by Quiring et al. (2011); Vmax and Rmax are in nautical miles. 
    Expression herein converted in km"""
    Vm= Vmax * 0.5399568
    Rmax = ((49.67 - 0.24 * Vm)) * 1.852 
    return Rmax

### 2.2 Holland shape parameter (B)
def B_P05(Vmax,Lat):
    """ Holland parameter B estimated with the statical regression formula 
    proposed by Powell et al (2005), eq.(12.2)"""
    b_shape = 0.886 + 0.0177 * Vmax - 0.0094 * Lat
    return b_shape

### 2.3 Wind speed

#The following formulas for the windfield calculation refers to the 
#Grey and Liu (2019) paper.

def wind_speed(Vmax, Rmax, r, B):
    """ cyclonic wind speed calculation according to 
    according to Grey and Liu (2019); 
    Note: the formula has been devided in x and y to ease the computation"""
    x = 1 -((Rmax / r) ** B)
    y = (Rmax / r) ** B
    Vc = Vmax * (y * np.exp(x)) ** 0.5
    return Vc

def b(r, Rmax):
    """ Inflow angle of the cyclonic wind fields direction 
    according to Grey and Liu (2019)"""
    if r < Rmax:
        b = 10 * r / Rmax
    elif r >= 1.2*Rmax:
        b = 10
    else:
        b = (75 * r / Rmax) - 65
    return b

### 2.4 Wind categorization 
def get_wind_category(row):
    """ define the Hurricane category according to the Saffir- Simpson 
    Hurrican Scale.
    Note: the scale is for a [m/s] wind velocity"""
    
    if row.V <= 33:
        return 0
    elif 33 < row.V <= 42:
        return 1
    elif 42 < row.V <= 49: 
        return 2
    elif 49 < row.V <= 58: 
        return 3
    elif 58 < row.V <= 70: 
        return 4
    else:
        return 5
    
### 2.5 Exposure grid

def exposure_grid(buffered_gdf, height, width,epsg):
    """ Function to create a grid around the windtrack for the exposure 
    estimation; it use the following nested functions to estimate the buffer
    and the grid;
    Note:
        - R is the buffer lenght (= MAX(RMW) of the gdf )
        - width and hight corresponds to the resolution of the grid
        - coordinate system adopt is WGS84
    """
    #buffered_gdf = buffer_grid(gdf, R)
    grid = get_grid(buffered_gdf, height, width)
    
    #set the coordinate system
    #grid = grid.set_crs(epsg)  
    #grid = grid.to_crs(epsg)
    
    #define lat and lon
    grid['Lon_centroid'] = grid.centroid.x
    grid['Lat_centroid'] = grid.centroid.y
    
    #set the coordinate system
    grid = grid.set_crs(epsg)
    
    return grid

def buffer_grid(gdf, R):
    """ Function to buffered the windtrack points according to the Rmax"""
    buffer_length_in_meters = R * 1000 
    cpr_gdf= gdf.to_crs(epsg=32621) # projected coordinate system
    cpr_gdf['geometry'] = cpr_gdf.geometry.buffer(buffer_length_in_meters)
    return cpr_gdf

def get_grid(buffered_gdf, height, width):
    
    """ Function to devide the bundaries area of the windtrack into a 
    equaly space grid with a defined height and width resolution"""
    #shapefile boundaries
    xmin,ymin,xmax,ymax = buffered_gdf.total_bounds
    
    #count number of rows/column 
    rows = abs(int(np.ceil((ymax-ymin) /  height)))
    cols = abs(int(np.ceil((xmax-xmin) / width)))

    #divide the area_peril in a grid
    XleftOrigin = xmin
    XrightOrigin = xmin + width
    YtopOrigin = ymax
    YbottomOrigin = ymax- height
    polygons = []
    for i in range(cols):
        Ytop = YtopOrigin
        Ybottom =YbottomOrigin
        for j in range(rows):
            polygons.append(Polygon([(XleftOrigin, Ytop), 
                                     (XrightOrigin, Ytop), 
                                     (XrightOrigin, Ybottom), 
                                     (XleftOrigin, Ybottom)])) 
            Ytop = Ytop - height
            Ybottom = Ybottom - height
        XleftOrigin = XleftOrigin + width
        XrightOrigin = XrightOrigin + width

    #create a geodataframe of the grid
    grid = gpd.GeoDataFrame({'geometry':polygons})
    #set the ID for each polygon 
    grid['AreaPeril_ID'] = grid.index
    return grid

#################################### MAIN #####################################
def main():
    # load data data
    #df_event = load_data_IBTrACS(fp)
    #df_event = load_multiple_STORM(path)
    df_exposure = load_data_exposure(fp1)
    
    #print(df_event.info())
    #print(df_event.head())
    print(df_exposure.head())
    gdf_exposure = geolocalization(df_exposure)
    
    gdf_exposure.to_file(exposure_points)
     


##############################################################################
if __name__ == "__main__":
    main()