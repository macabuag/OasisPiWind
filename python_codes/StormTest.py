# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 09:59:18 2020

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
import storm_database as sd

########################### SET REFERENCES ####################################
#print(os.getcwd()) #this is my current working directory

#PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname('__file__'), '..'))
#print("Current working directory is:" + PACKAGE_DIR)
#sys.path.append(PACKAGE_DIR)
#FILE_DIR = os.path.dirname('__file__')
#FILE_DATA_DIR = os.path.join(TEST_DIR, 'inputs/')
########################  FILE DIRECTORIES ####################################

# input files
admin = Path.cwd() / 'inputs'/ 'gadm36_GLP_shp' / 'gadm36_GLP_0.shp'
gar = Path.cwd() / 'inputs'/ 'lac_glp' / 'lac_glp.shp'
#ibtracs = Path.cwd() / 'inputs'/ 'IBTrACS' / 'Georges1998.csv'
#storm = Path.cwd() / 'inputs' / 'STORM' / 'storm370805.csv' #this is a single event
storm = Path.cwd() / 'inputs' / 'STORM' / 'STORM_DATA_IBTRACS_NA.csv'

intensity_dic = Path.cwd() / 'inputs'/ 'Oasis_original' / 'intensity_bin_dict_modified.csv'

vulnerability = Path.cwd() / 'inputs'/ 'vulnerability' / '37.csv'
damage_dic = Path.cwd() / 'inputs'/ 'Oasis_original' / 'damage_bin_dict.csv'

#outputs
fp_grid = Path.cwd() / 'outputs' / 'shapefiles' / 'grid.shp'
fp_centroids = Path.cwd() / 'outputs' / 'shapefiles' / 'centroids.shp'
fp_nodes = Path.cwd() / 'outputs' / 'shapefiles' / 'nodes.shp'
fp_areaperil = Path.cwd() / 'outputs' / 'csv' / 'areaperil_dict.csv'
fp_gar = Path.cwd() / 'outputs' / 'shapefiles' / 'gar.shp'
#fp_exp1 = Path.cwd() / 'shapefiles' / 'exposure.shp'
fp_exp2 = Path.cwd() / 'outputs' / 'csv' / 'SourceLocOED.csv'
fp_wind = Path.cwd() / 'outputs' / 'shapefiles' / 'wind_dataset.shp'
fp_combine_1 = Path.cwd() / 'outputs' / 'csv' / 'combined_storm_exposure.csv'
fp_groupby_1 = Path.cwd() / 'outputs' / 'csv' / 'groupby_storm_exposure.csv'
fp_groupby_2 = Path.cwd() / 'outputs' / 'csv' / 'groupby_storm_exposure_intensity.csv'
fp_footprint = Path.cwd() / 'outputs' / 'csv' / 'footprint.csv'

fp_vf = Path.cwd() / 'outputs' / 'csv' / 'vulnerability_function.csv'
fp_vf_dic = Path.cwd() / 'outputs' / 'csv' / 'vulnerability_dict.csv'
fp_groupby_vf1 = Path.cwd() / 'outputs' / 'csv' / 'combine_vulnerability_damage.csv'
fp_vulnerability = Path.cwd() / 'outputs' / 'csv' / 'vulnerability.csv'
fp_events = Path.cwd() / 'outputs' / 'csv' / 'events_p.csv'
fp_occurrence = Path.cwd() / 'outputs' / 'csv' / 'occurrence_lt.csv'

fp_test = Path.cwd() / 'outputs' / 'csv' / 'test.csv'

path_vf = 'outputs/vulnerability_functions/{}.csv'
path_vf_png = Path.cwd() / 'outputs' / 'vulnerability_functions'

path_in = 'outputs/vulnerability_functions/'
path_out_df = 'outputs/vulnerability_functions/vulnerability.csv'


############################## FUNCTIONS ######################################

### 0. Unit conversion

def get_mph(velocity):
    """
       
    Returns
    -------
    convert m/s to miles per hour [mph].

    """
    velocity = velocity * 3600 /1852
    
    return velocity



### 1. Area_peril_id

def get_admin(admin, projected_crs):
    """ 
    Inputs:
    -------
        - shapefile of admin bundaries;
    Description:
    -----------
        - load the shapefile and convert it to a geodatagrame;
        - projected the geodataframe to the selected CRS;
        - plot the geodataframe;
    Returns:
    -------
        - geodataframe 
    """
    #load the shp
    gdf_admin = gpd.read_file(admin)
    #print(gdf_admin.crs)
    gdf_admin = gdf_admin.to_crs(projected_crs)
    
    #plot the country
    gdf_admin.plot(color = 'white', edgecolor ='k', linewidth=0.5)
    plt.axis("off")
    plt.show()
    
    return gdf_admin

def get_extent(gdf_admin):
    """
    Inputs:
    -------
        - admin geodataframe;
        - projected crs (i.e. Guadaloupe is 2970; check on: https://epsg.io/2970)
    
    Description:
    -----------
        This function take the admin geodataframe and create a grid over it
    
    Returns:
    -------
        - extension bundaries of the admin geodataframe

    """
    data = gdf_admin.total_bounds
     
    # create a dataframe: lat = y and lon = x
    column_names = ['min_x', 'min_y', 'max_x', 'max_y']
    
    df_extent = pd.DataFrame(data=data, index = column_names)
    df_extent = df_extent.transpose()
    
    return df_extent

def buffer_grid(gdf_admin, radius):
    
    data = gdf_admin.total_bounds
    #print(data)
    box_data = box(*data)
    #print(box_data)
    buffer = box_data.buffer(radius)
    #print(buffer)
    bounds_extent = buffer.bounds
    #print(bounds_extent)
    return bounds_extent

def get_grid(bounds_extent, height, width, projected_crs):
    """ 
    Inputs:
    -------
        - admin geodataframe;
        - height and width of each square of the grid (grid resolution);
        - projected crs (i.e. Guadaloupe is 2970; check on: https://epsg.io/2970)
    
    Description:
    -----------
        This function take the admin geodataframe and create a grid over it; 
        save the grid as shapefile.
    
    Returns:
    -------
        - grid
    """
    #shapefile boundaries
    xmin,ymin,xmax,ymax = bounds_extent
    
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
    grid['AREA_PERIL_ID'] = grid.index + 1
    grid = grid.set_crs(projected_crs)
    grid.to_crs("EPSG:4326")
    
    grid.to_file(fp_grid)
    return grid

def get_grid_vertices(grid, projected_crs): #https://stackoverflow.com/questions/58844463/how-to-get-a-list-of-every-point-inside-a-multipolygon-using-shapely
    """ 
    Inputs:
    -------
        - grid;
        - projected crs 
    
    Description:
    -----------
        Convert the vertices of the grid in nodes and extract the correspondent
        latitude and longitude; save the nodes as shapefile
    
    Outputs:
    -------
        - geodataframe of the nodes of the grid
    """    
    col = grid.columns.tolist()
    #print(col)
    nodes = gpd.GeoDataFrame(columns=col)
    for index, row in grid.iterrows():
        for pt in list(row['geometry'].exterior.coords[:-1]): 
            nodes = nodes.append({'AREA_PERIL_ID':row['AREA_PERIL_ID'], 'geometry':Point(pt) },ignore_index=True)
    
    nodes = nodes.set_crs(projected_crs)
    nodes = nodes.to_crs("EPSG:4326")
    #print(nodes.crs)
    nodes['lon'] = nodes['geometry'].x
    nodes['lat'] = nodes['geometry'].y
    nodes.to_file(fp_nodes)
    #print(nodes.head())
    return nodes
    
def get_grid_centroid(grid, projected_crs):
    
    centroids = grid.copy()
    #centro = centro.to_crs(projected_crs)
    #print(centroids.crs)
    centroids.geometry = centroids['geometry'].centroid
    #print(centroids.info())
    centroids = centroids.to_crs("EPSG:4326")
    #print(centroids.crs)
    centroids['Centro_Lon'] = centroids.geometry.x
    centroids['Centro_Lat'] = centroids.geometry.y
    centroids.to_file(fp_centroids)
    return centroids
    
def convert_to_oasis_areaperil(nodes, peril, coverage):
    """ 
    Inputs:
    -------
        - nodes dataframe;
        - type of peril;
        - type of coverage
    
    Description:
    -----------
        Convert the nodes dataframe to the oasis "area_id" dataframe format;
        save the new dataframe as csv file;
    
    Returns:
    -------
        - dataframe "area_peril_id"
    """     
    
    column_names = ['PERIL_ID', 'COVERAGE_TYPE', 
                     'LON1', 'LAT1',
                     'LON2', 'LAT2',
                     'LON3', 'LAT3',
                     'LON4', 'LAT4',
                     'AREA_PERIL_ID'
                     ]
    
    nodes.to_crs("EPSG:4326")
    #print(nodes.crs)
    
    nodes['id'] = nodes.groupby('AREA_PERIL_ID').cumcount()
    
    point3 = nodes.loc[nodes['id'] == 0].copy()
    point3.rename(columns={'lon': 'LON1', 'lat': 'LAT1'}, inplace=True)  
    point1 = nodes.loc[nodes['id'] == 1].copy()
    point1.rename(columns={'lon': 'LON2', 'lat': 'LAT2'}, inplace=True)
    point2 = nodes.loc[nodes['id'] == 2].copy()
    point2.rename(columns={'lon': 'LON3', 'lat': 'LAT3'}, inplace=True)
    point4 = nodes.loc[nodes['id'] == 3].copy()
    point4.rename(columns={'lon': 'LON4', 'lat': 'LAT4'}, inplace=True)

    dmerge = pd.merge(point3, point1, on='AREA_PERIL_ID', how='inner')
    dmerge = pd.merge(dmerge, point2, on='AREA_PERIL_ID', how='inner')
    dmerge = pd.merge(dmerge, point4, on='AREA_PERIL_ID', how='inner')
    dmerge = dmerge.drop(['id_x', 'id_y', 'geometry_x', 'geometry_y'], axis = 1)
    
    dmerge['PERIL_ID'] = peril
    dmerge['COVERAGE_TYPE'] = coverage
    
    #rename and reorder columns
    
    df_areaperil = dmerge[column_names]
    
    #print(df_areaperil.info())
    #print(df_areaperil.head())  
    return df_areaperil

def concatenate_df_areaperil(frames_areaperil):
    df_areaperil = pd.concat(frames_areaperil)
    
    #save as csv
    df_areaperil.to_csv(fp_areaperil, index=False)
    return df_areaperil

### 2. Exposure

def get_exposure(gar):
    """ 
    Inputs:
    -------
        - exposure shapefile from the GAR database;
    
    Description:
    -----------
        Convert the shapefile to a geodataframe;
        add the WGS84 coordinate system and save again as shapefile.
    
    Returns:
    -------
        - geodataframe of exposure points
    """ 
    #load the shp
    gdf_exposure = gpd.read_file(gar)
    gdf_exposure = gdf_exposure.set_crs("EPSG:4326")
    gdf_exposure.to_file(fp_gar)
    return gdf_exposure
    
def get_location(exposure):
    """ 
    Inputs:
    -------
        - exposure geodatagrame; 
    
    Description:
    -----------
        Select the column from the dataframe that are in common with the oasis
        format; add latitude and longitude of each point in WGS84 format;
    
    Returns:
    -------
        - exposure geodataframe with latitude and longitude;
    """ 
    #print(exposure.head())
    select_columns = ['ID_5X5', 'COUNTRY', 'CPX', 'USE_SECTOR',
                      'BS_TYPE', 'RES_TYPE', 'USE_CLASS',
                      'VALFIS', 'VALHUM', 'geometry']
    
    location = exposure[select_columns].copy()
    location['Latitude'] = location['geometry'].y
    location['Longitude'] = location['geometry'].x
    return location

 

def convert_to_oasis_exposure(location, PortNumber, AccNumber, IsTenant, 
                              BuildingID, OccupancyCode, ConstructionCode,
                              LocPerilsCovered, ContentsTIV, BITIV, CondNumber):
    """ 
    Inputs:
    -------
        - location dataframe;
        - missin values for the exposure dafarame according to the oasis format;
    
    Description:
    -----------
        convert the exposure dataframe into the oasis format and 
        save it as csv file
    
    Returns:
    -------
        - exposure geodataframe 
    """ 
    #list of column from the oasis dataframe
    list_columns = [
        'PortNumber',
        'AccNumber',
        'LocNumber',
        'IsTenant', 
        'BuildingID',
        'CountryCode',
        'Latitude',
        'Longitude',
        'OccupancyCode',
        'ConstructionCode',
        'LocPerilsCovered',
        'BuildingTIV',
        'ContentsTIV',
        'BITIV',
        'CondNumber',
        'PortNumber',
        'LocCurrency']
    
    
    LOC = location.copy()
    #print(LOC.info())
    LOC.rename(columns={
        #'ID_5X5': 'LocNumber', 
        'COUNTRY': 'CountryCode',
        'VALFIS': 'BuildingTIV'
        }, inplace=True)
    
    #add the missing values to the correspondent column
    LOC['LocNumber'] = LOC.index #I set equal to the index as I had many locaton with the same value; they need to be unique!
    LOC['PortNumber'] = PortNumber
    LOC['AccNumber'] = AccNumber
    LOC['IsTenant'] = AccNumber
    LOC['BuildingID'] = BuildingID
    LOC['OccupancyCode'] = OccupancyCode
    LOC['ConstructionCode'] = ConstructionCode
    LOC['LocPerilsCovered'] = LocPerilsCovered
    LOC['ContentsTIV'] = ContentsTIV
    LOC['BITIV'] = BITIV
    LOC['CondNumber'] = CondNumber
    LOC['LocCurrency'] = ' ' #set empty value

    
    exposure = LOC[list_columns].copy()
    #exposure.to_file(fp_exp1)
    # save the dataframe as csv file
    exposure.to_csv(fp_exp2, index=False)
    #print(exposure.info())
    
    return exposure

# 3. Intensity dictionary

def load_intensity(intensity_dic):
    #load the csv
    df_intensity = pd.read_csv(intensity_dic, sep=',')#, index_col='bin_index')
    #print(df_intensity)
    return df_intensity
    
# 4. Footprint

def merge_wind_w_grid(wind_points, grid_centro):
    """ 
    Inputs:
    -------
        - windtrack dataframe;
        - grid dataframe;
    
    Description:
    -----------
        - Merge the two dataframe;
        - Estimate the distance between each windtrack point and the centre of 
        of each centroid of the grid;
        - estimate the wind velocity at the centroid [Vc in m/s];
        - convert the velocity to nautical miles per hour;
        -export the combined dataframe as csv
    Returns:
    -------
        - combined dataframe 
    """
    
    df = (wind_points.assign(key=1)
          .merge(grid_centro.assign(key=1), on="key")
          .drop("key", axis=1))
    #estimate the distance in km
    df['r_km'] = wnf.haversine(df['Lon'], df['Lat'], df['Centro_Lon'], df['Centro_Lat'])
    
    #estimate the velocity
    df['Vc_ms'] = wnf.wind_speed(df['Vmax_ms'], df['RMW_km'], df['r_km'] ,df['B'])
    
    df ['Vc_mph'] = get_mph(df['Vc_ms']) 
    #save the combine dataframe as csv file
    #df.to_csv(fp_combine_1, index=False)
    return df

def groupby_grid(df):
    df_max = df.groupby(['AREA_PERIL_ID', 'Centro_Lon', 'Centro_Lat', 'event_id' ], as_index=False).agg({
    'Vc_mph': max, 
    'r_km': min, 
    'RMW_km': max,
    'Vmax_ms':max,
    'B': 'mean'}).reset_index().copy()
    #print(df_max.info())
    #print(df_max.head())
    #df_max.to_csv(fp_groupby_1, index=False)
    return df_max

def intensity_lookup(df_intensity, df_max):
    
    #create a key called 'event_id' for merging the dataframe
    df_max = df_max.assign(key=1).merge(
        df_intensity.assign(key=1), on='key', how='outer')
    
    #merge the two database
    df_merge = df_max[(df_max['bin_from'] < df_max['Vc_mph']) 
                      &  (df_max['bin_to'] > df_max['Vc_mph'])]
    #export as csv
    #df_merge.to_csv(fp_groupby_2, index=False)
    #print(df_merge.head())
    return df_merge
    
def get_footprint(df_merge):
    
    #list of columns in the footrpint dataframe
    footprint_columns = ['event_id', 'AREA_PERIL_ID', 'bin_index', 'probability']
    
    # create the probability column and assign the value 1 #### NOTE THIS MIGHT BE MODIFIED IN THE FUTURE TO CHANGE PROBABILITIES
    
    df_merge = df_merge.assign(probability = 1)
    
    #clean the dataframe
    df_footprint = df_merge[footprint_columns].copy()
    
    #rename columns
    df_footprint = df_footprint.rename(columns = {
        'AREA_PERIL_ID':'areaperil_id',
        'bin_index':'intensity_bin_id'})
    df_footprint =  df_footprint.sort_values(by=['event_id', 'areaperil_id'])
    #save as csv
    df_footprint.to_csv(fp_footprint, index=False)
    return df_footprint
    
### Vulnerability

# Load the function from the vulnerability database

def load_vf(vulnerability):
    #load the csv
    df_vf = pd.read_csv(vulnerability, sep=',', header=0)
    #print(df_vf.info())
    return df_vf

def save_as_single_vf(df):
    #path = '/path/to/output/files/{}.csv'

    df.groupby(level=0).apply(
        lambda x: x.to_csv(path_vf.format(x.name), index=False))
    return df

def plot_multi_vf(path_in, path_out_df):
    files = glob.glob(path_in + '*.csv') 
    
    for file in files:
        #print(file)
        df1=pd.read_csv(file, header=0, sep=',')
        #print(df1)
        df_val_vf = get_values_vf(df1)
        df_val_vf.to_csv( file , index=False)
        
        #create and save figure
        y = df_val_vf.Damage
        x = df_val_vf.IM_mph
        fig, ax = plt.subplots()
        ax.plot(x,y, color = 'darkred')
        ax.grid(ls = ':')
        ax.set_xlabel('Intensity [mph]')
        ax.set_ylabel('Damage [%]')
        plt.tight_layout()
        plt.savefig(file + '.png', format="PNG")
        plt.close()
    
    combined_csv = pd.concat([pd.read_csv(f) for f in files])
    combined_csv.to_csv( path_out_df, index=False, encoding='utf-8-sig')
    return combined_csv

def get_values_vf(df_vf):
    # indicate the colum values (list)
    y = df_vf.Y_vals 
    x = df_vf.IM_c #intensity measure in m/s
    
    vulnerability_id = int(df_vf.ID_set)
    
    # convert list to a new dataframe
    df_y = pd.DataFrame([sub.split(",") for sub in y])
    df_y = df_y.astype(float)
    df_y = df_y / 100 
    df_x = pd.DataFrame([sub.split(",") for sub in x])
    df_x = df_x.astype(float)
    
    #transpose columns to row
    df_y_transposed = df_y.transpose()
    df_x_transposed = df_x.transpose()
    #print(df_x_transposed.info())
    #print(df_y_transposed.tail())
    
    #concatenate database along column
    df_val_vf = pd.concat([df_y_transposed,df_x_transposed], axis=1).copy()
    df_val_vf.columns = ['Damage', 'IM_c']
    #print(df_val_vf)
    #print(df_val_vf.info())
    #print(df_val_vf.head())
    df_val_vf.columns = ['Damage', 'IM_c']
    
    df_val_vf['IM_mph'] = get_mph(df_val_vf['IM_c'])
    df_val_vf.loc[:, 'VULNERABILITY_ID'] = vulnerability_id
    
    #print(df_val_vf.tail())
    #save as csv
    #df_val_vf.to_csv(fp_vf, index=False)
    return df_val_vf

def plot_values_vf(df_val_vf):
    y = df_val_vf.Damage
    x = df_val_vf.IM_mph
    fig, ax = plt.subplots()
    ax.plot(x,y, color = 'darkred')
    ax.grid(ls = ':')
    
    ax.set_xlabel('Intensity [mph]')
    ax.set_ylabel('Damage [%]')
    plt.tight_layout()
    
    # for saving multiple figures
    #plt.savefig(os.getcwd()+ file +'.pdf',figsize=(5,5),dpi=600)
    plt.show()
    
def vulnerability_dic(df_val_vf):# OccupancyCode):
    
    select_columns = ['ID_set', 'Name', 'Coverage', 'Hazard']
    df = df_val_vf[select_columns].copy()
    print(df.info())
    
    #set peril values as in oasis
    df.loc[(df.Hazard == 'Wind'), 'PERIL_ID'] = 'WTC'
    
    #set the vulnerability_id
    df.rename(columns = {'ID_set':'VULNERABILITY_ID'}, inplace = True)
    
    #set coverage
    # coverage = 0 : No deductible / limit
    # coverage = 1 : Building
    # coverage = 2 : Other (typically appurtenant structures)
    # coverage = 3 : Contents
    # coverage = 4 : Business Interruption (BI)
    # coverage = 5 : Property Damage (PD: Building + Other + Contents)
    # coverage = 6 : All (PD + BI)
    
    df.loc[(df.Coverage == 'Buildings'), 'COVERAGE_TYPE'] = 1
    #df.loc[(df.Coverage == 'Buildings'), 'COVERAGE_TYPE'] = 3 #this should be a different coverage in the vulnerability dataset
       
    #occupancy code
    #df.loc[(df.Name == 'C1M L EDU PRIVATE'), 'OCCUPANCYCODE'] = 1102 #change occupancy code
    df['OCCUPANCYCODE'] = 1102
    
    # columns as in oasis
    columns_oasis = ['PERIL_ID', 'COVERAGE_TYPE', 'OCCUPANCYCODE', 'VULNERABILITY_ID']
    df_vf_dic = df[columns_oasis].copy()
    #save as csv
    df_vf_dic.to_csv(fp_vf_dic, index=False)
    return df_vf_dic
    
    
# damage dictionary 
def load_damage_dic(damage_dic):
    #load the csv
    df_damage_dic = pd.read_csv(damage_dic, sep=',', header=0)
   # print(df_damage_dic.info())
    return df_damage_dic

# combine vulnerability with damage

def damage_intensity_lookup(df_val_vf, df_damage_dic, df_intensity):
    
    #create a key called 'key' for merging the dataframe
    df = df_val_vf.assign(key=1).merge(
        df_damage_dic.assign(key=1), on='key', how='outer')
    
    #print(df.info())
    
    #merge the two database
    df_merge = df[(df['bin_from'] < df['Damage']) 
                      &  (df['bin_to'] > df['Damage'])]
    
    #print(df_merge.head())
    
    # create intensity key
    df2 = df_merge.merge(df_intensity.assign(key=1), on='key', how='outer')
    
    print(df2.info())
    
    #merge the two database
    df_merge2 = df2[(df2['bin_from_y'] < df2['IM_mph']) 
                      &  (df2['bin_to_y'] > df2['IM_mph'])]
    
    #set probabilities
    #df_merge2.loc[:,'probability'] = 1
    
    #export as csv
    #df_merge2.to_csv(fp_groupby_vf1, index=False)
    
    return df_merge2

def get_vulnerability(df_merge_2):
    df = df_merge_2.copy()
    
    #set probabilities
    df.loc[:,'probability'] = 1
    
    df.rename(columns = {'VULNERABILITY_ID':'vulnerability_id',
                         'bin_index_y':'intensity_bin_id',
                         'bin_index_x':'damage_bin_id',
                         }, inplace = True)
    
    
    
    columns_oasis_vf = ['vulnerability_id','intensity_bin_id', 'damage_bin_id',
                        'probability']
    
    # dropping ALL duplicte values 
    df.sort_values("intensity_bin_id", inplace = True)
    df.drop_duplicates( subset = ["vulnerability_id","intensity_bin_id"],
                       keep = 'last', inplace = True)
    
    df_osasis_vulnerability = df[columns_oasis_vf]
    
    #export as csv
    df_osasis_vulnerability.to_csv(fp_vulnerability, index=False)
    
    return df_osasis_vulnerability
    
### Events

def unique_events(df):
    """ Create a dataframe with the unique event_id and save as csv file 
    as input files for oasis"""
    
    events = df.filter(['ID_event', 'Year', 'Month'], axis=1)
    events = events.drop_duplicates()  

    df_events = pd.DataFrame(events)
    df_events = df_events.reset_index()
    df_events['event_id'] = df_events.index + 1
    
    df_events = df_events.drop(['index'], axis=1)
    
    #export as csv
    #df_events.to_csv(fp_test, index=False)
    
    return df_events

def merge_unique_events_windtrack(df_events, windtrack):
    """ add the oasis event_id number to the windtrack dataframe"""
    
    df = pd.merge(windtrack,df_events, on="ID_event")
    
    #export as csv
    df.to_csv(fp_test, index=False)
    return df
    

def oasis_events(df_events):
    df_osasis_events = df_events.copy()
    df_osasis_events = df_osasis_events.drop(['ID_event', 'Year', 'Month'], axis=1)
    print(df_osasis_events.info())
       
    #save as csv
    df_osasis_events.to_csv(fp_events, index=False)
    
def oasis_occurence(df_events):
    """ create a dataframe with the unique event_id and correspondent occurence;
    save it as csv file for oasis lmf"""
    
    occurence = df_events.copy()
    occurence = occurence.drop(['ID_event'], axis=1)
    
   
    df_occurence = pd.DataFrame(occurence)
    
    #column names as in oeasis
    df_occurence['period_no'] = df_occurence.Year
    df_occurence['occ_year'] = df_occurence.Year
    df_occurence['occ_month'] = df_occurence.Month
    df_occurence['occ_day'] = df_occurence.Month
    
    #select columns
    columns_selected = ['event_id', 'period_no', 'occ_year', 'occ_month', 'occ_day' ]
    
    df_osasis_occurence = df_occurence[columns_selected]
    
    print(df_osasis_occurence.info())
  
    #save as csv
    df_osasis_occurence.to_csv(fp_occurrence, index=False)
    
    return df_osasis_occurence


################################ MAIN #########################################

def main():
    
    ### Area Peril id
    # load admin boundaries
    gdf_1 = get_admin(admin, projected_crs = "EPSG:2970" ) #"EPSG:2970"
    
    ## export as area peril dataframe
    #gdf_1_extents = get_extent(gdf_1)
    gdf_1_buffer = buffer_grid(gdf_1, radius = 10000)
    
    gdf_1_grid = get_grid(gdf_1_buffer, height=10000, width=10000, projected_crs = "EPSG:2970")
    
    gdf_centroids = get_grid_centroid(gdf_1_grid, projected_crs = "EPSG:2970")
    
    gdf_1_vertices = get_grid_vertices(gdf_1_grid,projected_crs = "EPSG:2970" )
    
    area_peril_df1 = convert_to_oasis_areaperil(gdf_1_vertices, peril ='WTC', coverage = 1)
    area_peril_df2 = convert_to_oasis_areaperil(gdf_1_vertices, peril ='WTC', coverage = 3)
    
    frames_areaperil = [area_peril_df1, area_peril_df2] #add here all the areaperils dataframes (we have different df as we have different coverages and perils)
    area_peril = concatenate_df_areaperil(frames_areaperil)
    
    ###  Exposure
    
    gdf_gar = get_exposure(gar)
    

    gdf_locantion = get_location(gdf_gar)
    gdf_exposure = convert_to_oasis_exposure(gdf_locantion, 
                                             PortNumber=1, 
                                             AccNumber= 'A11111',
                                             IsTenant = 1,
                                             BuildingID = 1,
                                             OccupancyCode = 1102, #change this based on the vulnerability functions, I think it should be a list
                                             ConstructionCode = 5000,
                                             LocPerilsCovered = 'WTC',
                                             ContentsTIV = 0,
                                             BITIV = 0,
                                             CondNumber = 0)
    
    ### Windfield
    
    # load the data
    df_windtrack = sd.load_STORM_analysis(storm)   
    
    list_events = unique_events(df_windtrack)
    df_oasis_events = oasis_events(list_events)
    df_oasis_occurence = oasis_occurence(list_events)
    df_wind = merge_unique_events_windtrack(list_events, df_windtrack)

    #trasform it to a geodataframe
    gdf_event = wnf.geolocalization(df_wind)

    # save the windtrack as shapefile
    gdf_event.to_file(fp_wind, driver='ESRI Shapefile')
    #print('CHECK: The windtrack  as been saved as shapefile in /outputs/shapefiles !')
    
    ### Intensity dictionary
    
    df_intensity = load_intensity(intensity_dic)
    
    ### Footprint
    
    df_merge = merge_wind_w_grid(gdf_event, gdf_centroids)
    #print(df_merge.info())
    df_groupby = groupby_grid(df_merge)
    
    df_lookup = intensity_lookup(df_intensity, df_groupby)
    
    df_footprint = get_footprint(df_lookup)
    
    ### vulnerability
    glp_vf_1 = load_vf(vulnerability)
      
    vulnerability_dictionary = vulnerability_dic(glp_vf_1)#, 1102 ) #change this accordinly to the exposure
    
    df_glp_vf_1 = get_values_vf(glp_vf_1)
    #plot = plot_values_vf(df_glp_vf_1) # awful graphic, but the plot is correct
    
    #damage dictionary
    df_damage =  load_damage_dic(damage_dic)
    #print(df_damage.info())
    
    # combine damage and vulnerability
    df_lookup_2 = damage_intensity_lookup(df_glp_vf_1, df_damage, df_intensity)

    
    vulnerability_oasis = get_vulnerability(df_lookup_2)
    
##############################################################################
if __name__ == "__main__":
    main()