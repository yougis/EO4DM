# -*- coding: utf-8 -*-
"""
##############################################################################

GEOSTATS functions for processing spatial statistics of drought indicators

##############################################################################
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import features
import rasterio.mask
from rasterio.warp import reproject, Resampling, aligned_target
from zipfile import BadZipFile, ZipFile

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def extractAreasMask(ref_rast, areas_shp, areas_key=None):
    """
    Extract areas mask according to :
        - a reference raster grid (target output grid)
        - a shapefile containing geometries/areas to identify
        - a key identifying areas name in input shp (default 'nom')
    Note :  mask pixels are filled with the ObjectId values of the shapefile geometries
            A look-up table (dataframe) is generated to associate any ObjectId to the corresponding geometry (name)
    """

    if (areas_key is None) or (areas_key=='') : areas_key='nom'
    
    zip_in = 0
    if 'zip' in ref_rast:
        zip_in = 1
        zf_path = ref_rast.split('!')[0][6:]
        tf_path = ref_rast.split('!')[1]
        with ZipFile(zf_path) as zf:
            ref_rast = zf.open(tf_path)

    with rasterio.open(ref_rast, GEOREF_SOURCES='INTERNAL') as d_ds:
        mskarea_gdf = gpd.read_file(areas_shp)
        mskarea_gdf = mskarea_gdf.to_crs(d_ds.crs)
        mskarea_gdf = mskarea_gdf.sort_values(by=[areas_key]).reset_index(drop=True)

        if 'OBJECTID' not in mskarea_gdf.keys():
            mskarea_gdf['OBJECTID'] = mskarea_gdf.index.array + 1
        if areas_key not in mskarea_gdf.keys():
            logging.critical(f'Wrong areas column name in input geometry shp : missing key {areas_key}')
            raise Exception(('Wrong areas column name in input geometry'))
        
        maskarea_data = features.rasterize(
            [(mskarea_gdf['geometry'][i],int(mskarea_gdf['OBJECTID'][i])) for i in range(len(mskarea_gdf))],
            out_shape=d_ds.shape,
            transform=d_ds.profile['transform'])
        
        profile_out = d_ds.profile
        profile_out.update(dtype=rasterio.int8, nodata=0)
        if profile_out['count']>1:
            profile_out.update(count=1)
            
    area_df = mskarea_gdf.drop(columns=['geometry'])
    if areas_key != 'nom':
        area_df = area_df.rename(columns={areas_key: 'nom'})
    
    if zip_in==1: ref_rast.close()
    
    return maskarea_data, profile_out, area_df



def extractLandCoverMask(ref_rast, landcover_rast, classlist2mask):
    """
    Extract landcover mask according to :
        - a reference raster grid (target output grid)
        - a landcover raster (any .tif landcover product)
        - a classlist containing pixel values (correspond to those in landcover) to mask
    """
    
    zip_in = 0
    if 'zip' in ref_rast:
        zip_in = 1
        zf_path = ref_rast.split('!')[0][6:]
        tf_path = ref_rast.split('!')[1]
        with ZipFile(zf_path) as zf:
            ref_rast = zf.open(tf_path)

    with rasterio.open(landcover_rast) as landc_ds, \
        rasterio.open(ref_rast, GEOREF_SOURCES='INTERNAL') as d_ds:
        src_LANDC = landc_ds.read(1)
        src_profile = landc_ds.profile
        dst_profile = d_ds.profile
        dst_LANDC = np.zeros((dst_profile['height'],dst_profile['width']), src_LANDC.dtype)

    if zip_in==1: ref_rast.close()

    reproject(src_LANDC,
              dst_LANDC,
              src_transform=src_profile['transform'],
              src_crs=src_profile['crs'],
              dst_transform=dst_profile['transform'],
              dst_crs=dst_profile['crs'],
              resampling=Resampling.nearest)
    
    masklandcover_data = np.ones(dst_LANDC.shape)
    for class_val in classlist2mask:
        masklandcover_data[(dst_LANDC==class_val)] = 0
    
    profile_out = dst_profile
    profile_out.update(dtype=rasterio.int8, nodata=0)
    if profile_out['count']>1:
        profile_out.update(count=1)
    
    return masklandcover_data, profile_out



def prepareGeoStatsMasks(file_like, file_areas, outdir_masks, suffix=None, file_landcover=None, areas_key=None):
    """
    Prepare (sub)areas land masks and look-up table
        -> used to extract spatial statistics and confidence index
    """

    # --- Generate mask considering (sub)areas (any landcover) ---
    if suffix is None:
        mask_filename = 'mask_Areas.tif'
        lut_filename = 'ID_Name_Areas-lookup.csv'
    else:
        mask_filename = f'mask_Areas_{suffix}.tif'
        lut_filename = f'ID_Name_Areas-lookup_{suffix}.csv'
    
    maskAREA, profile_out, com_df = extractAreasMask(file_like, file_areas, areas_key)
    out_mask_filename = os.path.join(outdir_masks, mask_filename)
    with rasterio.open(out_mask_filename,'w',**profile_out) as out_ds:
        out_ds.write(maskAREA,1)
    del mask_filename, out_mask_filename

    out_lut_filename = os.path.join(outdir_masks, lut_filename)
    com_df.to_csv(out_lut_filename,
                  index = False,
                  float_format = '%.0f',
                  sep=';')   

    # --- Generate mask removing Dense Vegetation and Build Areas (_NoTrees_NoBuild) ---
    if (file_landcover is not None) and (file_landcover != []):
        if suffix is None: mask_filename = 'mask_Areas_NOTrees_NOBuild.tif'
        else: mask_filename = f'mask_Areas_NOTrees_NOBuild_{suffix}.tif'
        mask_NOTrees_NOBuild, profile_out = extractLandCoverMask(file_like, file_landcover, [2,7])
        mask_NOTrees_NOBuild = (maskAREA!=0) & (mask_NOTrees_NOBuild==1)
        out_mask_filename = os.path.join(outdir_masks, mask_filename)
        with rasterio.open(out_mask_filename,'w',**profile_out) as out_ds:
            out_ds.write(mask_NOTrees_NOBuild,1)
        del mask_filename, out_mask_filename

        # --- Generate mask keeping only Dense Vegetation (_Trees) ---
        if suffix is None: mask_filename = 'mask_Areas_Trees.tif'
        else: mask_filename = f'mask_Areas_Trees_{suffix}.tif'
        mask_Trees, profile_out = extractLandCoverMask(file_like, file_landcover, [0,1,3,4,5,6,7,8,9,10])
        out_mask_filename = os.path.join(outdir_masks, mask_filename)
        with rasterio.open(out_mask_filename,'w',**profile_out) as out_ds:
            out_ds.write(mask_Trees,1)
        del mask_filename, out_mask_filename



def extractGeoStats(data, qscore, date, mask, maskAREA, area_df, territory, mask_NOTrees_NOBuild=None, mask_Trees=None):
    '''
    Spatial statistics are estimated on the entire Territory and on every Sub-Areas (predefined in input masks).
    Statistics are min, max and std of indices.
    The quality score is also given, corresponding to the pourcentage of good data used per pixel.
    Three dataframes are given :
        - considering any landcover
        - removing Dense Vegetation and Build Areas (_NoTrees_NoBuild)
        - keeping only Dense Vegetation (_Trees)

    Note : (sub)-areas on which to compute drought statistics are defined in spatial masks (.tif)
           The pixels of spatial masks are filled with specific values that are associated to areas names in the input look-up table (area_df)
           Masks and look-pu table are generated from 
    '''
    
    GeoStats_df = pd.DataFrame(columns=['LOCATION','DATE','MEAN','MIN','MAX','STD','QSCORE'])
    GeoStats_df_NOTrees_NOBuild = GeoStats_df.copy()
    GeoStats_df_Trees = GeoStats_df.copy()

    # --- Statistics over TERRITORY ---
    data_mean = np.nanmean(data[mask==1])
    data_min = np.nanmin(data[mask==1])
    data_max = np.nanmax(data[mask==1])
    data_std = np.nanstd(data[mask==1])
    data_qscore = np.mean(qscore[mask==1])
    d = {'LOCATION':territory,'DATE':date,'MEAN':data_mean,'MIN':data_min,
            'MAX':data_max,'STD':data_std,'QSCORE':data_qscore}
    GeoStats_df = pd.concat([GeoStats_df, pd.DataFrame([d])], ignore_index=True)
    del data_mean,data_min,data_max,data_std,data_qscore,d
    
    if mask_NOTrees_NOBuild is not None:
        data_mean = np.nanmean(data[mask_NOTrees_NOBuild==1])
        data_min = np.nanmin(data[mask_NOTrees_NOBuild==1])
        data_max = np.nanmax(data[mask_NOTrees_NOBuild==1])
        data_std = np.nanstd(data[mask_NOTrees_NOBuild==1])
        data_qscore = np.mean(qscore[mask_NOTrees_NOBuild==1])
        d = {'LOCATION':territory,'DATE':date,'MEAN':data_mean,'MIN':data_min,
                'MAX':data_max,'STD':data_std,'QSCORE':data_qscore}
        GeoStats_df_NOTrees_NOBuild = pd.concat([GeoStats_df_NOTrees_NOBuild, pd.DataFrame([d])], ignore_index=True)
        del data_mean,data_min,data_max,data_std,data_qscore,d
    
    if mask_Trees is not None:
        data_mean = np.nanmean(data[mask_Trees==1])
        data_min = np.nanmin(data[mask_Trees==1])
        data_max = np.nanmax(data[mask_Trees==1])
        data_std = np.nanstd(data[mask_Trees==1])
        data_qscore = np.mean(qscore[mask_Trees==1])
        d = {'LOCATION':territory,'DATE':date,'MEAN':data_mean,'MIN':data_min,
                'MAX':data_max,'STD':data_std,'QSCORE':data_qscore}
        GeoStats_df_Trees = pd.concat([GeoStats_df_Trees, pd.DataFrame([d])], ignore_index=True)
        del data_mean,data_min,data_max,data_std,data_qscore,d
    
    # --- Statistics over SUB-AREAS ---
    for a in area_df['nom']:
        
        ind_a = area_df.loc[area_df['nom']==a, 'OBJECTID']
        mask_a = (maskAREA == ind_a.values[0])
        Nb_ALLDATA_a = np.sum(mask_a)
        
        # If no pixels on commune: go to next iteration
        if Nb_ALLDATA_a==0:
            continue
        
        if Nb_ALLDATA_a!=0:
            data_mean_a = np.nanmean(data[mask_a==1])
            data_min_a = np.nanmin(data[mask_a==1])
            data_max_a = np.nanmax(data[mask_a==1])
            data_std_a = np.nanstd(data[mask_a==1])
            data_qscore_a = np.mean(qscore[mask_a==1])
            d_a = {'LOCATION':a,'DATE':date,'MEAN':data_mean_a,'MIN':data_min_a,
                    'MAX':data_max_a,'STD':data_std_a,'QSCORE':data_qscore_a}
            GeoStats_df = pd.concat([GeoStats_df, pd.DataFrame([d_a])], ignore_index=True)
            del data_mean_a,data_min_a,data_max_a,data_std_a,data_qscore_a,d_a
        
        if mask_NOTrees_NOBuild is not None:
            mask_a_NOTrees_NOBuild = (mask_a==1) & (mask_NOTrees_NOBuild==1)
            Nb_ALLDATA_a_NOTrees_NOBuild = np.sum(mask_a_NOTrees_NOBuild)

            if Nb_ALLDATA_a_NOTrees_NOBuild!=0:                   
                data_mean_a = np.nanmean(data[mask_a_NOTrees_NOBuild==1])
                data_min_a = np.nanmin(data[mask_a_NOTrees_NOBuild==1])
                data_max_a = np.nanmax(data[mask_a_NOTrees_NOBuild==1])
                data_std_a = np.nanstd(data[mask_a_NOTrees_NOBuild==1])
                data_qscore_a = np.mean(qscore[mask_a_NOTrees_NOBuild==1])
                d_a = {'LOCATION':a,'DATE':date,'MEAN':data_mean_a,'MIN':data_min_a,
                        'MAX':data_max_a,'STD':data_std_a,'QSCORE':data_qscore_a}
                GeoStats_df_NOTrees_NOBuild = pd.concat([GeoStats_df_NOTrees_NOBuild, pd.DataFrame([d_a])], ignore_index=True)
                del data_mean_a,data_min_a,data_max_a,data_std_a,data_qscore_a,d_a
        
        if mask_Trees is not None:
            mask_a_Trees = (mask_a==1) & (mask_Trees==1)
            Nb_ALLDATA_a_Trees = np.sum(mask_a_Trees)

            if Nb_ALLDATA_a_Trees!=0:                   
                data_mean_a = np.nanmean(data[mask_a_Trees==1])
                data_min_a = np.nanmin(data[mask_a_Trees==1])
                data_max_a = np.nanmax(data[mask_a_Trees==1])
                data_std_a = np.nanstd(data[mask_a_Trees==1])
                data_qscore_a = np.mean(qscore[mask_a_Trees==1])
                d_a = {'LOCATION':a,'DATE':date,'MEAN':data_mean_a,'MIN':data_min_a,
                        'MAX':data_max_a,'STD':data_std_a,'QSCORE':data_qscore_a}
                GeoStats_df_Trees = pd.concat([GeoStats_df_Trees, pd.DataFrame([d_a])], ignore_index=True)
                del data_mean_a,data_min_a,data_max_a,data_std_a,data_qscore_a,d_a
        
    return GeoStats_df, GeoStats_df_NOTrees_NOBuild, GeoStats_df_Trees


