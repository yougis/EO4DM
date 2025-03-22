# -*- coding: utf-8 -*-
"""
##############################################################################

GEE Main functions (check products, initialize, run)

##############################################################################
"""

import os
import sys
import ee
import time
import glob
import pandas as pd
from tqdm import tqdm
import dmpipeline.GEE_Processing.gee_accounts as geeauth
import dmpipeline.GEE_Processing.GEE_generic_functions as geegen
import dmpipeline.GEE_Processing.GEE_preprocessing_functions as geeprep
import dmpipeline.GEE_Processing.GEE_compositing_functions as geecomp

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def check_GEEProducts_MODIS(CONFIG):
    """
    Function to check MODIS products availability before launching drought processing :    
        - Extract collections on ROI/PERIOD (from config) and check if products are available
        - If max waiting time is exceeded (NB_WAIT_MAX) -> GO PROCESSING (last decade not controled)
        - If not, extracts products on LAST DECADE and compare number of available products to number of expected products
        - If last decade is not full :
            -> case of multiple decades : remove last decade and GO PROCESSING
            -> case of single decade :  stop and wait more data
    """

    logging.info('\n\n --- CHECKING GEE MODIS PRODUCTS AVAILABILITY ---\n')

    NB_WAIT_MAX = 15 # Waits max 15 days

    # --- GEE Authentification ---
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id

    # --- Read input CONFIG ---
    PERIOD_START = CONFIG['PERIOD_START'].split(',')[0]
    PERIOD_END = CONFIG['PERIOD_END'].split(',')[0]
    TERRITORY = CONFIG['TERRITORY'].replace('"', '')
    roi = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.stringContains('COUNTRY_NA', TERRITORY))

    # --- Extracts collections on ROI/PERIOD and CHECK IF PRODUCTS ARE AVALABLE ---
    LSTaqua_dataset = ee.ImageCollection('MODIS/061/MYD11A1').filterBounds(roi)
    LSTterra_dataset = ee.ImageCollection('MODIS/061/MOD11A1').filterBounds(roi)
    REFLECTterra_dataset = ee.ImageCollection('MODIS/061/MOD09GA').filterBounds(roi)
    
    if PERIOD_START=='' or PERIOD_START=='First product':
        ALL_dataset = LSTaqua_dataset.merge(LSTterra_dataset).merge(REFLECTterra_dataset)
        PERIOD_START = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
    
    if PERIOD_END=='':
        ALL_dataset = LSTaqua_dataset.merge(LSTterra_dataset).merge(REFLECTterra_dataset)
        PERIOD_END = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', False).first().date().advance(1, 'day').format('YYYY-MM-dd'), path2key)
            
    LSTaqua_dataset = LSTaqua_dataset.filterDate(PERIOD_START, PERIOD_END)
    LSTterra_dataset = LSTterra_dataset.filterDate(PERIOD_START, PERIOD_END)
    REFLECTterra_dataset = REFLECTterra_dataset.filterDate(PERIOD_START, PERIOD_END)

    Nb_LSTaqua = geegen.googleErrorsControl(LSTaqua_dataset.size(), path2key)
    Nb_LSTterra = geegen.googleErrorsControl(LSTterra_dataset.size(), path2key)
    Nb_LST_images = Nb_LSTaqua + Nb_LSTterra
    Nb_REFLECT_images = geegen.googleErrorsControl(REFLECTterra_dataset.size(), path2key)
    if Nb_LST_images==0 and Nb_REFLECT_images==0:
        logging.info('NO PRODUCTS -> STOP PROCESSING')
        go_modis = 0
        go_month = 0
        return go_modis, go_month, PERIOD_START, PERIOD_END, TERRITORY

    period_start = pd.to_datetime(PERIOD_START)
    period_end = pd.to_datetime(PERIOD_END) + pd.DateOffset(days=-1)        # -1 day since input period_end is exclusive
    NB_DEC = round((period_end - period_start).days/10)
    
    # --- Extracts products on LAST DECADE and COMPARE number of available products to number of expected products ---
    start_D = [1, 11, 21]
    end_D = [10, 20, 31]
    d_end = period_end.day
    m_end = period_end.month
    y_end = period_end.year
    for d in range(3):
        if d_end>=start_D[d] and d_end<=end_D[d]:
            break
    period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d]))
    period_end_d =  pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d])) + pd.DateOffset(days=10)
    if (d==2) and period_end_d.month==period_start_d.month:
        period_end_d = period_end_d + pd.DateOffset(days=1)
    
    # -> If waiting time exceeded : Do not control last decade
    date_now = pd.Timestamp.now()
    Nb_wait = date_now - period_end_d
    if Nb_wait.days > NB_WAIT_MAX:
        logging.info(f'WAITING TIME {Nb_wait.days} days greater than {NB_WAIT_MAX} days -> GO PROCESSING (last decade not controled)')
        go_modis = 1
    
    else:
        # -> Control MODIS:
        LSTaqua_D = LSTaqua_dataset.filterDate(period_start_d, period_end_d)
        LSTterra_D = LSTterra_dataset.filterDate(period_start_d, period_end_d)
        REFLECTterra_D = REFLECTterra_dataset.filterDate(period_start_d, period_end_d)
            
        Nb_LSTaqua_D = geegen.googleErrorsControl(LSTaqua_D.size(), path2key)
        Nb_LSTterra_D = geegen.googleErrorsControl(LSTterra_D.size(), path2key)
        Nb_REFLECT_D = geegen.googleErrorsControl(REFLECTterra_D.size(), path2key)
        Nb_EXPECT_images = (period_end_d - period_start_d).days

        # Case WITH expected number of products -> PROCESS FULL PERIOD
        if Nb_LSTaqua_D==Nb_EXPECT_images and Nb_LSTterra_D==Nb_EXPECT_images and Nb_REFLECT_D==Nb_EXPECT_images:
            logging.info(f'LAST DECADE D{d+1} : EXPECTED NUMBER OF PRODUCTS -> GO PROCESSING FOR FULL PERIOD')
            go_modis = 1

        else:
            # Case WITHOUT expected number of products and SINGLE DECADE to process -> STOP and WAIT
            if NB_DEC==1:
                logging.info(f'SINGLE DECADE D{d+1} : MISSING PRODUCTS -> STOP PROCESSING AND WAIT')
                go_modis= 0

            # Case WITHOUT expected number of products and MULTIPLE DECADES to process -> SETS PERIOD END TO PREVIOUS DECADE and GO PROCESSING
            elif NB_DEC>1:
                d_prev = d-1
                if d_prev>=0:
                    period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d_prev]))
                    period_end_d = period_start_d + pd.DateOffset(days=10)
                    period_end_d_inclu = period_end_d + pd.DateOffset(days=-1)
                else:
                    d_prev=2
                    period_end_d_inclu = period_end_d + pd.DateOffset(days=-(d_end+1))
                    period_end_d = period_end_d_inclu + pd.DateOffset(days=+1)
                    y_prev = period_end_d_inclu.year
                    m_prev = period_end_d_inclu.month
                    period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_prev,m_prev,start_D[d_prev]))
                PERIOD_START = period_start.strftime('%Y-%m-%d')
                PERIOD_END = period_end_d.strftime('%Y-%m-%d')
                PERIOD_END_STR = period_end_d_inclu.strftime('%Y-%m-%d')                
                d = d_prev
                go_modis = 1

                logging.info(f'LAST DECADE : MISSING PRODUCTS -> REDUCE PERIOD FROM {PERIOD_START} TO {PERIOD_END_STR} and GO PROCESSING\n')

    # --- SET MONTH PROCESSING (go_month=1) ---

    # -> If last decade is 3rd and go_modis :
    if go_modis==1 and d==2:
        go_month=1

    # -> If last decade is 2nd and more than 2 decades (and go_modis) :
    elif go_modis==1 and d==1 and NB_DEC>2:
        go_month=1
    
    # -> If last decade is 1rst and more than 1 decade (and go_modis) :
    elif go_modis==1 and d==0 and NB_DEC>1:
        go_month=1
        
    else:
        go_month=0
    
    return go_modis, go_month, PERIOD_START, PERIOD_END, TERRITORY



def check_GEEProducts_LANDSAT_S2(CONFIG):
    """
    Function to check LANDSAT and SENTINEL-2 products availability before launching drought processing :
        - Extract collections on ROI/PERIOD (from config) and check if products are available
        - If max waiting time is exceeded (NB_WAIT_MAX) -> GO PROCESSING (last decade not controled)
        - If not, extracts products on LAST DECADE and compare number of available products to number of expected products
        - If last decade is not full :
            -> case of multiple decades : remove last decade and GO PROCESSING
            -> case of single decade :  stop and wait more data
    """

    logging.info('\n\n --- CHECKING GEE LANDSAT/SENTINEL-2 PRODUCTS AVAILABILITY ---\n')
    
    NB_WAIT_MAX = 15        # Waits maximum 15 days
    NB_EXPECT_LANDSAT = 2   # Expects 2 landsats products per tile/decade
    NB_EXPECT_S2 = 5        # Expects 5 sentinel-2 products per tile/decade

    # --- GEE Authentification ---
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id
    gee_workdir = 'projects'+'/'+project_id+'/'+'assets'

    # --- Read input CONFIG ---
    PERIOD_START = CONFIG['PERIOD_START'].split(',')[0]
    PERIOD_END = CONFIG['PERIOD_END'].split(',')[0]
    TERRITORY = CONFIG['TERRITORY'].replace('"', '')

    if (CONFIG['LANDMASK_ROI'] is None) or (CONFIG['LANDMASK_ROI']==''): LANDMASK_ROI = 0
    else: LANDMASK_ROI = int(CONFIG['LANDMASK_ROI'])

    if (CONFIG['TILES_L'] is None) or (CONFIG['TILES_L']==''): TILES_L=['']
    elif type(CONFIG['TILES_L']) is list: TILES_L = CONFIG['TILES_L']
    else: TILES_L = CONFIG['TILES_L'].split(',')

    if (CONFIG['TILES_S2'] is None) or (CONFIG['TILES_S2']==''): TILES_S2=['']
    elif type(CONFIG['TILES_S2']) is list: TILES_S2 = CONFIG['TILES_S2']
    else: TILES_S2 = CONFIG['TILES_S2'].split(',')
    
    if LANDMASK_ROI==1:
        roi = ee.FeatureCollection(gee_workdir+'/'+'Annex/Landmask_Grid_ROI')
    else:
        roi = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.stringContains('COUNTRY_NA', TERRITORY))
    
    # --- Extracts collections on ROI/PERIOD/TILES and CHECK IF PRODUCTS ARE AVALABLE ---
    L7_dataset = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
    L8_dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
    L9_dataset = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
    S2_sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(roi)
    
    if PERIOD_START=='' or PERIOD_START=='First product':
        ALL_dataset = L7_dataset.merge(L8_dataset).merge(L9_dataset).merge(S2_sr)
        PERIOD_START = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
    
    if PERIOD_END=='':
        ALL_dataset = L7_dataset.merge(L8_dataset).merge(L9_dataset).merge(S2_sr)
        PERIOD_END = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', False).first().date().advance(1, 'day').format('YYYY-MM-dd'), path2key)
            
    L7_dataset = L7_dataset.filterDate(PERIOD_START, PERIOD_END)
    L8_dataset = L8_dataset.filterDate(PERIOD_START, PERIOD_END)
    L9_dataset = L9_dataset.filterDate(PERIOD_START, PERIOD_END)
    S2_sr = S2_sr.filterDate(PERIOD_START, PERIOD_END)

    # -> Control Landsat tiles (if tile filter) :
    if TILES_L!=['']:
        L7_dataset_tiles = ee.ImageCollection([])
        L8_dataset_tiles = ee.ImageCollection([])
        L9_dataset_tiles = ee.ImageCollection([])

        for tile_L in tqdm(TILES_L, desc='CHECK LANDSAT TILES (FULL PERIOD)'):
            L7_tile = L7_dataset.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
            L8_tile = L8_dataset.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
            L9_tile = L9_dataset.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
            
            Nb_L7_tile = geegen.googleErrorsControl(L7_tile.size(), path2key)
            Nb_L8_tile = geegen.googleErrorsControl(L8_tile.size(), path2key)
            Nb_L9_tile = geegen.googleErrorsControl(L9_tile.size(), path2key)

            if Nb_L7_tile==0 and Nb_L8_tile==0 and Nb_L9_tile==0:
                logging.info(f'\nMissing ALL LANDSAT collections for tile {tile_L} -> STOP PROCESSING')
                go_landsat_s2 = 0
                go_month = 0
                return go_landsat_s2, go_month, PERIOD_START, PERIOD_END, TERRITORY, TILES_L, TILES_S2
            if Nb_L7_tile==0:
                logging.warning(f'\nMissing LANDSAT-7 collection for tile {tile_L}')
            if Nb_L8_tile==0:
                logging.warning(f'\nMissing LANDSAT-8 collection for tile {tile_L}')
            if Nb_L9_tile==0:
                logging.warning(f'\nMissing LANDSAT-9 collection for tile {tile_L}')
            
            L7_dataset_tiles = L7_dataset_tiles.merge(L7_tile)
            L8_dataset_tiles = L8_dataset_tiles.merge(L8_tile)
            L9_dataset_tiles = L9_dataset_tiles.merge(L9_tile)
            del L7_tile, L8_tile, L9_tile, Nb_L7_tile, Nb_L8_tile, Nb_L9_tile
    
        L7_dataset = L7_dataset_tiles
        L8_dataset = L8_dataset_tiles
        L9_dataset = L9_dataset_tiles

    # -> Control Sentinel-2 tiles (if tile filter) :
    if TILES_S2!=['']:
        S2_sr_tiles = ee.ImageCollection([])
        for tile_S2 in tqdm(TILES_S2, desc='CHECK S2 TILES (FULL PERIOD)'):
            S2_tile = S2_sr.filter(ee.Filter.stringContains('MGRS_TILE', tile_S2))

            Nb_S2_tile = geegen.googleErrorsControl(S2_tile.size(), path2key)

            if Nb_S2_tile==0:
                logging.warning(f'\nNo SENTINEL-2 products for S2 tile {tile_S2} -> CONTROL COLLECTION PERIOD or TILE NAME(S) ?')
            else:
                S2_sr_tiles = S2_sr_tiles.merge(S2_tile)
            
            del S2_tile, Nb_S2_tile
        
        S2_sr = S2_sr_tiles
        del S2_sr_tiles
    
    # -> Control All tiles :
    Nb_L7images_alltiles = geegen.googleErrorsControl(L7_dataset.size(), path2key)
    Nb_L8images_alltiles = geegen.googleErrorsControl(L8_dataset.size(), path2key)
    Nb_L9images_alltiles = geegen.googleErrorsControl(L9_dataset.size(), path2key)
    Nb_S2images_alltiles = geegen.googleErrorsControl(S2_sr.size(), path2key)

    if (Nb_L7images_alltiles==0) and (Nb_L8images_alltiles==0 or Nb_L9images_alltiles==0):
        logging.info('MISSING LANDSAT PRODUCTS -> STOP PROCESSING')
        go_landsat_s2 = 0
        go_month = 0
        return go_landsat_s2, go_month, PERIOD_START, PERIOD_END, TERRITORY, TILES_L, TILES_S2
    
    period_start = pd.to_datetime(PERIOD_START)
    period_end = pd.to_datetime(PERIOD_END) + pd.DateOffset(days=-1)        # -1 day since input period_end is exclusive
    NB_DEC = round((period_end - period_start).days/10)

    # --- Extracts products on LAST DECADE and COMPARE number of available products to number of expected products ---
    start_D = [1, 11, 21]
    end_D = [10, 20, 31]
    d_end = period_end.day
    m_end = period_end.month
    y_end = period_end.year
    for d in range(3):
        if d_end>=start_D[d] and d_end<=end_D[d]:
            break
    period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d]))
    period_end_d =  pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d])) + pd.DateOffset(days=10)
    if (d==2) and period_end_d.month==period_start_d.month:
        period_end_d = period_end_d + pd.DateOffset(days=1)

    # -> If waiting time exceeded : Do not control last decade
    date_now = pd.Timestamp.now()
    Nb_wait = date_now - period_end_d
    if Nb_wait.days > NB_WAIT_MAX:
        logging.info(f'WAITING TIME {Nb_wait.days} days greater than {NB_WAIT_MAX} days -> GO PROCESSING (last decade not controled)')
        go_landsat_s2 = 1
    
    else:
        # -> Control Landsat tiles (recent products L8/L9) :
        L8_D = L8_dataset.filterDate(period_start_d, period_end_d)
        L9_D = L9_dataset.filterDate(period_start_d, period_end_d)
        if TILES_L==['']:
            TILES_L = geegen.listLandsatTiles(L8_dataset.merge(L9_dataset), path2key)

        d_prev_passed = 0
        for tile_L in tqdm(TILES_L, desc='CHECK LANDSAT TILES (LAST DECADE)'):
            L8_tile_D = L8_D.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
            L9_tile_D = L9_D.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))

            Nb_L8tile_D = geegen.googleErrorsControl(L8_tile_D.size(), path2key)
            Nb_L9tile_D = geegen.googleErrorsControl(L9_tile_D.size(), path2key)

            # Case WITH expected number of products -> PROCESS FULL PERIOD
            if (Nb_L8tile_D+Nb_L9tile_D)>=NB_EXPECT_LANDSAT:
                go_landsat_s2 = 1

            else:
                # Case WITHOUT expected number of products and SINGLE DECADE to process -> STOP and WAIT
                if NB_DEC==1:
                    logging.info(f'\nSINGLE DECADE D{d+1} : MISSING LANDSAT PRODUCTS FOR AT LEAST TILE {tile_L} -> STOP PROCESSING AND WAIT')
                    go_landsat_s2 = 0
                    go_month = 0
                    return go_landsat_s2, go_month, PERIOD_START, PERIOD_END, TERRITORY, TILES_L, TILES_S2
                
                # Case WITHOUT expected number of products and MULTIPLE DECADES to process -> SETS PERIOD END TO PREVIOUS DECADE and GO PROCESSING
                elif NB_DEC>1:
                    if d_prev_passed==0:
                        d_prev = d-1
                        if d_prev>=0:
                            period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d_prev]))
                            period_end_d = period_start_d + pd.DateOffset(days=10)
                            period_end_d_inclu = period_end_d + pd.DateOffset(days=-1)
                        else:
                            d_prev=2
                            period_end_d_inclu = period_end_d + pd.DateOffset(days=-(d_end+1))
                            period_end_d = period_end_d_inclu + pd.DateOffset(days=+1)
                            y_prev = period_end_d_inclu.year
                            m_prev = period_end_d_inclu.month
                            period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_prev,m_prev,start_D[d_prev]))
                        PERIOD_START = period_start.strftime('%Y-%m-%d')
                        PERIOD_END = period_end_d.strftime('%Y-%m-%d')
                        PERIOD_END_STR = period_end_d_inclu.strftime('%Y-%m-%d')
                        d_prev_passed = 1
                        d = d_prev
                    
                    logging.info(f'\nLAST DECADE : MISSING LANDSAT PRODUCTS FOR TILE {tile_L} -> REDUCE PERIOD FROM {PERIOD_START} TO {PERIOD_END_STR} and GO PROCESSING\n')
                    go_landsat_s2 = 1
    
        # -> Control Sentinel-2 tile (if s2 products) :
        if Nb_S2images_alltiles!=0:
            S2_D = S2_sr.filterDate(period_start_d, period_end_d)
            if TILES_S2==['']:
                TILES_S2 = geegen.listSentinelTiles(S2_D, path2key)

            for tile_S2 in tqdm(TILES_S2, desc='CHECK S2 TILES (LAST DECADE)'):
                S2_tile_D = S2_D.filter(ee.Filter.stringContains('MGRS_TILE', tile_S2))

                Nb_S2_tile = geegen.googleErrorsControl(S2_tile_D.size(), path2key)

                # Case WITH expected number of products -> PROCESS FULL PERIOD
                if Nb_S2_tile>=NB_EXPECT_S2:
                    go_landsat_s2 = 1

                else:
                    # Case WITHOUT expected number of products and SINGLE DECADE to process -> STOP and WAIT
                    if NB_DEC==1:
                        logging.info(f'\nSINGLE DECADE D{d+1} : MISSING S2 PRODUCTS FOR AT LEAST TILE {tile_S2} -> STOP PROCESSING AND WAIT')
                        go_landsat_s2 = 0
                        go_month = 0
                        return go_landsat_s2, go_month, PERIOD_START, PERIOD_END, TERRITORY, TILES_L, TILES_S2
                    
                    # Case WITHOUT expected number of products and MULTIPLE DECADES to process -> SETS PERIOD END TO PREVIOUS DECADE and GO PROCESSING
                    elif NB_DEC>1:
                        if d_prev_passed==0:
                            d_prev = d-1
                            if d_prev>=0:
                                period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d_prev]))
                                period_end_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_end,m_end,start_D[d_prev])) + pd.DateOffset(days=10)
                                period_end_d_inclu = period_end_d + pd.DateOffset(days=-1)
                            else:
                                d_prev=2
                                period_end_d_inclu = period_end_d + pd.DateOffset(days=-(d_end+1))
                                period_end_d = period_end_d_inclu + pd.DateOffset(days=+1)
                                y_prev = period_end_d_inclu.year
                                m_prev = period_end_d_inclu.month
                                period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_prev,m_prev,start_D[d_prev]))
                            PERIOD_START = period_start.strftime('%Y-%m-%d')
                            PERIOD_END = period_end_d.strftime('%Y-%m-%d')
                            PERIOD_END_STR = period_end_d_inclu.strftime('%Y-%m-%d')
                            d_prev_passed = 1
                            d = d_prev
                        
                        logging.info(f'\nLAST DECADE : MISSING S2 PRODUCTS FOR TILE {tile_S2} -> REDUCE PERIOD FROM {PERIOD_START} TO {PERIOD_END_STR} and GO PROCESSING\n')
                        go_landsat_s2 = 1

    # --- If 3rd decade and go_landsat_s2, SET MONTH PROCESSING (go_month=1) ---
    if go_landsat_s2==1 and d==2:
        go_month=1
    else:
        go_month=0
    
    return go_landsat_s2, go_month, PERIOD_START, PERIOD_END, TERRITORY, TILES_L, TILES_S2 



def initialize_GEEProcessing_MODIS(CONFIG, EXPORT_FOLDER):
    """
    Initiliaze GEE processing of MODIS products:
        - extract collections
        - prepare assets folders
        - prepare global parameters
    """
    
    logging.info('\n\n --- INITIALIZING GEE MODIS ---\n')

    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']
    TERRITORY = CONFIG['TERRITORY']
    TERRITORY_str = TERRITORY.replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_GLOBAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']
    if (CONFIG['ASSET_EXPORT_MOD'] is None) or (CONFIG['ASSET_EXPORT_MOD']==''): ASSET_EXPORT_MOD = 0
    else: ASSET_EXPORT_MOD = int(CONFIG['ASSET_EXPORT_MOD'])
    if (CONFIG['CLEAN_GEEFOLDER'] is None) or (CONFIG['CLEAN_GEEFOLDER']==''): CLEAN_GEEFOLDER = 0
    else: CLEAN_GEEFOLDER = int(CONFIG['CLEAN_GEEFOLDER'])
    if (CONFIG['CLEAN_GEECOL'] is None) or (CONFIG['CLEAN_GEECOL']==''): CLEAN_GEECOL = 0
    else: CLEAN_GEECOL = int(CONFIG['CLEAN_GEECOL'])

    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id
    gee_workdir = 'projects'+'/'+project_id+'/'+'assets'

    # --- Load ESA Land Cover collection and extract permanent LAND MASK ---
    esa_landc_2020 = ee.ImageCollection("ESA/WorldCover/v100").first()
    esa_landc_2021 = ee.ImageCollection("ESA/WorldCover/v200").first()
    landmask_2020 = esa_landc_2020.neq(80)
    landmask_2021 = esa_landc_2021.neq(80)
    landmask = landmask_2020.And(landmask_2021)

    # --- Load and Filter MODIS Collections to the specific AREA (territory) and PERIOD (if needed) ---
    roi = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.stringContains('COUNTRY_NA', TERRITORY))
    LSTaqua_dataset = (ee.ImageCollection('MODIS/061/MYD11A1').filterBounds(roi))
    LSTterra_dataset = (ee.ImageCollection('MODIS/061/MOD11A1').filterBounds(roi))
    REFLECTterra_dataset = ee.ImageCollection('MODIS/061/MOD09GA').filterBounds(roi)

    if PERIOD_START!=[''] and PERIOD_END!=['']:
        LSTaqua_dataset = LSTaqua_dataset.filterDate(PERIOD_START, PERIOD_END)
        LSTterra_dataset = LSTterra_dataset.filterDate(PERIOD_START, PERIOD_END)
        REFLECTterra_dataset = REFLECTterra_dataset.filterDate(PERIOD_START, PERIOD_END)
        ALL_dataset = LSTaqua_dataset.merge(LSTterra_dataset).merge(REFLECTterra_dataset)
    else:
        ALL_dataset = LSTaqua_dataset.merge(LSTterra_dataset).merge(REFLECTterra_dataset)
        PERIOD_START = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        PERIOD_END = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', False).first().date().advance(1, 'day').format('YYYY-MM-dd'), path2key)

    Nb_LSTaqua = geegen.googleErrorsControl(LSTaqua_dataset.size(), path2key)
    Nb_LSTterra = geegen.googleErrorsControl(LSTterra_dataset.size(), path2key)
    Nb_LST_images = Nb_LSTaqua + Nb_LSTterra
    Nb_REFLECT_images = geegen.googleErrorsControl(REFLECTterra_dataset.size(), path2key)

    if Nb_LST_images==0 and Nb_REFLECT_images==0 :
        logging.info('\nNO MODIS PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
        sys.exit()
    if Nb_LST_images==0:
        logging.info(f'NO MODIS LST PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
        # answer = input(' -> Do you want to continue ? (y/n)')
        # if answer.lower() in "n": sys.exit()
    if Nb_REFLECT_images==0:
        logging.info(f'NO MODIS REFLECT PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
        # answer = input(' -> Do you want to continue ? (y/n)')
        # if answer.lower() in "n": sys.exit()

    # --- Prepare GEE Assets Folder and Collections ---
    date_start_str = geegen.googleErrorsControl(ee.Date(PERIOD_START).format('YYYYMMdd'), path2key)
    date_end_str = geegen.googleErrorsControl(ee.Date(PERIOD_END).advance(-1, 'day').format('YYYYMMdd'), path2key)   
    if ASSET_EXPORT_MOD==1:
        gee_folder = geegen.createAssetsFolder(f'PREPROC_GLOBAL_INDICES_{TERRITORY_str}', gee_workdir, CLEAN_GEEFOLDER)
        new_collection = geegen.createAssetsCollection(f'MODIS', gee_folder, CLEAN_GEECOL)
      
    else:
        new_collection=[]

    # --- Prepare table dicts and counts the number of products per collection
    Count_dict = []
    Count_filename = f'COUNTCollectionGEE_MODIS_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    Count_columns = ['PRODUCT','DATE START','DATE END','NUMBER OF IMAGES']
    QA_dict = []
    QAtable_filename = f'GEEPREPROC_MODIS_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    QAtable_columns = ['DATE','NAN SCORE LST','NAN SCORE NDWI','PREPROC TIME (sec)']
    COMP_dict = []
    COMPtable_filename = f'COUNTCompositeDecade_MODIS_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    COMPtable_columns = ['DATE','COMPOSITE','AquaTerra LST','Terra REFLECT','COMPOSITE NAN SCORE LST','COMPOSITE NAN SCORE NDWI','COMPOSITE TIME (sec)']

    if Nb_LSTaqua!=0:
        date_start_LSTaqua = geegen.googleErrorsControl(LSTaqua_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        date_end_LSTaqua = geegen.googleErrorsControl(LSTaqua_dataset.limit(1, 'system:time_start', False).first().date().format('YYYY-MM-dd'), path2key)
    else:
        date_start_LSTaqua=''
        date_end_LSTaqua=''
    Count_dict += [ee.Feature(None, {'PRODUCT': 'LST AQUA',
                                     'DATE START': date_start_LSTaqua,
                                     'DATE END': date_end_LSTaqua,
                                     'NUMBER OF IMAGES': Nb_LSTaqua})]
    if Nb_LSTterra!=0:
        date_start_LSTterra = geegen.googleErrorsControl(LSTterra_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        date_end_LSTterra = geegen.googleErrorsControl(LSTterra_dataset.limit(1, 'system:time_start', False).first().date().format('YYYY-MM-dd'), path2key)
    else:
        date_start_LSTterra=''
        date_end_LSTterra=''
    Count_dict += [ee.Feature(None, {'PRODUCT': 'LST TERRA',
                                     'DATE START': date_start_LSTterra,
                                     'DATE END': date_end_LSTterra,
                                     'NUMBER OF IMAGES': Nb_LSTterra})]
    if Nb_REFLECT_images!=0:
        date_start_REFLECT = geegen.googleErrorsControl(REFLECTterra_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        date_end_REFLECT = geegen.googleErrorsControl(REFLECTterra_dataset.limit(1, 'system:time_start', False).first().date().format('YYYY-MM-dd'), path2key)
    else:
        date_start_REFLECT=''
        date_end_REFLECT=''
    Count_dict += [ee.Feature(None, {'PRODUCT': 'REFLECT TERRA',
                                    'DATE START': date_start_REFLECT,
                                    'DATE END': date_end_REFLECT,
                                    'NUMBER OF IMAGES': Nb_REFLECT_images})]

    logging.info(f'\n - TERRITORY : {TERRITORY_str}\n - PERIOD : {date_start_str} -> {date_end_str}\n - NB LST AQUA : {Nb_LSTaqua}\n - NB LST TERRA : {Nb_LSTterra}\n - NB REFLECT TERRA : {Nb_REFLECT_images}')
    
    # --- Export Count Info into csv data frame ---
    geegen.exportDataFrame(DRIVE_FOLDER, Count_dict, Count_filename, Count_columns, EXPORT_FOLDER, path2key) 

    # --- Prepare outputs parameters ---
    CRS_OUT = 'EPSG:4326'
    GRID_OUT = ee.Geometry.BBox(float(CONFIG['lon_min_modis']), float(CONFIG['lat_min_modis']),
                                float(CONFIG['lon_max_modis']), float(CONFIG['lat_max_modis']))
    SCALE_REFLECT = REFLECTterra_dataset.first().select('sur_refl_b02').projection().nominalScale()
    SCALE_OUT = SCALE_REFLECT
    date_start = ALL_dataset.limit(1, 'system:time_start', True).first().date()
    date_end = ALL_dataset.limit(1, 'system:time_start', False).first().date()
    
    COLLECTIONS = (LSTaqua_dataset, LSTterra_dataset, REFLECTterra_dataset, new_collection, landmask)
    PARAM = (GRID_OUT, SCALE_OUT, CRS_OUT, date_start, date_end)
    DICT = (QA_dict, QAtable_filename, QAtable_columns, COMP_dict, COMPtable_filename, COMPtable_columns)

    return COLLECTIONS, PARAM, DICT



def initialize_GEEProcessing_LANDSAT_S2(CONFIG, EXPORT_FOLDER):
    """
    Initiliaze GEE processing of LANDSAT and SENTINEL-2 products:
        - extract collections
        - prepare assets folders
        - prepare global parameters
    """
    
    logging.info('\n\n --- INITIALIZING GEE LANDSAT/SENTINEL-2 ---\n')

    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']
    TILES_L = CONFIG['TILES_L']
    TILES_S2 = CONFIG['TILES_S2']
    TERRITORY = CONFIG['TERRITORY']
    TERRITORY_str = TERRITORY.replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_LOCAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']
    if (CONFIG['LANDMASK_ROI'] is None) or (CONFIG['LANDMASK_ROI']==''): LANDMASK_ROI = 0
    else: LANDMASK_ROI = int(CONFIG['LANDMASK_ROI'])
    if (CONFIG['ASSET_EXPORT_L'] is None) or (CONFIG['ASSET_EXPORT_L']==''): ASSET_EXPORT_L = 0
    else: ASSET_EXPORT_L = int(CONFIG['ASSET_EXPORT_L'])
    if (CONFIG['ASSET_EXPORT_S2'] is None) or (CONFIG['ASSET_EXPORT_S2']==''): ASSET_EXPORT_S2 = 0
    else: ASSET_EXPORT_S2 = int(CONFIG['ASSET_EXPORT_S2'])
    if (CONFIG['CLEAN_GEEFOLDER'] is None) or (CONFIG['CLEAN_GEEFOLDER']==''): CLEAN_GEEFOLDER = 0
    else: CLEAN_GEEFOLDER = int(CONFIG['CLEAN_GEEFOLDER'])
    if (CONFIG['CLEAN_GEECOL'] is None) or (CONFIG['CLEAN_GEECOL']==''): CLEAN_GEECOL = 0
    else: CLEAN_GEECOL = int(CONFIG['CLEAN_GEECOL'])

    # --- GEE Authentification ---
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id
    gee_workdir = 'projects'+'/'+project_id+'/'+'assets'

    # --- Load ESA Land Cover collection and extract permanent LAND MASK (water/land) ---
    esa_landc_2020 = ee.ImageCollection("ESA/WorldCover/v100").first()
    esa_landc_2021 = ee.ImageCollection("ESA/WorldCover/v200").first()
    landmask_2020 = esa_landc_2020.neq(80)
    landmask_2021 = esa_landc_2021.neq(80)
    landmask = landmask_2020.And(landmask_2021)

    # --- IF LANDMASK_ROI = 1, load ROI and update Land Mask ---
    if LANDMASK_ROI==1:
        roi = ee.FeatureCollection(gee_workdir+'/'+'Annex/Landmask_Grid_ROI')
        geegen.googleErrorsControl(roi, path2key)
        landmask = landmask.clip(roi)
    else:
        roi = ee.FeatureCollection("USDOS/LSIB/2017").filter(ee.Filter.stringContains('COUNTRY_NA', TERRITORY))

    # --- Load and Filter LANDSAT/S2 Collections to the specific AREA (territory or landmask roi) and PERIOD (if needed) ---
    L7_fulldataset = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(roi)
    L8_fulldataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(roi)
    L9_fulldataset = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(roi)
    S2_fullsr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(roi)

    if PERIOD_START!=[''] and PERIOD_END!=['']:
        L7_dataset = L7_fulldataset.filterDate(PERIOD_START, PERIOD_END)
        L8_dataset = L8_fulldataset.filterDate(PERIOD_START, PERIOD_END)
        L9_dataset = L9_fulldataset.filterDate(PERIOD_START, PERIOD_END)
        S2_sr = S2_fullsr.filterDate(PERIOD_START, PERIOD_END)
    else:
        L7_dataset = L7_fulldataset.copy()
        L8_dataset = L8_fulldataset.copy()
        L9_dataset = L9_fulldataset.copy()
        S2_sr = S2_fullsr.copy()
        ALL_dataset = L7_dataset.merge(L8_dataset).merge(L9_dataset).merge(S2_sr)
        PERIOD_START = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        PERIOD_END = geegen.googleErrorsControl(ALL_dataset.limit(1, 'system:time_start', False).first().date().advance(1, 'day').format('YYYY-MM-dd'), path2key)

    Nb_L7images_alltiles = geegen.googleErrorsControl(L7_dataset.size(), path2key)
    Nb_L8images_alltiles = geegen.googleErrorsControl(L8_dataset.size(), path2key)
    Nb_L9images_alltiles = geegen.googleErrorsControl(L9_dataset.size(), path2key)
    Nb_S2images_alltiles = geegen.googleErrorsControl(S2_sr.size(), path2key)

    if Nb_L7images_alltiles==0 and Nb_L8images_alltiles==0 and Nb_L9images_alltiles==0 and Nb_S2images_alltiles==0:
        logging.info(f'\nNO PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
        sys.exit()
    if Nb_L7images_alltiles==0:
        logging.info(f'NO L7 PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
    if Nb_L8images_alltiles==0:
        logging.info(f'NO L8 PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
    if Nb_L9images_alltiles==0:
        logging.info(f'NO L9 PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
    if Nb_S2images_alltiles==0:
        logging.info('NO S2 PRODUCTS FOUND ON THE SELECTED AREA/PERIOD\n')
    
    # --- Load and Filter the S2_CLOUD_PROBABILITY Collection, and JOIN to the S2_SR Collection ---
    S2_cloudless= (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
                    .filterBounds(roi)
                    .filterDate(PERIOD_START, PERIOD_END))
    S2_dataset = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': S2_sr,
        'secondary': S2_cloudless,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'})
        }))

    # --- Prepare GEE Assets Folder ---
    date_start_str = geegen.googleErrorsControl(ee.Date(PERIOD_START).format('YYYYMMdd'), path2key)
    date_end_str = geegen.googleErrorsControl(ee.Date(PERIOD_END).advance(-1, 'day').format('YYYYMMdd'), path2key)   
    if ASSET_EXPORT_L==1 or ASSET_EXPORT_S2==1:
        gee_folder = geegen.createAssetsFolder(f'PREPROC_LOCAL_INDICES_{TERRITORY_str}', gee_workdir, CLEAN_GEEFOLDER)
    
    if TILES_L==['']:
        TILES_L = geegen.listLandsatTiles(L7_dataset.merge(L8_dataset).merge(L9_dataset), path2key)

    logging.info(f'\n - TERRITORY : {TERRITORY_str}\n - PERIOD : {date_start_str} -> {date_end_str}\n - TILES LANDSAT : {[t for t in TILES_L]}')
    
    # --- Prepare table dicts and counts the number of products per collection
    Count_dict = []
    Count_filename = f'COUNTCollectionGEE_LANDSAT_SENTINEL2_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    Count_columns = ['TILE','PRODUCT','DATE START','DATE END','NUMBER OF IMAGES']
    QA_dict_L = []
    QAtable_filename_L = f'GEEPREPROC_LANDSAT_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    QAtable_columns_L = ['TILE','FILE NAME','DATE','CLOUD LAND','IMAGE QUALITY','SLC MODE',
                         'NAN SCORE NDWI','NAN SCORE NDVI','PREPROC TIME (sec)']
    QA_dict_S2 = []
    QAtable_filename_S2 = f'GEEPREPROC_SENTINEL2_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    QAtable_columns_S2 = ['TILE','TILE S2','FILE NAME','DATE','CLOUD','PREPROC TIME (sec)']
    COMP_dict = []
    COMPtable_filename = f'COUNTComposite_LANDSAT_SENTINEL2_{TERRITORY_str}_{date_start_str}_{date_end_str}'
    COMPtable_columns = ['TILE','DATE','COMPOSITE','L7','L8','L9','S2','TOTAL','COMPOSITE NAN SCORE NDWI','COMPOSITE NAN SCORE NDVI','COMPOSITE TIME (sec)']

    # --- Extract input collections by tile and Prepare GEE ouput collections ---
    LANDSAT_COLLECTIONS = {}
    S2_COLLECTIONS = {}
    NEW_LANDSATCOLLECTIONS = {}
    NEW_S2COLLECTIONS = {}
    LANDSAT_GRIDS = {}
    LANDSAT_PROJ = {}
    DATES = {}

    for tile_L in TILES_L:

        LANDSAT_GRIDS[tile_L] = ee.FeatureCollection(gee_workdir+'/'+'Annex/Landsat_Grid_World').filter(ee.Filter.eq('PR',int(tile_L))).geometry()

        
        (L7_tile, L8_tile, L9_tile, S2_alltiles, TILES_S2_L,
         proj_tile_landsat, date_start, date_end, Count_dict) = geegen.extractCollectionsTile(
                                                                            L7_dataset, L8_dataset, L9_dataset, S2_dataset,
                                                                            tile_L, TILES_S2, LANDSAT_GRIDS[tile_L], Count_dict,
                                                                            L7_fulldataset, L8_fulldataset, L9_fulldataset, path2key)

        LANDSAT_COLLECTIONS[tile_L] = (L7_tile, L8_tile, L9_tile)
        S2_COLLECTIONS[tile_L] = S2_alltiles
        LANDSAT_PROJ[tile_L] = proj_tile_landsat
        DATES[tile_L] = (date_start, date_end)

        geegen.exportDataFrame(DRIVE_FOLDER, Count_dict, Count_filename, Count_columns, EXPORT_FOLDER, path2key) 

        if ASSET_EXPORT_L==1:
            NEW_LANDSATCOLLECTIONS[tile_L] = geegen.createAssetsCollection(f'LANDSAT_{tile_L}', gee_folder, CLEAN_GEECOL)
        else: 
            NEW_LANDSATCOLLECTIONS[tile_L]=[]
        
        if ASSET_EXPORT_S2==1:
            new_s2collections = {}
            for tile_S2 in TILES_S2_L:
                new_s2collections[tile_S2] = geegen.createAssetsCollection(f'SENTINEL2_{tile_S2}', gee_folder, CLEAN_GEECOL)
            NEW_S2COLLECTIONS[tile_L] = new_s2collections
        else:
            NEW_S2COLLECTIONS[tile_L]=[]

        del L7_tile, L8_tile, L9_tile, S2_alltiles, proj_tile_landsat, date_start, date_end, TILES_S2_L

    # --- Prepare outputs parameters ---
    SCALE_OUT = 10

    COLLECTIONS = (LANDSAT_COLLECTIONS, S2_COLLECTIONS,
                   NEW_LANDSATCOLLECTIONS, NEW_S2COLLECTIONS,
                   LANDSAT_GRIDS, landmask)
    PARAM = (SCALE_OUT, LANDSAT_PROJ, DATES)
    DICT = (QA_dict_L, QAtable_filename_L, QAtable_columns_L,
            QA_dict_S2, QAtable_filename_S2, QAtable_columns_S2,
            COMP_dict, COMPtable_filename, COMPtable_columns)

    return COLLECTIONS, PARAM, DICT



def run_GEEProcessing_MODIS(CONFIG, OUTDIR_PATHS, COLLECTIONS, PARAM, DICT, go_month):
    """
    Run GEE processing of MODIS products:
        - compute indices
        - temporal compositing
        - export composites (decades, month)
    """
    
    logging.info('\n\n --- START GEE PROCESSING MODIS ---\n')
    
    # --- GEE Authentification ---
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id
    
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_GLOBAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']
    if (CONFIG['ASSET_EXPORT_MOD'] is None) or (CONFIG['ASSET_EXPORT_MOD']==''): ASSET_EXPORT_MOD = 0
    else: ASSET_EXPORT_MOD = int(CONFIG['ASSET_EXPORT_MOD'])
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'MODIS')
    TRAN_OUT = CONFIG['TRAN_OUT']

    (LSTaqua_dataset, LSTterra_dataset, REFLECTterra_dataset, new_collection, landmask) = COLLECTIONS
    (GRID_OUT, SCALE_OUT, CRS_OUT, date_start, date_end) = PARAM
    (QA_dict, QAtable_filename, QAtable_columns, COMP_dict, COMPtable_filename, COMPtable_columns) = DICT
    

    # ========================================== LOOP OVER YEARS/MONTHS =================================

    years = geegen.googleErrorsControl(ee.List.sequence(date_start.get('year'), date_end.get('year')), path2key)
    months = geegen.googleErrorsControl(ee.List.sequence(1, 12), path2key)

    for y in years:
        for m in months:
        
            month2find = '{}{:02d}'.format(y,m)
            logging.info(f'MONTH TO FIND : {month2find}')
            
            LSTaqua_month = (LSTaqua_dataset
                             .filter(ee.Filter.calendarRange(y,y,'year'))
                             .filter(ee.Filter.calendarRange(m,m,'month')))
            
            LSTterra_month = (LSTterra_dataset
                              .filter(ee.Filter.calendarRange(y,y,'year'))
                              .filter(ee.Filter.calendarRange(m,m,'month')))
            REFLECTterra_month = (REFLECTterra_dataset
                                  .filter(ee.Filter.calendarRange(y,y,'year'))
                                  .filter(ee.Filter.calendarRange(m,m,'month')))
            NbLSTaqua_month = geegen.googleErrorsControl(LSTaqua_month.size(), path2key)
            NbLSTterra_month = geegen.googleErrorsControl(LSTterra_month.size(), path2key)
            NbLST_month = NbLSTaqua_month + NbLSTterra_month
            NbLREFLECT_month = geegen.googleErrorsControl(REFLECTterra_month.size(), path2key)
            
            # --- NO PRODUCTS FOUND ---
            if NbLST_month==0 and NbLREFLECT_month==0:
                logging.info('NO PRODUCTS')
                continue

            # --- COMPOSITING DECADES/MONTH ---
            (LSTcompD_col, NDWIcompD_col, LSTcompM_col, NDWIcompM_col, QA_dict, COMP_TYPES,
             COMP_NANSCORES, NBLST_DATES, NBREFLECT_DATES) = geecomp.processComposite_MODIS(LSTaqua_month, LSTterra_month, REFLECTterra_month, 
                                                                                            landmask, new_collection, QA_dict, CRS_OUT, GRID_OUT,
                                                                                            SCALE_OUT, go_month, path2key, ASSET_EXPORT_MOD)

            # --- EXPORTING DECADES ---
            N_dec = geegen.googleErrorsControl(LSTcompD_col.size(), path2key)
            LSTcompD_list = LSTcompD_col.toList(N_dec)
            NDWIcompD_list = NDWIcompD_col.toList(N_dec)

            for d in range(N_dec):

                LST_comp = ee.Image(LSTcompD_list.get(d))
                NDWI_comp = ee.Image(NDWIcompD_list.get(d))
                
                comp_lst_filename = f'MODIS_LST_{month2find}_{COMP_TYPES[d]}'
                comp_ndwi_filename = f'MODIS_NDWI_{month2find}_{COMP_TYPES[d]}'

                start_time = time.time()
                geegen.exportImage(DRIVE_FOLDER, LST_comp, comp_lst_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=CRS_OUT, data_transform=TRAN_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)
                geegen.exportImage(DRIVE_FOLDER, NDWI_comp, comp_ndwi_filename,  export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=CRS_OUT, data_transform=TRAN_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)
                elapsed_time = round(time.time() - start_time)

                COMP_dict += [ee.Feature(None, {'DATE': month2find,
                                                'COMPOSITE': COMP_TYPES[d],
                                                'AquaTerra LST': NBLST_DATES[d],
                                                'Terra REFLECT': NBREFLECT_DATES[d],
                                                'COMPOSITE NAN SCORE LST': COMP_NANSCORES[d]['LST'],
                                                'COMPOSITE NAN SCORE NDWI': COMP_NANSCORES[d]['NDWI'],
                                                'COMPOSITE TIME (sec)': elapsed_time})]
                del start_time, elapsed_time, comp_lst_filename, comp_ndwi_filename, LST_comp, NDWI_comp
            del COMP_TYPES, COMP_NANSCORES, NBLST_DATES, NBREFLECT_DATES

            # --- SAVE TABLES (DATA FRAMES) INTO CSV FILES (LOCAL MACHINE) ---
            geegen.exportDataFrame(DRIVE_FOLDER, QA_dict, QAtable_filename, QAtable_columns, OUTDIR_PATHS[3], path2key)
            geegen.exportDataFrame(DRIVE_FOLDER, COMP_dict, COMPtable_filename, COMPtable_columns, OUTDIR_PATHS[3], path2key)

            # --- EXPORTING MONTH (IF go_month=1) ---
            if go_month==1:
                comp_lst_filename = f'MODIS_LST_{month2find}_COMPM'
                comp_ndwi_filename = f'MODIS_NDWI_{month2find}_COMPM'
                geegen.exportImage(DRIVE_FOLDER, LSTcompM_col, comp_lst_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=CRS_OUT, data_transform=TRAN_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)
                geegen.exportImage(DRIVE_FOLDER, NDWIcompM_col, comp_ndwi_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=CRS_OUT, data_transform=TRAN_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)
                del comp_lst_filename, comp_ndwi_filename
            del LSTcompM_col, NDWIcompM_col

            # --- COPYING TO DATA_HISTO (update historic ref dir) ---
            new_files_d = glob.glob(os.path.join(OUTDIR_PATHS[1], f'*{month2find}*.tif'))
            for fd in new_files_d:
                geegen.copyfile_Errorscontrol(fd, os.path.join(DATA_HISTO, 'DECADE', os.path.basename(fd)))
            del new_files_d
            
            if go_month==1:
                new_files_m = glob.glob(os.path.join(OUTDIR_PATHS[2], f'*{month2find}*.tif'))
                for fm in new_files_m:
                    geegen.copyfile_Errorscontrol(fm, os.path.join(DATA_HISTO, 'MONTH', os.path.basename(fm)))
                del new_files_m



def run_GEEProcessing_LANDSAT_S2(CONFIG, OUTDIR_PATHS, COLLECTIONS, PARAM, DICT):
    """
    Run GEE processing of LANDSAT and SENTINEL-2 products:
        - compute indices
        - temporal compositing
        - export composites (decades)
    """
    
    logging.info('\n\n --- START GEE PROCESSING LANDSAT/SENTINEL-2 ---\n')

    # --- GEE Authentification ---
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    global project_globalId
    project_globalId = project_id

    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_LOCAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']
    if (CONFIG['ASSET_EXPORT_L'] is None) or (CONFIG['ASSET_EXPORT_L']==''): ASSET_EXPORT_L = 0
    else: ASSET_EXPORT_L = int(CONFIG['ASSET_EXPORT_L'])
    if (CONFIG['ASSET_EXPORT_S2'] is None) or (CONFIG['ASSET_EXPORT_S2']==''): ASSET_EXPORT_S2 = 0
    else: ASSET_EXPORT_S2 = int(CONFIG['ASSET_EXPORT_S2'])

    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'LANDSAT_SENTINEL2')

    (LANDSAT_COLLECTIONS, S2_COLLECTIONS,
     NEW_LANDSATCOLLECTIONS, NEW_S2COLLECTIONS, LANDSAT_GRIDS, landmask) = COLLECTIONS
    (SCALE_OUT, LANDSAT_PROJ, DATES) = PARAM
    (QA_dict_L, QAtable_filename_L, QAtable_columns_L,
     QA_dict_S2, QAtable_filename_S2, QAtable_columns_S2,
     COMP_dict, COMPtable_filename, COMPtable_columns) = DICT
    TILES_L = list(LANDSAT_COLLECTIONS.keys())

    months = geegen.googleErrorsControl(ee.List.sequence(1, 12), path2key)
    

    # ===================================== LOOP OVER LANDSAT TILES =================================

    for tile_L in TILES_L:
        logging.info(f'TILE LANDSAT {tile_L} :')

        (L7_tile, L8_tile, L9_tile) = LANDSAT_COLLECTIONS[tile_L]
        S2_alltiles = S2_COLLECTIONS[tile_L]

        proj_out = LANDSAT_PROJ[tile_L]
        (date_start, date_end) = DATES[tile_L]
        years = geegen.googleErrorsControl(ee.List.sequence(date_start.get('year'), date_end.get('year')), path2key)


    # ===================================== LOOP OVER YEARS/MONTHS =================================

        for y in years:
            for m in months:
            
                month2find = '{}{:02d}'.format(y,m)
                logging.info(f'MONTH TO FIND : {month2find}')

                COMP_NANSCORES={}

                L7_month = (L7_tile
                            .filter(ee.Filter.calendarRange(y,y,'year'))
                            .filter(ee.Filter.calendarRange(m,m,'month')))
                L8_month = (L8_tile
                            .filter(ee.Filter.calendarRange(y,y,'year'))
                            .filter(ee.Filter.calendarRange(m,m,'month')))
                L9_month = (L9_tile
                            .filter(ee.Filter.calendarRange(y,y,'year'))
                            .filter(ee.Filter.calendarRange(m,m,'month')))
                S2_month = (S2_alltiles
                            .filter(ee.Filter.calendarRange(y,y,'year'))
                            .filter(ee.Filter.calendarRange(m,m,'month')))
                
                NbL7images_month = geegen.googleErrorsControl(L7_month.size(), path2key)
                NbL8images_month = geegen.googleErrorsControl(L8_month.size(), path2key)
                NbL9images_month = geegen.googleErrorsControl(L9_month.size(), path2key)
                NbS2images_month = geegen.googleErrorsControl(S2_month.size(), path2key)

                Nbimages_month = NbL7images_month + NbL8images_month + NbL9images_month + NbS2images_month

                # --- NO PRODUCTS FOUND ---
                if Nbimages_month==0:
                    logging.info('NO PRODUCTS')
                    continue

                # --- SINGLE PRODUCT FOUND ---
                elif Nbimages_month==1:
                    logging.info('SINGLE PRODUCT')

                    if NbL7images_month!=0:
                        if ASSET_EXPORT_L==1: new_landsatcollection_tmp = NEW_LANDSATCOLLECTIONS[tile_L]
                        else: new_landsatcollection_tmp = ''
                        indices_preproc, QA_dict_L = geeprep.preprocessingL7(L7_month.first(), landmask, LANDSAT_GRIDS[tile_L], tile_L, QA_dict_L, path2key, new_landsatcollection_tmp, ASSET_EXPORT_L)
                        if indices_preproc=='not considered': continue
                        indices_preproc = geegen.calibrateData(indices_preproc, 'L7')
                    
                    if NbL8images_month!=0 :
                        if ASSET_EXPORT_L==1: new_landsatcollection_tmp = NEW_LANDSATCOLLECTIONS[tile_L]
                        else: new_landsatcollection_tmp = ''                        
                        indices_preproc, QA_dict_L = geeprep.preprocessingL8L9(L8_month.first(), landmask, LANDSAT_GRIDS[tile_L], tile_L, QA_dict_L, path2key, new_landsatcollection_tmp, ASSET_EXPORT_L)
                        if indices_preproc=='not considered': continue
                    
                    if NbL9images_month!=0:
                        if ASSET_EXPORT_L==1: new_landsatcollection_tmp = NEW_LANDSATCOLLECTIONS[tile_L]
                        else: new_landsatcollection_tmp = ''                        
                        indices_preproc, QA_dict_L = geeprep.preprocessingL8L9(L9_month.first(), landmask, LANDSAT_GRIDS[tile_L], tile_L, QA_dict_L, path2key, new_landsatcollection_tmp, ASSET_EXPORT_L)
                        if indices_preproc=='not considered': continue
                        indices_preproc = geegen.calibrateData(indices_preproc, 'L9')

                    if NbS2images_month!=0:
                        if ASSET_EXPORT_S2==1: new_s2collection_tile = NEW_S2COLLECTIONS[tile_L][geegen.googleErrorsControl(S2_month.first().get('MGRS_TILE'), path2key)]
                        else: new_s2collection_tile = ''
                        indices_preproc, QA_dict_S2_D = geeprep.preprocessingS2(S2_month.first(), landmask, tile_L, QA_dict_S2_D, path2key, new_s2collection_tile, ASSET_EXPORT_S2)
                        if indices_preproc=='not considered': continue
                        indices_preproc = geegen.calibrateData(indices_preproc, 'S2')
                    
                    INDICES_comp = {}
                    INDICES_comp['NDWI'] = indices_preproc.select('NDWI')
                    INDICES_comp['NDVI'] = indices_preproc.select('NDVI')
                    comp_type = 'DAY'
                    comp_ndwi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDWI_{comp_type}'
                    comp_ndvi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDVI_{comp_type}'
                    COMP_NANSCORES['NDWI']=''
                    COMP_NANSCORES['NDVI']=''

                    # --- Exporting --- 
                    comp_ndwi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDWI_{comp_type}'
                    comp_ndvi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDVI_{comp_type}'

                    start_time = time.time()
                    geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDWI'], comp_ndwi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                    geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDVI'], comp_ndvi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                    elapsed_time = round(time.time() - start_time)

                    COMP_dict += [ee.Feature(None, {'TILE': tile_L,
                                                    'DATE': month2find,
                                                    'COMPOSITE': comp_type,
                                                    'L7': NbL7images_month,
                                                    'L8': NbL8images_month,
                                                    'L9': NbL9images_month,
                                                    'S2': NbS2images_month,
                                                    'TOTAL': Nbimages_month,
                                                    'COMPOSITE NAN SCORE NDWI': COMP_NANSCORES['NDWI'],
                                                    'COMPOSITE NAN SCORE NDVI': COMP_NANSCORES['NDVI'],
                                                    'COMPOSITE TIME (sec)': elapsed_time})]
                    del start_time, elapsed_time, INDICES_comp, COMP_NANSCORES, comp_type

                # --- SEVERAL PRODUCTS FOUND -> MONTH COMPOSITING (period before 2018/2019 : with only L7 and/or L8) ---
                elif (Nbimages_month>1 and NbS2images_month==0 and NbL9images_month==0):
                    logging.info('SEVERAL PRODUCTS (with L7/L8)')
                    indicesMonth_preproc = []

                    if NbL7images_month!=0:
                        L7_list = L7_month.toList(NbL7images_month)
                        for i in range(NbL7images_month):
                            l7 = ee.Image(L7_list.get(i))
                            if ASSET_EXPORT_L==1: new_landsatcollection_tmp = NEW_LANDSATCOLLECTIONS[tile_L]
                            else: new_landsatcollection_tmp = ''                                                
                            indices_preproc, QA_dict_L = geeprep.preprocessingL7(l7, landmask, LANDSAT_GRIDS[tile_L], tile_L, QA_dict_L, path2key, new_landsatcollection_tmp, ASSET_EXPORT_L)
                            if indices_preproc=='not considered': continue
                            indicesMonth_preproc += [geegen.calibrateData(indices_preproc, 'L7')]
                            del l7, indices_preproc
                        del L7_list
                    
                    if NbL8images_month!=0:
                        L8_list = L8_month.toList(NbL8images_month)
                        for i in range(NbL8images_month):
                            l8 = ee.Image(L8_list.get(i))
                            if ASSET_EXPORT_L==1: new_landsatcollection_tmp = NEW_LANDSATCOLLECTIONS[tile_L]
                            else: new_landsatcollection_tmp = ''  
                            indices_preproc, QA_dict_L = geeprep.preprocessingL8L9(l8, landmask, LANDSAT_GRIDS[tile_L], tile_L, QA_dict_L, path2key, new_landsatcollection_tmp, ASSET_EXPORT_L)
                            if indices_preproc=='not considered': continue
                            indicesMonth_preproc += [indices_preproc]
                            del l8, indices_preproc
                        del L8_list
                    
                    indicesMonth_preproc = ee.ImageCollection(indicesMonth_preproc)
                    Nbgoodimages_month = geegen.googleErrorsControl(indicesMonth_preproc.size(), path2key)
                    
                    if Nbgoodimages_month==0:
                        continue
                    elif Nbgoodimages_month==1:
                        comp_type = 'DAY'
                        INDICES_comp = {}
                        INDICES_comp['NDWI'] = indicesMonth_preproc.first().select('NDWI')
                        INDICES_comp['NDVI'] = indicesMonth_preproc.first().select('NDVI')
                        COMP_NANSCORES['NDWI'] = ''
                        COMP_NANSCORES['NDVI'] = ''
                    elif Nbgoodimages_month>1:
                        comp_type = 'COMPM'
                        INDICES_comp = geecomp.extractComposite(indicesMonth_preproc, path2key, LANDSAT_GRIDS[tile_L])
                        COMP_NANSCORES = geegen.computeNanScores(INDICES_comp['NDWI'].addBands(INDICES_comp['NDVI'], overwrite=True).select(['NDWI', 'NDVI']), path2key, landmask, data_geom=LANDSAT_GRIDS[tile_L], scale=SCALE_OUT)
                    del indicesMonth_preproc, Nbgoodimages_month

                    # --- Exporting --- 
                    comp_ndwi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDWI_{comp_type}'
                    comp_ndvi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDVI_{comp_type}'

                    start_time = time.time()
                    geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDWI'], comp_ndwi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                    geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDVI'], comp_ndvi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                    elapsed_time = round(time.time() - start_time)

                    COMP_dict += [ee.Feature(None, {'TILE': tile_L,
                                                    'DATE': month2find,
                                                    'COMPOSITE': comp_type,
                                                    'L7': NbL7images_month,
                                                    'L8': NbL8images_month,
                                                    'L9': NbL9images_month,
                                                    'S2': NbS2images_month,
                                                    'TOTAL': Nbimages_month,
                                                    'COMPOSITE NAN SCORE NDWI': COMP_NANSCORES['NDWI'],
                                                    'COMPOSITE NAN SCORE NDVI': COMP_NANSCORES['NDVI'],
                                                    'COMPOSITE TIME (sec)': elapsed_time})]
                    del start_time, elapsed_time, INDICES_comp, COMP_NANSCORES, comp_type

                # --- SEVERAL PRODUCTS FOUND -> DECADE COMPOSITING (period from 2018/2019 : with S2, L7, L8, and/or L9) ---
                elif (Nbimages_month>1 and NbS2images_month!=0):
                    logging.info('SEVERAL PRODUCTS (with S2/L7/L8/L9)')

                    (NDWIcomp_col, NDVIcomp_col, QA_dict_L, QA_dict_S2, COMP_TYPES, COMP_NANSCORES,
                     NBL7_DATES, NBL8_DATES, NBL9_DATES, NBS2_DATES) = geecomp.processCompositeDecade_LANDSAT_S2(
                                                                                                  L7_tile, L8_tile, L9_tile, S2_alltiles,
                                                                                                  y, m, date_start, date_end,
                                                                                                  landmask, LANDSAT_GRIDS[tile_L], tile_L, SCALE_OUT,
                                                                                                  NEW_LANDSATCOLLECTIONS[tile_L], NEW_S2COLLECTIONS[tile_L],
                                                                                                  QA_dict_L, QA_dict_S2, path2key,
                                                                                                  ASSET_EXPORT_L, ASSET_EXPORT_S2)
                    
                    N_dec_ndwi = geegen.googleErrorsControl(NDWIcomp_col.size(), path2key)
                    N_dec_ndvi = geegen.googleErrorsControl(NDVIcomp_col.size(), path2key)
                    
                    if (N_dec_ndwi==0) or (N_dec_ndvi==0):
                        continue
                    
                    # --- Exporting --
                    NDWIcomp_list = NDWIcomp_col.toList(N_dec_ndwi)
                    NDVIcomp_list = NDVIcomp_col.toList(N_dec_ndwi)

                    for d in range(N_dec_ndwi):
                        INDICES_comp = {}
                        INDICES_comp['NDWI'] = ee.Image(NDWIcomp_list.get(d))
                        INDICES_comp['NDVI'] = ee.Image(NDVIcomp_list.get(d))
                        
                        comp_ndwi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDWI_{COMP_TYPES[d]}'
                        comp_ndvi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDVI_{COMP_TYPES[d]}'

                        start_time = time.time()
                        geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDWI'], comp_ndwi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                        geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDVI'], comp_ndvi_filename, export_folder=OUTDIR_PATHS[1], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])  
                        elapsed_time = round(time.time() - start_time)

                        COMP_dict += [ee.Feature(None, {'TILE': tile_L,
                                                        'DATE': month2find,
                                                        'COMPOSITE': COMP_TYPES[d],
                                                        'L7': NBL7_DATES[d],
                                                        'L8': NBL8_DATES[d],
                                                        'L9': NBL9_DATES[d],
                                                        'S2': NBS2_DATES[d],
                                                        'TOTAL': Nbimages_month,
                                                        'COMPOSITE NAN SCORE NDWI': COMP_NANSCORES[d]['NDWI'],
                                                        'COMPOSITE NAN SCORE NDVI': COMP_NANSCORES[d]['NDVI'],
                                                        'COMPOSITE TIME (sec)': elapsed_time})]
                        del start_time, elapsed_time, INDICES_comp, comp_ndwi_filename, comp_ndvi_filename

                    del NDWIcomp_col, NDVIcomp_col, COMP_TYPES, COMP_NANSCORES, NBL7_DATES, NBL8_DATES, NBL9_DATES, NBS2_DATES

                # --- SAVE TABLES (DATA FRAMES) INTO CSV FILES (LOCAL MACHINE) ---
                geegen.exportDataFrame(DRIVE_FOLDER, QA_dict_L, QAtable_filename_L, QAtable_columns_L, OUTDIR_PATHS[3], path2key)
                geegen.exportDataFrame(DRIVE_FOLDER, QA_dict_S2, QAtable_filename_S2, QAtable_columns_S2, OUTDIR_PATHS[3], path2key)
                geegen.exportDataFrame(DRIVE_FOLDER, COMP_dict, COMPtable_filename, COMPtable_columns, OUTDIR_PATHS[3], path2key)
                
                del L7_month, L8_month, L9_month, S2_month, NbL8images_month, NbL9images_month, NbS2images_month, Nbimages_month

                # --- COPYING TO DATA_HISTO (update historic ref dir) ---
                new_files_d = glob.glob(os.path.join(OUTDIR_PATHS[1], f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}*.tif'))
                for fd in new_files_d:
                    geegen.copyfile_Errorscontrol(fd, os.path.join(DATA_HISTO, 'DECADE', os.path.basename(fd)))
                del new_files_d, month2find



def run_GEECompositingMonth_MODIS(CONFIG, OUTDIR_PATHS, PARAM):
    """
    Month compositing of already processed MODIS LST and NDWI single-date indices.
    This function is to apply once indices have been exported to GEE assets (not implemented in pipeline).
    
    Note : The number of good data used in the month period is also given per pixel.
    """

    logging.info('\n\n --- GEE MONTH COMPOSITING MODIS ---\n')

    (GRID_OUT, SCALE_OUT, CRS_OUT, date_start, date_end) = PARAM
    
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_GLOBAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']

    # --- Localize gee Collection (assets) with indices corresponding to the territory ---
    folderName = f'PREPROC_GLOBAL_INDICES_{TERRITORY_str}'
    gee_workdir = 'projects'+'/'+project_id+'/'+'assets'
    gee_folder = gee_workdir+'/'+folderName
    modis_collection = ee.ImageCollection(gee_folder+'/'+'MODIS').filterDate(date_start, date_end)
    

    # ========================================== LOOP OVER YEARS/MONTHS =================================
    
    years = geegen.googleErrorsControl(ee.List.sequence(date_start.get('year'), date_end.get('year')), path2key)
    months = geegen.googleErrorsControl(ee.List.sequence(1, 12), path2key)

    for y in years:
        for m in months:
            month2find = '{}{:02d}'.format(y,m)
            logging.info(f'MONTH : {month2find}')

            MODIS_month = (modis_collection
                           .filter(ee.Filter.calendarRange(y,y,'year'))
                           .filter(ee.Filter.calendarRange(m,m,'month')))
            
            Nb_month = geegen.googleErrorsControl(MODIS_month.size(), path2key)
            
            # --- NO PRODUCTS FOUND ---
            if Nb_month==0:
                continue

            # --- MONTH COMPOSITING ---
            INDICES_comp = geecomp.extractComposite(MODIS_month)

            # --- EXPORTING ---                
            comp_lst_filename = f'MODIS_LST_{month2find}_COMPM'
            comp_ndwi_filename = f'MODIS_NDWI_{month2find}_COMPM'
            geegen.exportImage(DRIVE_FOLDER, INDICES_comp['LST'], comp_lst_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=CRS_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)
            geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDWI'], comp_ndwi_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=CRS_OUT, data_scale=SCALE_OUT, data_region=GRID_OUT)

            del INDICES_comp, Nb_month, MODIS_month, comp_lst_filename, comp_ndwi_filename



def run_GEECompositingMonth_LANDSAT_S2(CONFIG, OUTDIR_PATHS, COLLECTIONS, PARAM):
    """
    Month compositing of already processed LANDSAT/S2 NDWI and NDVI single-date indices.
    This function is to apply once indices have been exported to GEE assets (not implemented in pipeline).
    
    Note : The number of good data used in the month period is also given per pixel.
    """

    logging.info('\n\n --- GEE MONTH COMPOSITING LANDSAT/SENTINEL-2 ---\n')
    
    path2key = os.path.dirname(geeauth.__file__)
    project_id = geegen.googleAuthentification(path2key)
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    if (CONFIG['DRIVE_FOLDER'] is None) or (CONFIG['DRIVE_FOLDER']==''): DRIVE_FOLDER = 'EO4DM_EXPORT_'+TERRITORY_str+'_LOCAL'
    else: DRIVE_FOLDER = CONFIG['DRIVE_FOLDER']

    (LANDSAT_COLLECTIONS, S2_COLLECTIONS,
     NEW_LANDSATCOLLECTIONS, NEW_S2COLLECTIONS, LANDSAT_GRIDS, landmask) = COLLECTIONS
    (SCALE_OUT, LANDSAT_PROJ, DATES) = PARAM
    TILES_L = list(NEW_LANDSATCOLLECTIONS.keys())


    # ===================================== LOOP OVER LANDSAT TILES =================================

    for tile_L in TILES_L:
        logging.info(f'TILE LANDSAT : {tile_L}')

        # --- Extract gee Collection(s) (assets) corresponding to the Landsat Tile ---
        proj_out = LANDSAT_PROJ[tile_L]
        (date_start, date_end) = DATES[tile_L]
        landsats2_collection = ee.ImageCollection(NEW_LANDSATCOLLECTIONS[tile_L]).filterDate(date_start, date_end)
        proj_out = geegen.googleErrorsControl(landsats2_collection.first().select(1).projection(), path2key)

        for tile_S2 in NEW_S2COLLECTIONS[tile_L]:
            landsats2_collection = ee.merge(landsats2_collection,
                                            ee.ImageCollection(NEW_S2COLLECTIONS[tile_L][tile_S2]).filterDate(date_start, date_end))


        # ========================================== LOOP OVER YEARS/MONTHS =================================
        
        years = geegen.googleErrorsControl(ee.List.sequence(date_start.get('year'), date_end.get('year')), path2key)
        months = geegen.googleErrorsControl(ee.List.sequence(1, 12), path2key)

        for y in years:
            for m in months:
            
                month2find = '{}{:02d}'.format(y,m)
                logging.info(f'MONTH : {month2find}')

                landsats2_month = (landsats2_collection
                                    .filter(ee.Filter.calendarRange(y,y,'year'))
                                    .filter(ee.Filter.calendarRange(m,m,'month')))
                
                Nb_month = geegen.googleErrorsControl(landsats2_month.size(), path2key)
                
                # --- NO PRODUCTS FOUND ---
                if Nb_month==0:
                    continue

                # --- MONTH COMPOSITING ---
                INDICES_comp = geecomp.extractComposite(landsats2_month, path2key, LANDSAT_GRIDS[tile_L])

                # --- EXPORTING ---                
                comp_ndwi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDWI_COMPM'
                comp_ndvi_filename = f'LANDSAT_SENTINEL2_0{tile_L}_{month2find}_NDVI_COMPM'
                geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDWI'], comp_ndwi_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                geegen.exportImage(DRIVE_FOLDER, INDICES_comp['NDVI'], comp_ndvi_filename, export_folder=OUTDIR_PATHS[2], path2key=path2key, data_crs=proj_out['crs'], data_scale=SCALE_OUT, data_region=LANDSAT_GRIDS[tile_L])
                del INDICES_comp, Nb_month, landsats2_month, comp_ndwi_filename, comp_ndvi_filename
        
        del landsats2_collection, proj_out, date_start, date_end

