# -*- coding: utf-8 -*-
"""
##############################################################################

DROUGHT functions for processing Local drought indicator

##############################################################################
"""

import os
import glob
import shutil
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio.windows import Window
from rasterio.warp import reproject, Resampling
from tqdm import tqdm

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def prepare_DATAFolders(CONFIG):
    """
    Create folders tree for historic data (DATA_HISTO) and controls annex directory (ANNEX_DIR) :
        - DATA_HISTO will contain all indices computed and needed for the calculation of drought indicators
        - ANNEX_DIR must be created and filled by user (Landsat_Grids used in VAI computation to clip tiles)
    """

    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    
    # --- Set current directory to WRK_DIR ---

    WRK_DIR = os.path.normpath(CONFIG['WRK_DIR'])
    try:
        os.chdir(WRK_DIR)
    except FileNotFoundError as e:
        logging.critical(f'Error WRK directory : {e}\n-> Maybe due to different path naming convention on Windows (C:/...) or Linux (/C/...)\n')
        raise Exception('Error WRK directory')

    # --- DATA HISTO directory ---
    
    HISTO_DIR = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str)
    TREE_DIR = [os.path.join(HISTO_DIR, '0_INDICES', 'LANDSAT_SENTINEL2', 'DECADE'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'LOCAL', 'DECADE')]
    os.umask(0) # used to reset the directories permissions

    try:
        # Create the main dir if not exist
        os.makedirs(HISTO_DIR, exist_ok=True)
        os.chmod(HISTO_DIR, 0o777)

        # Browse folders list and create new ones if not already exist
        for dir in TREE_DIR:
            os.makedirs(dir, exist_ok=True)

    except OSError as e:
        logging.critical(f'Error when creating DATA HISTO folders: {e}')
        raise Exception('Error DATA HISTO directory')
    
    del TREE_DIR
    
     # --- ANNEX directory ---
    
    ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'])
    TREE_DIR = [os.path.join(ANNEX_DIR, 'Landsat_Grid_World')]

    # Control that directories were well created and filled
    if not os.path.exists(ANNEX_DIR):
        logging.critical(f'ANNEX directory was not found : need to be filled for LOCAL DROUGHT computation !\n')
        raise Exception('Error ANNEX directory')
    
    # Browse and control folders list
    for dir in TREE_DIR:
        if not os.path.exists(dir):
            logging.critical(f'In ANNEX, {os.path.basename(dir)} directory was not found : needed for LOCAL DROUGHT computation !\n')
            raise Exception('Error ANNEX directory')
        elif len(glob.glob(os.path.join(dir, '*')))==0:
            logging.critical(f'In ANNEX, {os.path.basename(dir)} directory is empty : needed for LOCAL DROUGHT computation !\n')
            raise Exception('Error ANNEX directory')
        
    return WRK_DIR



def prepare_RUNFolders(CONFIG):
    """
    Prepare folders that will contain output data (rasters, csv, etc.) for a specific run
    Verify if the main folder exists :
            If no -> create folder (and sub-folders)
            If yes -> refers to "CONFIG['CLEAN_RUNFOLDER']" to delete or keep it
    """
    
    # --- Prepare main folder name ---
    PERIOD_START = CONFIG['PERIOD_START'].split(',')[0]
    PERIOD_END = CONFIG['PERIOD_END'].split(',')[0]
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    date_start_str = PERIOD_START.replace('-', '')
    date_end = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)
    date_end_str = date_end.strftime('%Y%m%d')

    if (CONFIG['TILES_L'] is None) or (CONFIG['TILES_L']==''): TILES_str=''
    else:
        if type(CONFIG['TILES_L']) is list: TILES_L = CONFIG['TILES_L']
        else: TILES_L = CONFIG['TILES_L'].split(',')
        TILES_str = ''
        if len(TILES_L)<=4:
            for t in TILES_L: TILES_str = f'{TILES_str}_{t}'

    MODE = CONFIG['MODE']
    if (CONFIG['CLEAN_RUNFOLDER'] is None) or (CONFIG['CLEAN_RUNFOLDER']==''): CLEAN_RUNFOLDER = 0
    else: CLEAN_RUNFOLDER = int(CONFIG['CLEAN_RUNFOLDER'])

    if 'AUTO' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_AUTO_LOCAL_{TERRITORY_str}{TILES_str}_{date_start_str}_{date_end_str}')
    elif 'MANUAL' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_MANUAL_LOCAL_{TERRITORY_str}{TILES_str}_{date_start_str}_{date_end_str}')
    elif 'INDICES' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_INDICES_LOCAL_{TERRITORY_str}{TILES_str}_{date_start_str}_{date_end_str}')
    elif 'DROUGHT' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_DROUGHT_LOCAL_{TERRITORY_str}{TILES_str}_{date_start_str}_{date_end_str}')
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')
    
    # --- Generate directory ---
    os.umask(0) # used to reset the directories permission
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        os.chmod(outdir, 0o777)
    elif CLEAN_RUNFOLDER==1:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
        os.chmod(outdir, 0o777)

    # --- Generate sub-directories ---
    outdir_comp = os.path.normpath(outdir + os.sep + '1_GEE_COMPOSITES/')
    outdir_compdecade = os.path.normpath(outdir_comp + os.sep + 'DECADE/')
    outdir_compmonth = os.path.normpath(outdir_comp + os.sep + 'MONTH/')
    outdir_compstats = os.path.normpath(outdir_comp + os.sep + 'STATS/')
    outdir_drought = os.path.normpath(outdir + os.sep + '2_DROUGHT_INDICATORS/')
    outdir_vai = os.path.normpath(outdir_drought + os.sep + '1_VAI/')
    outdir_postproc = os.path.normpath(outdir_drought + os.sep + '2_POSTPROC/')
    outdir_droughtstats = os.path.normpath(outdir_drought + os.sep + 'STATS/')

    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):
        if not os.path.exists(outdir_comp): os.makedirs(outdir_comp)
        if not os.path.exists(outdir_compdecade): os.makedirs(outdir_compdecade)
        if not os.path.exists(outdir_compmonth): os.makedirs(outdir_compmonth)
        if not os.path.exists(outdir_compstats): os.makedirs(outdir_compstats)
    
    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('DROUGHT' in MODE):
        if not os.path.exists(outdir_drought): os.makedirs(outdir_drought)
        if not os.path.exists(outdir_vai): os.makedirs(outdir_vai)
        if not os.path.exists(outdir_postproc): os.makedirs(outdir_postproc)
        if not os.path.exists(outdir_droughtstats): os.makedirs(outdir_droughtstats)

    # --- Concenate dir paths into an output parameter ---
    OUTDIR_PATHS = (outdir_comp, outdir_compdecade, outdir_compmonth, outdir_compstats, outdir_vai, outdir_postproc, outdir_droughtstats) 

    return OUTDIR_PATHS



def auto_processingPERIOD(CONFIG):
    """
    Automatic settint of processing period :
        - extract period of already processed data (DATA_HISTO)
        - find last decade to be possibly ready for processing (the decade before the current decade)
        - give new period to process
        - if last decade was already processed => NO PROCESSING (go_process=0)

    Note : if mode (in config) is not set to 'AUTO' -> go out of function and return go_process=1
    """

    MODE = CONFIG['MODE']

    if 'AUTO' in MODE:
        logging.info(f'\n\n --- AUTOMATIC SETTING OF PROCESSING PERIOD (MODE = {MODE}) ---\n')

    elif ('MANUAL' in MODE) or ('INDICES' in MODE):
        PERIOD_START = CONFIG['PERIOD_START']
        PERIOD_END = CONFIG['PERIOD_END']
        if (PERIOD_START is None)  or (PERIOD_START==''):
            logging.critical(f'PERIOD START variable is None : need to be filled in {MODE} mode !\n')
            raise Exception('Wrong PERIOD variable')
        if (PERIOD_END is None) or (PERIOD_END==''):
            logging.critical(f'PERIOD END variable is None : need to be filled in {MODE} mode !\n')
            raise Exception('Wrong PERIOD variable')

        logging.info(f'\n\n --- MANUAL SETTING OF PROCESSING PERIOD (MODE = {MODE}) ---\n\nGO PROCESSING : from {PERIOD_START} to {PERIOD_END}\n')
        go_process = 1
        return go_process, PERIOD_START, PERIOD_END
    
    elif 'DROUGHT' in MODE:
        PERIOD_START = CONFIG['PERIOD_START']
        PERIOD_END = CONFIG['PERIOD_END']

        if ((PERIOD_START is None) or (PERIOD_START=='')) and ((PERIOD_END is None) or (PERIOD_END=='')):
            logging.info(f'\n\n --- AUTOMATIC SETTING OF PROCESSING PERIOD (MODE = {MODE}) ---\n')
            go_process = 1
        elif ((PERIOD_START is not None) or (PERIOD_START!='')) and ((PERIOD_END is not None) or (PERIOD_END!='')):
            logging.info(f'\n\n --- MANUAL SETTING OF PROCESSING PERIOD (MODE = {MODE}) ---\n')
            go_process = 1
            return go_process, PERIOD_START, PERIOD_END
        else:
            logging.critical(f'One of PERIOD START or PERIOD END variable is None : both need to be filled in {MODE} mode !\n')
            raise Exception('Wrong PERIOD variable')

    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')


    start_D = [1, 11, 21]
    end_D = [10, 20, 31]
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'LANDSAT_SENTINEL2', 'DECADE')
    
    # Controls input landsat TILES for filtering (or not) already processed products
    if (CONFIG['TILES_L'] is None) or (CONFIG['TILES_L']==''):
        logging.warning(f'No input landsat tiles : the last historic product will be the earliest detected on any tile !\n -> Fill TILES_L if a specific tile is needed \n')
        histo_files = glob.glob(os.path.join(DATA_HISTO, '*.tif'))

    else:
        if type(CONFIG['TILES_L']) is list: TILES_L = CONFIG['TILES_L']
        else: TILES_L = CONFIG['TILES_L'].split(',')

        if len(TILES_L)==1:
            histo_files = glob.glob(os.path.join(DATA_HISTO, f'*{TILES_L[0]}*.tif'))
            logging.info(f'Single input landsat tile {TILES_L[0]} : the last historic product will be detected considering this specific tile !\n')
        else:
            histo_files = glob.glob(os.path.join(DATA_HISTO, '*.tif'))
            logging.warning(f'Several input landsat tile {TILES_L} : the last historic product detected will be the earliest on any tile !\n -> Fill TILES_L if a specific tile is needed \n')

    # Set new PERIOD END = end of the last full decade (previous one)
    date_now = pd.Timestamp.now()
    d_now = date_now.day
    m_now = date_now.month
    y_now = date_now.year
    for d in range(3):
        if d_now>=start_D[d] and d_now<=end_D[d]:
            d_prev = d-1
            if d_prev>=0:
                period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_now,m_now,start_D[d_prev]))
                period_end_d =  pd.to_datetime('{}-{:02d}-{:02d}'.format(y_now,m_now,start_D[d_prev])) + pd.DateOffset(days=10)
            else:
                d_prev=2
                period_end_d = date_now + pd.DateOffset(days=-d_now+1)
                y_prev = period_end_d.year
                m_prev = period_end_d.month
                period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_prev,m_prev,start_D[d_prev]))
            break
    period_end = pd.Timestamp(period_end_d)
    PERIOD_END = period_end.strftime('%Y-%m-%d')
    
    # IF HISTO PRODUCTS : Set new PERIOD START according to last histo product processed
    if len(histo_files)>0:
        dates_list = [os.path.basename(f).split('_')[3] for f in histo_files]
        dates = np.unique(dates_list)
        dates = pd.to_datetime(dates, format='%Y%m')
        period_start_histo = pd.Series(dates).min()
        period_end_histo = pd.Series(dates).max()
        PERIOD_START_HISTO = period_start_histo.strftime('%Y-%m-%d')
        PERIOD_END_HISTO = period_end_histo.strftime('%Y%m')

        if len(TILES_L)==1:
            list_end_histo =  glob.glob(os.path.join(DATA_HISTO, f'*{TILES_L[0]}_{PERIOD_END_HISTO}*.tif'))
        else:
            list_end_histo =  glob.glob(os.path.join(DATA_HISTO, f'*{PERIOD_END_HISTO}*.tif'))
        D_END_HISTO = list(np.unique([f.split('_')[-1][:-4] for f in list_end_histo]))
        d_number_list = []
        for f in D_END_HISTO:
            for d in f:
                if d.isdigit(): d_number_list.append(int(d))
        if d_number_list==[]:
            # Case with 'COMPM' or 'DAY' -> decade_number is 3rd decade
            d_number = 2
        else:
            # Other cases -> decade_number is max number-1
            d_number =  np.max(d_number_list)-1
        
        d_end_histo = end_D[d_number]
        if d_number==2:
            # 3rd decade -> add number of days in the month (cause can vary)
            d_offset = period_end_histo.days_in_month
        else :
            # Other cases -> add d_en_histo
            d_offset = d_end_histo
        
        period_end_histo = period_end_histo + pd.DateOffset(days=d_offset-1)
        PERIOD_END_HISTO = period_end_histo.strftime('%Y-%m-%d')
        logging.info(f'Histo data were already processed from {PERIOD_START_HISTO} to {PERIOD_END_HISTO}')

        if 'DROUGHT' in MODE:
            period_end_histo = period_end_histo + pd.DateOffset(days=1)  # period end is exclusive
            PERIOD_END_HISTO =  period_end_histo.strftime('%Y-%m-%d')
            go_process = 1
            return go_process, PERIOD_START_HISTO, PERIOD_END_HISTO
        
        period_start = period_end_histo + pd.DateOffset(days=1)
        PERIOD_START = period_start.strftime('%Y-%m-%d')

    # IF NO HISTO PRODUCTS : Set PERIOD START as empty (will be set later by the first product available)
    else:
        if 'DROUGHT' in MODE:
            logging.critical(f'No histo data processed : needed in {MODE} mode !\n')
            raise Exception('No HISTO DATA')
        else:
            logging.info('No histo data processed')
            PERIOD_START = 'First product'
            period_start = period_end + pd.DateOffset(days=-1)

    # CONTROLS if last decade was already processed or not :
    if period_start<period_end:
        logging.info(f'GO PROCESSING : from {PERIOD_START} to {PERIOD_END}')
        go_process = 1
    else:
        PERIOD_START_D = period_start_d.strftime('%Y-%m-%d')
        logging.info(f'NO PROCESSING : Last decade D{d_prev+1} ({PERIOD_START_D} -> {PERIOD_END}) already processed')
        go_process = 0

    return go_process, PERIOD_START, PERIOD_END



def copyfile_Errorscontrol(src, dst):
    """
    Procedure to copy files and control in case of permission errors.
    """

    try:
        shutil.copyfile(src, dst)
    except PermissionError:
        if os.path.exists(dst):
            logging.warning(f'File already exists: {os.path.basename(src)} is replaced by new version')
            os.remove(dst)
            shutil.copyfile(src, dst)
        else:
            logging.critical(f'Copy PermissionError : {os.path.basename(src)} impossible to paste')
            raise Exception('Copy PermissionError : impossible to paste')



def extractProfilesLike_Tiles(data_dir, tiles_list):
    """
    Extract the profiles of historical products in data_dir, for each landsat tile listed in tiles_list
    """

    PROFILE_LIKE = {}
    
    for tile in tiles_list:
        file_like = glob.glob(os.path.join(data_dir, f'*{tile}*.tif'))[0]
        with rasterio.open(file_like) as f_ds:
            profile_tile = f_ds.profile
            profile_tile.update(count=2)
        
        PROFILE_LIKE[tile] = profile_tile
        del file_like, profile_tile
    
    return PROFILE_LIKE



def extractGridGeo(lgrid_shp, data_dir, tiles_list):
    """
    Extract Grids Geometries of landsat tiles from input shapefile containing landsat grids geometries.
    The output grids geometries are reprojected to the crs of an input reference raster.
    """

    GRIDS_GEO = {}

    if os.path.exists(lgrid_shp):
        for tile in tiles_list:
            file_like = glob.glob(os.path.join(data_dir, f'*{tile}*.tif'))[0]
            with rasterio.open(file_like) as f_ds:
                lgrid_ds = gpd.read_file(lgrid_shp)
                lgrid_ds = lgrid_ds.to_crs(f_ds.crs)
                lgrid_tile = lgrid_ds[lgrid_ds["PR"]==int(tile[1:])]
                lgrid_geo = [feature for feature in lgrid_tile["geometry"]]
            
            GRIDS_GEO[tile] = lgrid_geo
            del file_like, lgrid_geo
    
    else:
        logging.warning('Missing input shapefile containing geometries to repoject/clip rasters if new indices do not have same sizes as historic')

    return GRIDS_GEO



def reprojectData(data, src_profile, dst_profile, resampling_algo=Resampling.bilinear):
    """
    Function for reprojecting data matrix from input src_profile, dst_profile and resampling algorithm
    """

    data_rpj = np.full((dst_profile['height'],dst_profile['width']),
                       src_profile['nodata'],
                       dtype=src_profile['dtype'])
    reproject(data,
              data_rpj,
              src_transform=src_profile['transform'],
              src_crs=src_profile['crs'],
              src_nodata=src_profile['nodata'],
              dst_transform=dst_profile['transform'],
              dst_crs=dst_profile['crs'],
              dst_nodata=src_profile['nodata'],
              resampling=resampling_algo)
    
    return data_rpj



def clipData(b_ds, geo):
    """
    Function for clipping input matrix in a data source (read from rasterio) to an input geometry
    """

    b_out, transform_out = rasterio.mask.mask(b_ds, geo, nodata=np.nan, crop=True)
    
    profile_out = b_ds.profile
    profile_out.update(height=b_out.shape[1],
                       width=b_out.shape[2],
                       transform=transform_out,
                       count=b_out.shape[0])
    
    return b_out, profile_out



def process_ReprojectionClipping_Rasters(in_ds, dst_profile, lgrid_geo, out_path):
    """
    Processing the reprojection and clipping of rasters in case of different sizes of indices (new VS histo).
    This function is needed to then compute VAI : indices matrices must be of exact same size to compute per pixel historical anomalies.
    """

    dst_profile_count = dst_profile.copy()
    dst_profile_count.update(nodata=0, dtype='uint8')

    # READ AND REPROJECT BAND 1 (index)
    src_profile = in_ds.profile
    src_profile.update(nodata=np.nan)
    dst_profile.update(count=src_profile['count'])
    index_in = in_ds.read(1)
    index_rpj = reprojectData(index_in, src_profile, dst_profile)
    del index_in
    
    # READ AND REPROJECT BAND 2 (count if needed)
    if src_profile['count']==2:
        count_in = in_ds.read(2)
        src_profile_count = src_profile.copy()
        src_profile_count.update(nodata=0, dtype='uint8')
        dst_profile_count.update(count=src_profile_count['count'])
        count_rpj = reprojectData(count_in, src_profile_count, dst_profile_count, resampling_algo=Resampling.nearest)
        del count_in, src_profile_count

    # SAVE REPROJECTED RASTER
    with rasterio.open(out_path, 'w', **dst_profile) as out_ds:
        out_ds.write(index_rpj, 1)
        if src_profile['count']==2:
            out_ds.write(count_rpj, 2)
            del count_rpj
    
    del src_profile, index_rpj

    # CLIP RASTER
    with rasterio.open(out_path) as in_ds:
        data_out, profile_out = clipData(in_ds, lgrid_geo)
        
    # UPDATE type to float32
    profile_out.update(dtype='float32')
    
    # SAVE CLIPPED RASTER
    with rasterio.open(out_path, 'w', **profile_out) as out_ds:
        out_ds.write(data_out[0,:,:], 1)
        if profile_out['count']==2:
            out_ds.write(data_out[1,:,:], 2)



def removeCurrentYear_DATAList(DATAfiles):
    """
    In list of DATA files, remove file corresponding to the current year.
    Used here to remove the current ndwi landsat/s2 indices before computing historical mean/std for VAI calculation.
    """

    curr_year = pd.Timestamp.now().year
    DATAfiles_histo = []
    DATAfiles_current = []

    for f in DATAfiles:

        date_str = os.path.basename(f).split('_')[3]
        try:
            date_year = pd.to_datetime(date_str, format='%Y%m').year
        except ValueError:
            logging.info(f'Extracting date in Data files list : Wrong date format for {date_str}, then not considered')
            continue
        if date_year < curr_year:
            DATAfiles_histo.append(f)
        else:
            DATAfiles_current.append(f)

    return DATAfiles_histo, DATAfiles_current



def extractVAI(files, profile_like, ndays_max, outdir, period_indic=''):
    """
    Computing Vegetation Anomaly Index from Amri et al.(2011) and Peters et al.(2002):
        - VAI = (NDWIi - NDWImean) / NDWIstd)
        - QSCORE_comp = count_i/ndays_max
        - Qscore = mean(Qscore_i,Qscore_histo)
    
    Note 1: The composite Quality Score (QSCORE_comp) is computed for each specific date,
            from the ratio of valid pixels counted per indice (count_i) divided by the maximum
            number of days possibly observed on the period (here 1 image/day).
            A FINAL VAI Quality Score (QSCORE_vai) is computed taking into account the composite score
            QSCORE_comp and the seasonal QSCORES mean on the period observed (QSCORE_histo).
    
    Note 2 : NDWI of current year are not included in the computation of historical NDWImean,NDWIstd,QSCORE_histo.
    
    Note 3 : Rasters are divided into 4 blocks to avoid memory overload
    """

    # Extracting ndwi of current year
    files_histo, files_current = removeCurrentYear_DATAList(files)
    N_files = len(files)
    N_files_histo = len(files_histo)

    # Dividing into 4 blocks
    W_b1 = int(profile_like['width']/2)
    H_b1 = int(profile_like['height']/2)
    W_b2 = int(np.ceil(profile_like['width']/2))
    H_b2 = int(np.ceil(profile_like['height']/2))

    col_off = [0, W_b1, 0,    W_b1]
    row_off = [0, 0,    H_b1, H_b1]
    W = [W_b1, W_b2, W_b1, W_b2]
    H = [H_b1, H_b1, H_b2, H_b2]

    profile_in = profile_like.copy()

    # --- RASTERS ARE DIVIDED INTO 4 BLOCKS (to avoid memory overload) ---
    for b in range(4):

        logging.info(f'Block {b+1}/4')
    
        # --- COMPUTING HISTORICAL MEAN/STD and QSCORE ---
        DATA = np.full((H[b],W[b],N_files), np.nan, 'float32')
        QSCORE_year = np.full_like(DATA, np.nan, 'float32')

        logging.info('Computing historical mean/std and qscore')
        for i in tqdm(range(N_files_histo)):
            with rasterio.open(files_histo[i]) as d_ds:
                DATA[:,:,i] = d_ds.read(1, window=Window(col_off[b], row_off[b], W[b], H[b]))
                if d_ds.count==2:
                    COUNT = d_ds.read(2, window=Window(col_off[b], row_off[b], W[b], H[b]))
                elif d_ds.count==1:
                    COUNT = ~np.isnan(DATA[:,:,i])
            QSCORE_year[:,:,i] = COUNT/ndays_max
            del COUNT

        MEAN_histo = np.nanmean(DATA,axis=2)
        STD_histo = np.nanstd(DATA,axis=2)
        QSCORE_histo = np.mean(QSCORE_year[:,:,:i], axis=2)

        # Here adding ndwi of current year
        for f in files_current:
            i+=1
            with rasterio.open(f) as d_ds:
                DATA[:,:,i] = d_ds.read(1, window=Window(col_off[b], row_off[b], W[b], H[b]))
                if d_ds.count==2:
                    COUNT = d_ds.read(2, window=Window(col_off[b], row_off[b], W[b], H[b]))
                elif d_ds.count==1:
                    COUNT = ~np.isnan(DATA[:,:,i])
            QSCORE_year[:,:,i] = COUNT/ndays_max
            del COUNT

        # --- APPLY VAI EQUATION TO EACH IMAGE ON THE PERIOD ---
        logging.info('Apply VAI equation to each image on the period')
        for i in tqdm(range(N_files)):

            #  Compute VAI
            VAI = (DATA[:,:,i] - MEAN_histo)/STD_histo
            VAI[(STD_histo == 0)] = np.nan
            VAI_QSCORE = np.mean(np.stack((QSCORE_histo,QSCORE_year[:,:,i]),axis=2),axis=2)
            VAI_QSCORE[np.isnan(DATA[:,:,i])==1] = 0
            
            # Save file
            in_file_name = os.path.basename(files[i]).split('.tif')[0]
            tile = in_file_name.split('_')[2]
            date = in_file_name.split('_')[3]
            out_file_name = os.path.join(outdir, f'VAI_{tile}_{date}{period_indic}.tif')
            if b==0:
                # Writing rasters (block 1)
                mode_rio='w'
            elif b==1:
                # Updating rasters (blocks 2,3,4)
                mode_rio='r+'
                # Remove multiple keys from dictionary (avoid warnings)
                rem_list = ['blockxsize','blockysize','interleave','compress','tiled']
                for key in rem_list:
                    try: profile_in.pop(key)
                    except KeyError: pass
            
            with rasterio.open(out_file_name, mode=mode_rio, **profile_in) as out_ds:
                out_ds.write(VAI, window=Window(col_off[b], row_off[b], W[b], H[b]), indexes=1)
                out_ds.write(VAI_QSCORE, window=Window(col_off[b], row_off[b], W[b], H[b]), indexes=2)
            del VAI, VAI_QSCORE, in_file_name, out_file_name, tile, date
    
        del DATA, QSCORE_year, MEAN_histo, STD_histo, QSCORE_histo



def process_LocalDrought_VAI(CONFIG, OUTDIR_PATHS):
    """
    Vegetation Anomaly Index (VAI) is estimated from NDWI anomalies.
    A quality score is given (Qscore) according to the number valid pixels,
    and the maximum number of days possibly oberved on the compositing period (Ndays_max).
    
    NDAYS_MAX_D (decade) is set to 10.
    We admit that the best case is a product available each day of the decade.
    This is the case for MODIS, but not for Landat/S2 (maximum is ~6 per decade).
    This is set in order to permit the comparison of the quality of both indicators (VHI vs VAI)
    
    Final products are saved (updated) into date_histo directory.
    
    Note :  vai is given at S2 10m resolution
            2 bands are saved b1=VAI, b2=Qscore
    """

    logging.info('\n\n--- PROCESSING LOCAL DROUGHT INDICATOR (VAI, 10 m) ---\n')

    NDAYS_MAX_D = 10
    (outdir_comp, outdir_compdecade, outdir_compmonth, outdir_compstats,
     outdir_vai, outdir_postproc, outdir_droughtstats) = OUTDIR_PATHS
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'LANDSAT_SENTINEL2')
    DATA_DROUGHT = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'LOCAL')
    LANDSAT_GRID = os.path.join(CONFIG['ANNEX_DIR'], 'Landsat_Grid_World', 'Landsat_Grid_World.shp')
    MODE = CONFIG['MODE']
    
    # --- AUTOMATIC/MANUAL/INDICES MODE -> process VAI to corresponding new indices months/decades ---
    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):
        
        logging.info(f'{MODE} MODE : process VAI to corresponding new indices months/decades')
        
        new_files = glob.glob(os.path.join(outdir_compdecade, '*.tif'))

        m_list = [os.path.basename(f).split('_')[3][-2:] for f in new_files]
        m_str = np.unique(m_list)
        tiles_list = [os.path.basename(f).split('_')[2] for f in new_files]
        TILES_L = np.unique(tiles_list)
        PROFILES_L = extractProfilesLike_Tiles(os.path.join(DATA_HISTO, 'DECADE'), TILES_L)
        GRIDS_L = extractGridGeo(LANDSAT_GRID, os.path.join(DATA_HISTO, 'DECADE'), TILES_L)
        
        # Control that new indices sizes are like histo indices -> if not, reprojection and clipping according to histo profile
        for f in new_files:
            tile_L = os.path.basename(f).split('_')[2]
            profile_like = PROFILES_L[tile_L]
            lgrid_geo = GRIDS_L[tile_L]
            with rasterio.open(f) as f_ds:
                if f_ds.profile['height']!=profile_like['height'] \
                    or f_ds.profile['width']!=profile_like['width']:
                    process_ReprojectionClipping_Rasters(f_ds, profile_like, lgrid_geo, os.path.join(DATA_HISTO, 'DECADE', os.path.basename(f)))
            del tile_L, profile_like
    
    # --- DROUGHT MODE -> compute/update vai products from data historic directory ---
    elif 'DROUGHT' in MODE:
        if (CONFIG['TILES_L'] is None) or (CONFIG['TILES_L']==''):
            logging.info(f'{MODE} MODE : recompute vai products from data historic directory')
            new_files = glob.glob(os.path.join(DATA_HISTO, 'DECADE', '*.tif'))
            tiles_list = [os.path.basename(f).split('_')[2] for f in new_files]
            TILES_L = np.unique(tiles_list)
        else:
            if type(CONFIG['TILES_L']) is list: TILES_L = CONFIG['TILES_L']
            else: TILES_L = CONFIG['TILES_L'].split(',')
            if TILES_L !=['']:
                logging.info(f'{MODE} MODE : recompute vai products from data historic directory for tiles {TILES_L}')
            else:
                logging.info(f'{MODE} MODE : recompute vai products from data historic directory')
                new_files = glob.glob(os.path.join(DATA_HISTO, 'DECADE', '*.tif'))
                tiles_list = [os.path.basename(f).split('_')[2] for f in new_files]
                TILES_L = np.unique(tiles_list)
        PERIOD_START = CONFIG['PERIOD_START'].split(',')[0]
        PERIOD_END = CONFIG['PERIOD_END'].split(',')[0]
        if PERIOD_START!='' and PERIOD_END!='':
            period_start = pd.to_datetime(PERIOD_START, format='%Y-%m-%d')
            period_end_inclusive = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)
            if period_start.month==period_end_inclusive.month:
                m_dt = pd.Timestamp(period_start)
            else:
                m_dt = pd.date_range(start=PERIOD_START, end=period_end_inclusive.strftime('%Y-%m-%d'), freq='M')
        else:
            m_dt = pd.date_range(start='2020-01-01', end='2021-01-01', freq='M')
        m_str = np.unique(m_dt.strftime("%m"))
        PROFILES_L = extractProfilesLike_Tiles(os.path.join(DATA_HISTO, 'DECADE'), TILES_L)
    
    # --- WRONG MODE ---
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')
    
    
    # ===================================== LOOP OVER LANDSAT TILES =================================

    for tile_L in tqdm(TILES_L, desc='TILES'):
        # logging.info(f'TILE LANDSAT : {tile_L}')    


        # ========================================== LOOP OVER MONTHS =================================
        
        for month in tqdm(m_str, desc='MONTHS'):
            # logging.info(f'MONTH : {month}')

            # Here, month vai is not estimated in the local proc chain to improve time processing (can be added if necessary)


        # ========================================== LOOP OVER DECADES =================================

            for d in range(3):
                logging.info(f'DECADE : {d+1}')

                ndwi_files_d = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*{tile_L}_*{month}_NDWI_COMPD{d+1}*.tif'))
                ndwi_files_dayd = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*{tile_L}_*{month}_NDWI_DAY{d+1}*.tif'))
                ndwi_files_day = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*{tile_L}_*{month}_NDWI_DAY.tif'))
                ndwi_files_m = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*{tile_L}_*{month}_NDWI_COMPM.tif'))
                
                ndwi_files = ndwi_files_d + ndwi_files_dayd + ndwi_files_day + ndwi_files_m

                # --- MISSING INDICES ---
                if ndwi_files==[]:
                    logging.warning('NO PRODUCTS -> PASS')
                    del ndwi_files_d, ndwi_files_dayd, ndwi_files_day, ndwi_files_m, ndwi_files
                    continue

                # --- SINGLE YEAR ---
                if len(ndwi_files)==1:
                    logging.warning('SINGLE YEAR -> PASS')
                    del ndwi_files_d, ndwi_files_dayd, ndwi_files_day, ndwi_files_m, ndwi_files
                    continue

                # --- COMPUTING DECADE VAI ---
                extractVAI(ndwi_files, PROFILES_L[tile_L], NDAYS_MAX_D, outdir_vai, f'D{d+1}')
                del ndwi_files_d, ndwi_files_dayd, ndwi_files_day, ndwi_files_m, ndwi_files

                # --- COPYING TO DATA_HISTO ---
                new_vai_d = glob.glob(os.path.join(outdir_vai,f'VAI_*{tile_L}_*{month}D{d+1}.tif'))
                for fd in tqdm(new_vai_d):
                    copyfile_Errorscontrol(fd, os.path.join(DATA_DROUGHT, 'DECADE', os.path.basename(fd)))
                del new_vai_d
 

