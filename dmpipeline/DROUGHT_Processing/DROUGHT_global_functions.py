# -*- coding: utf-8 -*-
"""
##############################################################################

DROUGHT functions for processing Global drought indicator

##############################################################################
"""

import os
import glob
import shutil
import numpy as np
import pandas as pd
import rasterio
from tqdm import tqdm
import dmpipeline.GEOSTATS_Processing.GEOSTATS_processing_functions as geostats

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def prepare_DATAFolders(CONFIG):
    """
    Create folders tree for historic data (DATA_HISTO) and controls annex directory (ANNEX_DIR) :
        - DATA_HISTO will contain all indices computed and needed for the calculation of drought indicators
        - ANNEX_DIR must be created and filled by user if geostatistics are intended to be computed (DROUGHT_STATS=1)
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
    TREE_DIR = [os.path.join(HISTO_DIR, '0_INDICES', 'MODIS', 'DECADE'),
                os.path.join(HISTO_DIR, '0_INDICES', 'MODIS', 'MONTH'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'GLOBAL', 'DECADE'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'GLOBAL', 'MONTH'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'GLOBAL', 'STATS')]
    os.umask(0) # used to reset the directories permission

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

    # --- ANNEX directory (IF DROUGHT_STATS) ---

    if (CONFIG['DROUGHT_STATS'] is None) or (CONFIG['DROUGHT_STATS']==''): go_stats = 0
    else: go_stats = int(CONFIG['DROUGHT_STATS'])

    if go_stats==1:
    
        ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
        TREE_DIR = [os.path.join(ANNEX_DIR, 'Areas'),
                    os.path.join(ANNEX_DIR, 'Landcover')]

        # Control that directories were well created and filled
        if not os.path.exists(ANNEX_DIR):
            list_dir = glob.glob(os.path.join(CONFIG['ANNEX_DIR'], '*'))
            
            if len(list_dir)==0:
                logging.critical(f'ANNEX directory was not found for {TERRITORY_str} : need to be filled for drought stats computation !\n-> Create/fill correctly or set DROUGHT_STATS=0 to avoid error\n')
                raise Exception('Error ANNEX directory')
            
            n_dir=0
            for dir in list_dir:
                dir_name = os.path.basename(dir)
                if dir_name in TERRITORY_str:
                    try:
                        os.rename(dir, ANNEX_DIR)
                        break
                    except OSError as e:
                        logging.critical(f'Error when renaming ANNEX folder: {e}')
                        raise Exception('Error ANNEX directory')
                n_dir+=1
                del dir_name

            if n_dir==len(list_dir):
                logging.critical(f'ANNEX directory was not found for {TERRITORY_str} : need to be filled for drought stats computation !\n-> Create/fill correctly or set DROUGHT_STATS=0 to avoid error\n')
                raise Exception('Error ANNEX directory')          
        
        # Browse and control folders list
        for dir in TREE_DIR:
            dir_name = os.path.basename(dir)

            # -> Mandatory 'Areas' folder (critical if not detected)
            if 'Areas' in dir_name :
                if not os.path.exists(dir):
                    logging.critical(f'In ANNEX, {dir_name} directory was not found : needed for drought stats computation !\n-> Create/fill correctly or set DROUGHT_STATS=0 to avoid error\n')
                    raise Exception('Error ANNEX directory')
                elif len(glob.glob(os.path.join(dir, '*')))==0:
                    logging.critical(f'In ANNEX, {dir_name} directory is empty : needed for drought stats computation !\n-> Create/fill correctly or set DROUGHT_STATS=0 to avoid error\n')
                    raise Exception('Error ANNEX directory')
                
            # -> Optional 'Landcover' folder (info if not detected)
            elif 'Landcover' in dir_name:
                if not os.path.exists(dir):
                    logging.info(f'In ANNEX, {dir_name} directory was not found : no landcover mask will be applied for drought stats\n')
                elif len(glob.glob(os.path.join(dir, '*')))==0:
                    logging.info(f'In ANNEX, {dir_name} directory is empty : no landcover mask will be applied for drought stats\n')
            
            # -> Optional other folders (info if not detected)
            else:
                if not os.path.exists(dir):
                    logging.info(f'In ANNEX, {dir_name} directory was not found\n')
                elif len(glob.glob(os.path.join(dir, '*')))==0:
                    logging.info(f'In ANNEX, {dir_name} directory is empty\n')
            
            del dir_name

    return WRK_DIR



def prepare_RUNFolders(CONFIG):
    """
    Prepare folders that will contain output data (rasters, csv, etc.) for a specific run
    Verify if the main folder exists :
            If no -> create folder (and sub-folders)
            If yes -> refers to CONFIG['CLEAN_RUNFOLDER'] to delete or keep it
    """
    
    # --- Prepare main folder name ---
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    date_start_str = CONFIG['PERIOD_START'].replace('-', '')
    date_end = pd.to_datetime(CONFIG['PERIOD_END'], format='%Y-%m-%d') + pd.DateOffset(days=-1)
    date_end_str = date_end.strftime('%Y%m%d')

    MODE = CONFIG['MODE']
    if (CONFIG['CLEAN_RUNFOLDER'] is None) or (CONFIG['CLEAN_RUNFOLDER']==''): CLEAN_RUNFOLDER = 0
    else: CLEAN_RUNFOLDER = int(CONFIG['CLEAN_RUNFOLDER'])
    
    if 'AUTO' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_AUTO_GLOBAL_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'MANUAL' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_MANUAL_GLOBAL_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'INDICES' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_INDICES_GLOBAL_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'DROUGHT' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_DROUGHT_GLOBAL_{TERRITORY_str}_{date_start_str}_{date_end_str}')
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
    outdir_tci = os.path.normpath(outdir_drought + os.sep + '0_TCI/')
    outdir_vci = os.path.normpath(outdir_drought + os.sep + '1_VCI/')
    outdir_vhi = os.path.normpath(outdir_drought + os.sep + '2_VHI/')
    outdir_droughtstats = os.path.normpath(outdir_drought + os.sep + 'STATS/')

    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):
        if not os.path.exists(outdir_comp): os.makedirs(outdir_comp)
        if not os.path.exists(outdir_compdecade): os.makedirs(outdir_compdecade)
        if not os.path.exists(outdir_compmonth): os.makedirs(outdir_compmonth)
        if not os.path.exists(outdir_compstats): os.makedirs(outdir_compstats)

    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('DROUGHT' in MODE):
        if not os.path.exists(outdir_drought): os.makedirs(outdir_drought)
        if not os.path.exists(outdir_tci): os.makedirs(outdir_tci)
        if not os.path.exists(outdir_vci): os.makedirs(outdir_vci)
        if not os.path.exists(outdir_vhi): os.makedirs(outdir_vhi)
        if not os.path.exists(outdir_droughtstats): os.makedirs(outdir_droughtstats)

    # --- Concenate dir paths ---
    OUTDIR_PATHS = (outdir_comp, outdir_compdecade, outdir_compmonth, outdir_compstats, outdir_tci, outdir_vci, outdir_vhi, outdir_droughtstats) 

    return OUTDIR_PATHS



def auto_processingPERIOD(CONFIG):
    """
    Automatic setting of processing period :
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
        if (PERIOD_START is None) or (PERIOD_START==''):
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
            logging.critical(
                f'One of PERIOD START or PERIOD END variable is None : both need to be filled in {MODE} mode !\n')
            raise Exception('Wrong PERIOD variable')
    
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')


    start_D = [1, 11, 21]
    end_D = [10, 20, 31]
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'MODIS', 'DECADE')
    histo_files = glob.glob(os.path.join(DATA_HISTO, '*.tif'))

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
                period_end_d_inclu = period_end_d + pd.DateOffset(days=-1)
                y_prev = period_end_d_inclu.year
                m_prev = period_end_d_inclu.month
                period_start_d = pd.to_datetime('{}-{:02d}-{:02d}'.format(y_prev,m_prev,start_D[d_prev]))
            break
    period_end = pd.Timestamp(period_end_d).floor('D')
    PERIOD_END = period_end.strftime('%Y-%m-%d')
    
    # IF HISTO PRODUCTS : Set new PERIOD START according to last histo product processed
    if len(histo_files)>0:
        dates_list = [os.path.basename(f).split('_')[2] for f in histo_files]
        dates = np.unique(dates_list)
        dates = pd.to_datetime(dates, format='%Y%m')
        period_start_histo = pd.Series(dates).min()
        period_end_histo = pd.Series(dates).max()
        PERIOD_START_HISTO = period_start_histo.strftime('%Y-%m-%d')
        PERIOD_END_HISTO = period_end_histo.strftime('%Y%m')
        
        list_end_histo =  glob.glob(os.path.join(DATA_HISTO, f'*{PERIOD_END_HISTO}*.tif'))
        D_END_HISTO = np.nanmax(np.unique([int(f.split('_COMPD')[1][0]) for f in list_end_histo]))
        d_number = int(D_END_HISTO)-1
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



def auto_processingBBOX(CONFIG):
    """
    Automatic setting of processing Bounding Box (lon_min, lat_min, lon_max, lat_max).
    IF one of lat/lon input var (config) is None or mode is 'AUTO' :
        - extract bbox of already processed data (DATA_HISTO)
        - if no histo data, critical message for manually setting the bbox (config)
    IF NOT : read lat/lon input var (manual setting) 
    """

    MODE = CONFIG['MODE']

    if (('AUTO' in MODE) or (CONFIG['lon_min_modis'] is None)  or (CONFIG['lon_min_modis']=='')
        or (CONFIG['lon_max_modis'] is None)  or (CONFIG['lon_max_modis']=='')
        or (CONFIG['lat_min_modis'] is None) or (CONFIG['lat_min_modis']=='')
        or (CONFIG['lat_max_modis'] is None) or (CONFIG['lat_max_modis']=='')): 
        
        TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
        DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'MODIS', 'DECADE')
        histo_files = glob.glob(os.path.join(DATA_HISTO, '*.tif'))

        if len(histo_files)>0:
            logging.info(f'\n\n --- AUTOMATIC SETTING OF BOUNDING-BOX ---\n')
            with rasterio.open(histo_files[0]) as d_ds :
                (CONFIG['lon_min_modis'], CONFIG['lat_min_modis'],
                 CONFIG['lon_max_modis'], CONFIG['lat_max_modis']) = d_ds.bounds
                transform_out = d_ds.transform
                transform_out = list(transform_out)[0:6]
        else:
            if ((CONFIG['lon_min_modis'] is None)  or (CONFIG['lon_min_modis']=='')
                or (CONFIG['lon_max_modis'] is None)  or (CONFIG['lon_max_modis']=='')
                or (CONFIG['lat_min_modis'] is None) or (CONFIG['lat_min_modis']=='')
                or (CONFIG['lat_max_modis'] is None) or (CONFIG['lat_max_modis']=='')):
                logging.critical(f'BOUNDING_BOX can not be automaticly defined : NO HISTO DATA and one of LON/LAT MODIS variables is None !\n -> Control data histo OR Fill correctly lat/lon to avoid error\n')
                raise Exception('Wrong PERIOD variable')
            else:
                logging.info(f'\n\n --- MANUAL SETTING OF BOUNDING-BOX ---\n')
                transform_out=None

    else:
        logging.info(f'\n\n --- MANUAL SETTING OF BOUNDING-BOX ---\n')
        transform_out=None
    
    return CONFIG['lon_min_modis'], CONFIG['lat_min_modis'], CONFIG['lon_max_modis'], CONFIG['lat_max_modis'], transform_out



def removeCurrentYear_DATAList(DATAfiles):
    """
    In list of DATA files, remove file corresponding to the current year.
    Used here to remove the current lst and ndwi modis indices before computing historical min/max for VHI calculation.
    """

    curr_year = pd.Timestamp.now().year
    DATAfiles_histo = []
    DATAfiles_current = []

    for f in DATAfiles:

        date_str = os.path.basename(f).split('_')[2]
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



def extractTCI(files, profile_like, ndays_max, outdir, period_indic=''):
    """
    Computing Temperature Condition Index from Kogan (1995) :
        - TCI = (LSTmax - LSTi) / (LSTmax - LSTmin)
        - QSCORE_comp = count_i/ndays_max 
        - QSCORE_tci = mean(QSCORE_comp,QSCORE_histo)
    
    Note 1 :  The composite Quality Score (QSCORE_comp) is computed for each specific date,
              from the ratio of valid pixels counted per indice (count_i) divided by the maximum
              number of days possibly observed on the period (here 1 image/day).
              A FINAL TCI Quality Score (QSCORE_tci) is computed taking into account the composite score
              QSCORE_comp and the seasonal QSCORES mean on the period observed (QSCORE_histo).

    Note 2 : LST of current year are not included in the computation of historical LSTmax,LSTmin,QSCORE_histo.
    """

    logging.info('Apply TCI equation')

    # --- COMPUTING HISTORICAL MIN/MAX and QSCORE ---

    # Extracting lst of current year
    files_histo, files_current = removeCurrentYear_DATAList(files)
    N_files = len(files)
    N_files_histo = len(files_histo)
    LST_MEAN = np.full((profile_like['height'],profile_like['width'],N_files), np.nan, 'float32')
    LST_COUNT = np.full_like(LST_MEAN, np.nan, 'float32')

    for i in range(N_files_histo):
        with rasterio.open(files_histo[i]) as d_ds:
            LST_MEAN[:,:,i] = d_ds.read(1)
            if d_ds.count==2:
                LST_COUNT[:,:,i] = d_ds.read(2)
            elif d_ds.count==1:
                LST_COUNT[:,:,i] = ~np.isnan(LST_MEAN[:,:,i])

    LST_QSCORE = LST_COUNT[:,:,:i]/ndays_max
    MIN_histo = np.nanmin(LST_MEAN,axis=2)
    MAX_histo = np.nanmax(LST_MEAN,axis=2)
    QSCORE_histo = np.mean(LST_QSCORE,axis=2)

    # Here adding lst of current year
    for f in files_current:
        i+=1
        with rasterio.open(f) as d_ds:
            LST_MEAN[:,:,i] = d_ds.read(1)
            if d_ds.count==2:
                LST_COUNT[:,:,i] = d_ds.read(2)
            elif d_ds.count==1:
                LST_COUNT[:,:,i] = ~np.isnan(LST_MEAN[:,:,i])
    LST_QSCORE = LST_COUNT/ndays_max

    # --- APPLY TCI EQUATION TO EACH IMAGE ON THE PERIOD ---
    for i in tqdm(range(N_files)):
        
        # Compute TCI
        TCI = (MAX_histo - LST_MEAN[:,:,i])/(MAX_histo - MIN_histo)
        TCI[(MAX_histo == MIN_histo)] = np.nan
        TCI_QSCORE = np.mean(np.stack((QSCORE_histo,LST_QSCORE[:,:,i]),axis=2),axis=2)
        TCI_QSCORE[LST_QSCORE[:,:,i]==0] = 0
        
        # Save file
        in_file_name = os.path.basename(files[i]).split('.tif')[0]
        date = in_file_name.split('_')[2]
        out_file_name = os.path.join(outdir, f'TCI_{date}{period_indic}.tif')
        with rasterio.open(out_file_name, 'w', **profile_like) as out_ds:
            out_ds.write(TCI, 1)
            out_ds.write(TCI_QSCORE, 2)
        
        del TCI, TCI_QSCORE, out_file_name



def extractVCI(files, profile_like, ndays_max, outdir, period_indic=''):
    """
    Computing Vegetation Condition Index from Kogan (1995, 2000) :
        - VCI = (NDWIi - NDWImin) / (NDWImax - NDWImin)
        - QSCORE_comp = count_i/ndays_max 
        - QSCORE_vci = mean(QSCORE_comp,QSCORE_histo)
    
    Note 1 :  The composite Quality Score (QSCORE_comp) is computed for each specific date,
              from the ratio of valid pixels counted per indice (count_i) divided by the maximum
              number of days possibly observed on the period (here 1 image/day).
              A FINAL VCI Quality Score (QSCORE_vci) is computed taking into account the composite score
              QSCORE_comp and the seasonal QSCORES mean on the period observed (QSCORE_histo).

    Note 2 : NDWI of current year are not included in the computation of historical NDWImax,NDWImin,QSCORE_histo.
    """

    logging.info('Apply VCI equation')

    # --- COMPUTING HISTORICAL MIN/MAX and QSCORE ---

    # Extracting ndwi of current year
    files_histo, files_current = removeCurrentYear_DATAList(files)
    N_files = len(files)
    N_files_histo = len(files_histo)
    NDWI_MEAN = np.full((profile_like['height'],profile_like['width'],N_files), np.nan, 'float32')
    NDWI_COUNT = np.full_like(NDWI_MEAN, np.nan, 'float32')

    for i in range(N_files_histo):
        with rasterio.open(files_histo[i]) as d_ds:
            NDWI_MEAN[:,:,i] = d_ds.read(1)
            if d_ds.count==2:
                NDWI_COUNT[:,:,i] = d_ds.read(2)
            elif d_ds.count==1:
                NDWI_COUNT[:,:,i] = ~np.isnan(NDWI_MEAN[:,:,i])

    NDWI_QSCORE = NDWI_COUNT[:,:,:i]/ndays_max
    MIN_histo = np.nanmin(NDWI_MEAN,axis=2)
    MAX_histo = np.nanmax(NDWI_MEAN,axis=2)
    QSCORE_histo = np.mean(NDWI_QSCORE,axis=2)

    # Here adding ndwi of current year
    for f in files_current:
        i+=1
        with rasterio.open(f) as d_ds:
            NDWI_MEAN[:,:,i] = d_ds.read(1)
            if d_ds.count==2:
                NDWI_COUNT[:,:,i] = d_ds.read(2)
            elif d_ds.count==1:
                NDWI_COUNT[:,:,i] = ~np.isnan(NDWI_MEAN[:,:,i])
    NDWI_QSCORE = NDWI_COUNT/ndays_max

    # --- APPLY VCI EQUATION TO EACH IMAGE ON THE PERIOD ---
    for i in tqdm(range(N_files)):
        
        # Compute VCI
        VCI = (NDWI_MEAN[:,:,i] - MIN_histo)/(MAX_histo - MIN_histo)
        VCI[(MAX_histo == MIN_histo)] = np.nan
        VCI_QSCORE = np.mean(np.stack((QSCORE_histo,NDWI_QSCORE[:,:,i]),axis=2),axis=2)
        VCI_QSCORE[NDWI_QSCORE[:,:,i]==0] = 0
        
        # Save file
        in_file_name = os.path.basename(files[i]).split('.tif')[0]
        date = in_file_name.split('_')[2]
        out_file_name = os.path.join(outdir, f'VCI_{date}{period_indic}.tif')
        with rasterio.open(out_file_name, 'w', **profile_like) as out_ds:
            out_ds.write(VCI, 1)
            out_ds.write(VCI_QSCORE, 2)
        
        del VCI, VCI_QSCORE, out_file_name



def extractVHI(period_date, period_indic, profile_like, outdir_tci, outdir_vci, outdir_vhi, go_stats=0, outdir_droughtstats=None, annex_dir=None, territory=None, areas_key=None):
    """
    Computing Vegetation Health Index from Kogan (1997):
        - VHI = alpha * VCI + (1-alpha) * TCI
        - Qscore = mean(Qscore_vci,Qscore_tci)
    alpha value was set to 0.5 according to literature.
    IF go_stats=1, spatial statistics are also estimated on all territory and each sub-area (predefined in input masks).
    """
  
    ALPHA = 0.5

    logging.info('Apply VHI equation')

    vci_files = glob.glob(os.path.join(outdir_vci,f'VCI_*{period_date}{period_indic}.tif'))
    tci_files = glob.glob(os.path.join(outdir_tci,f'TCI_*{period_date}{period_indic}.tif'))
    if len(vci_files)==0:
        logging.critical(f'No input VCI products')
        raise Exception('No input VCI products')
    if len(tci_files)==0:
        logging.critical(f'No input TCI products')
        raise Exception('No input TCI products')

    # --- IF go_stats=1, Prepare input masks and look-up table (for estimating geostats) ---
    if go_stats==1:

        mask_areas_ok = len(glob.glob(os.path.join(outdir_droughtstats,'mask_Areas.tif')))>0
        if mask_areas_ok==0:
            file_areas = glob.glob(os.path.join(annex_dir, 'Areas', '*.shp'))
            file_areas_ok = len(file_areas)
            if file_areas_ok==0:
                logging.warning('Drought spatial stats will not be estimated : missing input shapefile containing geometries/areas to identify')
                go_stats=0
            else:
                file_areas = file_areas[0]
        
        mask_landcover_ok = len(glob.glob(os.path.join(outdir_droughtstats,'mask_Areas_*.tif')))>0    
        if mask_landcover_ok==0:
            file_landcover = glob.glob(os.path.join(annex_dir, 'Landcover', '*.tif'))
            file_landcover_ok = len(file_landcover)
            if file_landcover_ok==0:
                logging.info('Landcover masks will not be applied : missing input landcover file containing classes to mask')
                go_landcover=0
            else:
                file_landcover = file_landcover[0]
                go_landcover=1
        else:
            file_landcover=[]
            go_landcover=1
        
        if mask_areas_ok==0 and file_areas_ok==1:
            geostats.prepareGeoStatsMasks(vci_files[0], file_areas, outdir_droughtstats, file_landcover=file_landcover, areas_key=areas_key)
    
    else:
        logging.info('Drought spatial stats will not be estimated')
    
    if go_stats==1:
        logging.info('Drought spatial stats will be estimated')

        with rasterio.open(glob.glob(os.path.join(outdir_droughtstats,'mask_Areas.tif'))[0]) as area_ds:
            maskAREA = area_ds.read(1)
            mask = (maskAREA != 0)
        area_lut = pd.read_csv(glob.glob(os.path.join(outdir_droughtstats,'*.csv'))[0], sep=';')
        if go_landcover==1:
            with rasterio.open(glob.glob(os.path.join(outdir_droughtstats,'mask_Areas_NOTrees_NOBuild.tif'))[0]) as NOTrees_ds, \
                rasterio.open(glob.glob(os.path.join(outdir_droughtstats,'mask_Areas_Trees.tif'))[0]) as Trees_ds :
                mask_NOTrees_NOBuild = NOTrees_ds.read(1)
                mask_Trees = Trees_ds.read(1)

        # --- AND Verify if output stats dataframes already exist (yes -> do not add hearder in csv file) ---
        if period_indic=='M':
            stats_ok = len(glob.glob(os.path.join(outdir_droughtstats,'VHI_STATS_M*.csv')))>0
        else:
            stats_ok = len(glob.glob(os.path.join(outdir_droughtstats,'VHI_STATS_D*.csv')))>0
    
    
    count_head = 0

    for vci_f in tqdm(vci_files):

        # --- Read VCI/TCI files ---
        in_file_name = os.path.basename(vci_f).split('.tif')[0]
        full_date = in_file_name.split('_')[1]
        tci_f = glob.glob(os.path.join(outdir_tci,f'TCI_*{full_date}.tif'))
        
        if tci_f==[]:
            logging.warning(f'Missing TCI file for {full_date}')
            continue
        else:
            tci_f = tci_f[0]
        
        with rasterio.open(vci_f) as vci_ds, \
            rasterio.open(tci_f) as tci_ds:
            
            VCI = vci_ds.read(1)
            VCI_QSCORE = vci_ds.read(2)
            TCI = tci_ds.read(1)
            TCI_QSCORE = tci_ds.read(2)
        
        # --- Compute VHI ---
        VHI = ALPHA * VCI + (1-ALPHA) * TCI
        del VCI, TCI
        
        # --- Compute combined Qscore ---
        VHI_QSCORE = np.mean(np.stack((VCI_QSCORE,TCI_QSCORE),axis=2),axis=2)
        VHI_QSCORE[(VCI_QSCORE==0) | (TCI_QSCORE==0)] = 0
        del VCI_QSCORE, TCI_QSCORE
        
        # --- Save file ---
        out_file_name = os.path.join(outdir_vhi, f'VHI_{full_date}.tif')
        with rasterio.open(out_file_name, 'w', **profile_like) as out_ds:
            out_ds.write(VHI, 1)
            out_ds.write(VHI_QSCORE, 2)

        # --- IF go_stats=1, Estimate spatial stats ---
        if go_stats==1:
            if period_indic=='M': 
                date_df = pd.to_datetime(full_date[:-1], format='%Y%m')
                period_df = period_indic
            elif period_indic=='D1':
                date_df = pd.to_datetime(full_date[:-2]+'01', format='%Y%m%d')
                period_df = 'D'
            elif period_indic=='D2':
                date_df = pd.to_datetime(full_date[:-2]+'11', format='%Y%m%d')
                period_df = 'D'
            elif period_indic=='D3':
                date_df = pd.to_datetime(full_date[:-2]+'21', format='%Y%m%d')
                period_df = 'D'

            if go_landcover==1:
                GeoStats_df, GeoStats_df_NOTrees_NOBuild, GeoStats_df_Trees = geostats.extractGeoStats(VHI, VHI_QSCORE, date_df, mask, maskAREA,
                                                                                                       area_lut, territory, mask_NOTrees_NOBuild, mask_Trees)
                GeoStats_df_NOTrees_NOBuild = GeoStats_df_NOTrees_NOBuild.sort_values(by=['LOCATION','DATE'])
                GeoStats_df_Trees = GeoStats_df_Trees.sort_values(by=['LOCATION','DATE'])

                GeoStats_df_NOTrees_NOBuild.to_csv(
                    os.path.join(outdir_droughtstats, f'VHI_STATS_{period_df}_NoTrees_NoBuild.csv'),
                    index = False,
                    float_format='%.2f',
                    decimal = '.',
                    sep = ';',
                    mode='a',
                    header = (count_head==0 and stats_ok==0))
                GeoStats_df_Trees.to_csv(
                    os.path.join(outdir_droughtstats, f'VHI_STATS_{period_df}_Trees.csv'),
                    index = False,
                    float_format='%.2f',
                    decimal = '.',
                    sep = ';',
                    mode='a',
                    header = (count_head==0 and stats_ok==0))
                
                del GeoStats_df_NOTrees_NOBuild, GeoStats_df_Trees

            else:
                GeoStats_df, _, _ = geostats.extractGeoStats(VHI, VHI_QSCORE, date_df, mask, maskAREA, area_lut, territory)
            
            GeoStats_df = GeoStats_df.sort_values(by=['LOCATION','DATE'])

            GeoStats_df.to_csv(os.path.join(outdir_droughtstats, f'VHI_STATS_{period_df}.csv'),
                index = False,
                float_format='%.2f',
                decimal = '.',
                sep = ';',
                mode='a',
                header = (count_head==0 and stats_ok==0))
            
            count_head += 1
            
            del GeoStats_df, period_df, date_df
        
        del VHI, VHI_QSCORE, tci_f, in_file_name, full_date



def process_GlobalDrought_VHI(CONFIG, OUTDIR_PATHS):
    """
    Vegetation Health Index (VHI) is estimated from the combination of :
        - Temperature Condition Index TCI = f(LST)
        - Vegetation Condition Index VCI = f(NDWI)
    A quality score is given (Qscore) according to the number valid pixels,
    and the maximum number of days possibly oberved on the compositing period (Ndays_max) :
        - NDAYS_MAX_M (month) is set to 30
        - NDAYS_MAX_D (decade) is set to 10
    IF CONFIG['DROUGHT_STATS']=1, spatial statistics are also estimated on all territory and each sub-area.
    
    Final products and statistics are saved (updated) into data_histo directory.
    
    Note :  finals products are given at NDWI 500 m resolution
            2 bands are saved b1=VHI, b2=Qscore
    """

    logging.info('\n\n--- PROCESSING GLOBAL DROUGHT INDICATOR (VHI, 500 m) ---\n')

    NDAYS_MAX_M = 30
    NDAYS_MAX_D = 10
    (outdir_comp, outdir_compdecade, outdir_compmonth, outdir_compstats,
     outdir_tci, outdir_vci, outdir_vhi, outdir_droughtstats) = OUTDIR_PATHS
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'MODIS')
    DATA_DROUGHT = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'GLOBAL')
    MODE = CONFIG['MODE']
    if (CONFIG['DROUGHT_STATS'] is None) or (CONFIG['DROUGHT_STATS']==''): go_stats = 0
    else: go_stats = int(CONFIG['DROUGHT_STATS'])
    
    # --- AUTOMATIC/MANUAL/INDICES MODE -> process VHI to corresponding new indices months/decades ---
    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):
        
        logging.info(f'{MODE} MODE : process VHI to corresponding new indices months/decades')
        
        new_files_m = glob.glob(os.path.join(outdir_compmonth, '*.tif'))
        new_files_d = glob.glob(os.path.join(outdir_compdecade, '*.tif'))
        go_month = len(new_files_m)>0

        m_list = [os.path.basename(f).split('_')[2][-2:] for f in new_files_d]
        m_str = np.unique(m_list)
        d_list = [os.path.basename(f).split('_')[3][-5:-4] for f in new_files_d]
        d_str = np.unique(d_list)

    # --- DROUGHT MODE -> compute/update all vhi products from data historic directory ---
    elif 'DROUGHT' in MODE:

        logging.info(f'{MODE} MODE : compute/update VHI products from data historic directory')

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
        d_str = ['1','2','3']
        go_month = 1
    
    # --- WRONG MODE ---
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')

    ndwi_like = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*NDWI*.tif'))[0]
    lst_like = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*LST*.tif'))[0]
    with rasterio.open(ndwi_like) as ndwi_ds, \
        rasterio.open(lst_like) as lst_ds:
        profile_ndwilike = ndwi_ds.profile
        profile_lstlike = lst_ds.profile


    # ========================================== LOOP OVER MONTHS =================================
    
    for month in m_str:
        logging.info(f'MONTH : {month}')

        if go_month==1:
            ndwi_files_m = glob.glob(os.path.join(DATA_HISTO, 'MONTH', f'*NDWI_*{month}_COMPM.tif'))
            lst_files_m = glob.glob(os.path.join(DATA_HISTO, 'MONTH', f'*LST_*{month}_COMPM.tif'))

            # --- MISSING INDICES ---
            if ndwi_files_m==[] or lst_files_m==[]:
                logging.warning('NO PRODUCTS -> PASS')
                continue

            # --- SINGLE YEAR ---
            elif len(ndwi_files_m)==1 or len(lst_files_m)==1:
                logging.warning('SINGLE YEAR -> PASS')
                continue
            
            # --- COMPUTING MONTH TCI and VCI ---
            extractTCI(lst_files_m, profile_lstlike, NDAYS_MAX_M, outdir_tci, 'M')
            extractVCI(ndwi_files_m, profile_ndwilike, NDAYS_MAX_M, outdir_vci, 'M')
            del ndwi_files_m, lst_files_m

            # --- COMPUTING MONTH VHI ---
            extractVHI(month, 'M', profile_ndwilike, outdir_tci, outdir_vci, outdir_vhi, go_stats, outdir_droughtstats, ANNEX_DIR, CONFIG['TERRITORY'], CONFIG['KEY_STATS'])

            # --- COPYING TO DATA_HISTO ---
            new_vhi_m = glob.glob(os.path.join(outdir_vhi,f'VHI_*{month}M.tif'))
            for fm in tqdm(new_vhi_m):
                copyfile_Errorscontrol(fm, os.path.join(DATA_DROUGHT, 'MONTH', os.path.basename(fm)))
            del new_vhi_m


        else : logging.info('Month VHI not processed')


    # ========================================== LOOP OVER DECADES =================================

        for d in d_str:
            logging.info(f'DECADE : {d}')

            ndwi_files_d = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*NDWI_*{month}_COMPD{d}.tif'))
            lst_files_d = glob.glob(os.path.join(DATA_HISTO, 'DECADE', f'*LST_*{month}_COMPD{d}.tif'))

            # --- MISSING INDICES ---
            if ndwi_files_d==[] or lst_files_d==[]:
                logging.warning('NO PRODUCTS -> PASS')
                continue
                
            # --- SINGLE YEAR ---
            elif len(ndwi_files_d)==1 or len(lst_files_d)==1:
                logging.warning('SINGLE YEAR -> PASS')
                continue
            
            # --- COMPUTING DECADE TCI and VCI ---
            extractTCI(lst_files_d, profile_lstlike, NDAYS_MAX_D, outdir_tci, f'D{d}')
            extractVCI(ndwi_files_d, profile_ndwilike, NDAYS_MAX_D, outdir_vci, f'D{d}')
            del ndwi_files_d, lst_files_d

            # --- COMPUTING DECADE VHI ---
            extractVHI(month, f'D{d}', profile_ndwilike, outdir_tci, outdir_vci, outdir_vhi, go_stats, outdir_droughtstats, ANNEX_DIR, CONFIG['TERRITORY'], CONFIG['KEY_STATS'])

            # --- COPYING TO DATA_HISTO ---
            new_vhi_d = glob.glob(os.path.join(outdir_vhi,f'VHI_*{month}D{d}.tif'))
            for fd in tqdm(new_vhi_d):
                copyfile_Errorscontrol(fd, os.path.join(DATA_DROUGHT, 'DECADE', os.path.basename(fd)))
            del new_vhi_d
    

    # --- Copy/Update VHI Geo Statistics to data histo directory ---
    if go_stats==1:
        logging.info('Copy/Update VHI Geo Statistics to data histo directory')
        DATA_STATS = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'GLOBAL', 'STATS')
        GEOSTATS_DF = glob.glob(os.path.join(outdir_droughtstats,'VHI_STATS*.csv'))

        for f in GEOSTATS_DF:
            logging.info(f'Dataframe : {os.path.basename(f)}')
            f_datahisto = os.path.join(DATA_STATS, os.path.basename(f))

            if os.path.exists(f_datahisto):
                # Update : replace new stats to histo stats
                df_new = pd.read_csv(f, sep=';')
                try:
                    df_new['DATE'] = pd.to_datetime(df_new.DATE, format='%Y-%m-%d')
                except ValueError:
                    df_new['DATE'] = pd.to_datetime(df_new.DATE, format='%d/%m/%Y')
                df_histo = pd.read_csv(f_datahisto, sep=';')
                try:
                    df_histo['DATE'] = pd.to_datetime(df_histo.DATE, format='%Y-%m-%d')
                except ValueError:
                    df_histo['DATE'] = pd.to_datetime(df_histo.DATE, format='%d/%m/%Y')

                for i in tqdm(range(df_new.shape[0])):
                    location = df_new['LOCATION'][i]
                    date = df_new['DATE'][i]
                    mask_histo = (df_histo['LOCATION']==location) & (df_histo['DATE']==date)
                    mask_new = (df_new['LOCATION']==location) & (df_new['DATE']==date)
                    if mask_histo.any():
                        df_tmp = df_new.loc[mask_new].copy()
                        df_histo.loc[mask_histo] = df_tmp.values
                        del df_tmp
                    else:
                        df_histo = pd.concat([df_histo, df_new.loc[mask_new]], ignore_index=True)
                    del mask_histo, mask_new, location, date

                try:
                    df_histo.to_csv(f_datahisto,
                                    index = False,
                                    float_format='%.2f',
                                    decimal = '.',
                                    sep = ';')
                except PermissionError:
                    os.remove(f_datahisto)
                    df_histo.to_csv(f_datahisto,
                                    index = False,
                                    float_format='%.2f',
                                    decimal = '.',
                                    sep = ';')

                del df_histo, df_new

            else:
                # Copy stats to data histo :
                shutil.copyfile(f, f_datahisto)
            del f_datahisto
                   

