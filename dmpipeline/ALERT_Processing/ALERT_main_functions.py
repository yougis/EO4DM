# -*- coding: utf-8 -*-
"""
#############################################################################################

ALERT functions for processing and combining satellites and meteorological drought indicators

#############################################################################################
"""

import os
import shutil
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats
import rasterio
import rasterio.mask
import copy
from zipfile import BadZipFile, ZipFile
from pathlib import Path
import fnmatch
from tqdm import tqdm
import dmpipeline.GEOSTATS_Processing.GEOSTATS_processing_functions as geostats

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def filterPERIOD_DATA(DATAfiles, DATAType, period_start, period_end):
    """
    Filter list of DATA files acccording to input period start and period end
    """

    DATAfiles_filt = []
    dates_filt = []

    for f in DATAfiles:

        if DATAType=='ASCAT': date_str = os.path.basename(f).split('_')[1][:8]
        elif DATAType=='VHI': date_str = os.path.basename(f).split('_')[1][:-5]+'01'
        else: logging.critical('Wrong inptu data type in filtering period function')
        
        try:
            date = pd.to_datetime(date_str, format='%Y%m%d')
        except ValueError:
            logging.info(f'Filtering period DATA : Wrong date format for {date_str}, then not considered')
            continue

        if period_start <= date <= period_end:
            DATAfiles_filt.append(f)
            dates_filt.append(date)

    return DATAfiles_filt, dates_filt



def auto_processingPERIOD(CONFIG):
    """
    Automatic settint of processing period :
        - extract period of already processed data (DATA_HISTO)
        - find last month to be possibly ready for processing (the month before the current decade)
        - give new period to process
        - if last month was already processed => NO PROCESSING (go_process=0)
    
    Note : if mode (in config) is not set to 'AUTO' -> go out of function and return go_process=1.
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
            logging.critical(f'One of PERIOD START or PERIOD END variable is None : both need to be filled in {MODE} mode !\n')
            raise Exception('Wrong PERIOD variable')
    
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')


    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'ASCAT', 'MONTH')
    histo_files = glob.glob(os.path.join(DATA_HISTO, '*.tif'))

    # Set new PERIOD END = end of the previous month (PERIOD_END is exclusive)
    date_now = pd.Timestamp.now()
    period_end = date_now + pd.DateOffset(days=-date_now.day+1)
    period_end = period_end.normalize()
    PERIOD_END = period_end.strftime('%Y-%m-%d')

    # IF HISTO PRODUCTS : Set new PERIOD START according to last histo product processed
    if len(histo_files)>0:
        dates_list = [os.path.basename(f).split('_')[2] for f in histo_files]
        dates = np.unique(dates_list)
        dates = pd.to_datetime(dates, format='%Y%m')
        period_start_histo = pd.Series(dates).min()
        period_end_histo = pd.Series(dates).max()
        PERIOD_START_HISTO = period_start_histo.strftime('%Y-%m')
        PERIOD_END_HISTO = period_end_histo.strftime('%Y-%m')
        logging.info(f'Histo data were already processed from {PERIOD_START_HISTO} to {PERIOD_END_HISTO}')

        if 'DROUGHT' in MODE:
            period_end_histo = period_end_histo + pd.DateOffset(months=1)  # period end is exclusive
            PERIOD_END_HISTO =  period_end_histo.strftime('%Y-%m-%d')
            go_process = 1
            return go_process, PERIOD_START_HISTO, PERIOD_END_HISTO
        
        period_start = period_end_histo + pd.DateOffset(months=1)
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

    # CONTROLS if last month was already processed or not :
    if period_start<period_end:
        logging.info(f'GO PROCESSING : from {PERIOD_START} to {PERIOD_END}')
        go_process = 1
    else:
        logging.info(f'NO PROCESSING : Previous month ({PERIOD_END_HISTO}) already processed')
        go_process = 0

    return go_process, PERIOD_START, PERIOD_END



def Check_Products_ASCAT(CONFIG):
    """
    Function to check ASCAT products availability before launching drought processing :    
        - Extract collections on PERIOD (from config) and check if products are available
        - If max waiting time is exceeded (NB_WAIT_MAX) -> GO PROCESSING (last month not controled)
        - If not, extracts products on LAST MONTH and compare number of available products to number of expected products
        - If last month is not full :
            -> case of multiple months : remove last month and GO PROCESSING
            -> case of single month :  stop and wait more data
    """
    
    logging.info('\n\n --- CHECKING ASCAT PRODUCTS AVAILABILITY ---\n')

    NB_WAIT_MAX = 15 # Waits max 15 days

    # --- Read input CONFIG ---
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_RAW = os.path.join(CONFIG['DATA_HISTO'],TERRITORY_str,'0_INDICES','ASCAT','0_RAW')
    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']

    # --- Extracts ASCAT on the specific PERIOD and CHECK IF PRODUCTS ARE AVAILABLE ---
    ASCATdataset = glob.glob(os.path.join(DATA_RAW,'SWI_*'))
    if len(ASCATdataset)==0:
        logging.critical('\nNO ASCAT PRODUCTS AVAILABLE\n')
        raise Exception('NO ASCAT PRODUCTS AVAILABLE')
    
    if PERIOD_START==[''] or PERIOD_START=='First product':
        period_start = pd.to_datetime('2007-01-01', format='%Y-%m-%d')
    else:
        period_start = pd.to_datetime(PERIOD_START, format='%Y-%m-%d')

    if PERIOD_END==['']:
        period_end = pd.Timestamp.now()
    else:
        period_end = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)

    ASCATdataset_filt, dates_filt = filterPERIOD_DATA(ASCATdataset, 'ASCAT', period_start, period_end)
    
    Nb_ASCAT = len(ASCATdataset_filt)
    if Nb_ASCAT==0:
        logging.info('NO PRODUCTS -> STOP PROCESSING')
        go_ascat = 0
        DATES = ('','')
        return go_ascat, ASCATdataset_filt, '', ''
    
    period_start = pd.Series(dates_filt).min()
    period_end = pd.Series(dates_filt).max()
    period_start_str = period_start.strftime('%Y-%m-%d')
    period_end_str_excl = (period_end+pd.DateOffset(days=1)).strftime('%Y-%m-%d')
    NB_MONTH = len(np.unique(pd.Series(dates_filt).dt.to_period('M')))

    # --- Extracts products on LAST MONTH and COMPARE number of available products to number of expected products ---
    # -> If waiting time exceeded : Do not control last decade
    date_now = pd.Timestamp.now()
    Nb_wait = date_now - period_end
    if Nb_wait.days > NB_WAIT_MAX:
        logging.info(f'WAITING TIME {Nb_wait.days} days greater than {NB_WAIT_MAX} days -> GO PROCESSING (last month not controled)')
        go_ascat = 1
    
    else:
        # -> Control ASCAT:
        period_start_m = period_end + pd.DateOffset(days=-period_end.day+1)
        period_end_m = period_end
        ASCATdataset_filt_m, dates_filt_m = filterPERIOD_DATA(ASCATdataset_filt, 'ASCAT', period_start_m, period_end_m)

        Nb_ASCAT_m = len(ASCATdataset_filt_m)
        Nb_EXPECT_images = period_end.days_in_month
        period_end_m_str = period_end_m.strftime('%Y-%m')

        # Case WITH expected number of products -> PROCESS FULL PERIOD
        if Nb_ASCAT_m==Nb_EXPECT_images:
            logging.info(f'LAST MONTH {period_end_m_str} : EXPECTED NUMBER OF PRODUCTS -> GO PROCESSING FOR FULL PERIOD')
            go_ascat = 1

        else:
            # Case WITHOUT expected number of products and SINGLE MONTH to process -> STOP and WAIT
            if NB_MONTH==1:
                logging.info(f'SINGLE MONTH {period_end_m_str} : MISSING PRODUCTS -> STOP PROCESSING AND WAIT')
                go_ascat = 0

            # Case WITHOUT expected number of products and MULTIPLE MONTHS to process -> SETS PERIOD END TO PREVIOUS MONTH and GO PROCESSING
            elif NB_MONTH>1:
                period_start_m = period_start_m + pd.DateOffset(months=-1)
                period_end_m = period_start_m + pd.DateOffset(days=period_start_m.days_in_month-1)
                period_end_str = period_end_m.strftime('%Y-%m-%d')
                period_end_str_excl = (period_end_m+pd.DateOffset(days=1)).strftime('%Y-%m-%d')
                logging.info(f'LAST MONTH {period_end_m_str} : MISSING PRODUCTS -> REDUCE PERIOD FROM {period_start_str} TO {period_end_str} and GO PROCESSING')
                ASCATdataset_filt, dates_filt = filterPERIOD_DATA(ASCATdataset_filt, 'ASCAT', period_start, period_end_m)
                go_ascat = 1

    return go_ascat, ASCATdataset_filt, period_start_str, period_end_str_excl



def Control_Data_Meteo(DATAType, DATA_HISTO_METEO, DATA_ANNEX, PERIOD_START, PERIOD_END, NB_WAIT_MAX):
    """
    Function to check METEO (SPI or SPEI) products availability before launching drought processing :    
        - Extract monthly data on PERIOD and check if products are available
        - Controls that all expected stations are available
        - Per station, controls that number of products is equal to number of months expected in the period
        - If expected number -> GO PROCESSING
        - If not :
            -> case where max waiting time IS NOT exceeded (NB_WAIT_MAX): stop and waits more data
            -> case where max waiting time IS EXCEEDED (NB_WAIT_MAX) : Warning message, and GO PROCESSING
    """

    if DATAType=='SPI': data_val = 'SPI3_MENS'
    elif DATAType=='SPEI': data_val = 'SPEI_3'
    else: logging.critical('Wrong inptu data type in check meteo data function')

    # --- Extracts METEO on the specific PERIOD and CHECK IF PRODUCTS ARE AVAILABLE ---
    data_csv = os.path.join(DATA_HISTO_METEO, f'{DATAType}_ref_1991_2020.csv')
    if data_csv==[]:
        logging.critical(f'\nNO {DATAType} PRODUCTS AVAILABLE\n')
        raise Exception(f'NO {DATAType} PRODUCTS AVAILABLE')
    data_df = pd.read_csv(data_csv,sep=';',decimal=',')
    data_df = data_df[['NOM','DATE',data_val]]
    data_df['DATE'] = pd.to_datetime(data_df.DATE, format='%Y%m')

    if PERIOD_START==[''] or PERIOD_START=='First product':
        period_start = pd.to_datetime('2007-01-01', format='%Y-%m-%d')
    else:
        period_start = pd.to_datetime(PERIOD_START, format='%Y-%m-%d')
    if PERIOD_END==['']:
        period_end = pd.Timestamp.now()
    else:
        period_end = pd.to_datetime(PERIOD_END, format='%Y-%m-%d')

    data_filt_df = data_df[(data_df['DATE']>=period_start) & (data_df['DATE']<=period_end)].reset_index(drop=True)
    
    # --- Read expected stations list ---
    stations_csv = os.path.join(DATA_ANNEX, 'Stations', f'{DATAType}_communes_stations.csv')
    if stations_csv==[]:
        logging.critical(f'\nNO METEO STATIONS LIST ({DATAType}) AVAILABLE\n')
        raise Exception('NO METEO STATIONS LIST AVAILABLE')
    stations_df = pd.read_csv(stations_csv,sep=';')

    # --- HANDLE SPEI ERROR for Station HOUAILOU P -> replace by "HOUAILOU P." ---
    if DATAType=='SPEI':
        data_filt_df.loc[data_filt_df['NOM']=='HOUAILOU P', 'NOM'] = 'HOUAILOU P.'

    # --- COMPARE number of available products to number of expected products ---
    # -> Control all expected stations from input stations list:
    stations_expect = pd.concat([stations_df[colonne].dropna() for colonne in stations_df.columns[1:]], axis=0, ignore_index=True).unique()
    stations_filt = data_filt_df['NOM'].unique()
    miss_stations = []
    for station in stations_expect:
        if station not in stations_filt:
            miss_stations.append(station)

    if miss_stations!=[]:
        logging.critical(f'\nMISSING {DATAType} STATION(S) : {miss_stations}\n')
        raise Exception('MISSING METEO STATION(S)')
    
    else:
        # -> Per station, control number of months :
        Nb_EXPECT_months = (period_end.year-period_start.year)*12 + (period_end.month-period_start.month+1)
        MISS_MONTHS_STATIONS = [len(np.unique(data_filt_df.loc[data_filt_df['NOM']==s, 'DATE'].dt.to_period('M')))<Nb_EXPECT_months for s in stations_filt]

        # Case WITH expected number of months -> GO PROCESSING
        if np.sum(MISS_MONTHS_STATIONS)==0:
            logging.info(f'{DATAType} : EXPECTED NUMBER OF MONTHS PER STATION -> GO PROCESSING')
            go_meteo = 1

        elif np.sum(MISS_MONTHS_STATIONS)>0:
            # Case WITHOUT expected number of products -> STOP PROCESSING and WAIT
            missmonth_stations = stations_filt[MISS_MONTHS_STATIONS]
            date_now = pd.Timestamp.now()
            Nb_wait = date_now - period_end
            logging.info(f'\nMISSING {DATAType} MONTH(S) ON STATION(S) : {missmonth_stations}\n')
            if Nb_wait.days > NB_WAIT_MAX:
                logging.warning(f'WAITING TIME {Nb_wait.days} days greater than {NB_WAIT_MAX} days -> GO PROCESSING')
                go_meteo = 1
            else:
                logging.warning(f'-> STOP PROCESSING')
                go_meteo = 0

    return go_meteo



def Check_Products_METEO(CONFIG):
    """
    Function to check METEO products availability before launching drought processing :    
        - Controls SPI : expected stations, expected months per station
        - Controls SPEI : expected stations, expected months per station
    """
    
    logging.info('\n\n --- CHECKING METEO PRODUCTS AVAILABILITY ---\n')

    NB_WAIT_MAX = 15 # Waits max 15 days

    # --- Read input CONFIG ---
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO_METEO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'ALERT', 'METEO')
    DATA_ANNEX = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']

    # --- CHECK SPI PRODUCTS ---
    go_spi = Control_Data_Meteo('SPI', DATA_HISTO_METEO, DATA_ANNEX, PERIOD_START, PERIOD_END, NB_WAIT_MAX)

    # --- CHECK SPEI PRODUCTS ---
    go_spei = Control_Data_Meteo('SPEI', DATA_HISTO_METEO, DATA_ANNEX, PERIOD_START, PERIOD_END, NB_WAIT_MAX)

    return go_spi and go_spei



def Check_Products_VHI(CONFIG):
    """
    Function to check VHI products availability before launching drought processing :    
        - Extract monthly collection on PERIOD (from config) and check if products are available
        - Controls that number of products is equal to number of months expected in the period
        - If expected number -> go_processing
        - If not :
            If missing first month or month(s) inside the period -> GO PROCESSING (but warning message)
            If missing last month -> STOP PROCESSING AND WAIT
    """
    
    logging.info('\n\n --- CHECKING VHI PRODUCTS AVAILABILITY ---\n')

    NB_WAIT_MAX = 15 # Waits max 15 days

    # --- Read input CONFIG ---
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO_VHI = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'GLOBAL', 'MONTH')
    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']

    # --- Extracts VHI on the specific PERIOD and CHECK IF PRODUCTS ARE AVAILABLE ---
    VHIdataset = glob.glob(os.path.join(DATA_HISTO_VHI,'VHI_*'))
    if len(VHIdataset)==0:
        logging.critical('\nNO VHI PRODUCTS AVAILABLE\n')
        raise Exception('NO VHI PRODUCTS AVAILABLE')
    
    if PERIOD_START==[''] or PERIOD_START=='First product':
        period_start = pd.to_datetime('2007-01-01', format='%Y-%m-%d')
    else:
        period_start = pd.to_datetime(PERIOD_START, format='%Y-%m-%d')

    if PERIOD_END==['']:
        period_end = pd.Timestamp.now()
    else:
        period_end = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)

    VHIdataset_filt, dates_filt = filterPERIOD_DATA(VHIdataset, 'VHI', period_start, period_end)
    
    Nb_VHI = len(VHIdataset_filt)
    if Nb_VHI==0:
        logging.info('NO PRODUCTS -> STOP PROCESSING')
        go_vhi = 0
        return go_vhi
    
    period_start_vhi = pd.Series(dates_filt).min()
    period_end_vhi = pd.Series(dates_filt).max()
    Nb_EXPECT_images = (period_end.year-period_start.year)*12 + (period_end.month-period_start.month+1)

    # Case WITH expected number of products -> GO PROCESSING
    if Nb_VHI==Nb_EXPECT_images:
        logging.info('EXPECTED NUMBER OF PRODUCTS -> GO PROCESSING')
        go_vhi = 1
    
    elif Nb_VHI<Nb_EXPECT_images:
        # Case WITHOUT expected number of products on PERIOD START -> STOP PROCESSING
        if period_start_vhi!=period_start:
            logging.warning('MISSING PRODUCTS ON PERIOD START -> CONTROL COLLECTION but GO PROCESSING')
            go_vhi = 1
            return go_vhi
        
        # Case WITHOUT expected number of products on PERIOD END -> STOP PROCESSING and WAIT
        if period_end_vhi!=period_end:
            date_now = pd.Timestamp.now()
            Nb_wait = date_now - period_end
            logging.info('MISSING PRODUCTS ON PERIOD END -> STOP PROCESSING AND WAIT')
            if Nb_wait.days > NB_WAIT_MAX:
                logging.warning(f'WAITING TIME {Nb_wait.days} days greater than {NB_WAIT_MAX} days')
            go_vhi = 0
            return go_vhi

        # Case WITHOUT expected number of products INSIDE PERIOD -> GO PROCESSING
        logging.warning(f'MISSING PRODUCTS (INSIDE PERIOD) -> CONTROL COLLECTION but GO PROCESSING')
        go_vhi = 1
    
    return go_vhi



def CheckInit_Products_ALERT(CONFIG):
    """
    CHECKING/PREPARING input products for alert processing chain :
        - ASCAT
        - METEO (SPI/SPEI)
        - VHI
    The chain is launched only if all products are availble for the period to process (from config).
    If it is the case, output folders for current running are prepared, and a dataframe with ascat informations is exported. 
    """
    
    # --- Controls products availability and updates PERIOD START/END ---
    go_ascat, COLLECTION_ASCAT, PERIOD_START_ASCAT, PERIOD_END_ASCAT = Check_Products_ASCAT(CONFIG)
    PERIOD_START = CONFIG['PERIOD_START']
    PERIOD_END = CONFIG['PERIOD_END']
    go_meteo = Check_Products_METEO(CONFIG)
    go_vhi = Check_Products_VHI(CONFIG)

    # --- IF all products available ---
    if go_ascat==1 and go_meteo==1 and go_vhi==1:

        # -> Prepare output Folders :
        OUTDIR_PATHS = prepare_RUNFolders(CONFIG)

        # -> Export ascat csv (nb products, dates) ---
        PERIOD_START_str = PERIOD_START.replace('-','')
        period_end = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)
        PERIOD_END_str = period_end.strftime('%Y%m%d')

        period_end_ascat = pd.to_datetime(PERIOD_END_ASCAT, format='%Y-%m-%d') + pd.DateOffset(days=-1)
        period_end_ascat_str = period_end_ascat.strftime('%Y-%m-%d')

        TERRITORY = CONFIG['TERRITORY']
        TERRITORY_str = TERRITORY.replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
        Count_dict = []
        Count_filename = f'COUNTCollection_ASCAT_{TERRITORY_str}_{PERIOD_START_str}_{PERIOD_END_str}.csv'
        
        Count_dict = pd.DataFrame(data={'PRODUCT':'SWI ASCAT',
                                        'DATE START': PERIOD_START_ASCAT,
                                        'DATE END': period_end_ascat_str,
                                        'NUMBER OF IMAGES': len(COLLECTION_ASCAT)},
                                        index=[1])
        Count_dict.to_csv(
            os.path.join(OUTDIR_PATHS[3], Count_filename),
            index = False,
            decimal = '.',
            sep = ',')
        
        # -> Set go_products to 1 :
        go_products = 1

    else:
        go_products = 0
        OUTDIR_PATHS = ''
    
    
    return go_products, COLLECTION_ASCAT, OUTDIR_PATHS



def prepare_DATAFolders(CONFIG):
    """
    Create folders tree for historic data (DATA_HISTO) and controls annex directory (ANNEX_DIR) :
        - DATA_HISTO will contain all indices computed and needed for the calculation of alerts
        - ANNEX_DIR must be created and filled by user
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
    TREE_DIR = [os.path.join(HISTO_DIR, '0_INDICES', 'ASCAT', '0_RAW'),
                os.path.join(HISTO_DIR, '0_INDICES', 'ASCAT', 'DAY'),
                os.path.join(HISTO_DIR, '0_INDICES', 'ASCAT', 'MONTH'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'ALERT', 'MAI', 'MONTH'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'ALERT', 'MAI', 'STATS'),
                os.path.join(HISTO_DIR, '1_INDICATEURS', 'ALERT', 'METEO')]
    os.umask(0) # used to reset the directories permission

    try:
        # Create the main dir if not exist
        os.makedirs(HISTO_DIR, exist_ok=True)
        os.chmod(HISTO_DIR, 0o777)

        # Browse folders list and create new ones if not already exist
        for dir in TREE_DIR:
            os.makedirs(dir, exist_ok=True)

    except OSError as e:
        logging.info(f'Error when creating DATA FOLDERS: {e}')
    
    del TREE_DIR


    # --- ANNEX directory ---
    
    ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    TREE_DIR = [os.path.join(ANNEX_DIR, 'Areas'),
                os.path.join(ANNEX_DIR, 'Stations')]

    # Control that directories were well created and filled
    if not os.path.exists(ANNEX_DIR):
        logging.critical(f'ANNEX directory was not found for {TERRITORY_str} : need to be filled for ALERT computation !\n')
        raise Exception('Error ANNEX directory')
    
    # Browse and control folders list
    for dir in TREE_DIR:
        if not os.path.exists(dir):
            logging.critical(f'In ANNEX, {os.path.basename(dir)} directory was not found : needed for ALERT computation !\n')
            raise Exception('Error ANNEX directory')
        elif len(glob.glob(os.path.join(dir, '*')))==0:
                logging.critical(f'In ANNEX, {os.path.basename(dir)} directory is empty : needed for ALERT computation !\n')
                raise Exception('Error ANNEX directory')
    
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
    if 'AUTO' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_AUTO_ALERT_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'MANUAL' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_MANUAL_ALERT_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'INDICES' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_INDICES_ALERT_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    elif 'DROUGHT' in MODE:
        outdir = os.path.join(CONFIG['WRK_DIR'], f'RUN_DROUGHT_ALERT_{TERRITORY_str}_{date_start_str}_{date_end_str}')
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')
    
    # --- Generate directory ---
    os.umask(0) # used to reset the directories permission
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        os.chmod(outdir, 0o777)
    elif int(CONFIG['CLEAN_RUNFOLDER'])==1:
        shutil.rmtree(outdir)
        os.makedirs(outdir)
        os.chmod(outdir, 0o777)

    # --- Generate sub-directories ---
    outdir_moisture = os.path.normpath(outdir + os.sep + '1_MOISTURE/')
    outdirmoisture_preproc = os.path.normpath(outdir_moisture + os.sep + '1_PREPROC/')
    outdirmoisture_prepocswi = os.path.normpath(outdirmoisture_preproc + os.sep + 'SWI/')
    outdirmoisture_preprocstats = os.path.normpath(outdirmoisture_preproc + os.sep + 'STATS/')
    outdirmoisture_comp = os.path.normpath(outdir_moisture + os.sep + '2_COMPOSITE/')
    outdirmoisture_compmonth = os.path.normpath(outdirmoisture_comp + os.sep + 'MONTH/')
    outdirmoisture_drought = os.path.normpath(outdir_moisture + os.sep + '3_DROUGHT/')
    outdirmoisture_mai = os.path.normpath(outdirmoisture_drought + os.sep + f'MAI/')
    outdirmoisture_droughtstats = os.path.normpath(outdirmoisture_drought + os.sep + f'STATS/')
    outdir_alert = os.path.normpath(outdir + os.sep + '2_ALERT/')
    outdir_maskareas = os.path.normpath(outdir + os.sep + 'MASK_AREAS/')

    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):
        if not os.path.exists(outdir_moisture): os.makedirs(outdir_moisture)
        if not os.path.exists(outdirmoisture_preproc): os.makedirs(outdirmoisture_preproc)
        if not os.path.exists(outdirmoisture_prepocswi): os.makedirs(outdirmoisture_prepocswi)
        if not os.path.exists(outdirmoisture_preprocstats): os.makedirs(outdirmoisture_preprocstats)
        if not os.path.exists(outdirmoisture_comp): os.makedirs(outdirmoisture_comp)
        if not os.path.exists(outdirmoisture_compmonth): os.makedirs(outdirmoisture_compmonth)
        if not os.path.exists(outdirmoisture_drought): os.makedirs(outdirmoisture_drought)
        if not os.path.exists(outdirmoisture_mai): os.makedirs(outdirmoisture_mai)
        if not os.path.exists(outdirmoisture_droughtstats): os.makedirs(outdirmoisture_droughtstats)

    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('DROUGHT' in MODE):
        if not os.path.exists(outdir_alert): os.makedirs(outdir_alert)
    
    if not os.path.exists(outdir_maskareas): os.makedirs(outdir_maskareas)

    # --- Concenate dir paths into an output parameter ---
    OUTDIR_PATHS = (outdir_moisture, outdirmoisture_preproc, outdirmoisture_prepocswi, 
                    outdirmoisture_preprocstats, outdirmoisture_comp, outdirmoisture_compmonth,
                    outdirmoisture_drought, outdirmoisture_mai, outdirmoisture_droughtstats,
                    outdir_alert, outdir_maskareas)

    return OUTDIR_PATHS



def getSWI020FileFromInput(path):
    """
    From input zip folder, extract path to Soil Water Index file for time length T=20 (*SWI-020_*.tiff)
    """

    if isinstance(path, str):
        path = Path(path)
    if path.is_dir():
        file_list = [el.relative_to(path.parents[0]) for el in path.rglob("*")]
        prefix = f"{str(path.parents[0])}/"
    elif path.is_file() and path.suffix == ".zip":
        try:
            with ZipFile(path) as zf:
                file_list = [el for el in zf.namelist()]
                prefix =  f"zip://{path}!"
            
            swi_file = [el for el in file_list if fnmatch.fnmatch(el, "*_SWI-SWI-020*.tiff")]
            
            if len(swi_file) > 0:
                return {f"SWI-{s.split('_')[2][-3:]}": f"{prefix}{s}" for s in swi_file}
        
        except BadZipFile:
            logging.info(f"!! BadZipFile !! :\n  {path}\n")   
            return "BadZipFileError"
    else:
        logging.info(f"Error file {path} does not exists")
        return "FileNotFoundError"



def getMaskFilesFromInput(path):
    """
    From input zip folder, extract path to :
        - Quality Flag files (*QFLAG_*.tiff)
        - Surface State Flag file (*SSF_*.tiff)
    """

    if isinstance(path, str):
        path = Path(path)
    if path.is_dir():
        file_list = [el.relative_to(path.parents[0]) for el in path.rglob("*")]
        prefix = f"{str(path.parents[0])}/"
    elif path.is_file() and path.suffix == ".zip":
        try:
            with ZipFile(path) as zf:
                file_list = [el for el in zf.namelist()]
                prefix = f"zip://{path}!"
            
            qflag_file = [el for el in file_list if fnmatch.fnmatch(el, "*_SWI-QFLAG-020*.tiff")]
            ssf_file = [el for el in file_list if fnmatch.fnmatch(el, "*_SWI-SSF_*.tiff")]
            
            masks = {}

            if len(qflag_file) > 0:
                masks["QFLAG"] = f"{prefix}{qflag_file[0]}"

            if len(ssf_file) > 0:
                masks["SSF"] = f"{prefix}{ssf_file[0]}"
        
            return masks
                
        except BadZipFile:
            logging.info(f"!! BadZipFile !! :\n  {path}\n")   
   
    else:
        logging.info(f"Error file {path} does not exists")



def extractSWI(file):
    """
    Extracting ASCAT Soil Water Index :
        - rescaling
        - setting no data
    """

    zf_path = str(file['SWI-020']).split('!')[0][6:]
    tf_path = str(file['SWI-020']).split('!')[1]

    with ZipFile(zf_path) as zf:
        with zf.open(tf_path) as tf:
            with rasterio.open(tf, GEOREF_SOURCES='INTERNAL') as swi_ds:
                profile_out = swi_ds.profile
                
                # Rescaling
                swi_dn = swi_ds.read(1).astype(np.float32)
                swi_phys = swi_dn * 0.005
                    
                # Setting no data
                nodata_dn = 255
                nodata_phys = np.nan
                swi_phys[swi_dn==nodata_dn] = nodata_phys
                
                profile_out.update(dtype=rasterio.float32, count=1, nodata=nodata_phys)
    
    return swi_phys, profile_out



def extractMask(masks, maskName):
    """
    Extracting ASCAT Masks :
        - rescaling
        - setting no data
    """
     
    if maskName == "QFLAG":
        if masks is not None:
            qflag = masks['QFLAG']
            zf_path = qflag.split('!')[0][6:]
            tf_path = qflag.split('!')[1]
            with ZipFile(zf_path) as zf:
                with zf.open(tf_path) as tf:
                    with rasterio.open(tf) as mask_ds:
                        profile_out = mask_ds.profile
                        # Rescaling
                        m_dn = mask_ds.read(1).astype(np.float32)
                        m = m_dn * 0.005
                        # Setting no data
                        nodata = 255
                        m[m_dn==nodata] = nodata
        else:
            logging.info("QFLAG extraction needs mask param.")
            raise Exception()
    
    elif maskName == "SSF":
        if masks is not None:
            ssf = masks['SSF']
            zf_path = ssf.split('!')[0][6:]
            tf_path = ssf.split('!')[1]
            with ZipFile(zf_path) as zf:
                with zf.open(tf_path) as tf:
                    with rasterio.open(tf) as mask_ds:
                        profile_out = mask_ds.profile
                        m = mask_ds.read(1)
        else:
            logging.info("SSF extraction needs masks param.")
            raise Exception()
    
    else:
        logging.info(f"Mask name ({maskName}) not valid.")
        raise Exception()

    return m, profile_out



def extractComposite(mat_in):
    """
    Temporal compositing on the different images in mat_in:
        - average is computed (MEAN_comp)
        - count is also computed to count number of valid data used per pixel (COUNT_comp)
    """
    
    COUNT = np.zeros_like(mat_in, dtype='uint8')
    COUNT = ~np.isnan(mat_in)
    
    MEAN_comp = np.nanmean(mat_in,axis=2)
    COUNT_comp = np.sum(COUNT,axis=2)
    
    return MEAN_comp, COUNT_comp



def processingASCAT(CONFIG, OUTDIR_PATHS, COLLECTION):
    """
    Preprocessing and Compositing of ASCAT Soil Water Index products (SWI):
        - quality masking
        - estimating scores over land (qflag, ssflag)
        - month compositing
    
    IN: folder containing for each date :
            SWI files for different depths T (*SWI-001_*.tiff, ..., *SWI-100_*.tiff)
            Quality Flag files fo different depths T (*QFLAG_*.tiff)
            Surface State Flag file (*SSF_*.tiff)

    OUT: - Preprocessed daily ASCAT SWI products (*_PREPROC.tif)
         - Scores for each daily product (PREPROC_ASCAT_SWI_*.csv)
         - Composited ASCAT SWI product (*_COMPM.tif)
    
    Note: Output composite product have 2 bands (b1=average data, b2=count)
    """
    
    logging.info('\n\n--- PREPROCESSING and COMPOSITING ASCAT ---\n')
    
    date_start_str = CONFIG['PERIOD_START'].replace('-','')
    date_end = pd.to_datetime(CONFIG['PERIOD_END'], format='%Y-%m-%d') + pd.DateOffset(days=-1) 
    date_end_str = date_end.strftime('%Y%m%d')
    KEY_STATS = CONFIG['KEY_STATS']

    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    SWIFolders = COLLECTION

    QA_dict = pd.DataFrame(columns=['FILE NAME','DATE','NAN SCORE SWI-020','QFLAG-020','SSFLAG'])
    QAtable_filename = f'PREPROC_ASCAT_{TERRITORY_str}_{date_start_str}_{date_end_str}.csv'
    
    SWIFile_like = glob.glob(os.path.join(SWIFolders[0],'*.zip'))[0]
    swi_like = getSWI020FileFromInput(SWIFile_like)
    SWI_data, profile_like = extractSWI(swi_like)
    profile_out = copy.deepcopy(profile_like)
    profile_out.update(count=2)

    y_list = [os.path.basename(f).split('_')[1][:4] for f in SWIFolders]
    y_str = np.unique(y_list)
    m_list = [os.path.basename(f).split('_')[1][4:6] for f in SWIFolders]
    m_str = np.unique(m_list)

    # --- Prepare input mask (for estimating nanscores on land) ---
    mask_ok = len(glob.glob(os.path.join(OUTDIR_PATHS[-1],'mask_Areas_MAI.tif')))>0
    if mask_ok==0:
        file_areas = glob.glob(os.path.join(ANNEX_DIR, 'Areas', '*.shp'))
        if len(file_areas)==0:
            logging.critical('Missing input landmask shapefile for nanscores estimation')
            raise Exception ('Missing input landmask shapefile for nanscores estimation')
        file_areas = file_areas[0]
        geostats.prepareGeoStatsMasks(swi_like['SWI-020'], file_areas, OUTDIR_PATHS[-1], suffix='MAI', areas_key=KEY_STATS)

    with rasterio.open(glob.glob(os.path.join(OUTDIR_PATHS[-1],'mask_Areas_MAI.tif'))[0]) as area_ds:
        maskAREA = area_ds.read(1)
        mask = (maskAREA != 0)
    Nb_ALLDATA_land = np.sum(mask)


    # ========================================== LOOP OVER YEARS/MONTHS =================================

    for year in y_str:
        for month in m_str:
            month2find = f'{year}{month}'
            logging.info(f'MONTH TO FIND : {month2find}')

            SWIFolders_m = [el for el in SWIFolders if fnmatch.fnmatch(el, f'*SWI_{month2find}*')]
            NFolders_m = len(SWIFolders_m)

            if NFolders_m==0:
                logging.info('NO PRODUCTS')
                continue

            # Preparing matrices for compositing
            MEAN_swi = np.full((profile_like['height'],profile_like['width']), np.nan, 'float32')
            DATA = np.full((profile_like['height'],profile_like['width'],NFolders_m), np.nan,'float32')


            # ===================================== LOOP OVER DAILY PRODUCTS ============================

            for i in tqdm(range(NFolders_m)):
                
                SWIFile = glob.glob(os.path.join(SWIFolders_m[i],'*.zip'))[0]
                in_file_name = os.path.basename(SWIFile).split('.zip')[0]
                date = in_file_name.split('_')[3][:8]
                swi = getSWI020FileFromInput(SWIFile)
                
                # If error : Save file name and go to next iteration    
                if swi=="BadZipFileError":
                    with open(os.path.join(OUTDIR_PATHS[3],f'Badzipfiles_error_ASCAT_{date_start_str}_{date_end_str}.txt'), 'a') as f:
                        f.write(f'{in_file_name}\n')
                    continue
                if swi=="FileNotFoundError":
                    with open(os.path.join(OUTDIR_PATHS[3],f'FileNotFound_error_ASCAT_{date_start_str}_{date_end_str}.txt'), 'a') as f:
                        f.write(f'{in_file_name}\n')
                    continue
                Qmasks = getMaskFilesFromInput(SWIFile)
            
                # --- Reading SWI and Masks, and Estimating scores (land) ---
                SWI_data, profile_swi = extractSWI(swi)
                QF_data, _ = extractMask(Qmasks, "QFLAG")
                SSF_data, _ = extractMask(Qmasks, "SSF")
                
                # QFLAG mask : 50 % of threshold for ssm
                QFmask = (QF_data <= 0.5)
                Qf_score = np.sum(QFmask[mask==1]==1) / Nb_ALLDATA_land
                
                # SSF mask : Unknown=0, Frozen=2, Water=3, NotDetermined=255
                SSFmask = (SSF_data==0) | (SSF_data==2) | (SSF_data==3) | (SSF_data==255)
                ssf_score = np.sum(SSFmask[mask==1]==1) / Nb_ALLDATA_land
            
                # --- Applying Quality Masks and Land Mask on SWI ---
                QMask = (QFmask==1) | (SSFmask==1)
                SWI_Qmasked = SWI_data
                SWI_Qmasked[(QMask==1) | (mask==0)] = np.nan
                del SWI_data,QMask,QFmask,SSFmask
                
                # --- Estimating Percentage of No Data AFTER Quality Masking (only on LAND) ---
                Nb_NODATA = np.sum(np.isnan(SWI_Qmasked[mask==1]))
                NAN_score = Nb_NODATA / Nb_ALLDATA_land
                
                # --- Saving Preprocessed SWI and Quality scores IF NANSCORE < 100 % ---
                if NAN_score<1 :
                    
                    # Saving files
                    swi_file_name = os.path.basename(swi['SWI-020']).split('.tiff')[0]
                    date_old = swi_file_name.split('_')[3]
                    out_file_name = os.path.join(OUTDIR_PATHS[2], f'{swi_file_name.replace(date_old,date)}_PREPROC.tif')
                    with rasterio.open(out_file_name, 'w', **profile_swi) as out_ds:
                        out_ds.write(SWI_Qmasked, 1)
                    del out_file_name, swi_file_name, date_old
                    
                    d = {'FILE NAME':in_file_name,'DATE':pd.to_datetime(date, format='%Y%m%d'),
                         'NAN SCORE SWI-020':NAN_score,'QFLAG-020':Qf_score, 'SSFLAG':ssf_score}
                    QA_dict = pd.concat([QA_dict, pd.DataFrame([d])], ignore_index=True)
                    del d
                    
                    # Concatenating matrices
                    DATA[:,:,i] = SWI_Qmasked
                    del SWI_Qmasked
                
                del in_file_name, date, NAN_score, Qf_score, ssf_score
            
            QA_dict.to_csv(
                os.path.join(OUTDIR_PATHS[3], QAtable_filename),
                index = False,
                float_format='%.2f',
                decimal = '.',
                sep = ',')
        
            # --- Compositing SWI and computing Geostatistics ---
            MEAN_swi, COUNT_swi = extractComposite(DATA)
            outputcomp_swi_file = os.path.join(OUTDIR_PATHS[5], f'ASCAT_SWI_{month2find}_COMPM.tif')
            with rasterio.open(outputcomp_swi_file, 'w', **profile_out) as out_ds:
                out_ds.write(MEAN_swi, 1)
                out_ds.write(COUNT_swi, 2)



def extractMAI(files, profile_like, ndays_max, outdir_mai, period_indic='', go_stats=0, mask_dir=None, outdir_droughtstats=None, annex_dir=None, territory=None, areas_key=None):
    """
    Computing Moisture Anomaly Index from Amri (2012):
        - MAI = (SWIi - SWImean) / SWIstd
        - QSCORE_comp = count_i/ndays_max 
        - QSCORE_mai = mean(QSCORE_comp,QSCORE_histo)
    
    Note :  The composite Quality Score (QSCORE_comp) is computed for each specific date,
            from the ratio of valid pixels counted per indice (count_i) divided by the maximum
            number of days possibly observed on the period (here 1 image/day).
            A FINAL MAI Quality Score (QSCORE_mai) is computed taking into account the composite score
            QSCORE_comp and the seasonal QSCORES mean on the period observed (QSCORE_histo)
    """

    logging.info('Apply MAI equation')

    # --- IF go_stats=1, Check input masks and look-up table (for estimating geostats) ---
    if go_stats==1:
        mask_ok = len(glob.glob(os.path.join(mask_dir,'mask_Areas_MAI.tif')))>0
        if mask_ok==0:
            file_areas = glob.glob(os.path.join(annex_dir, 'Areas', '*.shp'))
            if len(file_areas)==0:
                logging.info('MAI spatial stats will not be estimated : missing input shapefile containing geometries/areas to identify')
                go_stats=0
            else:
                file_areas = file_areas[0]
                geostats.prepareGeoStatsMasks(files[0], file_areas, mask_dir, suffix='MAI', areas_key=areas_key)
    else:
        logging.info('MAI spatial stats will not be estimated')
    
    if go_stats==1:
        logging.info('MAI spatial stats will be estimated')

        with rasterio.open(glob.glob(os.path.join(mask_dir,'mask_Areas_MAI.tif'))[0]) as area_ds:
            maskAREA = area_ds.read(1)
            mask = (maskAREA != 0)
        area_lut = pd.read_csv(glob.glob(os.path.join(mask_dir,'ID_Name_Areas-lookup_MAI.csv'))[0], sep=';')

        # --- AND Verify if output stats dataframes already exist (yes -> do not add hearder in csv file) ---
        stats_ok = len(glob.glob(os.path.join(outdir_droughtstats,'MAI_STATS_M*.csv')))>0
    count_head = 0

    # --- COMPUTING HISTORICAL MEAN/STD and QSCORE ---
    N_files = len(files)
    SWI_MEAN = np.full((profile_like['height'],profile_like['width'],N_files), np.nan, 'float32')
    SWI_COUNT = np.full_like(SWI_MEAN, np.nan, 'float32')

    for i in range(N_files):
        with rasterio.open(files[i]) as d_ds:
            SWI_MEAN[:,:,i] = d_ds.read(1)
            SWI_COUNT[:,:,i] = d_ds.read(2)

    SWI_QSCORE = SWI_COUNT/ndays_max
    MEAN_histo = np.nanmean(SWI_MEAN,axis=2)
    STD_histo = np.nanstd(SWI_MEAN,axis=2)
    QSCORE_histo = np.mean(SWI_QSCORE,axis=2)

    # --- APPLY MAI EQUATION TO EACH IMAGE ON THE PERIOD ---
    for i in tqdm(range(N_files)):
        
        # Compute MAI
        MAI = (SWI_MEAN[:,:,i] - MEAN_histo)/STD_histo
        MAI[(STD_histo == 0)] = np.nan
        MAI_QSCORE = np.mean(np.stack((QSCORE_histo,SWI_QSCORE[:,:,i]),axis=2),axis=2)
        MAI_QSCORE[SWI_QSCORE[:,:,i]==0] = 0
        
        # Save file
        in_file_name = os.path.basename(files[i]).split('.tif')[0]
        full_date = in_file_name.split('_')[2]
        out_file_name = os.path.join(outdir_mai, f'MAI_{full_date}{period_indic}.tif')
        with rasterio.open(out_file_name, 'w', **profile_like) as out_ds:
            out_ds.write(MAI, 1)
            out_ds.write(MAI_QSCORE, 2)
        
        # IF go_stats=1, Estimate spatial stats
        if go_stats==1:
            date_df = pd.to_datetime(full_date, format='%Y%m')
            GeoStats_df, _, _ = geostats.extractGeoStats(MAI, MAI_QSCORE, date_df, mask, maskAREA, area_lut, territory)
            GeoStats_df = GeoStats_df.sort_values(by=['LOCATION','DATE'])
            
            GeoStats_df.to_csv(os.path.join(outdir_droughtstats, f'MAI_STATS_{period_indic}.csv'),
                index = False,
                float_format='%.2f',
                decimal = '.',
                sep = ';',
                mode='a',
                header = (count_head==0 and stats_ok==0))

            count_head += 1
            del GeoStats_df, date_df

        del MAI, MAI_QSCORE, out_file_name, full_date



def process_MoistureDeficit_MAI(CONFIG, OUTDIR_PATHS):
    '''
    Moisture Anomaly Index (MAI) is estimated for each month :
        - MAI = (SWIi - SWImean) / SWIstd
        - Qscore = mean(Qscore_i,Qscore_histo)
    A quality score is given (Qscore) according to the number valid pixels,
    and the maximum number of days possibly oberved on the compositing period (NDAYS_MAX_M=30)
        
    Final products and statistics are saved (updated) into date_histo directory.
    
    Note : 2 bands are saved b1=MAI, b2=Qscore
    '''
    
    logging.info('\n\n--- PROCESSING MOISTURE DEFICIT (MAI, 12 km) ---\n')

    NDAYS_MAX_M = 30
    (outdir_moisture, outdirmoisture_preproc, outdirmoisture_prepocswi, 
     outdirmoisture_preprocstats, outdirmoisture_comp, outdirmoisture_compmonth,
     outdirmoisture_drought, outdirmoisture_mai, outdirmoisture_droughtstats, outdir_alert, outdir_maskareas) = OUTDIR_PATHS
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    ANNEX_DIR = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    DATA_HISTO = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '0_INDICES', 'ASCAT')
    MODE = CONFIG['MODE']
    if (CONFIG['DROUGHT_STATS'] is None) or (CONFIG['DROUGHT_STATS']==''): go_stats = 0
    else: go_stats = int(CONFIG['DROUGHT_STATS'])
    
    # --- AUTOMATIC/MANUAL/INDICES MODE -> copy the new indices composite file(s) to data histo directory and extract the new months to process ---
    if ('AUTO' in MODE) or ('MANUAL' in MODE) or ('INDICES' in MODE):

        logging.info(f'{MODE} MODE : copy new indices to data historic directory and process corresponding MAI')
                
        new_files_m = glob.glob(os.path.join(outdirmoisture_compmonth, '*.tif'))
        go_month = len(new_files_m)>0
        if go_month==1:
            logging.info('Month Composites copying')
            for fm in tqdm(new_files_m):
                fm_datahisto = os.path.join(DATA_HISTO, 'MONTH', os.path.basename(fm))
                try:
                    shutil.copyfile(fm, fm_datahisto)
                except PermissionError:
                    if os.path.exists(fm_datahisto):
                        logging.info(f'File already exists: {os.path.basename(fm)} is replaced by new version')
                        os.remove(fm_datahisto)
                        shutil.copyfile(fm, fm_datahisto)
                    else:
                        logging.critical(f'Copy PermissionError : {os.path.basename(fm)} impossible to paste')
                        raise Exception('Copy PermissionError : impossible to paste')
                del fm_datahisto
        else:
            logging.critical('Month composites products not found')
            raise Exception ('Month composites products not found')
        m_list = [os.path.basename(f).split('_')[2][-2:] for f in new_files_m]
        m_str = np.unique(m_list)

        new_files_d = glob.glob(os.path.join(outdirmoisture_prepocswi, '*.tif'))
        logging.info('Day indices copying')
        for fd in tqdm(new_files_d):
            fd_datahisto = os.path.join(DATA_HISTO, 'DAY', os.path.basename(fd))
            try:
                shutil.copyfile(fd, fd_datahisto)
            except PermissionError:
                if os.path.exists(fd_datahisto):
                    logging.info(f'File already exists: {os.path.basename(fd)} is replaced by new version')
                    os.remove(fd_datahisto)
                    shutil.copyfile(fd, fd_datahisto)
                else:
                    logging.critical(f'Copy PermissionError : {os.path.basename(fd)} impossible to paste')
                    raise Exception('Copy PermissionError : impossible to paste')
            del fd_datahisto

    # --- DROUGHT MODE -> compute/update all mai products from data historic directory ---
    elif 'DROUGHT' in MODE:
        logging.info(f'{MODE} MODE : compute/update MAI products from data historic directory')
        PERIOD_START = CONFIG['PERIOD_START'].split(',')[0]
        PERIOD_END = CONFIG['PERIOD_END'].split(',')[0]
        if PERIOD_START!='' and PERIOD_END!='':
            period_end_inclusive = pd.to_datetime(PERIOD_END, format='%Y-%m-%d') + pd.DateOffset(days=-1)
            m_dt = pd.date_range(start=PERIOD_START, end=period_end_inclusive.strftime('%Y-%m-%d'), freq='M')
        else:
            m_dt = pd.date_range(start='2020-01-01', end='2021-01-01', freq='M')
        m_str = np.unique(m_dt.strftime("%m"))
    
    # --- WRONG MODE ---
    else:
        logging.critical(f'Wrong imput for processing mode : {MODE}')
        raise Exception('Wrong PROCESSING MODE')

    swi_like = glob.glob(os.path.join(DATA_HISTO, 'MONTH', f'*SWI*.tif'))[0]
    with rasterio.open(swi_like) as swi_ds:
        profile_like = swi_ds.profile


    # ========================================== LOOP OVER MONTHS =================================
    
    for month in m_str:
        logging.info(f'MONTH : {month}')

        swi_files_m = glob.glob(os.path.join(DATA_HISTO, 'MONTH', f'*SWI_*{month}_COMPM.tif'))

        # --- MISSING INDICES ---
        if swi_files_m==[]:
            logging.warning('NO PRODUCTS -> PASS')
            continue

        # --- SINGLE YEAR ---
        if len(swi_files_m)==1:
            logging.warning('SINGLE YEAR -> PASS')
            continue

        # --- COMPUTING MONTH MAI ---
        extractMAI(swi_files_m, profile_like, NDAYS_MAX_M, outdirmoisture_mai, 'M', go_stats, outdir_maskareas, outdirmoisture_droughtstats, ANNEX_DIR, CONFIG['TERRITORY'], CONFIG['KEY_STATS'])

        del swi_files_m
    
    
    # --- Copy new MAI file(s) to data histo directory ---
    DATA_DROUGHT = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'ALERT', 'MAI')
    new_mai_m = glob.glob(os.path.join(outdirmoisture_mai,f'MAI_*M.tif'))

    logging.info('Copy new month MAI file(s) to data histo directory')
    for fm in tqdm(new_mai_m):
        fm_datahisto = os.path.join(DATA_DROUGHT, 'MONTH', os.path.basename(fm))
        try:
            shutil.copyfile(fm, fm_datahisto)
        except PermissionError:
            if os.path.exists(fm_datahisto):
                logging.info(f'File already exists: {os.path.basename(fm)} is replaced by new version')
                os.remove(fm_datahisto)
                shutil.copyfile(fm, fm_datahisto)
            else:
                logging.warning(f'Copy PermissionError : {os.path.basename(fm)} impossible to paste')
        del fm_datahisto
    
    # --- Copy/Update MAI Geo Statistics to data histo directory ---
    if go_stats==1:
        logging.info('Copy/Update MAI Geo Statistics to data histo directory')
        DATA_STATS = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'ALERT', 'MAI', 'STATS')
        geostats_df = glob.glob(os.path.join(outdirmoisture_droughtstats,'MAI_STATS_M.csv'))[0]
        f_datahisto = os.path.join(DATA_STATS, os.path.basename(geostats_df))

        if os.path.exists(f_datahisto):
            # Update : replace new stats to histo stats
            df_new = pd.read_csv(geostats_df, sep=';')
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
            shutil.copyfile(geostats_df, f_datahisto)



def prepare_Month_ALERT(new_month, spiHISTO_df, speiHISTO_df, vhiHISTO_df, area_lut, dir_mai, dir_vhi):
    '''
    Prepare different variables (annex, outputs) to treat during Alert processing.
    Extract drought products (spi/spei, mai, vhi) for current (and before) months.
    '''
    
    # Extract month and month before
    m_date = pd.to_datetime(new_month+'01', format='%Y%m%d')
    m_date_bf = m_date - pd.DateOffset(months=1)
    month_bf = m_date_bf.strftime("%Y%m")
    
    # Filter historical data frames to month
    spiMDATES_df = spiHISTO_df[(spiHISTO_df['DATE']==m_date) | (spiHISTO_df['DATE']==m_date_bf)].reset_index(drop=True)
    speiMDATES_df = speiHISTO_df[(speiHISTO_df['DATE']==m_date) | (speiHISTO_df['DATE']==m_date_bf)].reset_index(drop=True)
    vhiMDATES_df = vhiHISTO_df.loc[(vhiHISTO_df['DATE']==m_date)].reset_index(drop=True)

    # Read satellite data
    mai_file = glob.glob(os.path.join(dir_mai, 'MONTH', f'MAI_{new_month}*.tif'))[0]
    vhi_file = glob.glob(os.path.join(dir_vhi, 'MONTH', f'VHI_{new_month}*.tif'))[0]
    with rasterio.open(mai_file) as mai_ds,\
        rasterio.open(vhi_file) as vhi_ds:
            mai_arr = mai_ds.read(1)
            vhi_arr = vhi_ds.read(1)

    mai_file_bf = glob.glob(os.path.join(dir_mai, 'MONTH', f'MAI_{month_bf}*.tif'))
    if mai_file_bf!=[]:
        with rasterio.open(mai_file_bf[0]) as maibf_ds:
            maibf_arr = maibf_ds.read(1)
    else:
        maibf_arr = np.full(mai_arr.shape,np.nan)
    
    # Prepare output data frames and dictionaries    
    drought_cat = pd.DataFrame({'DCAT':[-1,0,1],
                 'PRECIPITATION':['No precipitation data','No precipitation deficit','Precipitation deficit'],
                 'EVAPOTRANSPIRATION':['No precip/temperature data','No evapotranspiration deficit','Evapotranspiration deficit'],
                 'SOIL_MOISTURE':['No soil moisture data','No soil moisture deficit','Soil moisture deficit'],
                 'VEGETATION':['No vegetation data','No vegetation stress','Vegetation stress']})
    
    DroughtAlert_df = area_lut[['OBJECTID','nom']].copy()
    DroughtAlert_df = DroughtAlert_df.sort_values(by=['nom']).reset_index(drop=True)
    DroughtAlert_df = DroughtAlert_df.rename(columns={'nom':'LOCATION'})
    
    DroughtAlert_df['DATE'] = np.full((len(DroughtAlert_df)),m_date.strftime("%Y-%m-%d"))
    DroughtAlert_df['ALERT'] = np.full((len(DroughtAlert_df)),'Not processed')
    DroughtAlert_df['VEGETATION'] = np.full((len(DroughtAlert_df)),'Not processed')
    DroughtAlert_df['SOIL_MOISTURE'] = np.full((len(DroughtAlert_df)),'Not processed')
    DroughtAlert_df['EVAPOTRANSPIRATION'] = np.full((len(DroughtAlert_df)),'Not processed')
    DroughtAlert_df['PRECIPITATION'] = np.full((len(DroughtAlert_df)),'Not processed')
    DroughtAlert_df['VHI_MEAN'] = np.full((len(DroughtAlert_df)),np.nan)
    DroughtAlert_df['CONF_INDEX'] = np.full((len(DroughtAlert_df)),np.nan)
   
    return (spiMDATES_df, speiMDATES_df, mai_arr, maibf_arr,
            vhi_arr, vhiMDATES_df, drought_cat, DroughtAlert_df)



def detect_Deficits(spi_df, spei_df, mai_arr, maibf_arr, vhi_arr):
    '''
    Function that applies thresholds on different drought indicators to detect deficit.
    Threshold values used to detect deficit are -1 for spi/spei/mai and 0.3 for vhi
    Three categories/values are obtained :
    - No Data = -1
    - No Drought = 0
    - Drought = 1 (DEFICIT)
    The drought categories here are given per station (spi/spei) or per pixel (mai/vhi).
    '''

    THRESH_METEO_MAI = -1
    THRESH_VHI = 0.3
    
    # ---- Precipitation deficit (from SPI) ---
    # uses also the month before (cf. corr SPI vs VHI)
    spi_df['DROUGHT'] = np.full_like(spi_df['SPI3_MENS'], -1, dtype='int16')    # No Data
    spi_df.loc[spi_df['SPI3_MENS'] > THRESH_METEO_MAI, 'DROUGHT'] = 0           # No Drought
    spi_df.loc[spi_df['SPI3_MENS'] <= THRESH_METEO_MAI, 'DROUGHT'] = 1          # Drought
    
    # ---- Precipitation-Evapotranspiration deficit (from SPEI) ---
    # uses also the month before (cf. corr SPEI vs VHI)
    spei_df['DROUGHT'] = np.full_like(spei_df['SPEI_3'], -1, dtype='int16')     # No Data
    spei_df.loc[spei_df['SPEI_3'] > THRESH_METEO_MAI, 'DROUGHT'] = 0            # No Drought
    spei_df.loc[spei_df['SPEI_3'] <= THRESH_METEO_MAI, 'DROUGHT'] = 1           # Drought
    
    #---- Soil Moisture deficit (from MAI) ---
    # uses also the month before (cf. corr MAI vs VHI)
    MAI_t = np.full_like(mai_arr, -1, dtype='int16')        # No Data
    MAI_t[mai_arr > THRESH_METEO_MAI] = 0                   # No Drought
    MAI_t[(mai_arr <= THRESH_METEO_MAI)] = 1                # Drought
    MAIbf_t = np.full_like(maibf_arr, -1, dtype='int16')
    MAIbf_t[maibf_arr > THRESH_METEO_MAI] = 0
    MAIbf_t[(maibf_arr <= THRESH_METEO_MAI)] = 1
    
    # ---- Vegetation stress (from VHI) ---
    VHI_t = np.full_like(vhi_arr, -1, dtype='int16')        # No Data
    VHI_t[vhi_arr > THRESH_VHI] = 0                         # No Drought
    VHI_t[(vhi_arr <= THRESH_VHI)] = 1                      # Drought

    return (spi_df, spei_df, MAI_t, MAIbf_t, VHI_t)



def aggregStation(meteo_df, lut_df, a):
    '''
    Aggregating spi values from stations of the same sub-area.
    Stations to aggregate are selected from lut_df.

    Note NC :   If no station on the sub-area (communes in NC),
                the stations from the closest communes were selected and listed in lut_df
                (and according to climatic regions)
    '''
    
    s_list = list(lut_df[lut_df['LOCATION']==a].dropna(axis='columns').iloc[0,1:])
    
    if len(s_list) == 0:
        return pd.DataFrame([])

    COL_NAME = list(meteo_df)[2]

    meteo_list_df = [meteo_df[meteo_df['NOM']==s].drop(columns=['NOM']).reset_index(drop=True) for s in s_list]
    merged_df  = meteo_list_df[0]
    if len(meteo_list_df)>1:
        for i in range(len(meteo_list_df)-1):
            merged_df = pd.merge(merged_df, meteo_list_df[i+1], on='DATE', how='outer', suffixes=(f'_{str(i)}', f'_{str(i+1)}'))

    col2merge = []
    for i in range(len(list(merged_df))):
        if COL_NAME in list(merged_df)[i]: col2merge.append(list(merged_df)[i])

    merged_df[COL_NAME] = merged_df[col2merge].mean(axis=1)
    out_merged_df = merged_df[['DATE',COL_NAME]].copy()

    return out_merged_df



def extractStation(meteo_df, lut_df, a):
    '''
    Extracting spi values from stations of the same sub-area.
    Stations to extract are selected from lut_df.
    
    Note NC :   If no station on the sub-area (communes in NC),
                the stations from the closest communes were selected and listed in lut_df
                (and according to climatic regions)
    '''
    
    s_list = list(lut_df[lut_df['LOCATION']==a].dropna(axis='columns').iloc[0,1:])
    
    if len(s_list) == 0:
        return pd.DataFrame([])
    
    meteo_list = [meteo_df.loc[meteo_df['NOM']==s].reset_index(drop=True) for s in s_list]
    concat_df = pd.concat(meteo_list, axis=0).reset_index(drop=True)
        
    return concat_df



def crosscorr_pearson(datax, datay, lag=0):
    """
    Lag-N Pearson cross correlation. 
    Shifted data wraped with end values
    -> RETURN P-VALUE
    
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    
    Returns
    ----------
    crosscorr r : float
    p-value : float
    """
    
    if lag > 0:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        r, p = stats.pearsonr(datax,shiftedy)
        return r, p
    
    elif lag < 0:        
        shiftedy = datay.shift(lag)
        shiftedy.iloc[lag:] = datay.iloc[:-lag].values
        r, p = stats.pearsonr(datax,shiftedy)
        return r, p
    
    elif lag == 0:        
        r, p = stats.pearsonr(datax,datay)
        return r, p



def updateRScore(meteo_df, sat_df):
    '''
    Estimating Pearson correlation (Rscore) between meteo and sat time series.
    
    Note: Rscore is the maximum score obtained from lag cross correlation varying
          between -4 and +4 months
    '''
    
    lagmax = 4
    date_start = min(sat_df['DATE'])
    date_end = max(sat_df['DATE'])
    
    # --- Extracting data on the same period + Interpolation (if needed) ---
    
    # Meteo time series
    ind_end = np.where(meteo_df['DATE'] == date_end)
    
    if ind_end[0].size==0:
        logging.critical('\nMeteo data is not available for the last month')
        raise Exception('Meteo data is not available for the last month')
    
    meteo_df_val = meteo_df[(meteo_df['DATE']>=date_start) & (meteo_df['DATE']<=date_end)].reset_index(drop=True)
    
    if meteo_df_val.iloc[:,1].isnull().values.any():
        count_nan = meteo_df_val.iloc[:,1].isnull().sum()
        qscore_meteo = 1 - (count_nan/len(meteo_df_val))
        meteo_df.iloc[:,1] = meteo_df.iloc[:,1].interpolate()
        if meteo_df.iloc[:,1].isnull().values.any():meteo_df.iloc[:,1]=meteo_df.iloc[:,1].fillna(method='backfill')
        meteo_df_val = meteo_df[(meteo_df['DATE']>=date_start) & (meteo_df['DATE']<=date_end)].reset_index(drop=True)
        # logging.info('\nMETEO INTERPOLATION : {} NaN'.format(str(count_nan)))
    else:
        qscore_meteo = 1
    
    # Satellite time series    
    sat_df_mean = sat_df.loc[(sat_df['DATE']>=date_start) & (sat_df['DATE']<=date_end),'MEAN'].reset_index(drop=True)
    sat_df_qscore = sat_df.loc[(sat_df['DATE']>=date_start) & (sat_df['DATE']<=date_end),'QSCORE'].reset_index(drop=True)
    
    qscore_sat = np.mean(sat_df_qscore)
    
    if sat_df_mean.isnull().values.any():
        count_nan = sat_df_mean.isnull().sum()
        if count_nan==len(sat_df_mean):
            logging.info('\nNo sat data !')
            rmax=np.nan
            lagpeak=np.nan
            ppeak=np.nan
            qscore_sat=np.nan
            return rmax, lagpeak, ppeak, qscore_meteo, qscore_sat
        sat_df['MEAN'] = sat_df['MEAN'].interpolate()
        if sat_df['MEAN'].isnull().values.any():sat_df['MEAN']=sat_df['MEAN'].fillna(method='backfill')
        sat_df_mean = sat_df[(sat_df['DATE']>=date_start) & (sat_df['DATE']<=date_end)].reset_index(drop=True)
        # logging.info('\nSAT INTERPOLATION : {} NaN'.format(str(count_nan)))
    
    # --- Computing Pearson time lag correlation ---
    out_crosscorr = [crosscorr_pearson(meteo_df_val.iloc[:,1] , sat_df_mean, lag) for lag in range(-lagmax,lagmax+1)]
    rs = np.array(out_crosscorr)[:,0]
    ps = np.array(out_crosscorr)[:,1]
    lagpeak = int(np.argmax(rs) - np.floor(len(rs)/2))
    rmax = np.max(rs)
    ppeak = ps[np.argmax(rs)]
    # logging.info(f"rmax={rmax:.2f} at {lagpeak} month(s) for p<{ppeak}\n")
    
    return rmax, lagpeak, ppeak, qscore_meteo, qscore_sat



def extractAreasDroughtCat(spi_c, spei_c, MAI_t, MAIbf_t, VHI_t, maskMAI_c, maskVHI_c, drought_cat):
    '''
    Extracting Drought Category on sub-areas for each product :
        DCat = -1 -> No data
        DCat = 0 -> No drought
        Dcat = 1 -> Drought
    
    Note 1: precip (spi), evapo (spei) and soil moisture (mai) status (DCat_) are first
            considered on current month, and then on month before if no drought was detected
            
    Note 2: the actual status is saved (strCatnow_) considering only the current month
    '''
    
    # --- SPI ---
    if spi_c.empty:
        DCatnow_spi = -1
        DCat_spi = DCatnow_spi
        
    else:
        # Look stations status on current month in priority
        spinow_c = spi_c.loc[spi_c['DATE'] == spi_c['DATE'].max()]
        allCat_now = spinow_c['DROUGHT']
        if allCat_now.empty:
            DCatnow_spi = -1    
        elif allCat_now.isin([1]).any():
            # "Drought" if at least one station detects drought
            DCatnow_spi = 1
        elif allCat_now.isin([0]).any() and allCat_now.isin([-1]).any():
            # "No data" if at least one station has no data
            # and another doesn't detect drought
            DCatnow_spi = -1
        elif allCat_now.isin([-1]).all():
            # "No data" if none of the stations has data
            DCatnow_spi = -1
        elif allCat_now.isin([0]).all():
            # "No drought" if none of the stations detects drought
            DCatnow_spi = 0
    
        # If no "Drought" was detected on current month, check the month before
        if allCat_now.empty or allCat_now.isin([1]).any()==0:
            del allCat_now
            spibf_c = spi_c.loc[spi_c['DATE'] == spi_c['DATE'].min()]
            allCat_bf = spibf_c['DROUGHT']
            if allCat_bf.empty:
                  DCat_spi = -1
            elif allCat_bf.isin([1]).any():
                # "Drought" if at least one station detects drought
                DCat_spi = 1
            elif allCat_bf.isin([0]).any() and allCat_bf.isin([-1]).any():
                # "No data" if at least one station has no data
                # and another doesn't detect drought
                DCat_spi = -1
            elif allCat_bf.isin([-1]).all():
                # "No data" if none of the stations has data
                DCat_spi = -1
            elif allCat_bf.isin([0]).all():
                # "No drought" if none of the stations detects drought
                DCat_spi = 0
            del allCat_bf
        
        else:
            DCat_spi = DCatnow_spi
    
    strCatnow_spi = drought_cat.loc[drought_cat['DCAT']==DCatnow_spi,'PRECIPITATION']
    strCatnow_spi = strCatnow_spi.values[0]
    # logging.info(f'{strCatnow_spi} (cat={DCatnow_spi})')
    
    # --- SPEI ---
    if spei_c.empty:
        DCatnow_spei = -1
        DCat_spei = DCatnow_spei
    
    else:
        # Look stations status on current month in priority
        speinow_c = spei_c.loc[spei_c['DATE'] == spei_c['DATE'].max()]
        allCat_now = speinow_c['DROUGHT']
        if allCat_now.empty:
            DCatnow_spei = -1    
        elif allCat_now.isin([1]).any():
            # "Drought" if at least one station detects drought
            DCatnow_spei = 1
        elif allCat_now.isin([0]).any() and allCat_now.isin([-1]).any():
            # "No data" if at least one station has no data
            # and another doesn't detect drought
            DCatnow_spei = -1
        elif allCat_now.isin([-1]).all():
            # "No data" if none of the stations has data
            DCatnow_spei = -1
        elif allCat_now.isin([0]).all():
            # "No drought" if none of the stations detects drought
            DCatnow_spei = 0
        
        # If no "Drought" was detected on current month, check the month before
        if allCat_now.empty or allCat_now.isin([1]).any()==0:
            del allCat_now
            speibf_c = spei_c.loc[spei_c['DATE'] == spei_c['DATE'].min()]
            allCat_bf = speibf_c['DROUGHT']
            if allCat_bf.empty:
                  DCat_spei = -1
            elif allCat_bf.isin([1]).any():
                # "Drought" if at least one station detects drought
                DCat_spei = 1
            elif allCat_bf.isin([0]).any() and allCat_bf.isin([-1]).any():
                # "No data" if at least one station has no data
                # and another doesn't detect drought
                DCat_spei = -1
            elif allCat_bf.isin([-1]).all():
                # "No data" if none of the stations has data
                DCat_spei = -1
            elif allCat_bf.isin([0]).all():
                # "No drought" if none of the stations detects drought
                DCat_spei = 0
            del allCat_bf
        
        else:
            DCat_spei = DCatnow_spei
    
    strCatnow_spei = drought_cat.loc[drought_cat['DCAT']==DCatnow_spei,'EVAPOTRANSPIRATION']
    strCatnow_spei = strCatnow_spei.values[0]
    # logging.info(f'{strCatnow_spei} (cat={DCatnow_spei})')
    
    # --- MAI ---
    Nb_ALLDATA_c = np.sum(maskMAI_c)
    
    # Look at pixel values on the current month in priority
    MAI_tc = np.full(MAI_t.shape,np.nan)
    MAI_tc[maskMAI_c==1] = MAI_t[maskMAI_c==1]
    NANScore_c = np.sum(np.isnan(MAI_tc[maskMAI_c==1]))
    if NANScore_c != Nb_ALLDATA_c:
        allCatnow_mai = np.unique(MAI_tc[~np.isnan(MAI_tc)])
        allCatnow_mai = allCatnow_mai[::-1] # reverse to get drought cat first when applying np.max (if same nb pixels in 2 cat: see below)
        nCat = len(np.unique(allCatnow_mai))
        Nb_pxl_class_now, bin_edges = np.histogram(MAI_tc[~np.isnan(MAI_tc)], bins=nCat)
        Nb_pxl_class_now = Nb_pxl_class_now[::-1]
        DCatnow_mai = int(allCatnow_mai[Nb_pxl_class_now==np.max(Nb_pxl_class_now)][0])
        del Nb_pxl_class_now,nCat,bin_edges
    else :
        DCatnow_mai = -1
        allCatnow_mai = []
    
    # If no "Drought" was detected on current month, check the month before
    if DCatnow_mai!=1:
        MAIbf_tc = np.full(MAIbf_t.shape,np.nan)
        MAIbf_tc[maskMAI_c==1] = MAIbf_t[maskMAI_c==1]
        NANScorebf_c = np.sum(np.isnan(MAIbf_tc[maskMAI_c==1]))
        if NANScorebf_c != Nb_ALLDATA_c:
            allCatbf_mai = np.unique(MAIbf_tc[~np.isnan(MAIbf_tc)])
            allCatbf_mai = allCatbf_mai[::-1]
            nCat = len(np.unique(allCatbf_mai))
            Nb_pxl_class_bf, bin_edges = np.histogram(MAIbf_tc[~np.isnan(MAIbf_tc)], bins=nCat)
            Nb_pxl_class_bf = Nb_pxl_class_bf[::-1]
            DCatbf_mai = int(allCatbf_mai[Nb_pxl_class_bf==np.max(Nb_pxl_class_bf)][0])
            del Nb_pxl_class_bf,nCat,bin_edges
        else :
            DCatbf_mai = -1
            allCatbf_mai = []
        
        # If still no "Drought" detected : check the proportion P of "Drought" pxls on the two-months period
        if DCatbf_mai!=1 and (1 in allCatnow_mai or 1 in allCatbf_mai):
            MAI_MAT =  np.stack((MAI_tc, MAIbf_tc), axis=2)
            allCat_MAT = np.unique(MAI_MAT[~np.isnan(MAI_MAT)])
            allCat_MAT = allCat_MAT[::-1]
            nCat = len(np.unique(allCat_MAT))
            Nb_pxl_class_MAT, bin_edges = np.histogram(MAI_MAT[~np.isnan(MAI_MAT)], bins=nCat)
            Nb_pxl_class_MAT = Nb_pxl_class_MAT[::-1]
            DCat_mai = int(allCat_MAT[Nb_pxl_class_MAT==np.max(Nb_pxl_class_MAT)][0])
            
            Nb_pxl_drought = Nb_pxl_class_MAT[allCat_MAT==1][0]
            Nb_pxl_maj = Nb_pxl_class_MAT[allCat_MAT==DCat_mai][0]
            P = Nb_pxl_drought/Nb_pxl_maj
            
            # If P greater than 50%, set category to "Drought"
            if P>0.5:
                DCat_mai = 1
            del Nb_pxl_drought,Nb_pxl_maj,P,Nb_pxl_class_MAT,allCat_MAT,nCat,bin_edges
        else :
            DCat_mai = DCatbf_mai
        
        del DCatbf_mai,MAIbf_tc,allCatnow_mai,allCatbf_mai
        
    else:
        DCat_mai = DCatnow_mai
    
    strCatnow_mai = drought_cat.loc[drought_cat['DCAT']==DCatnow_mai,'SOIL_MOISTURE']
    strCatnow_mai = strCatnow_mai.values[0]
    # logging.info(f'\n{strCatnow_mai} (cat={DCatnow_mai})')
    
    # --- VHI ---
    Nb_ALLDATA_c = np.sum(maskVHI_c)
    VHI_tc = np.full(VHI_t.shape,np.nan)
    VHI_tc[maskVHI_c==1] = VHI_t[maskVHI_c==1]
    NANScore_c = np.sum(np.isnan(VHI_tc[maskVHI_c==1]))
    
    if NANScore_c != Nb_ALLDATA_c:
        allCat = np.unique(VHI_tc[~np.isnan(VHI_tc)])
        nCat = len(np.unique(allCat))
        Nb_pxl_class, bin_edges = np.histogram(VHI_tc[~np.isnan(VHI_tc)], bins=nCat)
        DCat_vhi = int(allCat[Nb_pxl_class==np.max(Nb_pxl_class)][0])
        
        if DCat_vhi!=1 and 1 in allCat:
            # If "No Drought" or "No Data": check the proportion P of "Drought" pxls
            Nb_pxl_drought = Nb_pxl_class[allCat==1][0]
            Nb_pxl_maj = Nb_pxl_class[allCat==DCat_vhi][0]
            P = Nb_pxl_drought/Nb_pxl_maj
            if P>0.5:
                # If P greater than 50%, set category to "Drought"
                DCat_vhi = 1
            del Nb_pxl_drought,Nb_pxl_maj,P
        del allCat,nCat,Nb_pxl_class,bin_edges
    else:
        DCat_vhi = -1
    
    strCatnow_vhi = drought_cat.loc[drought_cat['DCAT']==DCat_vhi,'VEGETATION']
    strCatnow_vhi = strCatnow_vhi.values[0]
    # logging.info(f'\n{strCatnow_vhi} (cat={DCat_vhi})')
    
    return (DCat_spi, DCat_spei, DCat_mai, DCat_vhi,
            strCatnow_spi, strCatnow_spei, strCatnow_mai, strCatnow_vhi)



def applyAlertClassif(DCat_spi, DCat_spei, DCat_mai, DCat_vhi):
    '''
    Applying Alert Classification according to Sepulcre-Canto et al. (2012)
    '''
    
    if DCat_spi==-1 and DCat_spei!=1 and DCat_mai!=1 and DCat_vhi!=1:
        AlertCat = 'No Data'
    
    if DCat_spi==0 and DCat_spei!=1 and DCat_mai!=1 and DCat_vhi!=1:
        AlertCat = 'No Alert'
    
    if DCat_spi==1 and DCat_spei!=1 and DCat_mai!=1 and DCat_vhi!=1:
        AlertCat = 'Watch'
    
    if (DCat_mai==1 or DCat_spei==1) and DCat_vhi!=1:
        AlertCat = 'Warning'
    
    if DCat_vhi==1:
        AlertCat = 'Alert'
    
    return AlertCat



def process_DroughtAlert(CONFIG, OUTDIR_PATHS):
    '''
    General function for processing drought alerts from several drought indicators :
        - SPI/SPEI meteorological indicators
        - MAI soil moisture deficit
        - VHI vegetation stress
    A classification of drought alert levels is given according to Sepulcre-Canto al. (2012)
    '''
    
    logging.info('\n\n--- PROCESSING DROUGHT ALERT (on Sub-Areas) ---\n')

    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    DATA_HISTO_ALERT = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'ALERT')
    DATA_HISTO_METEO = os.path.join(DATA_HISTO_ALERT, 'METEO')
    DATA_HISTO_MAI = os.path.join(DATA_HISTO_ALERT, 'MAI')
    DATA_HISTO_VHI = os.path.join(CONFIG['DATA_HISTO'], TERRITORY_str, '1_INDICATEURS', 'GLOBAL')
    DATA_ANNEX = os.path.join(CONFIG['ANNEX_DIR'], TERRITORY_str)
    KEY_STATS = CONFIG['KEY_STATS']

    (outdir_moisture, outdirmoisture_preproc, outdirmoisture_prepocswi, 
     outdirmoisture_preprocstats, outdirmoisture_comp, outdirmoisture_compmonth,
     outdirmoisture_drought, outdirmoisture_mai, outdirmoisture_droughtstats,
     outdir_alert, outdir_maskareas) = OUTDIR_PATHS
    
    date_start_str = CONFIG['PERIOD_START']
    date_end_str = CONFIG['PERIOD_END']

    period_start = pd.to_datetime(date_start_str, format='%Y-%m-%d')
    period_end_inclusive = pd.to_datetime(date_end_str, format='%Y-%m-%d') + pd.DateOffset(days=-1)
    if period_start.month==period_end_inclusive.month:
        date_dt = pd.Timestamp(period_start)
        date_vect = [date_dt.strftime("%Y%m")]
    else:
        date_dt = pd.date_range(start=date_start_str, end=period_end_inclusive.strftime('%Y-%m-%d'), freq='M')
        date_vect = date_dt.strftime("%Y%m")
    
    r_go = 1 # Parameter to compute rscore juste once

    # --- Prepare/Read input masks and look-up table (for sub-areas delimitation) ---
    mai_like = glob.glob(os.path.join(DATA_HISTO_MAI, 'MONTH', f'MAI*.tif'))[0]
    maskmai_name = 'mask_Areas_MAI.tif'
    maskmai_ok = len(glob.glob(os.path.join(outdir_maskareas, maskmai_name)))>0
    if maskmai_ok==0:
        file_areas = glob.glob(os.path.join(DATA_ANNEX, 'Areas', '*.shp'))
        if len(file_areas)==0:
            logging.critical('Missing input shapefile for sub-areas delimitation')
            raise Exception ('Missing input shapefile for sub-areas delimitation')
        file_areas = file_areas[0]
        geostats.prepareGeoStatsMasks(mai_like, file_areas, outdir_maskareas, suffix='MAI', areas_key=KEY_STATS)
    with rasterio.open(glob.glob(os.path.join(outdir_maskareas, maskmai_name))[0]) as area_ds:
        maskAREA_mai = area_ds.read(1)
    
    vhi_like = glob.glob(os.path.join(DATA_HISTO_VHI, 'MONTH', f'VHI*.tif'))[0]
    maskvhi_name = 'mask_Areas_VHI.tif'
    maskvhi_ok = len(glob.glob(os.path.join(outdir_maskareas, maskvhi_name)))>0
    if maskvhi_ok==0:
        file_areas = glob.glob(os.path.join(DATA_ANNEX, 'Areas', '*.shp'))
        if len(file_areas)==0:
            logging.critical('Missing input shapefile for sub-areas delimitation')
            raise Exception ('Missing input shapefile for sub-areas delimitation')
        file_areas = file_areas[0]
        geostats.prepareGeoStatsMasks(vhi_like, file_areas, outdir_maskareas, suffix='VHI', areas_key=KEY_STATS)
    with rasterio.open(glob.glob(os.path.join(outdir_maskareas, maskvhi_name))[0]) as area_ds:
        maskAREA_vhi = area_ds.read(1)

    area_lut = pd.read_csv(glob.glob(os.path.join(outdir_maskareas,'ID_Name_Areas-lookup_VHI.csv'))[0], sep=';')
    head_alert = 1 # the first time, add header to alert data frame

    # --- Prepare/Read input data frames ---
    spi_csv = os.path.join(DATA_HISTO_METEO, 'SPI_ref_1991_2020.csv')
    spi_df = pd.read_csv(spi_csv,sep=';',decimal=',')
    spiHISTO_df = spi_df[['NOM','DATE','SPI3_MENS']].copy()
    spiHISTO_df['DATE'] = pd.to_datetime(spiHISTO_df.DATE, format='%Y%m')
    
    spei_csv = os.path.join(DATA_HISTO_METEO, 'SPEI_ref_1991_2020.csv')
    spei_df = pd.read_csv(spei_csv,sep=';',decimal=',')
    speiHISTO_df = spei_df[['NOM','DATE','SPEI_3']].copy()
    speiHISTO_df['DATE'] = pd.to_datetime(speiHISTO_df.DATE, format='%Y%m')

    vhistats_csv = os.path.join(DATA_HISTO_VHI, 'STATS','VHI_STATS_M_NoTrees_NoBuild.csv')
    vhistats_df = pd.read_csv(vhistats_csv,sep=';')
    try:
        vhistats_df['DATE'] = pd.to_datetime(vhistats_df.DATE, format='%Y-%m-%d')
    except ValueError:
        vhistats_df['DATE'] = pd.to_datetime(vhistats_df.DATE, format='%d/%m/%Y')
    vhiHISTO_df = vhistats_df[['LOCATION','DATE','MEAN','QSCORE']].copy()

    stations_spi_csv = os.path.join(DATA_ANNEX, 'Stations', 'SPI_communes_stations.csv')
    stations_spi_df = pd.read_csv(stations_spi_csv,sep=';')
    stations_spei_csv = os.path.join(DATA_ANNEX, 'Stations', 'SPEI_communes_stations.csv')
    stations_spei_df = pd.read_csv(stations_spei_csv,sep=';')

    rscore_df = pd.DataFrame(columns=['LOCATION','RMAX','LAGMAX','PVMAX','QSCORE SPI','QSCORE VHI'])


    # ========================================== LOOP OVER MONTHS =================================

    for month in tqdm(date_vect, desc='MONTH'):

        # --- PREPARE FOR CURRENT MONTH ---
        (spiMDATES_df, speiMDATES_df, mai_arr, maibf_arr,
         vhi_arr, vhiMDATES_df, drought_cat, DroughtAlert_df) = prepare_Month_ALERT(month,
                                                                                    spiHISTO_df,
                                                                                    speiHISTO_df,
                                                                                    vhiHISTO_df,
                                                                                    area_lut,
                                                                                    DATA_HISTO_MAI,
                                                                                    DATA_HISTO_VHI)
        
        # --- DEFICITS DETECTION ---
        spi_df, spei_df, MAI_t, MAIbf_t, VHI_t = detect_Deficits(spiMDATES_df,
                                                                 speiMDATES_df,
                                                                 mai_arr,
                                                                 maibf_arr,
                                                                 vhi_arr)


    # ========================================== LOOP OVER SUB-AREAS =================================
        
        for a in tqdm(DroughtAlert_df['LOCATION'], desc='SUB-AREA'):
            
            # --- RSCORE FOR CURRENT SUB-AREA (juste ONCE) ---
            if r_go==1:
                # Extract time series on sub-area
                vhiHISTO_a = vhiHISTO_df.loc[vhiHISTO_df['LOCATION']==a].reset_index(drop=True)
                spiHISTO_a = aggregStation(spiHISTO_df, stations_spi_df, a)
                
                # Compute Pearson correlation (RSCORE)
                if spiHISTO_a.empty:
                    r_dict = {'LOCATION':a,'RMAX':np.nan,'LAGMAX':np.nan,'PVMAX':np.nan,
                              'QSCORE SPI':np.nan,'QSCORE VHI':np.nan}
                else:
                    vhiHISTO_a = vhiHISTO_a.sort_values(by=['DATE'])
                    spiHISTO_a = spiHISTO_a.sort_values(by=['DATE'])
                    rmax, lagpeak, ppeak, qscore_spi, qscore_vhi = updateRScore(spiHISTO_a, vhiHISTO_a)
                    r_dict = {'LOCATION':a,'RMAX':round(rmax,2),'LAGMAX':lagpeak,'PVMAX':ppeak,
                              'QSCORE SPI':round(qscore_spi,2),'QSCORE VHI':round(qscore_vhi,2)}
                    del rmax, lagpeak, ppeak, qscore_spi, qscore_vhi
                rscore_df = pd.concat([rscore_df, pd.DataFrame([r_dict])], ignore_index=True)
                del r_dict, vhiHISTO_a, spiHISTO_a

                            
            # --- FILL ALERT data frame with CONF_INDEX and VHI_MEAN
            vhi_mean_a = vhiMDATES_df.loc[vhiMDATES_df['LOCATION']==a,'MEAN']
            vhi_mean_a = vhi_mean_a.values[0]
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'VHI_MEAN'] = vhi_mean_a
            vhi_qscore_a = vhiMDATES_df.loc[vhiMDATES_df['LOCATION']==a,'QSCORE']
            vhi_rscore_a = rscore_df.loc[rscore_df['LOCATION']==a,'RMAX']
            conf_indx_a = np.mean([vhi_qscore_a,vhi_rscore_a])
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'CONF_INDEX'] = np.round(conf_indx_a,2)
            del vhi_mean_a,vhi_qscore_a,vhi_rscore_a,conf_indx_a
            
            # --- SUB-AREA SYNTHESIS : Counting the majority category ---
            spi_a = extractStation(spiMDATES_df, stations_spi_df, a)
            spei_a = extractStation(speiMDATES_df, stations_spei_df, a)
            maskMAI_a = (maskAREA_mai==int(DroughtAlert_df['OBJECTID'][DroughtAlert_df['LOCATION']==a]))
            maskVHI_a = (maskAREA_vhi==int(DroughtAlert_df['OBJECTID'][DroughtAlert_df['LOCATION']==a]))
            
            (DCat_spi, DCat_spei, DCat_mai, DCat_vhi, strCatnow_spi, strCatnow_spei,
             strCatnow_mai, strCatnow_vhi) = extractAreasDroughtCat(spi_a, spei_a,
                                                                    MAI_t, MAIbf_t, VHI_t,
                                                                    maskMAI_a, maskVHI_a,
                                                                    drought_cat)
            
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'PRECIPITATION'] = strCatnow_spi     
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'EVAPOTRANSPIRATION'] = strCatnow_spei
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'SOIL_MOISTURE'] = strCatnow_mai
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'VEGETATION'] = strCatnow_vhi
            
            del strCatnow_spi, strCatnow_spei, strCatnow_mai, strCatnow_vhi
            del spi_a, spei_a, maskMAI_a, maskVHI_a
            
            # --- ALERT CLASSIFICATION ---
            AlertCat = applyAlertClassif(DCat_spi, DCat_spei, DCat_mai, DCat_vhi)
            DroughtAlert_df.loc[DroughtAlert_df['LOCATION']==a,'ALERT'] = AlertCat
            del AlertCat, DCat_spi, DCat_spei, DCat_mai, DCat_vhi

        # --- Saving Alerts data frame ---
        DroughtAlert_df = DroughtAlert_df.drop(columns=['OBJECTID'])
        DroughtAlert_df = DroughtAlert_df.sort_values(by=['LOCATION','DATE'])
        DroughtAlert_df.to_csv(os.path.join(outdir_alert, 'ALERT_DROUGHT.csv'),
            index = False,
            float_format='%.2f',
            decimal = '.',
            sep=';',
            mode='a',
            header = head_alert)
        head_alert = 0 # remove header the next times

        # --- Savin Rscore data frame (just ONCE) ---
        if r_go==1:
            rscore_df.to_csv(os.path.join(outdir_alert,'VHI_SPI_RSCORE_QSCORES_NoTrees_NoBuild.csv'),
                index = False,
                decimal = '.',
                sep=';')
        r_go = 0 # no rscores for the next times

        del (spiMDATES_df, speiMDATES_df, mai_arr, maibf_arr, DroughtAlert_df,
             vhi_arr, vhiMDATES_df, drought_cat, spi_df, spei_df, MAI_t, MAIbf_t, VHI_t)
    
    # --- Copy/Update ALERT_DROUGHT Data frame to data histo directory ---
    logging.info(f'Copy/Update Dataframe ALERT_DROUGHT to data histo directory')
    f_datahisto = os.path.join(DATA_HISTO_ALERT, 'ALERT_DROUGHT.csv')
    f = os.path.join(outdir_alert, 'ALERT_DROUGHT.csv')
    if os.path.exists(f_datahisto):
        # Concatenate : add new alerts at the end of histo alerts
        df_histo = pd.read_csv(f_datahisto, sep=';')
        df_new = pd.read_csv(f, sep=';')
        try:
            df_histo['DATE'] = pd.to_datetime(df_histo.DATE, format='%Y-%m-%d')
        except ValueError:
            df_histo['DATE'] = pd.to_datetime(df_histo.DATE, format='%d/%m/%Y')
        try:
            df_new['DATE'] = pd.to_datetime(df_new.DATE, format='%Y-%m-%d')
        except ValueError:
            df_new['DATE'] = pd.to_datetime(df_new.DATE, format='%d/%m/%Y')

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
        # df_histo = pd.concat([df_histo, df_new], ignore_index=True)
        df_histo = df_histo.sort_values(by=['LOCATION','DATE']).reset_index(drop=True)

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
    del f_datahisto, f

    # --- Copy/Update VHI_SPI_RSCORE_QSCORES Data frame to data histo directory ---
    logging.info(f'Copy/Update Dataframe VHI_SPI_RSCORE_QSCORES to data histo directory')
    f_datahisto = os.path.join(DATA_HISTO_ALERT, 'VHI_SPI_RSCORE_QSCORES_NoTrees_NoBuild.csv')
    f = os.path.join(outdir_alert, 'VHI_SPI_RSCORE_QSCORES_NoTrees_NoBuild.csv')   
    try:
        shutil.copyfile(f, f_datahisto)
    except PermissionError:
        if os.path.exists(f_datahisto):
            logging.info(f'File already exists: {os.path.basename(f_datahisto)} is replaced by new version')
            os.remove(f_datahisto)
            shutil.copyfile(f, f_datahisto)
        else:
            logging.critical(f'Copy PermissionError : {os.path.basename(f)} impossible to paste')
            raise Exception('Copy PermissionError : impossible to paste')
    del f_datahisto


