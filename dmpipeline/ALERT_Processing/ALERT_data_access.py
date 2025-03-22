# -*- coding: utf-8 -*-
"""
###################################################################################

ALERT Functions for accessing/collecting complementary data involved in alert chain

###################################################################################
"""

import os
import glob
import json
from ftplib import FTP
from tqdm import tqdm
import dmpipeline.ALERT_Processing.ftp_accounts as alertauth

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def read_credentials(credentials_file):
    """
    # Function to read identification information from JSON file
    """
    with open(credentials_file, "r") as json_file:
        credentials = json.load(json_file)
    return credentials



def download_copernicus_ftp(credentials_file, local_coper_dir):
    """
    Function for downloading directories/files from copernicus-vito ftp server to local directory
    """
    # Identification to ftp server
    credentials = read_credentials(credentials_file)
    ftp = FTP(credentials["ftp_server"])
    ftp.login(user=credentials["ftp_user"], passwd=credentials["ftp_password"])
    remote_dir = credentials["ftp_directory"]

    # List directories on server
    dir_list = ftp.nlst()
    for dir in dir_list:
        # logging.info(dir)
        if dir!=remote_dir: continue
        else: break
    
    # Missing remote_dir -> STOP
    if dir!=remote_dir:
        logging.info(f'{remote_dir} was not found on server')
    
    # OR go to remote_dir and list dir -> CONTINUE
    else:
        os.makedirs(local_coper_dir, exist_ok=True)
        ftp.cwd(remote_dir)
        remotedir_path = ftp.pwd()
        new_dir_list = ftp.nlst()

        # Go to new_remote_dir list
        for new_dir in tqdm(new_dir_list, desc='COPERNICUS DOWNLOADING'):
            # logging.info(f'Directory {new_dir}')
            new_local_dir = os.path.join(local_coper_dir, new_dir)
            os.makedirs(new_local_dir, exist_ok=True)
            os.chdir(new_local_dir)
            ftp.cwd(new_dir)

            # Local files downloading
            files_list = ftp.nlst()
            for file in files_list:
                try:
                    with open(file, "wb") as f:
                        ftp.retrbinary("RETR " + file, f.write)
                        logging.info(f"{file} Downloaded with success")
                except Exception as e:
                    if '[Errno 13] Permission denied' in str(e):
                        logging.info(f"{file} File already exists : not downloaded")
                    else:
                        logging.critical(f"An unknown error has occured : {str(e)}")
            
            ftp.cwd(remotedir_path)
            os.chdir(local_coper_dir)
        
        ftp.quit()
        logging.info("Downloading Copernicus completed")



def download_meteo_ftp(credentials_file, local_meteo_dir):
    """
    Function for downloading (csv) files from meteo-france ftp server to local directory
    """

    # Identification to ftp server
    credentials = read_credentials(credentials_file)
    ftp = FTP(credentials["ftp_server"])
    ftp.login(user=credentials["ftp_user"], passwd=credentials["ftp_password"])

    # List files on server
    files_list = ftp.nlst()

    # Local downloading
    os.makedirs(local_meteo_dir, exist_ok=True)
    os.chdir(local_meteo_dir)
    for file in tqdm(files_list, desc='METEO DOWNLOADING'):
        try:
            with open(file, "wb") as f:
                ftp.retrbinary("RETR " + file, f.write)
                logging.info(f"{file} Downloaded with success")
        except Exception as e:
            if '[Errno 13] Permission denied' in str(e):
                os.remove(file)
                with open(file, "wb") as f:
                    ftp.retrbinary("RETR " + file, f.write)
                    logging.info(f'{file} File already exists: replaced by new version')
            else:
                logging.critical(f"An unknown error has occured : {str(e)}")
    
    ftp.quit()
    logging.info("Downloading Meteo completed")



def download_newProducts(CONFIG):
    """
    Downloading new Copernicus and Meteo products from FTP servers to local data historic directories :
        - Copernicus = last ASCAT raw daily products
        - Meteo = last SPI-SPEI month indicators
    """

    logging.info('\n\n --- COLLECTING NEW COPERNICUS/METEO PRODUCTS ---\n')

    path2keys = os.path.dirname(alertauth.__file__)
    TERRITORY_str = CONFIG['TERRITORY'].replace('"', '').replace(' ', '_').replace('(', '').replace(')', '')
    local_dir_coper = os.path.join(CONFIG['DATA_HISTO'],TERRITORY_str,'0_INDICES','ASCAT','0_RAW')  
    local_dir_meteo = os.path.join(CONFIG['DATA_HISTO'],TERRITORY_str,'1_INDICATEURS','ALERT','METEO')

    # --- Donwload COPERNICUS : ASCAT raw daily products ---
    keyfile_coper = glob.glob(os.path.join(path2keys,'*copernicus*.json'))[0]
    download_copernicus_ftp(keyfile_coper, local_dir_coper)

    # --- Donwload METEO-FRANCE : SPI-SPEI month indicators ---
    keyfile_meteo = glob.glob(os.path.join(path2keys,'*meteo*.json'))[0]
    download_meteo_ftp(keyfile_meteo, local_dir_meteo)


