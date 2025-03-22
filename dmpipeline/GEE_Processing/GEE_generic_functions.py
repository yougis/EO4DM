# -*- coding: utf-8 -*-
"""
##############################################################################

GEE Useful generic functions

##############################################################################
"""

import os
import sys
import ee
import time
import json
import glob
import shutil
from tqdm import tqdm
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
from oauth2client.service_account import ServiceAccountCredentials
import dmpipeline.GEE_Processing.GEE_main_functions as geemain

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def googleAuthentification(path2key=os.getcwd()):
    """
    Automatic authentification to GEE
    Reads key_file to access credentials of google service account (needs to be created)
    """
    
    key_file = glob.glob(os.path.join(path2key,'*.json'))[0]
    with open(key_file, 'r') as key_object:
        key_content = key_object.read()

    key_dict = json.loads(key_content)
    service_account = key_dict['client_email']
    project_id = key_dict['project_id']

    credentials = ee.ServiceAccountCredentials(service_account, key_file)
    ee.Initialize(credentials,project=project_id)

    return project_id



def googleErrorsControl(in_gee, path2key):
    """
    Control procedure used to track GEE Errors :
        - Computation timed out.
        - Image.load
        - Collection.loadTable
        - Other (Unknown, Python error)
    """

    MAX_GEETRY = 6
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:

        try:
            out_gee = in_gee.getInfo()
            google_error = 'no_google_error'

        except ee.EEException as e:
            if 'Computation timed out.' in str(e):
                if n_geetry<MAX_GEETRY:
                    logging.warning(f'Loose gee connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                    googleAuthentification(path2key)
                else: 
                    logging.critical(f'Error gee authentification ({e})')
                    raise Exception ('Error gee authentification')
            
            elif 'User memory limit exceeded.' in str(e):
                if n_geetry<MAX_GEETRY:
                    logging.warning(f'Too many gee inupts at the same time ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                    googleAuthentification(path2key)
                else: 
                    logging.critical(f'Error gee user memory ({e})')
                    raise Exception ('Error gee user memory')
            
            elif 'The service is currently unavailable.' in str(e):
                if n_geetry<MAX_GEETRY:
                    logging.warning(f'{e} -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                    googleAuthentification(path2key)
                else: 
                    logging.critical(f'{e}')
                    raise Exception ('Error gee unavailable service')
            
            elif 'Image.load' in str(e):
                logging.warning(f'No access to gee image ({e})\n -> Retry loading (If still not working, control that max size available on gee assets project was not reached)')
            
            elif 'Collection.loadTable' in str(e):
                logging.warning(f'No access to gee table ({e})\n -> Retry loading (If still not working, control collection table)')

            else:
                logging.critical(f'Unknown gee error ({e})')
                raise Exception ('Unknown gee error')
            
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)

        except Exception as e:
            if 'ConnectionError' in str(e):
                if n_geetry<MAX_GEETRY:
                    logging.warning(f'Loose google connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                    googleAuthentification(path2key)

            elif 'Connection aborted.' in str(e):
                if n_geetry<MAX_GEETRY:
                    logging.warning(f'Connection aborted ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                    googleAuthentification(path2key)
            
            else:
                logging.critical(f'Unknown Python error ({e})')
                raise Exception ('Unknown Python error')
            
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)
    
    return out_gee



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



def cleanAssets(asset):
    """
    Clean gee asset :
        delete folder, (sub)collection(s), (sub)image(s)
    """
    
    logging.info(f'\nCleaning asset(s) :')
    
    if asset['type']=='IMAGE':
        ee.data.deleteAsset(asset['id'])

    elif asset['type']=='IMAGE_COLLECTION':
        subimage = (ee.data.listAssets({'parent':asset['id']}))['assets']
        if subimage != []:
            for suba in tqdm(subimage): ee.data.deleteAsset(suba['id'])
        ee.data.deleteAsset(asset['id'])
        
    elif asset['type']=='FOLDER':
        subcoll = (ee.data.listAssets({'parent':asset['id']}))['assets']
        if subcoll != []:
            for suba in tqdm(subcoll):
                if suba['type']=='IMAGE_COLLECTION':
                    subimage = (ee.data.listAssets({'parent':suba['id']}))['assets']
                    if subimage != []:
                        for subsuba in subimage: ee.data.deleteAsset(subsuba['id'])
                ee.data.deleteAsset(suba['id'])

        ee.data.deleteAsset(asset['id'])



def createAssetsFolder(folderName, gee_workdir, clean_asset=None):
    '''
    Prepare gee folder
    Verify if the folder to process exists in the google assets products
            If no -> create preproc folder
            If yes -> refers to "clean_asset" to delete or keep it
    '''

    gee_folder = gee_workdir+'/'+folderName
    assets_list = (ee.data.listAssets({'parent':gee_workdir}))['assets']

    # --- Vérify if working folder doesn't exist, otherwise ask user to delete it ---
    assets_match = []
    for asset in assets_list:
        if gee_folder in asset["name"]:
            # logging.info(asset["name"])
            assets_match += [asset]

    if assets_match==[]:
        ee.data.createAsset({'type': 'FOLDER'}, opt_path=gee_folder) # Create new folder
    
    else:
        logging.info(f'\nFolder {gee_folder} already exists !')
        if clean_asset==1:
            cleanAssets(assets_match[0]) # Clean folder
            ee.data.createAsset({'type': 'FOLDER'}, opt_path=gee_folder) # Create new folder
    
    return gee_folder



def createAssetsCollection(collectionName, gee_folder, clean_asset=None):
    """
    Prepare gee collection
    Verify if the collection to process exist in the google assets products
            If no -> create preproc collection
            If yes -> refers to "clean_asset" to delete or keep it
    """

    gee_collection = gee_folder+'/'+collectionName

    col_list = (ee.data.listAssets({'parent':gee_folder}))['assets']
    col_match = []
    for asset in col_list:
        if gee_collection in asset["name"]:
            # logging.info(asset["name"])
            col_match += [asset]

    if col_match==[]:
        ee.data.createAsset({'type': 'ImageCollection'}, opt_path=gee_collection) # Create new collection
        pass 
    else:
        if clean_asset==1:                
            cleanAssets(col_match[0]) # Clean collection
            ee.data.createAsset({'type': 'ImageCollection'}, opt_path=gee_collection) # Create new collection
        # else:
        #     gee_collection = 'already processed'

    return gee_collection



def prepareFolders(folderName, collectionName, gee_workdir):
    """
    Prepare gee folders
    Verify if the folder/collection to process exist in the google assets products
            If no -> create preproc folder/collection
            If yes -> ask user to delete or leave it
    """

    gee_preproc = gee_workdir+'/'+folderName
    preproc_collection = gee_preproc+'/'+collectionName
    assets_list = (ee.data.listAssets({'parent':gee_workdir}))['assets']
    
    # --- Vérify if working folder doesn't exist, otherwise ask user to delete it ---
    assets_match = []
    for asset in assets_list:
        if gee_preproc in asset["name"]:
            # logging.info(asset["name"])
            assets_match += [asset]

    if assets_match==[]:
        ee.data.createAsset({'type': 'FOLDER'}, opt_path=gee_preproc) # Create new folder
        new_col = ee.data.createAsset({'type': 'ImageCollection'}, opt_path=preproc_collection) # Create new collection
    
    else:
        logging.info(f'\nFolder {gee_preproc} already exists !')
        answer = input(' -> Do you want to delete it ? (y/n) ')
        
        if answer.lower() in "y":
            cleanAssets(assets_match[0]) # Clean folder
            ee.data.createAsset({'type': 'FOLDER'}, opt_path=gee_preproc) # Create new folder
            new_col = ee.data.createAsset({'type': 'ImageCollection'}, opt_path=preproc_collection) # Create new collection

        elif answer.lower() in "n":

            # --- Vérify if new collection doesn't exist, otherwise ask user to delete it or stops execution ---
            col_list = (ee.data.listAssets({'parent':gee_preproc}))['assets']
            col_match = []
            for asset in col_list:
                if preproc_collection in asset["name"]:
                    # logging.info(asset["name"])
                    col_match += [asset]

            if col_match==[]:
                new_col = ee.data.createAsset({'type': 'ImageCollection'}, opt_path=preproc_collection) # Create new collection
                pass
            
            else:
                logging.info(f'\nCollection {preproc_collection} already exists !')
                answer = input(' -> Need to delete it, otherwise execution will stop, y/n ? ')
                
                if answer.lower() in "y":
                    cleanAssets(col_match[0]) # Clean collection
                    ee.data.createAsset({'type': 'ImageCollection'}, opt_path=preproc_collection) # Create new collection
                
                elif answer.lower() in "n":
                    sys.exit()

    return preproc_collection



def listLandsatTiles(landsat_collection, path2key):
    """
    List the disctinct tiles ('path0raw') in a landsat collection
    """

    dataset_distinct_path = landsat_collection.distinct('WRS_PATH')
    Nb_disctinct_path = googleErrorsControl(dataset_distinct_path.size(), path2key)
    list_distinct_path = dataset_distinct_path.toList(Nb_disctinct_path)
    path_list = [googleErrorsControl(ee.Image(list_distinct_path.get(i)).get('WRS_PATH'), path2key) for i in range(Nb_disctinct_path)]

    tiles_list = []
    for path in path_list:
        dataset_path = landsat_collection.filter(ee.Filter.eq('WRS_PATH', path))
        dataset_distinct_row = dataset_path.distinct('WRS_ROW')
        Nb_disctinct_row = googleErrorsControl(dataset_distinct_row.size(), path2key)
        list_distinct_row = dataset_distinct_row.toList(Nb_disctinct_row)
        row_list = [googleErrorsControl(ee.Image(list_distinct_row.get(i)).get('WRS_ROW'), path2key) for i in range(Nb_disctinct_row)]

        for row in row_list:
            tiles_list += ['{}{}{}'.format(path,0,row)]
        
        del dataset_path, dataset_distinct_row, Nb_disctinct_row, list_distinct_row, row_list

    return tiles_list



def listSentinelTiles(sentinel_collection, path2key):
    """
    List the disctinct tiles in a sentinel collection
    """

    dataset_distinct_tile = sentinel_collection.distinct('MGRS_TILE')
    Nb_disctinct_tile = googleErrorsControl(dataset_distinct_tile.size(), path2key)
    list_distinct_tile = dataset_distinct_tile.toList(Nb_disctinct_tile)
    tiles_list = [googleErrorsControl(ee.Image(list_distinct_tile.get(i)).get('MGRS_TILE'), path2key) for i in range(Nb_disctinct_tile)]

    return tiles_list



def extractCollectionsTile(L7collection, L8collection, L9collection, S2collection, tile_L, tiles_S2, landsat_grid, Count_dict,
                           L7fullcollection, L8fullcollection, L9fullcollection, path2key):
    """
    Extract the satellite collections on specified landsat tile_L.
    S2 tiles can be filtered according to the input parameter tiles_S2.
    Returns filtered collections, the names of S2 tiles found and it updates dictionnary with date start/end and number of images.
    """

    logging.info(f'\nTILE LANDSAT {tile_L} :')

     # --- Extract L7 collection for the specific tile ---
    L7_tile = L7collection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    L7fullcollection_tile = L7fullcollection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    Nb_L7images = googleErrorsControl(L7_tile.size(), path2key)
    if Nb_L7images != 0:
        date_start_L7 = googleErrorsControl(L7_tile.limit(1, 'system:time_start', True).first().get('DATE_ACQUIRED'), path2key)
        date_end_L7 = googleErrorsControl(L7_tile.limit(1, 'system:time_start', False).first().get('DATE_ACQUIRED'), path2key)
    else:
        date_start_L7 = ''
        date_end_L7 = ''
    Count_dict += [ee.Feature(None, {'TILE': tile_L,
                                     'PRODUCT': 'L7',
                                     'DATE START': date_start_L7,
                                     'DATE END': date_end_L7,
                                     'NUMBER OF IMAGES': Nb_L7images})]
    logging.info(f' - L7 = {Nb_L7images}')

    # --- Extract L8 collection for the specific landsat tile ---
    L8_tile = L8collection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    L8fullcollection_tile = L8fullcollection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    Nb_L8images = googleErrorsControl(L8_tile.size(), path2key)
    if Nb_L8images != 0:
        date_start_L8 = googleErrorsControl(L8_tile.limit(1, 'system:time_start', True).first().get('DATE_ACQUIRED'), path2key)
        date_end_L8 = googleErrorsControl(L8_tile.limit(1, 'system:time_start', False).first().get('DATE_ACQUIRED'), path2key)
    else:
        date_start_L8 = ''
        date_end_L8 = ''
    Count_dict += [ee.Feature(None, {'TILE': tile_L,
                                     'PRODUCT': 'L8',
                                     'DATE START': date_start_L8,
                                     'DATE END': date_end_L8,
                                     'NUMBER OF IMAGES': Nb_L8images})]
    logging.info(f' - L8 = {Nb_L8images}')

    # --- Extract L9 collection for the specific landsat tile ---
    L9_tile = L9collection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    L9fullcollection_tile = L9fullcollection.filter(ee.Filter.stringContains('LANDSAT_PRODUCT_ID', tile_L))
    Nb_L9images = googleErrorsControl(L9_tile.size(), path2key)
    if Nb_L9images != 0:
        date_start_L9 = googleErrorsControl(L9_tile.limit(1, 'system:time_start', True).first().get('DATE_ACQUIRED'), path2key)
        date_end_L9 = googleErrorsControl(L9_tile.limit(1, 'system:time_start', False).first().get('DATE_ACQUIRED'), path2key)
    else:
        date_start_L9 = ''
        date_end_L9 = ''
    Count_dict += [ee.Feature(None, {'TILE': tile_L,
                                     'PRODUCT': 'L9',
                                     'DATE START': date_start_L9,
                                     'DATE END': date_end_L9,
                                     'NUMBER OF IMAGES': Nb_L9images})]
    logging.info(f' - L9 = {Nb_L9images}')

    # --- Extract S2 collection(s) for the specific landsat tile (grid) ---
    S2_alltiles = S2collection.filterBounds(landsat_grid)
    tiles_S2_L = listSentinelTiles(S2_alltiles, path2key)
    
    if tiles_S2_L!=[]:

        if tiles_S2!=['']:
            S2_alltiles_new = ee.ImageCollection([])
            tiles_S2_new = []
            for t_s2 in tiles_S2:
                if t_s2 in tiles_S2_L:
                    S2_tile = S2_alltiles.filter(ee.Filter.stringContains('MGRS_TILE', t_s2))
                    S2_alltiles_new = S2_alltiles_new.merge(S2_tile)
                    tiles_S2_new.append(t_s2)
                    del S2_tile
            S2_alltiles = S2_alltiles_new
            tiles_S2_L = tiles_S2_new

        Nb_S2images =  googleErrorsControl(S2_alltiles.size(), path2key)
        date_start_S2 = googleErrorsControl(S2_alltiles.limit(1, 'system:time_start', True).first().date().format('YYYY-MM-dd'), path2key)
        date_end_S2 = googleErrorsControl(S2_alltiles.limit(1, 'system:time_start', False).first().date().format('YYYY-MM-dd'), path2key)
        
    else:
        Nb_S2images = 0
        date_start_S2 = ''
        date_end_S2 = ''
        
    Count_dict += [ee.Feature(None, {'TILE': tile_L,
                                     'PRODUCT': 'S2',
                                     'DATE START': date_start_S2,
                                     'DATE END': date_end_S2,
                                     'NUMBER OF IMAGES': Nb_S2images})]
    logging.info(f' - S2 = {Nb_S2images}')

    # --- Extract projection (landsat full) and period (landsat and s2) ---
    landsat_fulldataset = L7fullcollection_tile.merge(L8fullcollection_tile).merge(L9fullcollection_tile)
    proj_landsat = googleErrorsControl(landsat_fulldataset.first().select('SR_B4').projection(), path2key)
    all_dataset = L7_tile.merge(L8_tile).merge(L9_tile).merge(S2_alltiles)
    date_start = all_dataset.limit(1, 'system:time_start', True).first().date()
    date_end = all_dataset.limit(1, 'system:time_start', False).first().date()
            
    return (L7_tile, L8_tile, L9_tile, S2_alltiles, tiles_S2_L,
            proj_landsat, date_start, date_end, Count_dict)



def exportDataFrame(drive_folder, table_dict, filename, columns, export_folder=os.getcwd(), path2key=os.getcwd()):
    """
    Export dataframe to local machine (drive, then local)
    """
    
    # Converts dict to feature collection
    table = ee.FeatureCollection(table_dict)
    
    # Export to Drive
    task = ee.batch.Export.table.toDrive(**{
        'collection': table,
        'description': filename,
        'folder': drive_folder,
        'selectors': columns,
        'fileFormat': 'CSV'})
    
    # GOOGLE ERROR CONTROL PROCEDURE for GEE EXPORT TASK (in case of loosing google connection)
    MAX_GEETRY = 6
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:
        try:
            # Run task
            logging.info(f'\nExporting table to drive : {filename}')
            task.start()
            while task.active():
                print("-", end="", flush=True)
                time.sleep(0.5)
                print("\b", end="", flush=True)
                time.sleep(0.5)
            google_error = 'no_google_error'

        except Exception as e:
            # Reconnect and wait
            if n_geetry<MAX_GEETRY:
                logging.warning(f'Loose google connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                googleAuthentification(path2key)
            else:
                logging.critical(f'Error google authentification ({e})')
                raise Exception ('Error google authentification')
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)

    # Authenticate to Google Drive (of the Service account)
    gauth = GoogleAuth()
    scopes = ['https://www.googleapis.com/auth/drive']
    key_file = glob.glob(os.path.join(path2key,'*.json'))[0]
    gauth.credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file, scopes=scopes)
    drive = GoogleDrive(gauth)
    
    # GOOGLE ERROR CONTROL PROCEDURE for DRIVE EXPORT TASK (in case of loosing google connection)
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:
        try:
            # Extract google drive file list
            file_list = drive.ListFile({'q': "'root' in parents and trashed=false"}).GetList()
            
            # Retrieve the folder id - start searching from root
            # Choose starting point by inserting folder name
            folder_title = drive_folder
            folder_id = ''
            for file in file_list:
                if(file['title']==folder_title):
                    folder_id = file['id']
                    break
            
            # Build string dynamically (need to use escape characters to support single quote syntax)
            str_fold = "\'" + folder_id + "\'" + " in parents and trashed=false"    
            
            # Iterating over files and downloading
            file_list = drive.ListFile({'q': str_fold}).GetList()
            for file in file_list:
                filename = file['title']
                logging.info(f'\nLocal downloading : {filename}')
            
                # download file into working directory
                file.GetContentFile(os.path.join(export_folder, file['title']), mimetype="table")
            
                # delete file afterwards to keep the Drive empty
                file.Delete()

            google_error = 'no_google_error'

        except Exception as e:
            if n_geetry<MAX_GEETRY:
                logging.warning(f'Loose google drive connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                gauth = GoogleAuth()
                scopes = ['https://www.googleapis.com/auth/drive']
                key_file = glob.glob(os.path.join(path2key,'*.json'))[0]
                gauth.credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file, scopes=scopes)
                drive = GoogleDrive(gauth)
            else:
                logging.critical(f'Error google authentification ({e})')
                raise Exception ('Error google authentification')
    
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)



def exportImageToAsset(data, data_path, path2key=os.getcwd(), data_crs=None, data_transform=None, data_scale=None, data_region=None):
    """
    Export image to google asset
    """

    task = ee.batch.Export.image.toAsset(**{
                'image': data,
                'assetId': data_path,
                'crs': data_crs,
                'crsTransform':data_transform,
                'scale': data_scale,
                'region': data_region,
                'maxPixels': 1e10
                })

    # GOOGLE ERROR CONTROL PROCEDURE for GEE EXPORT TASK (in case of loosing google connection)
    MAX_GEETRY = 6
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:
        try:
            # Run task
            task.start()
            while task.active():
                print("-", end="", flush=True)
                time.sleep(0.5)
                print("\b", end="", flush=True)
                time.sleep(0.5)
            google_error = 'no_google_error'

        except Exception as e:
            # Reconnect and wait
            if n_geetry<MAX_GEETRY:
                logging.warning(f'Loose google connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                googleAuthentification(path2key)
            else:
                logging.critical(f'Error google authentification ({e})')
                raise Exception ('Error google authentification')
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)



def exportImage(drive_folder, data, filename, export_folder=os.getcwd(), path2key=os.getcwd(), data_crs=None, data_transform=None, data_scale=None, data_region=None):
    """
    Export image to local machine (drive, then local)
    """

    if data_transform is not None: data_scale=None

    task = ee.batch.Export.image.toDrive(**{
                'image': data,
                'description': filename,
                'folder': drive_folder,
                'crs': data_crs,
                'crsTransform': data_transform,
                'scale': data_scale,
                'fileFormat': 'GeoTIFF',
                'region': data_region,
                'maxPixels': 1e10
                })
    
    # GOOGLE ERROR CONTROL PROCEDURE for GEE EXPORT TASK (in case of loosing google connection)
    MAX_GEETRY = 6
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:
        try:
            # Run task
            logging.info(f'\nExporting image to drive : {filename}')
            task.start()
            while task.active():
                print("-", end="", flush=True)
                time.sleep(0.5)
                print("\b", end="", flush=True)
                time.sleep(0.5)
            google_error = 'no_google_error'
        
        except Exception as e:
            # Reconnect and wait
            if n_geetry<MAX_GEETRY:
                logging.warning(f'Loose google connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                googleAuthentification(path2key)
            else:
                logging.critical(f'Error google authentification ({e})')
                raise Exception ('Error google authentification')
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)
            
    # Authenticate to Google Drive (of the Service account) :
    gauth = GoogleAuth()
    scopes = ['https://www.googleapis.com/auth/drive']
    key_file = glob.glob(os.path.join(path2key,'*.json'))[0]
    gauth.credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file, scopes=scopes)
    drive = GoogleDrive(gauth)

    # GOOGLE ERROR CONTROL PROCEDURE for DRIVE EXPORT TASK (in case of loosing google connection)
    google_error = 'google_error'
    n_geetry = 1

    while google_error!='no_google_error' and n_geetry<=MAX_GEETRY:
        try:
            # Extract google drive file list
            file_list = drive.ListFile({'q': "'root' in parents and trashed=false"}).GetList()
            
            # Retrieve the folder id - start searching from root
            # Choose starting point by inserting folder name
            folder_title = drive_folder
            folder_id = ''
            for file in file_list:
                if(file['title']==folder_title):
                    folder_id = file['id']
                    break

            # Build string dynamically (need to use escape characters to support single quote syntax) :
            str_fold = "\'" + folder_id + "\'" + " in parents and trashed=false"    

            # Iterating over files and downloading :
            file_list = drive.ListFile({'q': str_fold}).GetList()

            for file in file_list:
                filename = file['title']
                logging.info(f'\nLocal downloading : {filename}')
                
                # download file into working directory (in this case a tiff-file)
                file.GetContentFile(os.path.join(export_folder, file['title']), mimetype="image/tiff")

                # delete file afterwards to keep the Drive empty
                file.Delete()

            google_error = 'no_google_error'

        except Exception as e:
            # Reconnect and wait
            if n_geetry<MAX_GEETRY:
                logging.warning(f'Loose google drive connection ({e}) -> Re-authentification number {n_geetry}/{MAX_GEETRY-1}')
                gauth = GoogleAuth()
                scopes = ['https://www.googleapis.com/auth/drive']
                key_file = glob.glob(os.path.join(path2key,'*.json'))[0]
                gauth.credentials = ServiceAccountCredentials.from_json_keyfile_name(key_file, scopes=scopes)
                drive = GoogleDrive(gauth)
            else:
                logging.critical(f'Error google authentification ({e})')
                raise Exception ('Error google authentification')
    
            google_error = 'google_error'
            n_geetry+=1
            time.sleep(5**n_geetry)



def removeOutliers(indices):
    """
    Remouving outlier to indices :
        - water bodies
        - values outside "valid" range
    """
    
    watermask = indices.select('NDVI').gte(0.03)
    outliersmask_ndvi = watermask.And(indices.select('NDVI').lte(1))
    outliersmask_ndwi = watermask.And(indices.select('NDWI').lte(1))
    
    indices_filt = indices.updateMask([outliersmask_ndwi,outliersmask_ndvi])

    return indices_filt



def computeNanScores(image, path2key, landmask, data_geom=None, scale=None, bands_list=None):
    """
    Compute the NaN scores of image bands according to specific landmask
    """
    
    if data_geom is None: data_geom = image.select(1).geometry()
    
    if scale is None: scale = image.select(1).projection().nominalScale()
    
    if bands_list is None: bands_list = googleErrorsControl(image.bandNames(), path2key)
    else: image = image.select(bands_list)

    # --- Extract number of land pixels ---
    Nb_ALLDATA_land = landmask.reduceRegion(**{
                                    'reducer': ee.Reducer.sum(),
                                    'geometry': data_geom,
                                    'scale': scale,
                                    'bestEffort': True})
    Nb_ALLDATA_land = googleErrorsControl(Nb_ALLDATA_land.getNumber(Nb_ALLDATA_land.keys().getString(0)), path2key)

    # --- Extract number of NaN pixels on land (for bands_list) ---
    Nb_NANDATA_land = googleErrorsControl(
                            image.mask()
                            .eq(0).And(landmask)
                            .reduceRegion(**{
                                'reducer': ee.Reducer.sum(),
                                'geometry': data_geom,
                                'scale': scale,
                                'bestEffort': True}), path2key)
    
    # --- Compute scores ---
    nanscores = Nb_NANDATA_land
    if Nb_ALLDATA_land!=0:
        for b in bands_list:
            nanscores[b] = round(Nb_NANDATA_land[b] / Nb_ALLDATA_land, 2)
    else:
        for b in bands_list:
            nanscores[b] = 1

    return nanscores



def calibrateData(indices, productName):
    """
    Radiometric calibration of satellite data to L8 radiometric values
    """

    # DEFINE CALIB PARAMETERS :
    if productName=='L7':
        a_ndwi = 0.97
        a_ndvi = 0.96
        b_ndwi = 0.04
        b_ndvi = 0.07
    elif productName=='L9':
        a_ndwi = 0.95
        a_ndvi = 0.97
        b_ndwi = 0.02
        b_ndvi = 0.03
    elif productName=='S2':
        a_ndwi = 0.98
        a_ndvi = 0.91
        b_ndwi = 0.07
        b_ndvi = 0.11
    
    # APPLY CALIB :
    NDWI_calib = indices.select('NDWI').multiply(a_ndwi).add(b_ndwi)
    NDVI_calib = indices.select('NDVI').multiply(a_ndvi).add(b_ndvi)
    indices_calib = ee.Image([NDWI_calib, NDVI_calib])

    # REMOVE OUTLIERS :
    indices_calib = removeOutliers(indices_calib)

    return indices_calib


