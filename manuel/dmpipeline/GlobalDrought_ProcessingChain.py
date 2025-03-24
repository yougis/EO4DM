# -*- coding: utf-8 -*-
"""
##############################################################################

Global Drought Processing chain

##############################################################################
Created on Mon June 26

@author: Mathis Neuhauser
"""

import sys
import dotenv
import dmpipeline.GEE_Processing.GEE_main_functions as geemain
import dmpipeline.DROUGHT_Processing.DROUGHT_global_functions as dglob

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def main():

    logging.info('\n\n\n =========== START GLOBAL SCRIPT =========== \n')
    
    CONFIG = dotenv.dotenv_values('Config_process.env')
    CONFIG['WRK_DIR'] = dglob.prepare_DATAFolders(CONFIG)
    go_modis = 1


    # ========================================= AUTOMATIC / MANUAL SETTING OF PROCESSING PERIOD =====================================

    go_process, CONFIG['PERIOD_START'], CONFIG['PERIOD_END'] = dglob.auto_processingPERIOD(CONFIG)

    if go_process==1:

        # =========================================== INDICES PROCESSING (NDWI, LST) ===============================================
        
        if ('AUTO' in CONFIG['MODE']) or ('MANUAL' in CONFIG['MODE']) or ('INDICES' in CONFIG['MODE']):
            (go_modis, go_month,
            CONFIG['PERIOD_START'], CONFIG['PERIOD_END'], CONFIG['TERRITORY']) = geemain.check_GEEProducts_MODIS(CONFIG)
        
            if go_modis==1:  
                OUTDIR_PATHS = dglob.prepare_RUNFolders(CONFIG)
                (CONFIG['lon_min_modis'], CONFIG['lat_min_modis'],
                 CONFIG['lon_max_modis'], CONFIG['lat_max_modis'], CONFIG['TRAN_OUT']) = dglob.auto_processingBBOX(CONFIG)
                COLLECTIONS, PARAM, DICT = geemain.initialize_GEEProcessing_MODIS(CONFIG, OUTDIR_PATHS[3])
                geemain.run_GEEProcessing_MODIS(CONFIG, OUTDIR_PATHS, COLLECTIONS, PARAM, DICT, go_month)

        # ================================================ DROUGHT PROCESSING (VHI) ================================================

        if ((go_modis==1) and ('AUTO' in CONFIG['MODE'] or 'MANUAL' in CONFIG['MODE'])) or ('DROUGHT' in CONFIG['MODE']):
            
            if ('DROUGHT' in CONFIG['MODE']):
                OUTDIR_PATHS = dglob.prepare_RUNFolders(CONFIG)
            dglob.process_GlobalDrought_VHI(CONFIG, OUTDIR_PATHS)
    
    
    logging.info('\n\n =========== END GLOBAL SCRIPT =========== \n')


if __name__ == "__main__":
    sys.exit(main())