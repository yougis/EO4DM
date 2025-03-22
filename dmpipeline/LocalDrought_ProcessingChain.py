# -*- coding: utf-8 -*-
"""
##############################################################################

Local Drought Processing chain

##############################################################################
Created on Mon June 26

@author: Mathis Neuhauser
"""

import sys
import dotenv
import dmpipeline.GEE_Processing.GEE_main_functions as geemain
import dmpipeline.DROUGHT_Processing.DROUGHT_local_functions as dloc

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def main():
    
    logging.info('\n\n\n =========== START LOCAL SCRIPT =========== \n')

    CONFIG = dotenv.dotenv_values('Config_process.env')
    CONFIG['WRK_DIR'] = dloc.prepare_DATAFolders(CONFIG)
    go_landsat_s2 = 1
    

    # ========================================= AUTOMATIC / MANUAL SETTING OF PROCESSING PERIOD =====================================

    go_process, CONFIG['PERIOD_START'], CONFIG['PERIOD_END'] = dloc.auto_processingPERIOD(CONFIG)

    if go_process==1 :

        # =========================================== INDICES PROCESSING (NDWI, NDVI) ===============================================
        
        if ('AUTO' in CONFIG['MODE']) or ('MANUAL' in CONFIG['MODE']) or ('INDICES' in CONFIG['MODE']):
            (go_landsat_s2, go_month,
            CONFIG['PERIOD_START'], CONFIG['PERIOD_END'],
            CONFIG['TERRITORY'], CONFIG['TILES_L'], CONFIG['TILES_S2']) = geemain.check_GEEProducts_LANDSAT_S2(CONFIG)

            if go_landsat_s2==1:
                OUTDIR_PATHS = dloc.prepare_RUNFolders(CONFIG)
                COLLECTIONS, PARAM, DICT = geemain.initialize_GEEProcessing_LANDSAT_S2(CONFIG, OUTDIR_PATHS[3])
                geemain.run_GEEProcessing_LANDSAT_S2(CONFIG, OUTDIR_PATHS, COLLECTIONS, PARAM, DICT)
        
        # ================================================= DROUGHT PROCESSING (VAI) ================================================

        if ((go_landsat_s2==1) and ('AUTO' in CONFIG['MODE'] or 'MANUAL' in CONFIG['MODE'])) or ('DROUGHT' in CONFIG['MODE']):

            if ('DROUGHT' in CONFIG['MODE']):
                OUTDIR_PATHS = dloc.prepare_RUNFolders(CONFIG)
            dloc.process_LocalDrought_VAI(CONFIG, OUTDIR_PATHS)
    
    
    logging.info('\n\n =========== END LOCAL SCRIPT =========== \n')


if __name__ == "__main__":
    sys.exit(main())