# -*- coding: utf-8 -*-
"""
##############################################################################

Alert Drought Processing chain

##############################################################################
Created on Tue Oct 31

@author: Mathis Neuhauser
"""

import sys
import dotenv
import dmpipeline.ALERT_Processing.ALERT_data_access as alertaccess
import dmpipeline.ALERT_Processing.ALERT_main_functions as alertmain

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def main():

    logging.info('\n\n\n =========== START ALERT SCRIPT =========== \n')
    
    CONFIG = dotenv.dotenv_values('Config_process.env')
    CONFIG['WRK_DIR']  = alertmain.prepare_DATAFolders(CONFIG)
    go_products = 1


    # ============================================ DOWNLOAD NEW COPERNICUS/METEO PRODUCTS ====================================
    
    alertaccess.download_newProducts(CONFIG)
    
    # ======================================== AUTOMATIC / MANUAL SETTING OF PROCESSING PERIOD ===============================

    go_process, CONFIG['PERIOD_START'], CONFIG['PERIOD_END'] = alertmain.auto_processingPERIOD(CONFIG)

    if go_process==1:

        # ===================================== INDICES (SWI) and MOISTURE DEFICIT (MAI) PROCESSING ===========================

        if ('AUTO' in CONFIG['MODE']) or ('MANUAL' in CONFIG['MODE']) or ('INDICES' in CONFIG['MODE']):
            (go_products, COLLECTION_ASCAT, OUTDIR_PATHS) = alertmain.CheckInit_Products_ALERT(CONFIG)
            
            if go_products==1:
                alertmain.processingASCAT(CONFIG, OUTDIR_PATHS, COLLECTION_ASCAT)
                alertmain.process_MoistureDeficit_MAI(CONFIG, OUTDIR_PATHS)

        # ================================================= DROUGHT ALERTS PROCESSING ==========================================

        if ((go_products==1) and ('AUTO' in CONFIG['MODE'] or 'MANUAL' in CONFIG['MODE'])) or ('DROUGHT' in CONFIG['MODE']):
            
            if ('DROUGHT' in CONFIG['MODE']):
                OUTDIR_PATHS = alertmain.prepare_RUNFolders(CONFIG)
            alertmain.process_DroughtAlert(CONFIG, OUTDIR_PATHS)
    
    
    logging.info('\n\n =========== END ALERT SCRIPT =========== \n')


if __name__ == "__main__":
    sys.exit(main())