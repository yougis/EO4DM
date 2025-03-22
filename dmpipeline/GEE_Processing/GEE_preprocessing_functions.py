# -*- coding: utf-8 -*-
"""
##############################################################################

GEE Functions to Preprocess satellite products

##############################################################################
"""


import os
import ee
import time
import dmpipeline.GEE_Processing.GEE_generic_functions as geegen

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def computeMODIS_LST(image, grid, landmask):
    """
    Computing MODIS LST :
        - clipping to grid
        - converting to °C
        - applying landmask
    """

    M = 0.02
    A = -273.15

    # CLIP TO GRID
    image = image.clip(grid)

    # CONVERTS TO °C
    LST = image.select('LST_Day_1km').multiply(M).add(A).rename('LST')

    # APPLY LANDMASK
    LST = LST.updateMask(landmask)  

    return LST



def computeMODISv21A1D_LST(image, grid, landmask):
    """
    Computing MODIS LST from products MYD21A1D (adapted to VIIRS):
        - clipping to grid
        - converting to °C
        - applying landmask
        - removing extreme high temperature outliers (no default filter in MYD21A1D) :
            -> silicaterock_mask (observed wrong emissity on these areas)
            -> threshold_mask (LST>60 removed) for remaining outliers
    """

    A = -273.15
    LST_THRESH = 60

    # CLIP TO GRID
    image = image.clip(grid)

    # QC MASK (silicaterock_mask)
    qc = image.select('QC')
    qc_bit10 = 1 << 10
    qc_bit11 = 1 << 11
    silicaterock_mask = ((qc.bitwiseAnd(qc_bit10).eq(0))
                         .And(qc.bitwiseAnd(qc_bit11).eq(0)))
    silicaterock_mask = silicaterock_mask.unmask().Not()
    
    # CONVERTS LST TO °C
    LST = image.select('LST_1KM').add(A).rename('LST')
    
    # THRESHOL EXTREME HIGH LST
    threshold_mask = LST.lte(LST_THRESH).clip(image.geometry())
    
    # COMBINE ALL MASKS (+ LANDMASK) :
    QMASK = (silicaterock_mask.And(threshold_mask).And(landmask).rename('QMASK'))

    # APPLY QMASK
    LST = LST.updateMask(QMASK)  

    return LST



def computeMODIS_NDWI(image, grid, landmask):
    """
    Computing MODIS NDWI :
        - clipping to grid
        - quality masking (clouds, shadows, cirrus, water/land)
        - computing ndwi
    """
    
    FILL_VALUE = -28672

    # CLIP TO GRID
    image = image.clip(grid)
    
    # FILL MASK
    fmask_bands = image.select(['sur_refl_b02','sur_refl_b06']).eq(FILL_VALUE).rename(['fmask_b02', 'fmask_b06'])
    fmask_indice = fmask_bands.select('fmask_b02').Or(fmask_bands.select('fmask_b06')).eq(0).rename('fmask')
    
    # QC MASK : include pixels that are only produced at ideal quality (bit0=0 and bit1=0)
    qc = image.select('QC_500m')
    qc_bit0 = 1 << 0
    qc_bit1 = 1 << 1
    globalmask = ((qc.bitwiseAnd(qc_bit0).eq(0))
                    .And(qc.bitwiseAnd(qc_bit1).eq(0)))
    
    # STATE_QA mask : clouds=bit0-1, shadows=bit2, cirrus=bit8-9
    state_qa = image.select('state_1km')
    cloudBitMask_0 = 1 << 0
    cloudBitMask_1 = 1 << 1
    shadowBitmask = 1 << 2
    cirrusBitmask_0 = 1 << 8
    cirrusBitmask_1 = 1 << 9
    cloudmask = ((state_qa.bitwiseAnd(cloudBitMask_0).eq(0))
                    .And(state_qa.bitwiseAnd(cloudBitMask_1).eq(0)))
    shadowmask = state_qa.bitwiseAnd(shadowBitmask).eq(0)
    cirrusmask = ((state_qa.bitwiseAnd(cirrusBitmask_0).eq(0))
                    .And(state_qa.bitwiseAnd(cirrusBitmask_1).eq(0)))

    # COMBINE ALL MASKS (+ LANDMASK) :
    QMASK = (fmask_indice.And(globalmask).And(cloudmask).And(shadowmask)
             .And(cirrusmask).And(landmask).rename('QMASK'))

    # COMPUTE NDWI
    ndwi = image.normalizedDifference(['sur_refl_b02','sur_refl_b06']).rename('NDWI')
    image_concat = image.addBands(ndwi)
    NDWI = image_concat.select('NDWI')

    # APPLY QMASK
    NDWI = NDWI.updateMask(QMASK)

    return NDWI



def preprocessingMODIS(lst_aqua, lst_terra, reflect, landmask, new_collection, QA_dict, crs_out, grid_out, scale_out, path2key=os.getcwd(), asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single MODIS product :
        1) computing indices
        2) quality masking (clouds, shadows, landmask)
        3) (IF asset_export==1 and nanscore<=95%) exporting preprocessed images to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)
    """  

    lst_aqua_ok = geegen.googleErrorsControl(lst_aqua, path2key)!=None
    lst_terra_ok = geegen.googleErrorsControl(lst_terra, path2key)!=None
    reflect_ok = geegen.googleErrorsControl(reflect, path2key)!=None

    if (lst_aqua_ok==0 and lst_terra_ok==0) or reflect_ok==0:
        return 'not considered', QA_dict
    
    product_date = geegen.googleErrorsControl(reflect.date().format('YYYY-MM-dd'), path2key)
    product_id = geegen.googleErrorsControl(reflect.id(), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE (compute nanscores) and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'MODIS {product_date} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])
                nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)
                QA_dict += [ee.Feature(None, {'DATE': product_date,
                                            'NAN SCORE LST': nanscores['LST'],
                                            'NAN SCORE NDWI' : nanscores['NDWI'],
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict
    
    # --- COMPUTE MASKED LST ---       
    if lst_aqua_ok==1:
        lst_aqua_preproc = computeMODIS_LST(lst_aqua, grid_out, landmask)
        lst_preproc = lst_aqua_preproc.float()
    if lst_terra_ok==1:
        lst_terra_preproc = computeMODIS_LST(lst_terra, grid_out, landmask)
        lst_preproc = lst_terra_preproc.float()

    # --- COMPOSITE LST (if Aqua AND Terra) ---
    if lst_aqua_ok==1 and lst_terra_ok==1:
        lst_preproc = ee.ImageCollection([lst_aqua_preproc, lst_terra_preproc]).mean().float().clip(grid_out)
    
    # --- COMPUTE MASKED NDWI ---
    ndwi_preproc = computeMODIS_NDWI(reflect, grid_out, landmask).float()
    
    # --- APPLY COMMON MASK TO BOTH INDICES
    qmask_lst = lst_preproc.mask()
    qmask_ndwi = ndwi_preproc.mask()
    qmask_all = qmask_lst.And(qmask_ndwi)
    lst_preproc = lst_preproc.updateMask(qmask_all)
    ndwi_preproc = ndwi_preproc.updateMask(qmask_all)

    # --- CONCATENATE INDICES and COMPUTE NANSCORES ---
    indices_preproc = ndwi_preproc.addBands(lst_preproc)
    nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)

    if (nanscores['LST']<=0.95 and nanscores['NDWI']<=0.95):
        
        logging.info(f'MODIS {product_date} : preprocessing product')

        # --- EXPORT TO ASSET ---
        if asset_export==1:
            MAX_GEETRY = 6
            gee_test=''
            n_geetry = 1

            while gee_test=='' and n_geetry<=MAX_GEETRY:
                start_time = time.time()
                geegen.exportImageToAsset(indices_preproc, preproc_filename, path2key, data_crs=crs_out, data_scale=scale_out, data_region=grid_out)
                elapsed_time = round(time.time() - start_time)

                # Read already exported products to improve computation time
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])

                # Test to ensure good asset export (if not, retry exporting)
                gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                if gee_test=='':
                    n_geetry+=1
                    time.sleep(5**n_geetry)

        else: elapsed_time = ''         
        
        # --- CONCATENATE SCORES INTO TABLE ---
        QA_dict += [ee.Feature(None, {'DATE': product_date,
                                      'NAN SCORE LST': nanscores['LST'],
                                      'NAN SCORE NDWI' : nanscores['NDWI'],
                                      'PREPROC TIME (sec)': elapsed_time})]
        
        return indices_preproc, QA_dict
        
    else:
        logging.info(f'MODIS {product_date} : not considered due to nanscores > 95%')
        return 'not considered', QA_dict



def preprocessingMODISv21A1D(lst, reflect, landmask, new_collection, QA_dict, crs_out, grid_out, scale_out, path2key=os.getcwd(), asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single MODIS product (version for combining with VIIRS products) :
        1) computing indices
        2) quality masking (clouds, shadows, landmask)
        3) (IF asset_export==1 and nanscore<=95%) exporting preprocessed images to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)

    Version v21A1D :
        - uses LST MYD21A1D products (Aqua) for combining with LST VIIRS VNP21A1D
        - calibrate MODIS on VIIRS (A COMPLETER)
    """

    lst_ok = geegen.googleErrorsControl(lst, path2key)!=None
    reflect_ok = geegen.googleErrorsControl(reflect, path2key)!=None

    if lst_ok==0 or reflect_ok==0:
        return 'not considered', QA_dict
    
    product_date = geegen.googleErrorsControl(reflect.date().format('YYYY-MM-dd'), path2key)
    product_id = geegen.googleErrorsControl(reflect.id(), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE (compute nanscores) and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'MODIS {product_date} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])
                nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)
                QA_dict += [ee.Feature(None, {'DATE': product_date,
                                            'NAN SCORE LST': nanscores['LST'],
                                            'NAN SCORE NDWI' : nanscores['NDWI'],
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict
    
    # --- COMPUTE MASKED LST ---       
    lst_preproc = computeMODISv21A1D_LST(lst, grid_out, landmask).float()

    # --- COMPUTE MASKED NDWI ---
    ndwi_preproc = computeMODIS_NDWI(reflect, grid_out, landmask).float()
    
    # --- APPLY COMMON MASK TO BOTH INDICES
    qmask_lst = lst_preproc.mask()
    qmask_ndwi = ndwi_preproc.mask()
    qmask_all = qmask_lst.And(qmask_ndwi)
    lst_preproc = lst_preproc.updateMask(qmask_all)
    ndwi_preproc = ndwi_preproc.updateMask(qmask_all)

    # --- CONCATENATE INDICES and COMPUTE NANSCORES ---
    indices_preproc = ndwi_preproc.addBands(lst_preproc)
    nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)

    if (nanscores['LST']<=0.95 and nanscores['NDWI']<=0.95):
        
        logging.info(f'MODIS {product_date} : preprocessing product')

        # --- EXPORT TO ASSET ---
        if asset_export==1:
            MAX_GEETRY = 6
            gee_test=''
            n_geetry = 1

            while gee_test=='' and n_geetry<=MAX_GEETRY:
                start_time = time.time()
                geegen.exportImageToAsset(indices_preproc, preproc_filename, path2key, data_crs=crs_out, data_scale=scale_out, data_region=grid_out)
                elapsed_time = round(time.time() - start_time)

                # Read already exported products to improve computation time
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])

                # Test to ensure good asset export (if not, retry exporting)
                gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                if gee_test=='':
                    n_geetry+=1
                    time.sleep(5**n_geetry)

        else: elapsed_time = ''         
        
        # --- CONCATENATE SCORES INTO TABLE ---
        QA_dict += [ee.Feature(None, {'DATE': product_date,
                                      'NAN SCORE LST': nanscores['LST'],
                                      'NAN SCORE NDWI' : nanscores['NDWI'],
                                      'PREPROC TIME (sec)': elapsed_time})]
        
        return indices_preproc, QA_dict
        
    else:
        logging.info(f'MODIS {product_date} : not considered due to nanscores > 95%')
        return 'not considered', QA_dict



def computeVIIRS_LST(image, grid, landmask):
    """
    Computing VIIRS LST :
        - clipping to grid
        - converting to °C
        - applying landmask
    """

    A = -273.15

    # CLIP TO GRID
    image = image.clip(grid)

    # CONVERTS TO °C
    LST = image.select('LST_1KM').add(A).rename('LST')

    # APPLY LANDMASK
    LST = LST.updateMask(landmask)  

    return LST



def computeVIIRS_NDWI(image, grid, landmask):
    """
    Computing VIIRS NDWI :
        - clipping to grid
        - quality masking (clouds, shadows, cirrus, water/land)
        - computing ndwi
    """
    
    FILL_VALUE = -28672

    # CLIP TO GRID
    image = image.clip(grid)
    
    # FILL MASK
    fmask_bands = image.select(['I2','I3']).eq(FILL_VALUE).rename(['fmask_I2', 'fmask_I3'])
    fmask_indice = fmask_bands.select('fmask_I2').Or(fmask_bands.select('fmask_I3')).eq(0).rename('fmask')
    
    # QF1 MASK
    qf1 = image.select('QF1')
    qf1_bit2 = 1 << 2
    qf1_bit3 = 1 << 3
    cloudmask_detect = ((qf1.bitwiseAnd(qf1_bit2).eq(0))
                          .And(qf1.bitwiseAnd(qf1_bit3).eq(0)))
    
    # QF2 MASK
    qf2 = image.select('QF2')
    qf2_bit3 = 1 << 3
    qf2_bit4 = 1 << 4
    qf2_bit6 = 1 << 6
    qf2_bit7 = 1 << 7
    shadowmask = qf2.bitwiseAnd(qf2_bit3).eq(0)
    aeromask = qf2.bitwiseAnd(qf2_bit4).eq(0)
    cirrusrefl_mask = qf2.bitwiseAnd(qf2_bit6).eq(0)
    cirrusemiss_mask = qf2.bitwiseAnd(qf2_bit7).eq(0)

    # QF4 MASK
    qf4 = image.select('QF4')
    qf4_bit2 = 1 << 2
    qf4_bit3 = 1 << 3
    badI2_mask = qf4.bitwiseAnd(qf4_bit2).eq(0)
    badI3_mask = qf4.bitwiseAnd(qf4_bit3).eq(0)
    
    # QF7 MASK
    qf7 = image.select('QF7')
    qf7_bit4 = 1 << 4
    cirrusflag_mask = qf7.bitwiseAnd(qf7_bit4).eq(0)

    # COMBINE ALL MASKS (+ LANDMASK) :
    QMASK = (fmask_indice.And(cloudmask_detect).And(shadowmask)
              .And(aeromask).And(cirrusrefl_mask).And(cirrusemiss_mask)
              .And(badI2_mask).And(badI3_mask).And(cirrusflag_mask).And(landmask).rename('QMASK'))

    # COMPUTE NDWI
    ndwi = image.normalizedDifference(['I2','I3']).rename('NDWI')
    image_concat = image.addBands(ndwi)
    NDWI = image_concat.select('NDWI')

    # APPLY QMASK
    NDWI = NDWI.updateMask(QMASK)

    return NDWI



def preprocessingVIIRS(lst, reflect, landmask, new_collection, QA_dict, crs_out, grid_out, scale_out, path2key=os.getcwd(), asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single VIIRS product :
        1) computing indices
        2) quality masking (clouds, shadows, landmask)
        3) (IF asset_export==1 and nanscore<=95%) exporting preprocessed images to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)
    """  

    lst_ok = geegen.googleErrorsControl(lst, path2key)!=None
    reflect_ok = geegen.googleErrorsControl(reflect, path2key)!=None

    if lst_ok==0 or reflect_ok==0:
        return 'not considered', QA_dict
    
    product_date = geegen.googleErrorsControl(reflect.date().format('YYYY-MM-dd'), path2key)
    product_id = geegen.googleErrorsControl(reflect.id(), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE (compute nanscores) and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'VIIRS {product_date} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])
                nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)
                QA_dict += [ee.Feature(None, {'DATE': product_date,
                                            'NAN SCORE LST': nanscores['LST'],
                                            'NAN SCORE NDWI' : nanscores['NDWI'],
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict
    
    # --- COMPUTE MASKED LST ---       
    lst_preproc = computeVIIRS_LST(lst, grid_out, landmask).float()

    # # --- COMPOSITE LST (if Aqua AND Terra) ---
    # if lst_aqua_ok==1 and lst_terra_ok==1:
    #     lst_preproc = ee.ImageCollection([lst_aqua_preproc, lst_terra_preproc]).mean().float().clip(grid_out)
    
    # --- COMPUTE MASKED NDWI ---
    ndwi_preproc = computeVIIRS_NDWI(reflect, grid_out, landmask).float()
    
    # --- APPLY COMMON MASK TO BOTH INDICES
    qmask_lst = lst_preproc.mask()
    qmask_ndwi = ndwi_preproc.mask()
    qmask_all = qmask_lst.And(qmask_ndwi)
    lst_preproc = lst_preproc.updateMask(qmask_all)
    ndwi_preproc = ndwi_preproc.updateMask(qmask_all)

    # --- CONCATENATE INDICES and COMPUTE NANSCORES ---
    indices_preproc = ndwi_preproc.addBands(lst_preproc)
    nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask, scale=scale_out)

    if (nanscores['LST']<=0.95 and nanscores['NDWI']<=0.95):
        
        logging.info(f'VIIRS {product_date} : preprocessing product')

        # --- EXPORT TO ASSET ---
        if asset_export==1:
            MAX_GEETRY = 6
            gee_test=''
            n_geetry = 1

            while gee_test=='' and n_geetry<=MAX_GEETRY:
                start_time = time.time()
                geegen.exportImageToAsset(indices_preproc, preproc_filename, path2key, data_crs=crs_out, data_scale=scale_out, data_region=grid_out)
                elapsed_time = round(time.time() - start_time)

                # Read already exported products to improve computation time
                indices_preproc = ee.Image(preproc_filename).select(['LST', 'NDWI'])

                # Test to ensure good asset export (if not, retry exporting)
                gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                if gee_test=='':
                    n_geetry+=1
                    time.sleep(5**n_geetry)

        else: elapsed_time = ''         
        
        # --- CONCATENATE SCORES INTO TABLE ---
        QA_dict += [ee.Feature(None, {'DATE': product_date,
                                      'NAN SCORE LST': nanscores['LST'],
                                      'NAN SCORE NDWI' : nanscores['NDWI'],
                                      'PREPROC TIME (sec)': elapsed_time})]
        
        return indices_preproc, QA_dict
        
    else:
        logging.info(f'VIIRS {product_date} : not considered due to nanscores > 95%')
        return 'not considered', QA_dict



def computeL7indices(image, landmask):
    """
    Computing L7 indices :
        - quality masking (clouds, opacity, saturation, water/land)
        - computing ndvi, ndwi
    """

    M = 2.75E-05 
    A = -0.2
    
    # LEVEL-1 CLOUD MASK : clouds=bit3, shadows=bit4, dilation=bit1
    qa_pxl = image.select('QA_PIXEL')
    
    cloudBitMask_L1 = 1 << 3
    shadowBitmask_L1 =  1 << 4
    dilationBitmask_L1 = 1 << 1
    
    cloudmask_L1 = ((qa_pxl.bitwiseAnd(cloudBitMask_L1).eq(0))
                    .And(qa_pxl.bitwiseAnd(shadowBitmask_L1).eq(0))
                    .And(qa_pxl.bitwiseAnd(dilationBitmask_L1).eq(0)))
    
    # LEVEL-2 CLOUD MASK : clouds=bit1, shadows=bit2, adjacent=bit3
    qa_sr = image.select('SR_CLOUD_QA')
    
    cloudBitMask_L2 = 1 << 1
    shadowBitmask_L2 =  1 << 2
    adjBitmask_L2 = 1 << 3
    
    cloudmask_L2 = (qa_sr.bitwiseAnd(cloudBitMask_L2)
                    .Or(qa_sr.bitwiseAnd(shadowBitmask_L2))
                    .Or(qa_sr.bitwiseAnd(adjBitmask_L2)))
    cloudmask_L2 = cloudmask_L2.unmask().Not()
    
    # COMBINE CLOUD MASKS :
    cloudmask = cloudmask_L1.And(cloudmask_L2)
    
    # LEVEL-1 WATER MASK : bit7
    waterBitMask_L1 = 1 << 7
    watermask_L1 = qa_pxl.bitwiseAnd(waterBitMask_L1).eq(0)
    
    # LEVEL-2 WATER MASK : bit5
    waterBitMask_L2 = 1 << 5
    watermask_L2 = qa_sr.bitwiseAnd(waterBitMask_L2)
    watermask_L2 = watermask_L2.unmask().Not()
    
    # COMBINE WATER MASKS :
    watermask = watermask_L1.And(watermask_L2)
    
    # OPACITY MASK : filter out pxls with opa>0.3
    qa_opa_dn = image.select('SR_ATMOS_OPACITY')
    qa_opa = qa_opa_dn.multiply(0.001)
    opamask = qa_opa.lte(0.3)
    
    # SATURATION MASK : band3=bit2, band4=bit3 , band5=bit4
    qa_sat = image.select('QA_RADSAT')
    satBitmask_b3 = 1 << 2
    satBitmask_b4 =  1 << 3
    satBitmask_b5 = 1 << 4
    satmask_b3 = (qa_sat.bitwiseAnd(satBitmask_b3).eq(0))
    satmask_b4 = (qa_sat.bitwiseAnd(satBitmask_b4).eq(0))
    satmask_b5 = (qa_sat.bitwiseAnd(satBitmask_b5).eq(0))

    # COMBINE ALL MASKS (+ LANDMASK) :
    QMASK_NDWI = (cloudmask.And(watermask).And(opamask).And(satmask_b4)
                  .And(satmask_b5).And(landmask)).rename('QMASK_NDWI')
    QMASK_NDVI = (cloudmask.And(watermask).And(opamask).And(satmask_b4)
                  .And(satmask_b3).And(landmask)).rename('QMASK_NDVI')
    QMASK_OUT = ee.Image([QMASK_NDWI, QMASK_NDVI])
    
    # VEGETATION INDICES : NDVI, NDWI
    opticalBands = image.select('SR_B.').multiply(M).add(A)
    image = image.addBands(opticalBands, None, True)
    ndwi = image.normalizedDifference(['SR_B4','SR_B5']).rename('NDWI')
    ndvi = image.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI')
    newBands = ee.Image([ndwi, ndvi])
    image_concat = image.addBands(newBands)
    INDICES_OUT = image_concat.select(['NDWI','NDVI'])
    
    # APPLY QMASK ON INDICES :
    INDICES_OUT = INDICES_OUT.updateMask(QMASK_OUT)
    
    # EXTRACT FILLMASK : SLC error mask (strip-lines) 
    Fmask = 1 << 0
    fmask = qa_pxl.bitwiseAnd(Fmask).eq(0).rename('FILLMASK')
    fmask = fmask.unmask(0)
    
    return INDICES_OUT, QMASK_OUT, fmask



def computeL8L9indices(image, landmask):
    """
    Computing L8 or L9 indices :
        - quality masking (clouds, opacity, saturation, water/land)
        - computing ndvi, ndwi
    """

    M = 2.75E-05 
    A = -0.2
    
    # CLOUD MASK : clouds=bit3, shadows=bit4, dilation=bit1, cirrus=bit2
    qa_pxl = image.select('QA_PIXEL')
    
    cloudBitMask = 1 << 3
    shadowBitmask =  1 << 4
    dilationBitmask = 1 << 1
    cirrusBitmask = 1 << 2
    
    cloudmask = ((qa_pxl.bitwiseAnd(cloudBitMask).eq(0))
                    .And(qa_pxl.bitwiseAnd(shadowBitmask).eq(0))
                    .And(qa_pxl.bitwiseAnd(dilationBitmask).eq(0))
                    .And(qa_pxl.bitwiseAnd(cirrusBitmask).eq(0)))
    
    # WATER MASK : bit7
    waterBitMask = 1 << 7
    watermask = qa_pxl.bitwiseAnd(waterBitMask).eq(0)
    
    # OPACITY MASK : high confidence level = bit6 & bit7
    qa_opa_dn = image.select('SR_QA_AEROSOL')

    opaBitmask6 = 1 << 6
    opaBitmask7 = 1 << 7
    
    opamask = (qa_opa_dn.bitwiseAnd(opaBitmask6)
               .And(qa_opa_dn.bitwiseAnd(opaBitmask7))).eq(0)

    # SATURATION MASK : band4=bit3, band5=bit4 , band6=bit5
    qa_sat = image.select('QA_RADSAT')
    satBitmask_b4 = 1 << 3
    satBitmask_b5 =  1 << 4
    satBitmask_b6 = 1 << 5
    satmask_b4 = (qa_sat.bitwiseAnd(satBitmask_b4).eq(0))
    satmask_b5= (qa_sat.bitwiseAnd(satBitmask_b5).eq(0))
    satmask_b6 = (qa_sat.bitwiseAnd(satBitmask_b6).eq(0))

    # TERRAIN OCCLUION MASK : bit11
    toBitmask = 1 << 11
    tomask = (qa_sat.bitwiseAnd(toBitmask).eq(0))

    # COMBINE ALL MASKS (+ LANDMASK) :
    QMASK_NDWI = (cloudmask.And(watermask).And(opamask).And(satmask_b5)
                  .And(satmask_b6).And(landmask).And(tomask)).rename('QMASK_NDWI')
    QMASK_NDVI = (cloudmask.And(watermask).And(opamask).And(satmask_b4)
                  .And(satmask_b5).And(landmask).And(tomask)).rename('QMASK_NDVI')
    QMASK = ee.Image([QMASK_NDWI, QMASK_NDVI])
    
    # VEGETATION INDICES : NDVI, NDWI
    opticalBands = image.select('SR_B.').multiply(M).add(A)
    image = image.addBands(opticalBands, None, True)
    ndwi = image.normalizedDifference(['SR_B5','SR_B6']).rename('NDWI')
    ndvi = image.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI')
    newBands = ee.Image([ndwi, ndvi])
    image_concat = image.addBands(newBands)
    INDICES_OUT = image_concat.select(['NDWI','NDVI'])
    
    # APPLY QMASK ON INDICES :
    INDICES_OUT = INDICES_OUT.updateMask(QMASK)
    
    return INDICES_OUT



def computeS2indices(image, landmask):
    """
    Computing S2 indices :
        - quality masking (clouds, opacity, saturation, water/land)
        - computing ndvi, ndwi
    """
    
    # Set preprocessing parameters:
    CLD_PRB_THRESH = 50
    NIR_DRK_THRESH = 0.15
    SR_BAND_SCALE = 1e4
    CLD_PRJ_DIST = 1
    BUFFER = 50
    OPA_THRESH = 0.3
    
    # 1) WATER Mask:
    watermask = image.select('SCL').neq(6)
    
    # 2.1) SEN2CLOUDLESS Cloud Probability Mask:
    cld_prb = ee.Image(image.get('s2cloudless')).select('probability')
    cloudmask = cld_prb.gt(CLD_PRB_THRESH).eq(0).clip(image.geometry())
    
    # 2.2) CLOUD SHADOW Mask:
    dark_pixels = (image.select('B8')
                   .lt(NIR_DRK_THRESH*SR_BAND_SCALE)
                   .multiply(watermask))
    shadow_azimuth = ee.Number(90).subtract(ee.Number(image.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
    cld_proj = (cloudmask.eq(0)
                .directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
                .reproject(**{'crs': image.select(0).projection(), 'scale': 100})
                .select('distance')
                .mask())
    shadowmask = cld_proj.multiply(dark_pixels).eq(0)
    cloudshdw_mask = cloudmask.And(shadowmask)
    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    cloudshdw_mask = (cloudshdw_mask.focalMax(2).focalMin(BUFFER*2/20)
                     .reproject(**{'crs': image.select([0]).projection(), 'scale': 20}))
    
    # 3) OPACITY Mask:
    qa_opa_dn = image.select('AOT')
    qa_opa = qa_opa_dn.divide(SR_BAND_SCALE)
    opamask = qa_opa.lte(OPA_THRESH)
    
    # 4) SATURATED/DEFECTIVE PXLS Mask:
    satmask = image.select('SCL').neq(1)
    
    # 5) COMBINE ALL MASKS :
    QMASK = (watermask.And(cloudshdw_mask).And(opamask)
             .And(satmask).And(landmask)).rename('QMASK')
    
    # 6) VEGETATION INDICES : NDWI, NDVI
    ndwi = image.normalizedDifference(['B8','B11']).rename('NDWI')
    ndvi = image.normalizedDifference(['B8','B4']).rename('NDVI')
    newBands = ee.Image([ndwi, ndvi])
    image_concat = image.addBands(newBands)
    INDICES_OUT = image_concat.select(['NDWI','NDVI'])
    
    # 7) APPLY QMASK ON INDICES :
    INDICES_OUT = INDICES_OUT.updateMask(QMASK)
    
    return INDICES_OUT



def gapfillingL7(indices, qmask, landmask):
    '''
    Gap-filling L7 indices :
        - fill gaps due to SLC error strip-lines, on 2003-2013
        - morphological filters applied on indices (and masks)
    '''

    RAD_INDICES = 1
    IT_INDICES = 25
    RAD_QMASK = 2
    IT_QMASK = 11

    # Gap-filling Indices with MEAN Morphological filter :
    indices_fill = indices.focal_mean(RAD_INDICES, 'square', 'pixels', IT_INDICES)
    indices_gapfill = indices.addBands(indices_fill.blend(indices),overwrite=True)
    
    # Gap-filling Qmask with MODE (dominant value) Morphological filter :
    qmask_fill = qmask.focal_mode(RAD_QMASK, 'square', 'pixels', IT_QMASK)
    qmask_gapfill = qmask_fill.blend(qmask)
    
    # Re-Apply LANDMASK to delete edge effects due to gap-filling :
    indices_gapfill = indices_gapfill.updateMask(landmask)
    qmask_gapfill = qmask_gapfill.updateMask(landmask)
    
    # Re-applying gap-filled QUALITY MASK to keep gap-filled values only on stripe lines :
    indices_gapfill = indices_gapfill.updateMask(qmask_gapfill)

    return indices_gapfill



def preprocessingL7(l7, landmask, landsat_grid, tile, QA_dict, path2key=os.getcwd(), new_collection=None, asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single L7 product :
        1) quality masking (clouds, opacity, saturation), computing indices
        2.1) gap-filling indices : method = focal_mean(rad=1, it=25)
        2.2) gap-filling quality masks : method = focal_mode(rad=2, it=11)
        3) post-processing : applying qmask (gap-filled)), removing last outliers, clipping
        4) (IF asset_export==1 and nanscore<=95%) exporting preprocessed image to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)
    """

    product_id = geegen.googleErrorsControl(l7.get('LANDSAT_PRODUCT_ID'), path2key)
    product_date = geegen.googleErrorsControl(l7.get('DATE_ACQUIRED'), path2key)
    product_proj = geegen.googleErrorsControl(l7.select('SR_B3').projection(), path2key)
    slc_sensor = geegen.googleErrorsControl(l7.get('SENSOR_MODE_SLC'), path2key)
    cloud_land = geegen.googleErrorsControl(l7.get('CLOUD_COVER_LAND'), path2key)
    im_quality = geegen.googleErrorsControl(l7.get('IMAGE_QUALITY'), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE (compute nanscores) and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'L7 {product_date} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])
                nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask)
                QA_dict += [ee.Feature(None, {'TILE':tile,
                                            'FILE NAME': product_id,
                                            'DATE': product_date,
                                            'CLOUD LAND': cloud_land,
                                            'IMAGE QUALITY' : im_quality,
                                            'SLC MODE': slc_sensor,
                                            'NAN SCORE NDWI' : nanscores['NDWI'],
                                            'NAN SCORE NDVI' : nanscores['NDVI'],
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict 
    
    if cloud_land<=90:
        
        # --- COMPUTE MASKED INDICES ---
        indices, qmask, fmask = computeL7indices(l7, landmask)
        
        # --- GAP-FILLING (IF SLC='OFF'=> only 2003-2013) ---
        if slc_sensor=='OFF':
            indices_gapfill = gapfillingL7(indices, qmask, landmask)
        else:
            indices_gapfill = indices

        # --- REMOUVING OUTLIERS and CLIPPING TO LANDSAT GRID ---
        indices_preproc = geegen.removeOutliers(indices_gapfill)
        indices_preproc = indices_preproc.clip(landsat_grid).float()
        nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask)

        if (nanscores['NDWI']<=0.95 and nanscores['NDVI']<=0.95):
            
            logging.info(f'L7 {product_date} : preprocessing product')

            # --- EXPORT TO ASSET ---
            if asset_export==1:
                MAX_GEETRY = 6
                gee_test=''
                n_geetry = 1

                while gee_test=='' and n_geetry<=MAX_GEETRY:
                    preproc = indices_preproc.addBands(fmask, overwrite=True)
                    start_time = time.time()
                    geegen.exportImageToAsset(preproc, preproc_filename, path2key, data_crs=product_proj['crs'], data_transform=product_proj['transform'])
                    elapsed_time = round(time.time() - start_time)

                    # Read already exported products to improve computation time
                    indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])

                    # Test to ensure good asset export (if not, retry exporting)
                    gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                    if gee_test=='':
                        n_geetry+=1
                        time.sleep(5**n_geetry)

            else: elapsed_time = ''

            # --- CONCATENATE SCORES INTO TABLE ---
            QA_dict += [ee.Feature(None, {'TILE':tile,
                                          'FILE NAME': product_id,
                                          'DATE': product_date,
                                          'CLOUD LAND': cloud_land,
                                          'IMAGE QUALITY' : im_quality,
                                          'SLC MODE': slc_sensor,
                                          'NAN SCORE NDWI' : nanscores['NDWI'],
                                          'NAN SCORE NDVI' : nanscores['NDVI'],
                                          'PREPROC TIME (sec)': elapsed_time})]
    
            return indices_preproc, QA_dict
        
        else:
            logging.info(f'L7 {product_date} : not considered due to nanscores > 95%')
            return 'not considered', QA_dict
    else:
        logging.info(f'L7 {product_date} : not considered due to {round(cloud_land,2)}% cloud cover')
        return 'not considered', QA_dict



def preprocessingL8L9(landsat, landmask, landsat_grid, tile, QA_dict, path2key=os.getcwd(), new_collection=None, asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single L8 or L9 product :
        1) computing indices
        2) quality masking (clouds, opacity, saturation)
        3) post-processing : removing last outliers and clipping
        4) (IF asset_export==1 and nanscore<=95%) exporting preprocessed image to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)
    """
    
    product_id = geegen.googleErrorsControl(landsat.get('LANDSAT_PRODUCT_ID'), path2key)
    product_date = geegen.googleErrorsControl(landsat.get('DATE_ACQUIRED'), path2key)
    product_proj = geegen.googleErrorsControl(landsat.select('SR_B4').projection(), path2key)
    cloud_land = geegen.googleErrorsControl(landsat.get('CLOUD_COVER_LAND'), path2key)
    im_quality = geegen.googleErrorsControl(landsat.get('IMAGE_QUALITY_OLI'), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE (compute nanscores) and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'L8/9 {product_date} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])
                nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask)
                QA_dict += [ee.Feature(None, {'TILE':tile,
                                            'FILE NAME': product_id,
                                            'DATE': product_date,
                                            'CLOUD LAND': cloud_land,
                                            'IMAGE QUALITY' : im_quality,
                                            'SLC MODE': '',
                                            'NAN SCORE NDWI' : nanscores['NDWI'],
                                            'NAN SCORE NDVI' : nanscores['NDVI'],
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict
    
    if cloud_land<=90:
        
        # --- COMPUTE MASKED INDICES ---       
        indices_preproc = computeL8L9indices(landsat, landmask)
         
        # --- REMOUVING OUTLIERS and CLIPPING TO LANDSAT GRID ---
        indices_preproc = geegen.removeOutliers(indices_preproc)
        indices_preproc = indices_preproc.clip(landsat_grid).float()
        nanscores = geegen.computeNanScores(indices_preproc, path2key, landmask)

        if (nanscores['NDWI']<=0.95 and nanscores['NDVI']<=0.95):
            
            logging.info(f'L8/9 {product_date} : preprocessing product')

            # --- EXPORT TO ASSET ---
            if asset_export==1:
                MAX_GEETRY = 6
                gee_test=''
                n_geetry = 1

                while gee_test=='' and n_geetry<=MAX_GEETRY:
                    start_time = time.time()
                    geegen.exportImageToAsset(indices_preproc, preproc_filename, path2key, data_crs=product_proj['crs'], data_transform=product_proj['transform'])
                    elapsed_time = round(time.time() - start_time)

                    # Read already exported products to improve computation time
                    indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])

                    # Test to ensure good asset export (if not, retry exporting)
                    gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                    if gee_test=='':
                        n_geetry+=1
                        time.sleep(5**n_geetry)

            else: elapsed_time = ''
            
            # --- CONCATENATE SCORES INTO TABLE ---
            QA_dict += [ee.Feature(None, {'TILE':tile,
                                          'FILE NAME': product_id,
                                          'DATE': product_date,
                                          'CLOUD LAND': cloud_land,
                                          'IMAGE QUALITY' : im_quality,
                                          'SLC MODE': '',
                                          'NAN SCORE NDWI' : nanscores['NDWI'],
                                          'NAN SCORE NDVI' : nanscores['NDVI'],
                                          'PREPROC TIME (sec)': elapsed_time})]

            return indices_preproc, QA_dict
        
        else:
            logging.info(f'L8/9 {product_date} : not considered due to nanscores > 95%')
            return 'not considered', QA_dict
    else:
        logging.info(f'L8/9 {product_date} : not considered due to {round(cloud_land,2)}% cloud cover')
        return 'not considered', QA_dict



def preprocessingS2(s2, landmask, tile_landsat, QA_dict, path2key=os.getcwd(), new_collection=None, asset_export=None):
    """
    Global function that calls sub-functions for preprocessing a single S2 product :
        1) computing indices
        2) quality masking (clouds, opacity, saturation)
        3) post-processing : removing last outliers
        4) (IF asset_export==1) exporting preprocessed image to gee asset
    
    Returns gee preproc indices and updated quality dictionnary (clouds, nan scores, etc.)
    Note : Condition "nanscore<=95%" is not used here to improve computation speed,
           Cloud condition is still used
    """
    
    product_id = geegen.googleErrorsControl(s2.get('PRODUCT_ID'), path2key)
    product_date = geegen.googleErrorsControl(s2.date().format('YYY-MM-dd'), path2key)
    product_proj = geegen.googleErrorsControl(s2.select('B4').projection(), path2key)
    cloud = geegen.googleErrorsControl(s2.get('CLOUDY_PIXEL_PERCENTAGE'), path2key)
    tile_s2 = geegen.googleErrorsControl(s2.get('MGRS_TILE'), path2key)

    # --- IF PRODUCT WAS ALREADY PREPROCESSED -> READ GEE, FILL TABLE and EXIT ---
    if asset_export==1:
        preproc_filename = new_collection+'/'+product_id
        prod_list = (ee.data.listAssets({'parent':new_collection}))['assets']
        for asset in prod_list:
            if preproc_filename in asset["name"]:
                logging.info(f'S2 {product_date} {tile_s2} : already preprocessed product')
                indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])
                QA_dict += [ee.Feature(None,{'TILE': tile_landsat,
                                            'TILE S2': tile_s2,
                                            'FILE NAME': product_id,
                                            'DATE': product_date,
                                            'CLOUD': cloud,
                                            'PREPROC TIME (sec)': ''})]
                return indices_preproc, QA_dict
    
    if cloud<=90:
        
        # --- COMPUTE MASKED INDICES ---       
        indices_preproc = computeS2indices(s2, landmask)
         
        # --- REMOUVING OUTLIERS ---
        indices_preproc = geegen.removeOutliers(indices_preproc).float()
        # nanscores = geegen.computeNanScores(s2indices_preproc, path2key, landmask)                

        # if (nanscores['NDWI']<=0.95 and nanscores['NDVI']<=0.95):
                
        logging.info(f'S2 {product_date} {tile_s2} : preprocessing product')

        # --- EXPORT TO ASSET ---
        if asset_export==1:
            MAX_GEETRY = 6
            gee_test=''
            n_geetry = 1

            while gee_test=='' and n_geetry<=MAX_GEETRY:
                start_time = time.time()
                geegen.exportImageToAsset(indices_preproc, preproc_filename, path2key, data_crs=product_proj['crs'], data_transform=product_proj['transform'])
                elapsed_time = round(time.time() - start_time)

                # Read already exported products to improve computation time
                indices_preproc = ee.Image(preproc_filename).select(['NDWI', 'NDVI'])

                # Test to ensure good asset export (if not, retry exporting)
                gee_test = geegen.googleErrorsControl(indices_preproc, path2key)
                if gee_test=='':
                    n_geetry+=1
                    time.sleep(5**n_geetry)

        else: elapsed_time = ''

        # --- CONCATENATE SCORES INTO TABLE ---        
        QA_dict += [ee.Feature(None,{'TILE': tile_landsat,
                                     'TILE S2': tile_s2,
                                     'FILE NAME': product_id,
                                     'DATE': product_date,
                                     'CLOUD': cloud,
                                    #  'NAN SCORE NDWI' : nanscore_ndwi,
                                    #  'NAN SCORE NDVI' : nanscore_ndvi,
                                     'PREPROC TIME (sec)': elapsed_time})]
    
        return indices_preproc, QA_dict
        
    else:
        logging.info(f'S2 {product_date} {tile_s2} : not considered due to {round(cloud,2)}% cloud cover')
        return 'not considered', QA_dict

