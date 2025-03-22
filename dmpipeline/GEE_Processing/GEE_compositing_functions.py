# -*- coding: utf-8 -*-
"""
##############################################################################

GEE Functions to Composite satellite products

##############################################################################
"""


import os
import ee
import dmpipeline.GEE_Processing.GEE_generic_functions as geegen
import dmpipeline.GEE_Processing.GEE_preprocessing_functions as geeprep

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %H:%M:%S', level=logging.INFO)



def extractComposite(collection, path2key, data_geom=None, data_scale=None, bands_list=None):
    """
    Extract indices composite per band :
        - average indices -> INDICES_mean
        - count number of valid pixels -> INDICES_count
    Add geometry with clip() (since compositing removes geometry information)
    """

    if data_geom is None: data_geom = collection.first().select(1).geometry()

    if data_scale is None: data_scale = collection.first().select(1).projection().nominalScale()

    if bands_list is None: bands_list = geegen.googleErrorsControl(collection.first().bandNames(), path2key)
    else: collection = collection.select(bands_list)

    bands_list_count = bands_list.copy()
    for i in range(len(bands_list)):
        bands_list_count[i] = geegen.googleErrorsControl(ee.String(bands_list[i]).cat('_count'), path2key)

    # INDICES : Average values and count number of good pixels
    INDICES_mean = collection.mean().float().clip(data_geom)
    INDICES_count = collection.count().unmask(0).float().clip(data_geom)
    INDICES_count = INDICES_count.rename(bands_list_count)
    
    # CONCATENATE INTO TWO COLLECTIONS (one per band) :
    INDICES_composite = {}
    for i in range(len(bands_list)):
        indice_composite = ee.Image([INDICES_mean.select(bands_list[i]), INDICES_count.select(bands_list_count[i])])
        INDICES_composite[bands_list[i]] = indice_composite
        del indice_composite
    
    return INDICES_composite



def processComposite_MODIS(LSTaqua_month, LSTterra_month, REFLECTterra_month, landmask, new_collection, QA_dict, crs_out, grid_out, scale_out, go_month, path2key=os.getcwd(), asset_export=None):
    """
    Decade and Month compositing of MODIS LST and NDWI products.
    The number of good data used in the period is also given per pixel.
    """

    COMP_NANSCORES = []
    COMP_TYPES = []
    NBLST_DATES = []
    NBREFLECT_DATES = []
    start_D = [1, 11, 21]
    end_D = [10, 20, 31]
    LSTcompD_col = ee.ImageCollection([])
    NDWIcompD_col = ee.ImageCollection([])
    LSTcompM_col = ee.ImageCollection([])
    NDWIcompM_col = ee.ImageCollection([])
    indicesM_preproc = []


    # ====== DECADE COMPOSITING ======

    for d in range(3):
        logging.info(f'DECADE : {d+1}')

        indicesD_preproc = []
        LSTaqua_D = LSTaqua_month.filter(ee.Filter.calendarRange(start_D[d],end_D[d],'DAY_OF_MONTH'))
        LSTterra_D = LSTterra_month.filter(ee.Filter.calendarRange(start_D[d],end_D[d],'DAY_OF_MONTH'))
        REFLECTterra_D = REFLECTterra_month.filter(ee.Filter.calendarRange(start_D[d],end_D[d],'DAY_OF_MONTH'))

        NbLSTaqua_D = geegen.googleErrorsControl(LSTaqua_D.size(), path2key)
        NbLSTterra_D = geegen.googleErrorsControl(LSTterra_D.size(), path2key)
        NbLST_D = NbLSTaqua_D + NbLSTterra_D
        NbLREFLECT_D = geegen.googleErrorsControl(REFLECTterra_D.size(), path2key)

        # --- MISSING PRODUCTS ---
        if NbLST_D==0 or NbLREFLECT_D==0:
            continue

        # --- LOOP OVER DAYS (in decade D) ---
        ALL_dataset_D = LSTaqua_D.merge(LSTterra_D).merge(REFLECTterra_D)
        date_start_D = ALL_dataset_D.limit(1, 'system:time_start', True).first().date()
        date_end_D = ALL_dataset_D.limit(1, 'system:time_start', False).first().date()
        days_D = geegen.googleErrorsControl(ee.List.sequence(date_start_D.get('day'), date_end_D.get('day')), path2key)

        for day in days_D:
            lst_aqua = LSTaqua_D.filter(ee.Filter.calendarRange(day,day,'DAY_OF_MONTH')).first()
            lst_terra = LSTterra_D.filter(ee.Filter.calendarRange(day,day,'DAY_OF_MONTH')).first()
            reflect = REFLECTterra_D.filter(ee.Filter.calendarRange(day,day,'DAY_OF_MONTH')).first()

            indices_preproc, QA_dict  = geeprep.preprocessingMODIS(lst_aqua, lst_terra, reflect, landmask, new_collection, QA_dict, crs_out, grid_out, scale_out, path2key, asset_export)
            if indices_preproc=='not considered': pass
            else:
                indicesD_preproc += [indices_preproc]
                if go_month==1: indicesM_preproc += [indices_preproc]
            del lst_aqua, lst_terra, reflect

        indicesD_preproc = ee.ImageCollection(indicesD_preproc)
        Nbgoodimages_D = geegen.googleErrorsControl(indicesD_preproc.size(), path2key)
        if Nbgoodimages_D==0:
            continue
        elif Nbgoodimages_D==1:
            indicesD_comp = {}
            indicesD_comp['LST'] =  indicesD_preproc.first().select('LST')
            indicesD_comp['NDWI'] =  indicesD_preproc.first().select('NDWI')
            nanscores_D = geegen.computeNanScores(indicesD_preproc.first(), path2key, landmask, scale=scale_out)
            comp_type_D = f'COMPD{d+1}'
        elif Nbgoodimages_D > 1:
            indicesD_comp = extractComposite(indicesD_preproc, path2key, grid_out, scale_out)
            nanscores_D = geegen.computeNanScores(indicesD_comp['LST'].addBands(indicesD_comp['NDWI'], overwrite=True).select(['LST', 'NDWI']), path2key, landmask, scale=scale_out)
            comp_type_D = f'COMPD{d+1}'
        del indicesD_preproc

        LSTcompD_col = LSTcompD_col.merge(indicesD_comp['LST'])
        NDWIcompD_col = NDWIcompD_col.merge(indicesD_comp['NDWI'])
        COMP_NANSCORES += [nanscores_D]
        COMP_TYPES += [comp_type_D]
        NBLST_DATES += [NbLST_D]
        NBREFLECT_DATES += [NbLREFLECT_D]

        del LSTaqua_D, LSTterra_D, REFLECTterra_D, indicesD_comp
    

    # ====== MONTH COMPOSITING (IF go_month=1) ======
        
    if go_month==1:

        indicesM_preproc = ee.ImageCollection(indicesM_preproc)
        Nbgoodimages_M = geegen.googleErrorsControl(indicesM_preproc.size(), path2key)
        if Nbgoodimages_M==0:
            pass
        elif Nbgoodimages_M==1:
            LSTcompM_col =  indicesM_preproc.first().select('LST')
            NDWIcompM_col =  indicesM_preproc.first().select('NDWI')
        elif Nbgoodimages_M > 1:
            indicesM_comp = extractComposite(indicesM_preproc, path2key, grid_out, scale_out)
            LSTcompM_col = indicesM_comp['LST']
            NDWIcompM_col = indicesM_comp['NDWI']
    

    return LSTcompD_col, NDWIcompD_col, LSTcompM_col, NDWIcompM_col, QA_dict, COMP_TYPES, COMP_NANSCORES, NBLST_DATES, NBREFLECT_DATES



def processCompositeDecade_LANDSAT_S2(L7collection, L8collection, L9collection, S2collection, y, m, date_start_collection, date_end_collection, landmask, landsat_grid, tile_landsat,
                                      scale_out, new_landsatcollection, new_s2collection, QA_dict_L, QA_dict_S2, path2key=os.getcwd(), asset_export_l=None, asset_export_s2=None):
    """
    Decade compositing of NDVI and NDWI products.
    The number of good data used in the period is also given per pixel.

    Rmq : ajout du pas de temps intermédiaire, et centre la recherche sur la décade étudiée (-5j +5j, puis -10j +10j)
          implique 4 autres types de composite (DAY1e -> un produit trouvé dans les 20j autour,
                                                DAY1m -> un produit trouvé dans les 30j autour,
                                                DAY1M -> un produit trouvé sur les jours du mois en cours (dernier recours),
                                                COMPD1e -> un composite représentant les 20j autour,
                                                COMPD1m -> un composite représentant dans les 30j autour
                                                           un peu différent du composite v0 qui se calculait sur les jours du mois en cours
                                                COMPD1M -> composite v0 calculé sur les jours du mois en cours (dernier recours)
                                                           si aucun produit n'est trouvé via étapes précédentes alors que des produits sont tout de même dans autre décade)
    """

    COMP_NAN_THRESH = 0.5
    COMP_NANSCORES = []
    COMP_TYPES = []
    NBL7_DATES = []
    NBL8_DATES = []
    NBL9_DATES = []
    NBS2_DATES = []
    QA_dict_L_D = QA_dict_L
    QA_dict_S2_D = QA_dict_S2
    start_D = [1, 11, 21]
    end_D = [10, 20, 31]

    date_start_collection = date_start_collection.update(hour=0, minute=0, second=0)
    date_end_collection = date_end_collection.update(hour=0, minute=0, second=0)
    y_start_collection = geegen.googleErrorsControl(date_start_collection.get('year'), path2key)
    y_end_collection = geegen.googleErrorsControl(date_end_collection.get('year'), path2key)
    m_start_collection = geegen.googleErrorsControl(date_start_collection.get('month'), path2key)
    m_end_collection = geegen.googleErrorsControl(date_end_collection.get('month'), path2key)
    
    NDWIcomp_col = ee.ImageCollection([])
    NDVIcomp_col = ee.ImageCollection([])

    for d in range(3):

        date_start_D = ee.Date('{}-{:02d}-{}'.format(y,m,start_D[d]))
        date_end_D = ee.Date('{}-{:02d}-{}'.format(y,m,end_D[d]))
        if d==2:
            day_range = date_start_D.getRange('month')
            date_end_D = day_range.end().advance(-1, 'day')
            del day_range

        # Do not process current decade if not inside valid collection period
        if y==y_start_collection and m==m_start_collection:
            date_diff = date_end_D.difference(date_start_collection, 'days')
            if geegen.googleErrorsControl(date_diff, path2key)<0:
                print(f'DECADE {d+1} IS BEFORE COLLECTION PERIOD')
                continue
            del date_diff
        if y==y_end_collection and m==m_end_collection:
            date_diff = date_start_D.difference(date_end_collection, 'day')
            if geegen.googleErrorsControl(date_diff, path2key)>0:
                print(f'DECADE {d+1} IS AFTER COLLECTION PERIOD')
                continue
            del date_diff

        logging.info(f'DECADE : {d+1}')
        comp_type_D = []
        nanscores_D = {'NDVI':1, 'NDWI':1}
                 
        L7col_D = L7collection.filterDate(date_start_D, date_end_D.advance(1, 'day'))   # +1 since date_end is exclusive here
        L8col_D = L8collection.filterDate(date_start_D, date_end_D.advance(1, 'day'))   # +1 since date_end is exclusive here
        L9col_D = L9collection.filterDate(date_start_D, date_end_D.advance(1, 'day'))   # +1 since date_end is exclusive here
        S2col_D = S2collection.filterDate(date_start_D, date_end_D.advance(1, 'day'))   # +1 since date_end is exclusive here
        S2col_D = S2col_D.map(lambda image: image.set('simpleDateMillis', ee.Date(image.date().format('YYYY-MM-dd')).millis()))
        
        NbL7images_D = geegen.googleErrorsControl(L7col_D.size(), path2key)
        NbL8images_D = geegen.googleErrorsControl(L8col_D.size(), path2key)
        NbL9images_D = geegen.googleErrorsControl(L9col_D.size(), path2key)
        NbS2images_D = geegen.googleErrorsControl(S2col_D.size(), path2key)
        Nbimages_D = NbL8images_D + NbL9images_D + NbS2images_D

        col_D = L7col_D.merge(L8col_D).merge(L9col_D).merge(S2col_D)
        
        # --- NO PRODUCTS ---
        if Nbimages_D==0: pass

        # --- PROCESS DECADE D ---
        elif Nbimages_D >= 1:
            indicesD_preproc = []

            if NbL7images_D != 0:
                L7_list = L7col_D.toList(NbL7images_D)
                for i in range(NbL7images_D):
                    l7 = ee.Image(L7_list.get(i))
                    if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                    else: new_landsatcollection_tmp = ''
                    indices_preproc, QA_dict_L_D = geeprep.preprocessingL7(l7, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                    if indices_preproc=='not considered': pass
                    else:
                        indicesD_preproc += [indices_preproc]
                    del l7, indices_preproc
                del L7_list
            
            if NbL8images_D != 0:
                L8_list = L8col_D.toList(NbL8images_D)
                for i in range(NbL8images_D):
                    l8 = ee.Image(L8_list.get(i))
                    if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                    else: new_landsatcollection_tmp = ''
                    indices_preproc, QA_dict_L_D = geeprep.preprocessingL8L9(l8, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                    if indices_preproc=='not considered': pass
                    else:
                        indicesD_preproc += [indices_preproc]
                    del l8, indices_preproc
                del L8_list
            
            if NbL9images_D != 0:
                L9_list = L9col_D.toList(NbL9images_D)
                for i in range(NbL9images_D):
                    l9 = ee.Image(L9_list.get(i))
                    if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                    else: new_landsatcollection_tmp = ''
                    indices_preproc, QA_dict_L_D = geeprep.preprocessingL8L9(l9, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                    if indices_preproc=='not considered': pass
                    else:
                        indicesD_preproc += [geegen.calibrateData(indices_preproc, 'L9')]
                    del l9, indices_preproc
                del L9_list
            
            if NbS2images_D != 0:
                S2_list = S2col_D.toList(NbS2images_D)
                for i in range(NbS2images_D):
                    s2 = ee.Image(S2_list.get(i))
                    if asset_export_s2==1: new_s2collection_tile = new_s2collection[geegen.googleErrorsControl(s2.get('MGRS_TILE'), path2key)]
                    else: new_s2collection_tile = ''
                    indices_preproc, QA_dict_S2_D = geeprep.preprocessingS2(s2, landmask, tile_landsat, QA_dict_S2_D, path2key, new_s2collection_tile, asset_export_s2)
                    
                    if indices_preproc=='not considered': pass
                    else:
                        indicesD_preproc += [geegen.calibrateData(indices_preproc, 'S2')]
                    del s2, indices_preproc
                del S2_list
            
            indicesD_preproc = ee.ImageCollection(indicesD_preproc)
            Nbgoodimages_D = geegen.googleErrorsControl(indicesD_preproc.size(), path2key)
            
            if Nbgoodimages_D==0:
                pass
            elif Nbgoodimages_D==1:
                indicesD_comp={}
                indicesD_comp['NDWI'] =  indicesD_preproc.first().select('NDWI')
                indicesD_comp['NDVI'] =  indicesD_preproc.first().select('NDVI')
                nanscores_D = geegen.computeNanScores(indicesD_preproc.first(), path2key, landmask, data_geom=landsat_grid, scale=scale_out)
                comp_type_D = f'COMPD{d+1}'
            elif Nbgoodimages_D > 1:
                indicesD_comp = extractComposite(indicesD_preproc, path2key, landsat_grid)
                nanscores_D = geegen.computeNanScores(indicesD_comp['NDWI'].addBands(indicesD_comp['NDVI'], overwrite=True).select(['NDWI', 'NDVI']), path2key, landmask, data_geom=landsat_grid, scale=scale_out)
                comp_type_D = f'COMPD{d+1}'
            del indicesD_preproc

        
        # --- IF NANSCORE > COMP_NAN_THRESH => EXTEND DECAD PERIOD [-5j +5j] then [-10j +10j]  ---
        # (Here indices are treated together)
        cpt_ext = 1
        date_start_D_EXT = date_start_D
        date_end_D_EXT = date_end_D.advance(1, 'day') # +1 since date_end is exclusive here

        while ((nanscores_D['NDWI']>COMP_NAN_THRESH) or (nanscores_D['NDVI']>COMP_NAN_THRESH)) and (cpt_ext<=3):
            logging.info(f'EXTEND DECAD PERIOD ({cpt_ext}/3)')
            QA_dict_L_D = QA_dict_L
            QA_dict_S2_D = QA_dict_S2
            
            if cpt_ext!=3:
                # Cases for (-5d +5d), then (-10d +10d)
                date_start_D_EXT = date_start_D_EXT.advance(-5, 'day')
                date_end_D_EXT = date_end_D_EXT.advance(5, 'day')

                L7col_D_EXT = L7collection.filterDate(date_start_D_EXT, date_end_D_EXT)
                L8col_D_EXT = L8collection.filterDate(date_start_D_EXT, date_end_D_EXT)
                L9col_D_EXT = L9collection.filterDate(date_start_D_EXT, date_end_D_EXT)
                S2col_D_EXT = S2collection.filterDate(date_start_D_EXT, date_end_D_EXT)
            
            else:
                # Case for all current month
                L7col_D_EXT = L7collection
                L8col_D_EXT = L8collection
                L9col_D_EXT = L9collection
                S2col_D_EXT = S2collection
            
            NbL7images_D_EXT = geegen.googleErrorsControl(L7col_D_EXT.size(), path2key)
            NbL8images_D_EXT = geegen.googleErrorsControl(L8col_D_EXT.size(), path2key)
            NbL9images_D_EXT = geegen.googleErrorsControl(L9col_D_EXT.size(), path2key)
            NbS2images_D_EXT = geegen.googleErrorsControl(S2col_D_EXT.size(), path2key)
            Nbimages_D_EXT = NbL8images_D_EXT + NbL9images_D_EXT + NbS2images_D_EXT
            
            col_D_EXT = L7col_D_EXT.merge(L8col_D_EXT).merge(L9col_D_EXT).merge(S2col_D_EXT)

            if Nbimages_D_EXT==Nbimages_D: pass
            
            elif Nbimages_D_EXT >= 1:
                indicesD_preproc = []

                if NbL7images_D_EXT != 0:
                    L7_list = L7col_D_EXT.toList(NbL7images_D_EXT)
                    for i in range(NbL7images_D_EXT):
                        l7 = ee.Image(L7_list.get(i))
                        if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                        else: new_landsatcollection_tmp = ''
                        indices_preproc, QA_dict_L_D = geeprep.preprocessingL7(l7, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                        if indices_preproc=='not considered': pass
                        else:
                            indicesD_preproc += [indices_preproc]
                        del l7, indices_preproc
                    del L7_list
                
                if NbL8images_D_EXT != 0:
                    L8_list = L8col_D_EXT.toList(NbL8images_D_EXT)
                    for i in range(NbL8images_D_EXT):
                        l8 = ee.Image(L8_list.get(i))
                        if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                        else: new_landsatcollection_tmp = ''
                        indices_preproc, QA_dict_L_D = geeprep.preprocessingL8L9(l8, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                        if indices_preproc=='not considered': pass
                        else:
                            indicesD_preproc += [indices_preproc]
                        del l8, indices_preproc
                    del L8_list
                
                if NbL9images_D_EXT != 0:
                    L9_list = L9col_D_EXT.toList(NbL9images_D_EXT)
                    for i in range(NbL9images_D_EXT):
                        l9 = ee.Image(L9_list.get(i))
                        if asset_export_l==1: new_landsatcollection_tmp = new_landsatcollection
                        else: new_landsatcollection_tmp = ''
                        indices_preproc, QA_dict_L_D = geeprep.preprocessingL8L9(l9, landmask, landsat_grid, tile_landsat, QA_dict_L_D, path2key, new_landsatcollection_tmp, asset_export_l)
                    
                        if indices_preproc=='not considered': pass
                        else:
                            indicesD_preproc += [geegen.calibrateData(indices_preproc, 'L9')]
                        del l9, indices_preproc
                    del L9_list
                
                if NbS2images_D_EXT != 0:
                    S2_list = S2col_D_EXT.toList(NbS2images_D_EXT)
                    for i in range(NbS2images_D_EXT):
                        s2 = ee.Image(S2_list.get(i))
                        if asset_export_s2==1: new_s2collection_tile = new_s2collection[geegen.googleErrorsControl(s2.get('MGRS_TILE'), path2key)]
                        else: new_s2collection_tile = ''
                        indices_preproc, QA_dict_S2_D = geeprep.preprocessingS2(s2, landmask, tile_landsat, QA_dict_S2_D, path2key, new_s2collection_tile, asset_export_s2)

                        if indices_preproc=='not considered': pass
                        else:
                            indicesD_preproc += [geegen.calibrateData(indices_preproc, 'S2')]
                        del s2, indices_preproc
                    del S2_list

                indicesD_preproc = ee.ImageCollection(indicesD_preproc)
                Nbgoodimages_D = geegen.googleErrorsControl(indicesD_preproc.size(), path2key)

                if Nbgoodimages_D==0:
                    pass
                elif Nbgoodimages_D==1:
                    indicesD_comp={}
                    indicesD_comp['NDWI'] = indicesD_preproc.first().select('NDWI')
                    indicesD_comp['NDVI'] = indicesD_preproc.first().select('NDVI')
                    nanscores_D = geegen.computeNanScores(indicesD_preproc.first(), path2key, landmask, data_geom=landsat_grid, scale=scale_out)
                    if cpt_ext==1: comp_type_D = f'DAY{d+1}e'
                    elif cpt_ext==2: comp_type_D = f'DAY{d+1}m'
                    elif cpt_ext==3: comp_type_D = f'DAY{d+1}M'
                elif Nbgoodimages_D>1:
                    indicesD_comp = extractComposite(indicesD_preproc, path2key, landsat_grid)
                    nanscores_D = geegen.computeNanScores(indicesD_comp['NDWI'].addBands(indicesD_comp['NDVI'], overwrite=True).select(['NDWI', 'NDVI']), path2key, landmask, data_geom=landsat_grid, scale=scale_out)
                    if cpt_ext==1: comp_type_D = f'COMPD{d+1}e'
                    elif cpt_ext==2: comp_type_D = f'COMPD{d+1}m'
                    elif cpt_ext==3: comp_type_D = f'COMPD{d+1}M'
                del indicesD_preproc

            col_D = col_D_EXT
            Nbimages_D = Nbimages_D_EXT
            S2col_D = S2col_D_EXT
            NbL7images_D = NbL7images_D_EXT
            NbL8images_D = NbL8images_D_EXT
            NbL9images_D = NbL9images_D_EXT
            del col_D_EXT, L7col_D_EXT, L8col_D_EXT, L9col_D_EXT, S2col_D_EXT, Nbimages_D_EXT, NbL7images_D_EXT, NbL8images_D_EXT, NbL9images_D_EXT, NbS2images_D_EXT
            cpt_ext += 1
        
        NDWIcomp_col = NDWIcomp_col.merge(indicesD_comp['NDWI'])
        NDVIcomp_col = NDVIcomp_col.merge(indicesD_comp['NDVI'])
        QA_dict_L = QA_dict_L_D
        QA_dict_S2 = QA_dict_S2_D
        COMP_NANSCORES += [nanscores_D]
        COMP_TYPES += [comp_type_D]
        NBL7_DATES += [NbL7images_D]
        NBL8_DATES += [NbL8images_D]
        NBL9_DATES += [NbL9images_D]
        NbS2_disctinct_dates_D = geegen.googleErrorsControl(S2col_D.aggregate_count_distinct('simpleDateMillis'), path2key)
        NBS2_DATES += [NbS2_disctinct_dates_D]

        del L7col_D, L8col_D, L9col_D, S2col_D, col_D

    return NDWIcomp_col, NDVIcomp_col, QA_dict_L, QA_dict_S2, COMP_TYPES, COMP_NANSCORES, NBL7_DATES, NBL8_DATES, NBL9_DATES, NBS2_DATES


