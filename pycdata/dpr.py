import numpy as np
from astroquery.mast import Observations
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from pycheops.core import load_config
import os

# Path of the current directory
p1 = os.getcwd()

# Path of the pycheops cache directory
config = load_config()
p3 = config['DEFAULT']['data_cache_path']

def tess_data(name, pdc=True, verbose=True):
    """
    Parameters:
    -----------
    name : str
        Name of the planet
    pdc : bool
        Whether to extract PDCSAP flux or not
        Default is True
    verbose : bool
        Boolean on whether to show print updates
        Default is True
    -----------
    return
        tgz file readable to
        pycheops
    """
    try:
        obt = Observations.query_object(name, radius=0.01)
    except:
        raise Exception('The name of the object does not seem to be correct.\nPlease try again...')
    # b contains indices of the timeseries observations from TESS
    b = np.array([])
    for j in range(len(obt['intentType'])):
        if obt['obs_collection'][j] == 'TESS' and obt['dataproduct_type'][j] == 'timeseries':
            b = np.hstack((b,j))
    if len(b) == 0:
        raise Exception('No TESS timeseries data available for this target.\nTry another target...')
    # To extract obs-id from the observation table
    sectors, ticids, pi_name, obsids, exptime = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for i in range(len(b)):
        data1 = obt['dataURL'][int(b[i])]
        if data1[-9:] == 's_lc.fits':
            fls = data1.split('-')
            for j in range(len(fls)):
                if len(fls[j]) == 5:
                    sec = fls[j]
                    tic = fls[j+1]
            sectors = np.hstack((sectors, sec))
            ticids = np.hstack((ticids, tic))
            obsids = np.hstack((obsids, obt['obsid'][int(b[i])]))
            pi_name = np.hstack((pi_name, obt['proposal_pi'][int(b[i])]))
            exptime = np.hstack((exptime, obt['t_exptime'][int(b[i])]))
    if verbose:
        print('Data products found over ' + str(len(sectors)) + ' sectors.')
        print('Downloading them...')
    for i in range(len(sectors)):
        dpr = Observations.get_product_list(obt[int(b[i])])
        cij = 0
        for j in range(len(dpr['obsID'])):
            if dpr['description'][j] == 'Light curves':
                cij = j
        tab = Observations.download_products(dpr[cij])
        lpt = tab['Local Path'][0][1:]
        # Reading fits
        hdul = fits.open(p1 + lpt)
        hdr = hdul[0].header
        dta = Table.read(hdul[1])
        # Available data products
        if pdc:
            fl = np.asarray(dta['PDCSAP_FLUX'])
            fle = np.asarray(dta['PDCSAP_FLUX_ERR'])
        else:
            fl = np.asarray(dta['SAP_FLUX'])
            fle = np.asarray(dta['SAP_FLUX_ERR'])
        mask = np.isfinite(fl)                                # Creating Mask to remove Nans
        bjd1 = np.asarray(dta['TIME'])[mask] + 2457000        # BJD without mask
        fl, fle = fl[mask], fle[mask]                         # Flux and Error in flux without Nans
        bjd2 = Time(bjd1, format='jd')                        # Astropy.Time object to store value in JD
        utc1 = bjd2.to_value('fits')                          # To store date and time in UTC
        cenx = np.asarray(dta['MOM_CENTR1'])[mask]             # Centroid-x
        ceny = np.asarray(dta['MOM_CENTR2'])[mask]             # Centroid-y
        if pdc:
            bg = np.zeros_like(fl)
        else:
            bg = np.asarray(dta['SAP_BKG'])[mask]             # Background for SAP Flux
        mjd1 = bjd1 - 2400000.
        status, event, dark, conta, conta_err, smear, smerr, roll, locx, locy =\
            np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl),\
            np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl)
        # Creating astropy Table for storing data
        tab = Table()
        tab['UTC_TIME'], tab['MJD_TIME'], tab['BJD_TIME'], tab['FLUX'], tab['FLUXERR'], tab['STATUS'],\
            tab['EVENT'], tab['DARK'], tab['BACKGROUND'], tab['CONTA_LC'], tab['CONTA_LC_ERR'], tab['SMEARING_LC'],\
            tab['SMEARING_LC_ERR'], tab['ROLL_ANGLE'], tab['LOCATION_X'], tab['LOCATION_Y'], tab['CENTROID_X'], tab['CENTROID_Y'] =\
            utc1, mjd1, bjd1, fl, fle, status, event, dark, bg, conta, conta_err, smear, smerr, roll, locx, locy, cenx, ceny
        # Writing table data to fits file
        tab.write('TIC_' + ticids[i] + '_' + sectors[i] + '.fits', format='fits')
        tb_fits = fits.open('TIC_' + ticids[i] + '_' + sectors[i] + '.fits')
        tb_fits_hdr = tb_fits[1].header
        hdr1 = hdr[8:]
        crds = hdr1.cards
        for j in range(len(crds)):
            tb_fits_hdr.append(crds[j])
        tb_fits_hdr.append(('PI_NAME', pi_name[i], 'Name of the PI'))                    # Adding PI Name to the header file
        tb_fits_hdr.append(('TARGNAME', name, 'Name of the target'))                     # Adding name of the target
        tb_fits_hdr.append(('OBSID', obsids[i], 'Observation ID'))                       # Adding Observation ID
        tb_fits_hdr.append(('NEXP', '1', 'Number of Images co-added'))                   # NEXP (Not important for TESS)
        tb_fits_hdr.append(('EXPTIME', exptime[i], '(sec) Exposure Time'))               # Adding Exposure time
        tb_fits_hdr.append(('TEXPTIME', exptime[i], '(sec) Total exposure time'))        # Adding Total Exposure time
        tb_fits_hdr.append(('SPECTYPE', str(hdr1['TEFF']), '(K) Temperature'))           # Adding Spectral type/not available, so adding Teff
        tb_fits_hdr.append(('PIPE_VER', '0.0.0', 'Pipeline'))                            # Adding Pipeline Version
        tb_fits_hdr.append(('AP_RADI', 11.0, '(px) Radius of Aperture'))                 # Aperture radius
        tb_fits_hdr.rename_keyword('RA_OBJ', 'RA_TARG')                                  # RA and DEC keyword name change
        tb_fits_hdr.rename_keyword('DEC_OBJ', 'DEC_TARG')
        os.system('rm ' + 'TIC_' + ticids[i] + '_' + sectors[i] + '.fits')
        os.system('rm -r mastDownload')
        name1 = 'CH_PR' + ticids[i][8:14] + '_TG' + ticids[i][14:17] + sectors[i][-4:] + '_V0000_TIC'
        name2 = 'CH_PR' + ticids[i][8:14] + '_TG' + ticids[i][14:17] + sectors[i][-4:] + '_SCI_COR_Lightcurve-DEFAULT_V0000_TIC'
        tb_fits.writeto(name2 + '.fits')
        os.system('tar -cvzf ' + name1 + '.tgz ' + name2 + '.fits')
        os.system('rm ' + name2 + '.fits')
        os.system('mv ' + name1 + '.tgz ' + p3 + '/' + name1 + '.tgz')
        if verbose:
            print('----------------------------------------------------------------------------------------')
            print('Name\t\tTIC-id\t\t\tSector\t\t\t.tgz name')
            print('----------------------------------------------------------------------------------------')
            print(name + '\t\t' + ticids[i][8:] + '\t\t' + sectors[i][-4:] + '\t\t' + name1)


def pipe_data(name, fileid, imgette=True, flag=0):
    """
    Parameters:
    -----------
    name : str
        name of the fits file
        containing the PIPE data
    fileid : str
        unique file id for the target visit
    imegette : bool
        Boolean on whether the data comes
        from imegette or sub-array
        Default is True.
    -----------
    return
        tgz file readable to
        pycheops
    """
    fileid = fileid[:-6]
    hdul = fits.open(name)
    hdr = hdul[1].header
    dta = Table.read(hdul[1])
    # Masking datasets
    flg = np.asarray(dta['FLAG'])                 # Flagged data
    msk = np.where(flg==0)[0]                     # Creating a mask to remove flg>0 values
    # Extracting data from file
    mjd1, bjd1 = np.asarray(dta['MJD_TIME'])[msk], np.asarray(dta['BJD_TIME'])[msk]                 # Various time formats stored
    bjd2 = Time(bjd1, format='jd')                                                                  # Converting to UTC time
    utc1 = bjd2.to_value('fits')
    # Storing Data
    fl, fle = np.asarray(dta['FLUX'])[msk], np.asarray(dta['FLUXERR'])[msk]
    roll, xc, yc, tft2 = np.asarray(dta['ROLL'])[msk], np.asarray(dta['XC'])[msk],\
        np.asarray(dta['YC'])[msk], np.asarray(dta['thermFront_2'])[msk]
    if imgette:
        bg = np.asarray(dta['BG'])[msk] * 50 * 50
    else:
        bg = np.asarray(dta['BG'])[msk] * 200 * 200
    # For relative weights of the first arbitrary PSF principal components
    Us_n = np.array([])      # To store the names of U
    Us = []                  # To store the values of U
    cols = dta.colnames      # Name of all columns
    for j in range(len(cols)):
        if cols[j][0] == 'U':
            Us_n = np.hstack((Us_n, cols[j]))
    for j in range(len(Us_n)):
        usn = np.asarray(dta[Us_n[j]])[msk]
        Us.append(usn)
    # Making an arrays of zero for other elements
    status, event, dark, conta, conta_err, smear, smerr, locx, locy =\
        np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl),\
        np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl), np.zeros_like(fl)
    # Creating an astropy table for storing the data
    tab = Table()
    tab['UTC_TIME'], tab['MJD_TIME'], tab['BJD_TIME'], tab['FLUX'], tab['FLUXERR'], tab['STATUS'],\
        tab['EVENT'], tab['DARK'], tab['BACKGROUND'], tab['CONTA_LC'], tab['CONTA_LC_ERR'], tab['SMEARING_LC'],\
        tab['SMEARING_LC_ERR'], tab['ROLL_ANGLE'], tab['LOCATION_X'], tab['LOCATION_Y'], tab['CENTROID_X'], tab['CENTROID_Y'] =\
        utc1, mjd1, bjd1, fl, fle, status, event, dark, bg, conta, conta_err, smear, smerr, roll, locx, locy, xc, yc
    for j in range(len(Us_n)):
        tab[Us_n[j]] = Us[j]
    # Saving the table first
    tab.write(fileid + '_V0000_PIPE.fits', format='fits')
    # Making changes to the saved file
    tb_fits = fits.open(fileid + '_V0000_PIPE.fits')
    tb_fits_hdr = tb_fits[1].header
    hdr1 = hdr[8:]
    hdr2 = hdr1[:-47]
    crds = hdr2.cards
    for j in range(len(crds)):
        tb_fits_hdr.append(crds[j])
    tb_fits_hdr.append(('AP_RADI', 0.0, '(px) Radius of Aperture/NA -- PSF Photometry'))                 # Aperture radius
    nn1 = tb_fits_hdr['TARGNAME']
    os.system('rm ' + fileid + '_V0000_PIPE.fits')
    tb_fits.writeto(fileid + '_SCI_COR_Lightcurve-DEFAULT_V0000_PIPE.fits')
    os.system('tar -cvzf ' + fileid + '_V0000_PIPE.tgz ' + fileid + '_SCI_COR_Lightcurve-DEFAULT_V0000_PIPE.fits')
    os.system('mv ' + fileid + '_V0000_PIPE.tgz ' + p3 + '/' + fileid + '_V0000_PIPE.tgz')
    os.system('rm ' + fileid + '_SCI_COR_Lightcurve-DEFAULT_V0000_PIPE.fits')
    print('-----------------------------------------------------------------------------')
    print('Name of the file\t\t\t\t.tgz file')
    print('-----------------------------------------------------------------------------')
    print(nn1 + '\t\t\t\t' + fileid + '_V0000_PIPE.tgz')
    # For meta files
    tab1 = Table()
    tab1['thermFront_2'] = tft2
    tab1.write(fileid + '-meta.fits')
    os.system('mv ' + fileid + '-meta.fits ' + p3 + '/' + fileid + '_V0000-meta.fits')
    print('meta\t\t\t\t\t' + fileid + '_V0000-meta.fits')