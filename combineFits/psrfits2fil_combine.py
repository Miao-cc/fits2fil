#!/home/mcc/software/py3env/bin/python

from __future__ import print_function
from builtins import range
import numpy as np
from presto import psrfits
from presto import filterbank
from presto import sigproc
import optparse
import sys
import os
import time

def translate_header(psrfits_file,**args):
    fits_hdr = psrfits_file.header
    subint_hdr = psrfits_file.fits['SUBINT'].header 
    subint_data = psrfits_file.fits['SUBINT'].data
    fil_header = {}

    if fits_hdr['TELESCOP'] in sigproc.telescope_ids:
        fil_header["telescope_id"] = \
                    sigproc.telescope_ids[fits_hdr['TELESCOP']]
    else:
        fil_header["telescope_id"] = -1
    if fits_hdr['BACKEND'] in sigproc.machine_ids:
        fil_header["machine_id"] = \
                    sigproc.machine_ids[fits_hdr['BACKEND']]
    else:
        fil_header["machine_id"] = 8

    fn = psrfits_file.filename
    #fil_header["rawdatafile"] = os.path.basename(fn.replace(".fits",".fil"))
    fil_header["rawdatafile"] = os.path.basename("")
    fil_header["data_type"] = 1 # filterbank
    fil_header["source_name"] = fits_hdr['SRC_NAME']
    fil_header["barycentric"] = 0 # always not barycentered?
    fil_header["pulsarcentric"] = 0 # whats pulsarcentric?
    fil_header["az_start"] = subint_data[0]['TEL_AZ']
    fil_header["za_start"] = subint_data[0]['TEL_ZEN']
    fil_header["src_raj"] = float(fits_hdr['RA'].replace(':',''))
    fil_header["src_dej"] = float(fits_hdr['DEC'].replace(':',''))
    fil_header["tstart"] = fits_hdr['STT_IMJD'] + \
                           fits_hdr['STT_SMJD']/86400.0 + \
                           fits_hdr['STT_OFFS']/86400.0
    fil_header["tsamp"] = subint_hdr['TBIN']

    fil_header["tobs"] = subint_hdr['TBIN'] * subint_hdr['NAXIS2'] * subint_hdr['NSBLK']
    fil_header["nsamples"] = subint_hdr['NAXIS2'] * subint_hdr['NSBLK']
    
    fil_header["nbits"] = None # set by user. Input should always be 4-bit.

    # first channel (fch1) in sigproc is the highest freq
    # foff is negative to signify this
    fil_header["fch1"] = fits_hdr['OBSFREQ'] + \
                         np.abs(fits_hdr['OBSBW'])/2.0 - \
                         np.abs(subint_hdr['CHAN_BW'])/2.0
    fil_header["foff"] = -1.0*np.abs(subint_hdr['CHAN_BW'])
    fil_header["nchans"] = subint_hdr['NCHAN']
    #fil_header["nifs"] = subint_hdr['NPOL']
    fil_header["nifs"] = 1

    for key, value in args.items():
        if value is None:
            continue
        elif key=="ra":
            fil_header["src_raj"] = float(value.replace(':',''))
        elif key=="dec":
            fil_header["src_dej"] = float(value.replace(':',''))
        else:
            fil_header[key] = value

    return fil_header

def main(fits_fn, outfn, gaindiff, baselinediff, nbits, \
            apply_weights, apply_scales, apply_offsets, ra, dec):

    with open(fits_fn) as fw:
        fits_list = fw.readlines()
        for num, filename in enumerate(fits_list):
            filename = filename.replace("\n","")
            if num==0:
                start = time.time()
                psrfits_file = psrfits.PsrfitsFile(filename)

                fil_header = translate_header(psrfits_file,ra=ra, dec=dec)
                fil_header['nbits'] = nbits
                outfil = filterbank.create_filterbank_file(outfn, fil_header, \
                                                    nbits=nbits, mode='write')

                # if frequency channels are in ascending order
                # band will need to be flipped
                if psrfits_file.fits['SUBINT'].header['CHAN_BW'] > 0:
                    flip_band=True
                    print("\nFits file frequencies in ascending order.")
                    print("\tFlipping frequency band.\n")
                else:
                    flip_band=False

                # check nbits for input

                print("\nCalculating statistics on first subintegration...")
                subint0 = psrfits_file.read_subint(0, \
                                apply_weights, apply_scales, apply_offsets, total_intensity=False)
                subint0 = (subint0[:,0,:] + subint0[:,1,:]*gaindiff + baselinediff)/2.
                new_max = 3 * np.median(subint0)
                print("\t3*median =",new_max)
                if new_max > 2.0**nbits:
                    scale = True
                    scale_fac = new_max / ( 2.0**nbits )
                    print("\tScaling data by",1/scale_fac)
                    print("\tValues larger than",new_max,"(pre-scaling) "\
                              "will be set to",2.0**nbits - 1,"\n")
                else:
                    scale_fac = 1

                print("Writing data...")
                sys.stdout.flush()
                oldpcnt = ""

            psrfits_file = psrfits.PsrfitsFile(filename)
            for isub in range(int(psrfits_file.nsubints)):
                subint = psrfits_file.read_subint(isub, apply_weights, apply_scales, apply_offsets, total_intensity=False)
                subint = (subint[:,0,:] + subint[:,1,:]*gaindiff + baselinediff)/2.
                if flip_band:
                    subint = np.fliplr(subint)
                subint /= scale_fac
                outfil.append_spectra(subint)

                pcnt = "%d" % (isub*100.0/psrfits_file.nsubints)
                if pcnt != oldpcnt:
                    sys.stdout.write("% 4s%% complete(finish %s, total %s)\r" % (pcnt, num, len(fits_list)))
                    sys.stdout.flush()
        print("Done                              ")
        outfil.close()

    print("Runtime:",time.time() - start)


if __name__=='__main__':
    parser = optparse.OptionParser(prog='psrfits2fil.py',
                    version="v0.2 Paul Scholz, Patrick Lazarus (Sept 2012)",
                    usage = "usage: %prog [options] input_fits")
    parser.add_option("-n",dest='nbits', action='store',
                      default=8, type='int',
                      help="The number of bits in the output .fil file. " +
                           "Default=8")
    parser.add_option("-o",dest='outfn',action='store',
                      default=None, type='string',
                      help="The filename of the output filterbank file. " +
                           "Default: same as .fits input but with .fil extn")
    parser.add_option("-g",dest='gaindiff',action='store',
                      default=None, type='float',
                      help="gain difference between AA and BB")
    parser.add_option("-b",dest='baselinediff',action='store',
                      default=None, type='float',
                      help="baseline difference between AA and BB")

    parser.add_option("--ra",dest='ra',action='store',
                      default=None, type='string',
                      help="ra (J2000)")
    parser.add_option("--dec",dest='dec',action='store',
                      default=None, type='string',
                      help="dec (J2000)")
    parser.add_option("--noweights", dest='apply_weights',
                      default=True, action="store_false",
                      help="Do not apply weights when converting data.")
    parser.add_option("--noscales", dest='apply_scales',
                      default=True, action="store_false",
                      help="Do not apply scales when converting data.")
    parser.add_option("--nooffsets", dest='apply_offsets',
                      default=True, action="store_false",
                      help="Do not apply offsets when converting data.")
    (options, args) = parser.parse_args()

    fits_fn = args[0]

    if options.outfn is None:
        print("A output filename is needed.")
    else:
        if (options.outfn).endswith('.fil'):
            outfn = options.outfn
        else:
            outfn = options.outfn + '.fil'

    if options.gaindiff:
        gaindiff = options.gaindiff
        baselinediff = options.baselinediff
    else:
        gaindiff = 1
        baselinediff = 0

    main(fits_fn, outfn, gaindiff, baselinediff, options.nbits, options.apply_weights,
            options.apply_scales, options.apply_offsets,options.ra,options.dec)
