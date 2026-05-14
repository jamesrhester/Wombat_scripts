# A series of routines for outputting data
import math,copy

# Output CIF data, including metadata
def write_cif_block(ds, sink):
    """Write the dataset in CIF format"""
    from CifFile import CifFile, CifLoopBlock
    from datetime import datetime
    block_name = sanitise(ds.title[0:17])
# Create a block name from dataset name and current time
    current_time =  datetime.now().isoformat()
    block_name = str(block_name) + str(current_time)
    metadata_store = ds.harvest_metadata("CIF")
    alldata.NewBlock(block_name,blockcontents=metadata_store)
    # Create a unique block id
    username = '?'
    try:
        username = sanitise(str(ds['user']))
    except:
        pass
    metadata_store["_pd_block_id"] = "%s|%s|%s" % (block_name,current_time,username) 
    metadata_store["_audit_creation_date"] = current_time
    metadata_store["_audit_creation_method"] = "Automatically generated from raw NeXuS data file by Gumtree routines"
    metadata_store.AddCifItem((
        (("_audit_conform_dict_name", "_audit_conform_dict_version", "_audit_conform_dict_location"),),
        ((("cif_core.dic","cif_pd.dic"),("2.3.1","1.0.1"),
         ("ftp://ftp.iucr.org/pub/cifdics/cif_core_2.3.1.dic","ftp://ftp.iucr.org/pub/cifdics/cif_pd_1.0.1.dic")),))
)   
    import time
    angles = map(lambda a:("%.5f" % a),ds.axes[0])
    ints = map(lambda a,b:format_esd(a,b),ds,ds.var)
    esds = map(lambda a:"%.5f" % math.sqrt(a),ds.var)
    metadata_store.AddCifItem((
            (("_pd_proc_2theta_corrected", "_pd_proc_intensity_net", "_pd_proc_intensity_net_esd"),),
  #          ((("%10.5f" % ds.axes[0][0],),(format_esd(ds[0],ds.var[0]),),("%15.5f" % math.sqrt(ds.var[0]),)),))
            ((angles,ints,esds),))

            )

def write_cif_data(ds, filename):

    from CifFile import CifFile

    alldata = CifFile()
    write_cif_block(ds, alldata)
    if not filename[-3:]=='cif':
        filename = filename+'.cif'
    fh = open(filename,"w")
    fh.write(str(alldata))
    fh.close()

def write_esg_data(segment_data, filename):
    """
    Write data for texture software. One file should have top/middle/bottom data
    blocks in CIF format, containing necessary metadata.
    """
    from CifFile import CifFile
    alldata = CifFile()
    for one_segment in segment_data:
        write_cif_block(segment_data[one_segment], alldata)

    if not filename[-3:]=='esg':
        filename = filename+'.esg'
    fh = open(filename,"w")
    fh.write(str(alldata))
    fh.close()
    
def sanitise(badstring):
    """Remove dodgy characters from username"""
    badstring = badstring.replace(' ','-')
    return badstring.replace('_','-')

# A fairly primitive formatter that will break on pathological strings,
# for example semicolons as the first character in a line in a multi-line string
#
def prepare_string(bad_string):
    """Format a string for a CIF file"""
    if "\r" in bad_string or "\n" in bad_string:
        #Wrap the text nicely
        good_string = "\n;\n"+bad_string+"\n;\n";
    else: 
       	# Note that this is not perfect; we may have a string containing whitespace
       	# immediately following a single or double quote, in which case we must use
       	# the semicolon-delimited variation.  We don't check for this here, but a
       	# simple regular expression should catch it
        if " " in bad_string or "\t" in bad_string:
            if "\"" in bad_string:
                good_string = "'" + bad_string + "'";
            else: 
                good_string = "\"" + bad_string + "\"";
        else: good_string = bad_string;
    return good_string;

#
def format_esd(number,var):
    """Format the esd as intensity(nn), where nn follows the ISO rules"""
    # The idea is to multiply the error until we have a number between
    # 2 and 20, and the number of significant figures after the
    # decimal point is adjusted according to the size of the error and
    # and the number of multiplications
    err = math.sqrt(var)
    if err<0.000000001: return "%.5f(0)" % number
    err_as_int = err
    outsigfigs = 0
    while err_as_int < 1.95: 
        outsigfigs+=1
        err_as_int *=10
    if outsigfigs <= 0:  #Error somewhere to left of decimal point
        # should format in scientific notation but too lazy
        flt_format = "%.0f(%d)"
    else:
        flt_format = "%%.%df(%%d)" % outsigfigs
    # print 'Format: %s, number %f, err %f, err_as_int %d' % (flt_format,number,err,int(round(err_as_int)))
    return flt_format % (number, int(round(err_as_int)))

###################################################
# Output xyd data (three column ASCII)
# Change comment char and extension for other variants
###################################################
def write_xyd_data(ds,filename,comment="#",extension="xyd"):
    from datetime import datetime
    current_time =  datetime.now().isoformat()
    angles = map(lambda a:("%.5f" % a),ds.axes[0])
    ints = map(lambda a:"%.2f" % a,ds)
    esds = map(lambda a:"%.5f" % math.sqrt(a),ds.var)
    if not filename[-3:]==extension[-3:]:
        filename = filename+'.' + extension
    fh = open(filename,"w")
    fh.write(comment+" Data from file %s, written %s\n" % (ds.title[0:17],str(current_time)))
    fh.write(comment+" %10s %10s %10s\n" % ("Angle","Intensity","Error"))
    for point in zip(angles,ints,esds):
        fh.write("  %10s %10s %10s\n" % (point[0],point[1],point[2]))
    fh.close()

###################################################
#  FXYE format (New GSAS)                         #
###################################################
def write_fxye_data(ds,filename):
    """
 *  This class implements a GSAS FXYE format exporter for powder data.  The specifications were
 *  obtained from Brian Toby.  Only one histogram is supported, although the format allows
 *  multiple histograms. 
 *  
 *  For full specification of the format, see GSAS Technical Manual, 26/09/2004 Edition p 217.
 *  Salient points for our output here: 
 *  (1) comments begin with '#' as the first character and each line should be padded
 *  to 80 characters.
 *  (2) Comments must not appear after the 'BANK' line
 *  (3) The file consists of one or more comments lines followed by data block lines
 *  (4) A data block starts with a 'BANK' line
 *  (5) BANK line format is 'BANK' IBANK NCHAN NREC BINTYP TYPE TYPE TYPE TYPE
 *  (6) IBANK is a unique ID, NCHAN is number of points, NREC is number of lines, TYPE is
 *  'FXYE' in this case
 *  (7) Data are then presented as one point per line: two theta in centidegrees, 
 *  intensities, errors: 3(F8.0,F7.4,F5.4) in Fortran notation
 *  (8) BINTYPE needs to be CONS even though it is not relevant for FXYE
 *  @author jrh"""

    from datetime import datetime
    current_time =  datetime.now().isoformat()
    angles = map(lambda a:("%8.0f" % a),ds.axes[0]*100.0)
    thlen = len(angles)
    avstep = (ds.axes[0][-1]-ds.axes[0][0])/(thlen - 1)
    ints = map(lambda a:"%7.4f" % a,ds)
    esds = map(lambda a:"%5.4f" % math.sqrt(a),ds.var)
    if not filename[-3:]=='xye':
        filename = filename+'.xye'
    fh = open(filename,"w")
    fh.write("%-80s\n" % "# Data from file %s, written %s. GSAS FXYE format" % (ds.title[0:17],str(current_time)))
    fh.write("BANK 01 %5d %5d CONS %s %8.2f 0 0 FXYE\n" % (thlen,thlen,angles[0],avstep*100))
    for ang,intensity,esd in zip(angles,ints,esds):
        fh.write("%s %s %s\n" % (ang,intensity,esd))
    fh.close()
