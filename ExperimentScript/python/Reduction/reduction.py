'''
@author:        davidm,jamesh
@organization:  ANSTO

@version:  1.7.0.1
@date:     14/02/2013
@source:   http://www.nbi.ansto.gov.au/echidna/scripts/reduction.py

'''

'''
Issues:
    !!! float_copy doesn't deep-copy meta data

    !!! Normalization doesn't effect monitor_data, detector_data, ...
    !!! Stitching - do I need to sum monitor_data, detector_data, etc ?
    
            rs.total_counts  += ds_frame.total_counts
            rs.monitor_data  += ds_frame.monitor_data
            rs.detector_time += ds_frame.detector_time
            rs.detector_data += ds_frame.detector_data

    following is assumed at the moment:
    
        bm1_counts ~ bm2_counts ~ bm3_counts ~ detector_time

'''

from WOM import WOM

from gumpy.nexus import *
from gumpy.commons.jutils import jintcopy

from math import *
import os.path as path
import AddCifMetadata

def copy_metadata_deep(ds, dfrom):
    ds.__copy_metadata__(dfrom, deep=True)
def copy_metadata_shallow(ds, dfrom):
    ds.__copy_metadata__(dfrom, deep=False)

# to supplement math-library 
def sqr(value):
    return value * value

class DefaultOkIter:
    def next(self):
        return 1 # always OK-pixel

epsilon = 1e-5 # has to be grater than zero
def _min(l, r):
        if l <= r:
            return l
        else:
            return r
def _max(l, r):
        if l >= r:
            return l
        else:
            return r
def getCenters(boundaries):
        # check dimensions
        if boundaries.ndim != 1:
            raise AttributeError('boundaries.ndim != 1')
        if boundaries.size < 2:
            raise AttributeError('boundaries.size < 2')
            
        # result
        rs = zeros(boundaries.size - 1) # result is one item shorter
        
        rs[:] = boundaries[0:-1]
        rs   += boundaries[1:]
        rs   *= 0.5

        print "Converted %f to %f (%d) to %f to %f (%d)" % (boundaries[0], boundaries[-1],
                                                            len(boundaries), rs[0], rs[-1],
                                                            len(rs))

        return rs
    
def applyNormalization(ds, reference, target=-1):
    """Normalise datasets ds by multiplying by target/reference.  Beam monitor counts, count time and total counts are
       all adjusted by this amount.  Reference is a string referring to a particular location in the dataset, and
       target is the target value to which they will be adjusted.  If target is not specified, the maximum value of
       the reference array is used and reported for further use. The variance of the values in the reference array
       is assumed to follow counting statistics.  We modify the input dataset rather than creating a new dataset
       as files on Wombat are so large. We need to make sure that the
       input dataset is converted to float, otherwise the multiplcation
       will be integer."""
    print 'normalization of', ds.title
    # Store reference name for later
    refname = str(reference)
    # Normalization
    reference = getattr(ds,reference)

    # check if reference/target is a number
    # TODO: gumpy doesn't allow us to handle a scalar with variance
    # for multiplying arrays, so we can't propagate variance at present
    numericReference = isinstance(reference, (int, long, float))
    
    # check arguments
    if not numericReference:
        if reference.ndim != 1:
            raise AttributeError('reference.ndim != 1')
        if reference.shape[0] != ds.shape[0]:
            raise AttributeError('reference.shape[0] != ds.shape[0] (%d != %d)' % (reference.shape[0],ds.shape[0]))

    def do_norm(rs, f, varf):
        # We propagate errors in the data, but not in
        # the ancillary values
        print 'In do_norm, given %f(%f)' % (f,varf)
        # Funny syntax below to make sure we write into the original area,
        # not assign a new value
        rs.var *= f * f
        rs.var += varf * rs * rs
        rs.storage       *= f
        try:             #These may be absent in some cases
            rs.bm1_counts    *= f
            rs.bm2_counts    *= f
            rs.bm3_counts    *= f
            rs.detector_time *= f
            rs.total_counts  *= f
        except AttributeError:
            pass
       
    # normalization
    if numericReference and target > 0:
        # We have a single number to refer to for normalisation, so
        # we are effectively scaling everything by a single number
        scale_factor = float(target)/reference
        variance = scale_factor * target/(reference*reference)
        do_norm(ds, scale_factor, variance)
        info_string = "Data multiplied by %f with variance %f" % (scale_factor,variance)
    elif not numericReference:
        # Each step has a different value, and we manually perform the
        # error propagation 
        reference = Data(reference)
        if target <= 0:
            target = reference.max()
        for i in xrange(ds.shape[0]):

            # If we have zero counts, we engineer a no-op
            
            if reference[i] == 0:
                reference[i] = float(target)
            f = float(target)/reference[i]
            v = f*target/(reference[i]*reference[i])
            # Funny syntax below to make sure we write into the original area,
            # not assign a new value
            tar_shape = [1,ds.shape[1],ds.shape[2]]
            tar_origin = [i,0,0]
            rss = ds.storage.get_section(tar_origin,tar_shape).get_reduced()
            rsv = ds.var.get_section(tar_origin,tar_shape).get_reduced()
            ds.var[i] = rsv*f * f
            ds.var[i] += v * rss * rss
            ds.storage[i] = ds.storage[i]*f
        info_string = "Data normalised to %f on %s with error propagation assuming counting statistics" % (float(target),refname)
        print info_string
    else:
        # interesting note - if we get here, we are passed a single reference number
        # and a negative target, meaning that we use the reference as the target and
        # end up multiplying by 1.0, so no need to do anything at all.
        target = reference
        info_string = "No normalisation applied to data."
        print info_string
    ds.add_metadata('_pd_proc_info_data_reduction',info_string, tag="CIF",append=True)
    print 'normalized:', ds.title
    return target

def getSummed(ds, applyStth=0.0, contribs = None, use_zeros = False):
    """ A faster version of Dataset.intg which discards a lot of the
    metadata handling that intg needs to do"""
    import time
    print 'summation of', ds.title

    print 'Applying stth of %d' % applyStth
    print 'Axis starts: %f %f %f' % (ds.axes[0][0], ds.axes[1][0], ds.axes[2][0])
    # check arguments
    if ds.ndim != 3:
        raise AttributeError('ds.ndim != 3')
    if applyStth and (ds.axes[2].title != 'x_pixel_angular_offset'):
        raise AttributeError('ds.axes[2].title != x_pixel_angular_offset')

    # sum first dimension of storage and variance
    frame_count = ds.shape[0]

    # make sure we are using pixel centres not boundaries

    if len(ds.axes[-1]) == ds.shape[-1] + 1:
        title = ds.axes[-1].title
        units = ds.axes[-1].units
        ds.axes[-1] = getCenters(ds.axes[-1])
        ds.axes[-1].title = title
        ds.axes[-1].units = units
        
    # detect single frames
    if frame_count == 1:
        rs = ds.get_reduced()
    else:
        base_data = zeros_like(ds.storage[0])
        base_var = zeros_like(ds.storage[0])
        for frame in xrange(0, frame_count):
            base_data          += ds.storage[frame]
            base_var           += ds.var[frame]

        # finalize result
        rs = Dataset(base_data)
        rs.var = base_var
        rs.axes[1] = ds.axes[2]
        rs.axes[0] = ds.axes[1]

    ok_map = ones(rs.shape)
    if contribs != None:
        ok_map = contribs[0]  #should be identical for all
    elif use_zeros:
        ok_map[rs <= 0] = 0

    rs.title = ds.title

    if applyStth:  #we check for identity
        rs.axes[1] += applyStth
        rs.axes[1].title = 'Two theta'

    print 'summed frames:', frame_count
    rs.copy_cif_metadata(ds)
    print 'Output dtype' + `rs.dtype`
    assert len(rs.axes[1]) == rs.shape[-1]
    return rs, ok_map

def boundaries_to_mdpts(ds, axno = 2):

    if len(ds.axes[axno]) == ds.shape[-1] + 1:
        title = ds.axes[-1].title
        units = ds.axes[-1].units
        ds.axes[axno] = getCenters(ds.axes[axno])
        ds.axes[axno].title = title
        ds.axes[axno].units = units

def getStepSummed(ds, contribs = None, use_zeros=False):
    """ As for `getSummed`, but additionally offsets each frame to model 2th
    movement between frames. `contribs`, if present, is non-zero for active
    pixels. If `use_zeros` is `True`, zero-valued pixels are ignored. `use_zeros`
    is ignored if `contribs` is not None."""
    import time
    print 'step summation of', ds.title

    # check arguments
    if ds.ndim != 3:
        raise AttributeError('ds.ndim != 3')
    if ds.axes[2].title != 'x_pixel_angular_offset':
        raise AttributeError('ds.axes[2].title != x_pixel_angular_offset, is %s' % ds.axes[2].title)

    # sum first dimension of storage and variance

    frame_count = ds.shape[0]

    # make sure we have pixel positions not boundaries

    if len(ds.axes[2]) == ds.shape[-1] + 1:
        title = ds.axes[-1].title
        units = ds.axes[-1].units
        ds.axes[2] = getCenters(ds.axes[2])
        ds.axes[2].title = title
        ds.axes[2].units = units
 
    # detect single frames

    if frame_count == 1:
        rs = ds.get_reduced()
        ok_map = ones(rs.shape)
        if use_zeros and contribs == None:
            ok_map[rs<=0]=0
        elif contribs != None:
            ok_map = contribs.get_reduced()
        new_axis = ds.axes[2] + ds.axes[0][0]
    else:
            # calculate step size as a proportion of x pixel step
    
        _, pixel_step, bin_size = get_wire_step(ds)

        if bin_size == 0:  #no significant detector movement
            return getSummed(ds, applyStth=ds.axes[0][0], contribs = contribs, use_zeros = use_zeros)
        
        if pixel_step < 0:     # step bigger than wire separation
            pixel_step = -1 * pixel_step
            extra_width = pixel_step * (frame_count - 1)
            extra_points =  ds.axes[2][-1] + bin_size*(arange(extra_width)+1) + ds.axes[0][0]# step by wire separation
        else:                  # step by stth step
            extra_width = pixel_step * (frame_count - 1)
            extra_points = ds.axes[2][-1] + ds.axes[0][1:]
        print "Pixel_step %d, bin size %f, extra %f" % (pixel_step, bin_size, extra_width)
        new_shape = ds.shape[1],ds.shape[2]+ int(extra_width)
        base_data = zeros(new_shape)
        base_var = zeros(new_shape)
        ok_map = zeros(new_shape)
        
        for frame in xrange(0, frame_count):
            start_index = frame*pixel_step
            finish_index = ds.shape[2]+frame*pixel_step
            base_data[:,start_index:finish_index]  += ds.storage[frame]
            base_var[:,start_index:finish_index]   += ds.var[frame]
            mask = ones(ds[frame].shape)
            if contribs != None:
                ok_map[:, start_index:finish_index] += contribs[frame]
            elif use_zeros:
                mask[ds.storage[frame]<=0]=0
                ok_map[:,start_index:finish_index] += mask
            else:
                ok_map[:,start_index:finish_index] += mask

        # finalize result
        rs = Dataset(base_data)
        rs.var = base_var
        new_axis = zeros(ds.axes[2].shape[0]+int(extra_width))
        print 'New axis length is %d + %d' % (ds.axes[2].shape[0], int(extra_width))
        new_axis[0:ds.axes[2].shape[0]] = ds.axes[2]+ds.axes[0][0]
        print 'Original 2theta: %s to %s' % (new_axis[0], new_axis[ds.axes[2].shape[0]-1])
        new_axis[ds.axes[2].shape[0]:] = extra_points
        print "Extra points " + repr(new_axis[ds.axes[2].shape[0]:])
        rs.axes[0] = ds.axes[1]

    rs.title = ds.title

    print 'summed frames with stth:', frame_count
    rs.copy_cif_metadata(ds)
    rs.set_axes([ds.axes[1],new_axis],anames=["Vertical offset","Two theta"],aunits=["mm","Degrees"])

    return rs, ok_map

def read_efficiency_cif(filename):
    """Return a dataset,variance stored in a CIF file as efficiency values"""
    import time
    from Formats import CifFile
    print 'Reading in %s as CIF at %s' % (filename,time.asctime())
    eff_cif = CifFile.CifFile(str(filename))
    print 'Finished reading in %s as CIF at %s' % (filename,time.asctime())
    eff_cif = eff_cif['efficiencies']
    eff_data = map(float,eff_cif['_[local]_efficiency_data']) 
    eff_var = map(float,eff_cif['_[local]_efficiency_variance']) 
    final_data = Dataset(Data(eff_data).reshape([128,128]))
    final_data.var = (Array(eff_var).reshape([128,128]))
    print 'Finished reading at %s' % time.asctime()
    return final_data,eff_cif

def get_collapsed(ds):
    """If multiple frames are present, reduce them to a single frame"""
    # This is currently a simple sum
    if ds.ndim <= 2:
        return ds
    rs = ds.sum(axis=0)
    print 'Collapsing axis 0:' + `ds.shape` + '-> ' + `rs.shape`
    return rs
    
def getVerticalIntegrated(ds, okmap=None, normalization=-1, axis=1,top=None,bottom=None):
    print 'vertical integration of', ds.title
    print 'Summing over axis %d in shape %s' % (axis, ds.shape)
    start_dim = ds.ndim

    # check shape
    if (okmap is not None) and (ds.shape != okmap.shape):
        raise AttributeError('ds.shape != okMap.shape')    

    if okmap is None:
        okmap = ones(ds.shape,dtype=int)
    # JRH strategy: we need to sum vertically, accumulating individual pixel
    # errors as we go, and counting the contributions.
    #
    # The okmap should give us contributions by summing vertically
    # Note that we are assuming all observed counts are positive
    
    import time
    if bottom is None or bottom < 0: bottom = 0
    if top is None or top >= ds.shape[axis]: top = ds.shape[axis]-1
    working_slice = ds.take(range(bottom,top), axis = axis)
    totals = working_slice.intg(axis=axis)
    contribs = okmap.intg(axis=axis)
    #
    # We have now reduced the scale of the problem by 100
    #
    # Normalise to the maximum number of contributors
    max_contribs = float(contribs.max())
    min_contribs = float(contribs.min())
    #
    print 'Maximum/minimum no of contributors %f/%f' % (max_contribs,min_contribs)
    contribs = contribs/max_contribs  #
    save_var = totals.var
    totals = totals / contribs        #Any way to avoid error propagation here?
    totals.var = save_var/contribs

    # finalize result
    totals.title = ds.title
    totals.copy_cif_metadata(ds)
    info_string = "Data were vertically integrated from pixels %d to %d (maximum number of contributors %d)." % (bottom,top,max_contribs)
    
    # normalize result if required
    if normalization > 0:
        totals *= (float(normalization) / totals.max())
        totals.title = totals.title
        info_string += "The maximum intensity was then normalised to %f counts." % float(normalization)
    # check if any axis needs to be converted from boundaries to centers
    new_axes = []
    for i in range(totals.ndim):
        if len(totals.axes[i]) == totals.shape[i] + 1:
            print "Correcting axis %d from boundaries to centers"
            new_axes.append(getCenters(totals.axes[i]))
        else:
            new_axes.append(totals.axes[i])
        print 'Axis %d: %s' % (i,totals.axes[i].title)
    old_names = map(lambda a:a.name,totals.axes)
    old_units = map(lambda a:a.units,totals.axes)
    old_names[-1] = 'Two theta'
    old_units[-1] = 'Degrees'
    totals.set_axes(new_axes,anames=old_names,aunits=old_units)
    return totals


def getEfficiencyCorrectionMap(van, bkg, okMap=None):
    print 'create efficiency correction map...'

    # [davidm] I do not recommend to automatically sum the input, because the result may not be normalized correctly
    # check if summation is required
    # (it does not matter if we subtract each bkg-frame from each van-frame and then sum the result or
    #  if we subtract the summed bkg-frames from the summed van-frames)
    #if van.ndim == 4:
    #    van = Reduction.getSummed(van, applyStth=False, show=False)
    #if bkg.ndim == 4:
    #    bkg = Reduction.getSummed(bkg, applyStth=False, show=False)

    # check dimensions
    if van.ndim != 2:
        raise AttributeError('van.ndim != 2')
    if bkg.ndim != 2:
        raise AttributeError('bkg.ndim != 2')
    if (okMap is not None) and (okMap.ndim != 2):
        raise AttributeError('okMap.ndim != 2')

    # check shape
    if van.shape != bkg.shape:
        raise AttributeError('van.shape != bkg.shape')
    if (okMap is not None) and (van.shape != okMap.shape):
        raise AttributeError('van.shape != okMap.shape')

    # result     
    rs  = van.float_copy()
    rs -= bkg

    # iterators
    rs_val_iter = rs.item_iter()
    rs_var_iter = rs.var.item_iter()
    # special check for okMap
    if okMap is not None:
        ok_iter = okMap.item_iter()
    else:
        ok_iter = DefaultOkIter()

    px = 0.0 # number of pixels used
    sm = 0.0 # sum of inverted values required for mean value

    try:
        while True:
            rs_val = rs_val_iter.next()
            rs_var_iter.next()

            if (ok_iter.next() > epsilon) and (rs_val > epsilon): # to compensate for floating-point error
                px += 1
                sm += rs_val
            else:
                # for bad-pixel set result to zero
                rs_val_iter.set_curr(0.0)
                rs_var_iter.set_curr(0.0)

    except StopIteration:
        pass

    # finalize result
    av = sm / px # average
    rs = av / rs # = 1 / (rs / av) # to obtain inverted result

    rs.title = 'Efficiency Map [based on %s]' % van.title

    return rs
 # incomplete

''' corrections '''

def getBackgroundCorrected(ds, bkg, norm_ref=None, norm_target=-1):
    """Subtract the background from the supplied dataset, after normalising the
    background to the specified counts on monitor given in reference"""
    print 'background correction of', ds.title

    # normalise
    if norm_ref:
            applyNormalization(bkg,norm_ref,norm_target)
    if ds.ndim == bkg.ndim:
        # check shape
        if ds.shape != bkg.shape:
            raise AttributeError('ds.shape != bkg.shape')

        # result
        rs = ds - bkg
        rs.copy_cif_metadata(ds)

        # ensure that result doesn't contain negative pixels
        rs[rs < 0] = 0

        print 'background corrected frames:', 1

    elif ds.ndim == 3 and bkg.ndim == 2:
        # check arguments
        if ds.axes[0].title != 'run_number':
            raise AttributeError('ds.axes[0].title != run_number')
        if ds.shape[1:] != bkg.shape:
            raise AttributeError('ds.shape[1:] != bkg.shape')

        # result
        rs = ds.__copy__()
        rs.copy_cif_metadata(ds)
        for frame in xrange(ds.shape[0]):
            rs[frame, 0] -= bkg

        # ensure that result doesn't contain negative pixels
        rs[rs < 0] = 0

        print 'background corrected frames:', ds.shape[0]

    else:
        raise AttributeError('ds.ndim != 2 or 3')
    # finalize result
    rs.title = ds.title
    info_string = 'Background subtracted using %s' % str(bkg.title)
    if norm_ref:
        info_string += 'after normalising to %f using monitor %s.' % (norm_target,norm_ref)
    else:
        info_string += 'with no normalisation of background.'
    rs.add_metadata("_pd_proc_info_data_reduction",info_string,tag="CIF",append=True)
    return rs

def getEfficiencyCorrected(ds, eff):
    print 'efficiency correction of', ds.title

    # check dimensions
    if eff.ndim != 2:
        raise AttributeError('eff.ndim != 2')

    if ds.ndim == 2:
        # check shape
        if ds.shape != eff.shape:
            raise AttributeError('ds.shape != eff.shape')

        # result
        rs = ds * eff

        print 'efficiency corrected frames:', 1

    elif ds.ndim == 3:
        # check arguments
        if ds.axes[1].title != 'y_pixel_offset':
            raise AttributeError('ds.axes[1].title != y_pixel_offset')
        if ds.shape[1:] != eff.shape:
            raise AttributeError('ds.shape[1:] != eff.shape')

        # result
        rs = zeros(ds.shape,dtype='float') #otherwise might be int array
        for frame in xrange(ds.shape[0]):
            rs[frame] = ds[frame] * eff
        print 'efficiency corrected frames:', rs.shape[0]
        rs.axes = ds.axes
    else:
        raise AttributeError('ds.ndim != 2 or 3')
    rs.title = ds.title
    rs.copy_cif_metadata(ds)
    # now include all the efficiency file metadata, except data reduction
    return rs

def sum_datasets(dslist):
    """Add the provided datasets together"""
    #Assume all same length, same axis values
    newds = zeros_like(dslist[0])
    AddCifMetadata.add_standard_metadata(newds)	
    title_info = ""
    proc_info = """This dataset was created by summing points from multiple datasets. Points were 
    assumed to coincide exactly. Data reduction information for the individual source datasets is as follows:"""
    for one_ds in dslist:
    		newds += one_ds
		title_info = title_info + one_ds.title + "+"
		proc_info += "\n\n===Dataset %s===\n" % str(one_ds.title) 
		try:
		    proc_info += one_ds.harvest_metadata("CIF")["_pd_proc_info_data_reduction"]
		except KeyError,AttributeError:
		    pass
    newds.title = title_info[:-1]  #chop off trailing '+'
    newds.axes[0] = dslist[0].axes[0]
    # Add some basic metadata based on metadata of first dataset
    newds.copy_cif_metadata(dslist[0])
    newds.add_metadata('_pd_proc_info_data_reduction',proc_info,"CIF")
    return newds

def convert_to_dspacing(ds):
    if ds.axes[0].name == 'd-spacing':
        return
    try:
        wavelength = float(ds.harvest_metadata("CIF")["_diffrn_radiation_wavelength"])
        print 'Wavelength for %s is %f' % (ds.title,wavelength)
    except KeyError:
        print 'Unable to find a wavelength, no conversion attempted'
        return   #Unable to convert anything
    # Funny call of sin below to avoid problems with sin being shadowed by the
    # standard maths library
    new_axis = wavelength/(2.0*(ds.axes[0]*3.14159/360.0).__sin__())
    ds.set_axes([new_axis],anames=['d-spacing'],aunits=['Angstroms'])
    return 'Changed'

def convert_to_twotheta(ds):
    if ds.axes[0].name == 'Two theta':
        return
    try:
        wavelength = float(ds.harvest_metadata("CIF")["_diffrn_radiation_wavelength"])
    except KeyError:
        print 'Unable to find a wavelength, no conversion attempted'
        return   #Unable to convert anything
    print 'Wavelength for %s is %f' % (ds.title,wavelength)
    new_axis = arcsin(wavelength/(2.0*ds.axes[0]))*360/3.14159
    ds.set_axes([new_axis],anames=['Two theta'],aunits=['Degrees'])
    return 'Changed'

#======== Following code based on Echidna code ============#

def get_wire_step(ds):
     # Determine horizontal pixels per vertical wire interval
    # If the step is larger than the wire separation, return
    # -1 * wires per step and bin size is then wire separation
    # If there is no step, bin size is zero
    wire_pos = ds.axes[-1]
    wire_sep = abs(wire_pos[0]-wire_pos[-1])/(len(wire_pos)-1)
    print "Wire sep %f for %d steps %f - %f" % (wire_sep,len(wire_pos),wire_pos[0],wire_pos[-1])
    det_steps = ds.axes[0]
    bin_size = abs(det_steps[0]-det_steps[-1])/(len(det_steps)-1)
    if bin_size < wire_sep/100: 
        print "Detector not stepped!"
        return wire_sep, 0, 0
    print "Bin size %f for %d steps %f - %f" % (bin_size,len(det_steps),det_steps[0],det_steps[-1])
    if abs(bin_size) < 0.0001:   #too small
        print "Movement too small (%d) ignored" % bin_size
        return det_steps, 0.0, 0.0
    pixel_step = wire_sep/bin_size
    if pixel_step > 0.9:
        pixel_step = int(round(pixel_step))
        bin_size = wire_sep/pixel_step
    else:   # detector step bigger than a single wire interval
        pixel_step = int(round(-1/pixel_step))
        bin_size = wire_sep
    print '%f wire separation, %d steps before overlap, ideal binsize %f' % (wire_sep,pixel_step,bin_size)
    return det_steps,pixel_step,bin_size

#
# Calculate adjusted gain based on matching intensities between overlapping
# sections of data from different detectors
def do_overlap(ds,iterno,algo="FordRollett",unit_weights=False,top=None,bottom=None,
               drop_frames='',drop_wires = '', use_gains = [],do_sum=False, dumpfile=None,
               fix_ignore=0):
    """Calculate rescaling factors for pixel columns based on overlapping data
    regions. Specifying unit weights
    = False will use the variances contained in the input dataset. Note that
    the output dataset will have been vertically integrated as part of the
    algorithm. The vertical integration limits are set by top and bottom, if
    None all points are included.  Drop_frames is a
    specially-formatted string giving a list of frames to be ignored. Drop_wires
    is a similarly-formatted string giving a list of vertical wires to be ignored. If 
    use_gains is not empty, these [val,esd] values will be used instead of those
    obtained from the iteration routine. Dumpfile, if set, will output
    starting values for use by other routines. do_sum will sum each
    detector step before refining. `fix_ignore` sets the number of wires that are
    ignored at the beginning and end of the scan and must match the number used
    for any pre-determined values in `use_gains`. `fix_ignore` is ignored if `use_gains`
    is empty."""
    import time
    from Reduction import overlap
    # Get sensible values
    if top is None: top = ds.shape[1]-1
    if bottom is None: bottom = 0
    # Dimensions are step,vertical,horizontal
    b = ds[:,bottom:top,:].intg(axis=1).get_reduced()
    det_steps,pixel_step,bin_size = get_wire_step(ds)
    if pixel_step <= 0:       #step is too large
        raise ValueError, "Detector step is larger than wire spacing"
    dropped_frames = parse_ignore_spec(drop_frames)
    dropped_cols = parse_ignore_spec(drop_wires)
    # Drop frames from the end as far as we can
    for empty_no in range(b.shape[0]-1,0,-1):
        print "Trying %d" % empty_no
        if empty_no not in dropped_frames:
            break
        dropped_frames.remove(empty_no)
    print "All frames after %d empty so dropped" % empty_no
    b = b[:empty_no+1]
    # Do we need to add dummy missing frames?
    extra_steps = b.shape[0]%pixel_step
    if extra_steps > 0:
        start_drop = b.shape[0]
        # gumpy has no resize
        new_b = zeros([((b.shape[0]/pixel_step)+1)*pixel_step,b.shape[1]])
        new_b[:b.shape[0]] = b
        b = new_b
        extra_dropped_frames = range(start_drop,b.shape[0])
        print "Filled out array from %d to %d with dummy frames" % (start_drop,b.shape[0])
        dropped_frames |= set(extra_dropped_frames)
    else:
        extra_dropped_frames = []
    # Zero out dropped frames
    print 'Dropped frames: ' + `dropped_frames`
    b_zeroed = copy(b)
    # Make a simple array to work out which sectors are missing frames
    frame_check = array.ones(b.shape[0])
    # Additionally zero out all matching steps
    all_zeroed = copy(b)
    region_starts = [a*pixel_step for a in range(b.shape[0]/pixel_step)]
    for frame_no in dropped_frames:
        b_zeroed[frame_no] = 0
        b_zeroed.var[frame_no] = 0
        dropped_step = frame_no%pixel_step
        ref_drop_steps = [r+dropped_step for r in region_starts]
        for drop_step in ref_drop_steps:
            frame_check[drop_step] = 0
            all_zeroed[drop_step] = 0
            all_zeroed.var[drop_step] = 0
    # Now drop out whole detectors
    for wire_no in dropped_cols:
        b_zeroed[:,wire_no] = 0
        b_zeroed.var[:,wire_no] = 0
        all_zeroed[:,wire_no] = 0
        all_zeroed.var[:,wire_no] = 0
    c = all_zeroed.reshape([b.shape[0]/pixel_step,pixel_step,b.shape[-1]])
    frame_check = frame_check.reshape([b.shape[0]/pixel_step,pixel_step])
    frame_sum = frame_check.intg(axis=1)
    print `b.shape` + "->" + `c.shape`
    print 'Relative no of frames: ' + `frame_sum`
    mult_fact = frame_sum[0]*len(frame_sum)
    print 'Multiplication factor at end: %d' % mult_fact
    # Output the starting data for external use
    if dumpfile is not None:
        dump_wire_intensities(dumpfile,raw=b_zeroed)
    if len(use_gains)==0:   #we have to calculate them
        if c.shape[0] == 1:   #can't be done, there is no overlap
            raise ValueError, "No overlapping scans present"
        if do_sum:
            # sum the individual unoverlapped sections. Reshape is required as the
            # intg function removes the dimension
            d = c.intg(axis=1).reshape([c.shape[0],1,c.shape[2]]) #array of [rangeno,stepno,tubeno]
            # normalise by the number of frames in each section
        else:
            d = c  #no op
        # Note gumpy can't do transposes of more than two axes at once
        e = d.transpose((2,0))  #array of [wireno,stepno,section]
        e = e.transpose((1,2))  #array of [wireno,section,stepno]
        print "Data shape: " + repr(e.shape)
        print "Check shape: " + repr(frame_sum.shape)
        # ignore initial wires completely
        ignore = 0
        for one_wire in range(len(e)):
           if not e[one_wire].any():   #all zero
               #print "Ignoring wire %d" % one_wire
               ignore+=1      #mask it out
           else:
               break
        ignore += len(frame_sum)  #to avoid any contamination
        print "Ignoring all wires within %d of edges" % ignore
        gain,dd,interim_result,residual_map,chisquared,oldesds,first_ave,weights = \
            iterate_data(e[ignore:-ignore],iter_no=iterno,unit_weights=unit_weights,pixel_mask=None)
    else:        #we have been provided with gains
        gain = use_gains
        ignore = fix_ignore
        chisquared=0.0
    # calculate errors based on full dataset
    # First get a full model
    reshape_ds = b_zeroed.reshape([b.shape[0]/pixel_step,pixel_step,b.shape[-1]])
    start_ds = reshape_ds.transpose((2,0))[ignore:-ignore] #array of [wireno,stepno,section]
    start_ds = start_ds.transpose((1,2))
    start_var = start_ds.var
    # Our new pixel mask has to have all of the steps in
    pixel_mask = array.ones_like(start_ds)
    for one_wire in range(len(start_ds)):
        if not start_ds[one_wire].any():   #all zero
            pixel_mask[one_wire] = 0      #mask it out
   # Normalise gains so that average is 1.0
    gain = gain*len(gain)/gain.sum()
    model,wd,model_var,esds = overlap.apply_gain(start_ds,1.0/start_var,gain,
                                            calc_var=True,bad_steps=dropped_frames,pixel_mask=pixel_mask)
    # model and model_var have shape wireno*pixel_step + no_steps (see shift_tube_add_new)
    print 'Have full model and errors at %f' % time.clock()
    # step size could be less than pixel_step if we have a short non-overlap scan
    real_step = pixel_step
    if len(det_steps)< pixel_step:
        real_step = len(det_steps)
        # and we have to prune the output data too
        holeless_model = zeros([real_step*start_ds.shape[0]])
        holeless_var = zeros_like(holeless_model)
        for wire_set in range(start_ds.shape[0]):
            holeless_model[wire_set*real_step:(wire_set+1)*real_step]=model[wire_set*pixel_step:(wire_set+1)*pixel_step]    
            holeless_var[wire_set*real_step:(wire_set+1)*real_step]=model_var[wire_set*pixel_step:(wire_set+1)*pixel_step]    
        model = holeless_model
        model_var = holeless_var
    model *= mult_fact
    model_var = model_var*mult_fact*mult_fact
    cs = Dataset(model)
    cs.var = model_var
    mult_string =  """Intensities were
            multiplied by %d to approximate total intensity measured at each point.""" % mult_fact
    # Now build up the important information
    cs.title = ds.title
    cs.copy_cif_metadata(ds)
    # construct the axes
    axis = arange(len(model))
    new_axis = axis*bin_size + ds.axes[0][0] + ignore*pixel_step*bin_size
    axis_string = """Following application of gain correction, two theta values were 
    calculated assuming a step size of %8.3f starting at %f.""" % (bin_size,ds.axes[0][0]+ignore*pixel_step*bin_size)
    cs.set_axes([new_axis],anames=['Two theta'],aunits=['Degrees'])
    print 'New axis goes from %f to %f in %d steps' % (new_axis[0],new_axis[-1],len(new_axis))
    print 'Total %d points in output data' % len(cs)
    # prepare info for CIF file
    import math
    detno = map(lambda a:"%d" % (a+ignore),range(len(gain)))
    gain_as_strings = map(lambda a:"%.4f" % a,gain)
    gain_esd = ["%.4f" % a for a in esds]
    cs.harvest_metadata("CIF").AddCifItem((
        (("_[local]_wire_number","_[local]_refined_gain","_[local]_refined_gain_esd"),),
        ((detno,gain_as_strings,gain_esd),))
        )
    if len(use_gains)==0:
        info_string = "After vertical integration between pixels %d and %d," % (bottom,top) + \
        """ individual wire gains were iteratively refined using the Ford/Rollett algorithm (Acta Cryst. (1968) B24,293). 
            Final gains are stored in the _[local]_refined_gain loop.""" + mult_string + axis_string
    else:
        info_string =  "After vertical integration between pixels %d and %d," % (bottom,top) + \
        " individual wire gains were corrected based on a previous iterative refinement using the Ford/Rollett algorithm. The gains used" + \
        "are stored in the _[local]_refined_gain loop." + mult_string + axis_string
    cs.add_metadata("_pd_proc_info_data_reduction",info_string,append=True)
    return cs,gain,esds,chisquared,c.shape[0],ignore

# Do an iterative refinement of the gain values. We calculate errors only when chisquared shift is
# small, and aim for a shift/esd of <0.1
def iterate_data(dataset,iter_no=5,pixel_mask=None,plot_clear=True,algo="FordRollett",unit_weights=False):
    """Iteratively refine the gain. The pixel_step is the number of steps
    a wire takes before it overlaps with the next wire. Parameter
    'dataset' is an n x m x s array, for m distinct angular regions of
    s steps each covered by wire number n.  Note that this layout is
    convenient as a reshape([m,s]) operation applied to the flat model
    array will put sequential values together. iter_no is the number
    of iterations, if negative the routine will iterate until
    chisquared does not change by more than 0.01 or abs(iter_no)
    steps, whichever comes first. Pixel_mask has a zero for any wire
    that should be excluded. Algo 'ford rollett' applies the algorithm
    of Ford and Rollet, Acta Cryst. (1968) B24, p293
    """
    import overlap
    start_gain = array.ones(len(dataset))
    if unit_weights is True:
        weights = array.ones_like(dataset)
    else:
        weights = 1.0/dataset.var
    # Use weights as the mask
    if pixel_mask is not None:
        weights = weights*pixel_mask
    if algo == "FordRollett":
        gain,first_ave,ar,esds,k = overlap.find_gain_fr(dataset,weights,start_gain,pixel_mask=pixel_mask)
    else:
        raise ValueError("No such algorithm: %s" % algo)
    chisquared,residual_map = overlap.get_statistics_fr(gain,first_ave,dataset,dataset.var,pixel_mask)
    old_result = first_ave    #store for later
    chisq_history = [chisquared]
    k_history = [k]
    if iter_no > 0: 
        no_iters = iter_no
    else:
        no_iters = abs(iter_no)
    for cycle_no in range(no_iters+1):
        esdflag = (cycle_no == no_iters)  # need esds as well, and flags the last cycle
        print 'Cycle %d' % cycle_no
        if cycle_no > 3 and iter_no < 0:
            esdflag = (esdflag or (abs(chisq_history[-2]-chisq_history[-1]))<0.005)
        if algo == "FordRollett":
            gain,interim_result,ar,esds,k = overlap.find_gain_fr(dataset,weights,gain,arminus1=ar,pixel_mask=pixel_mask,errors=esdflag)
        chisquared,residual_map = overlap.get_statistics_fr(gain,interim_result,dataset,dataset.var,pixel_mask)
        chisq_history.append(chisquared)
        k_history.append(k)
        if esdflag is True:
            break
    print 'Chisquared: ' + `chisq_history`
    print 'K: ' + `k_history`
    print 'Total cycles: %d' % cycle_no
    print 'Maximum shift/error: %f' % max(ar/esds)
    return gain,dataset,interim_result,residual_map,chisq_history,esds,first_ave,weights

def test_iterate_data():
    import random,math
    """Generate a simple dataset and check that it refines to reasonable values"""
    # Dimensions 5 x 2 x 10 i.e. 5 tubes, 20 steps overlapping after 10
    true_gains = rand(20) + 0.5   #20 gains between 0.5 and 1.5
    true_vals = ones([210])*100.0  #background is 100
    peak_1 = make_peak(4)
    true_vals[25:25+len(peak_1)] += 1000*peak_1
    true_vals[95:95+len(peak_1)] += 550*peak_1
    true_vals[143:143+len(peak_1)] += 2000*peak_1
    start_data = ones([20,20])     #20 tubes take 20 steps
    for i in range(len(start_data)):
        start_data[i] = true_vals[i*10:i*10+20]*true_gains[i]
        start_data[i] = [random.gauss(a,math.sqrt(a)) for a in start_data[i]] 
    start_data = Dataset(start_data).reshape([20,2,10])
    print 'Starting array, first 3: ' + repr(start_data.storage[:3])
    g,d,ir,rm,hist,esds,fa,wts = iterate_data(start_data,10,20)
    return g,true_gains,ir,true_vals
    
def load_regain_values(filename):
    """Load a list of gain values derived from a previous call to do_overlap"""
    gain_lines = open(filename,"r").readlines()
    gain_lines = [l.split() for l in gain_lines if len(l)>0 and l[0]!='#'] #remove comments and blanks
    ignored = int(gain_lines[0][1])
    tubes,gain_vals = zip(*[(int(l[0]),float(l[1])) for l in gain_lines[1:]])
    return Array(gain_vals),ignored

def store_regain_values(filename,gain_vals,gain_comment="",ignored=0):
    """Store values calculated from the do_overlap routine"""
    f = open(filename,"w")
    f.write("#Gain values calculated by Wombat reduction routine\n")
    f.write("#"+gain_comment+"\n")
    f.write("Ignored: %d\n" % ignored) 
    for pos,gain in zip(range(len(gain_vals)),gain_vals):
        f.write("%d   %8.3f\n" % (pos,gain))
    f.close()

def dump_gain_file(filename,raw=None,gain=None,model=None,chisq=[],name_prefix=""):
    """Dump data for external gain refinement. We are provided a [sector,steps,tubeno]
    array, where each tubes scan is split into sectors of steps length."""
    import os
    dirname,basename = os.path.split(filename)
    full_filename = os.path.join(dirname,name_prefix+basename)
    outfile = open(full_filename,"w")
    d = raw
    print 'Dumping gains to %s' % full_filename
    stepsize = d.shape[1]
    # Header
    outfile.write("#Wire Sectors Steps\n")
    outfile.write("%d %d %d\n" % (d.shape[-1],d.shape[0],d.shape[1]))
    # raw data
    for wire_no in range(d.shape[-1]):
        for sector_no in range(d.shape[0]):
            for step_no in range(stepsize):
                    outfile.write("%8.3f %7.3f\n" % (d.storage[sector_no,step_no,wire_no],d.var[sector_no,step_no,wire_no]))
    outfile.write("#intensities\n")
    for intensity in model:
        outfile.write("%8.3f\n" % intensity)
    outfile.write("#gains\n")
    for g in gain:
        outfile.write("%8.3f\n" % g)
    outfile.write("#chisq\n")
    for c in chisq:
        outfile.write("%8.3f\n" % c)
    outfile.close()

def parse_ignore_spec(ignore_string):
    """A helper function to parse a string of form a:b,c:d returning
    (a,b),(c,d)"""
    import re
    p = ignore_string.split(',')
    q = map(lambda a:a.split(':'),p)
    # q is now  a sequence of string values
    final = []
    for strval in q:
        try:
            ranges = map(lambda a:int(a.strip()),strval)
            if len(ranges) == 1:
                final.append(ranges[0])
            else:
                final = final + range(ranges[0],ranges[1]+1)
        except ValueError:
            pass
    return set(final)   # a set to avoid duplications
