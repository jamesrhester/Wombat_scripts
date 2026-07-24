# Script control setup area
__script__.title     = 'WOM Texture'
__script__.version   = '0.1'

import sys
# For direct access to the selected filenames
__datasource__ = __register__.getDataSourceViewer()

''' User Interface '''

# Output Folder
out_folder = Par('file')
out_folder.dtype = 'folder'
out_folder.title = 'Directory'
output_stem = Par('string','reduced')
output_stem.title = 'Append to filename'
output_restrict = Par('string','All')
output_restrict.title = 'These frames only (n:m)'
Group('Output File').add(output_stem, out_folder, output_restrict)

# Normalization
# We link the normalisation sources to actual dataset locations right here, right now
norm_table = {'Monitor 1':'bm1_counts','Monitor 2':'bm2_counts'
              ,'Detector time':'detector_time'}
norm_apply     = Par('bool', 'True')
norm_apply.title = 'Apply'
norm_reference = Par('string', 'Monitor 3', options = norm_table.keys())
norm_reference.title = 'Source'
norm_target = Par('int',-1)
norm_target.title = 'Normalise output datasets to (-1 for none):'
norm_plot = Act('plot_norm_proc()','Plot')
norm_plot_all = Act('plot_all_norm_proc()','Plot all')
Group('Normalization').add(norm_apply, norm_reference,norm_target,norm_plot_all,norm_plot)

# Efficiency Correction
eff_apply = Par('bool', 'True')
eff_apply.title = 'Apply'
eff_map   = Par('file', '')
eff_map.ext = '*.*'
eff_map.title = 'Efficiency File'
eff_show  = Act('eff_show_proc()', 'Show') 
Group('Efficiency Correction').add(eff_apply, eff_map, eff_show)

# Vertical Integration
vig_central_only = Par('bool', 'False')
vig_central_only.title = 'Centre only?'
Group('Vertical Integration').add(vig_central_only)

# ===== End of reduction settings ===== #

# Plot Helper
plh_from = Par('string', 'Plot 2', options = ['Plot 1', 'Plot 2', 'Plot 3'])
plh_from.title = 'From'
plh_to   = Par('string', 'Plot 3', options = ['Plot 1', 'Plot 2', 'Plot 3'])
plh_to.title= 'To'
plh_copy = Act('plh_copy_proc()', 'Copy')
Group('Copy 1D Datasets').add(plh_from, plh_to, plh_copy)

plh_plot    = Par('string', '', options = ['Plot 1', 'Plot 2', 'Plot 3'], command = 'plh_plot_changed()')
plh_plot.title = 'Plot Name'
plh_dataset = Par('string', '', options = ['All'])
plh_dataset.title = 'Dataset'
plh_delete  = Act('plh_delete_proc()', 'Delete')
Group('Delete 1D Datasets').add(plh_plot, plh_dataset, plh_delete)

''' Load Preferences '''

efficiency_file_uri     = __UI__.getPreference("au.gov.ansto.bragg.wombat.ui:efficiency_file_uri")
normalisation_reference = __UI__.getPreference("au.gov.ansto.bragg.wombat.ui:normalisation_reference")
user_output_dir         = __UI__.getPreference("au.gov.ansto.bragg.wombat.ui:user_output_dir")
#
# Set the optional values to those in the preferences file
#
if user_output_dir:
    out_folder.value = user_output_dir
if normalisation_reference:  #saved as location, need label instead
        vals = filter(lambda a:a[1]==normalisation_reference,norm_table.items())
        if vals: norm_reference.value = vals[0]
if efficiency_file_uri:
    eff_map.value = efficiency_file_uri

''' Button Actions '''

# Plot normalisation info
def plot_norm_proc():
    plot_norm_master()

def plot_all_norm_proc():
    """Plot all normalisation values found in file"""
    plot_norm_master(all_mons=True)

def plot_norm_master(all_mons = False):
    dss = __datasource__.getSelectedDatasets()
    if Plot2.ds:
        remove_list = copy(Plot2.ds)  #otherwise dynamically changes
        for ds in remove_list:
            Plot2.remove_dataset(ds)  #clear doesn't work
    for fn in dss:
        loc = fn.getLocation()
        dset = df[str(loc)]
        print 'Dataset %s' % os.path.basename(str(loc))
        for monitor_loc in norm_table.keys():
            if all_mons or monitor_loc == str(norm_reference.value):
                norm_source = norm_table[monitor_loc]
                plot_data = Dataset(getattr(dset,norm_source))
                if norm_apply.value or all_mons:
                    ave_val = plot_data.sum()/len(plot_data)
                    plot_data = plot_data/ave_val
                plot_data.title = os.path.basename(str(loc))+':' + str(monitor_loc) + '_'
                send_to_plot(plot_data,Plot2,add=True)
        
def show_helper(filename, plot, pre_title = ''):
    if filename:
        
        ds = Dataset(str(filename))
        
        if ds.ndim == 4:
            plot.set_dataset(ds[0, 0])
            plot.title = ds.title + " (first frame)"
        elif ds.ndim == 3:
            plot.set_dataset(ds[0])
            plot.title = ds.title + " (first frame)"
        else:
            plot.set_dataset(ds)
        
        if pre_title:
            plot.title = pre_title + plot.title
            
    else:
        print 'no valid filename was specified'

# show Efficiency Correction Map 
def eff_show_proc():
    global Plot1
    show_helper(eff_map.value, Plot1, "Efficiency Map:")

def plh_copy_proc():
    
    src = str(plh_from.value)
    dst = str(plh_to.value)
    
    plots = {'Plot 1': Plot1, 'Plot 2': Plot2, 'Plot 3': Plot3}

    if not src in plots:
        print 'specify source plot'
        return
    if not dst in plots:
        print 'specify target plot'
        return
    if src == dst:
        print 'specify a different target plot'
        return
        
    src_plot = plots[src]
    dst_plot = plots[dst]
    
    src_ds = src_plot.ds
    if type(src_ds) is not list:
        print 'source plot does not contain 1D datasets'
        return
    
    dst_ds = dst_plot.ds
    if type(dst_ds) is not list:
        dst_plot.clear()
        dst_ds = []
    
    dst_ds_ids = [id(ds) for ds in dst_ds]
    
    for ds in src_ds:
        if id(ds) not in dst_ds_ids:
            dst_plot.add_dataset(ds)

def plh_plot_changed():
    
    target = str(plh_plot.value)
    
    plots = {'Plot 1': Plot1, 'Plot 2': Plot2, 'Plot 3': Plot3}
    
    if not target in plots:
        print 'specify source plot'
        plh_dataset.options = []
        return
    
    target_plot = plots[target]
    target_ds   = target_plot.ds
    target_list = ['All']
    
    if (type(target_ds) is not list) or (len(target_ds) == 0):
        print 'target plot does not contain 1D datasets'
        plh_dataset.options = []
        return
    
    for ds in target_ds:
        target_list.append(ds.title)
    
    plh_dataset.options = target_list
    plh_dataset.value   = 'All'

def plh_delete_proc():
    
    target  = str(plh_plot.value)
    dataset = str(plh_dataset.value)
    
    plots = {'Plot 1': Plot1, 'Plot 2': Plot2, 'Plot 3': Plot3}
    
    if not target in plots:
        print 'specify source plot'
        plh_dataset.options = []
        return
    
    target_plot = plots[target]
    target_ds   = copy(target_plot.ds)
    
    if (type(target_ds) is not list) or (len(target_ds) == 0):
        print 'target plot does not contain 1D datasets'
        plh_dataset.options = []
        return
    
    if dataset == 'All':
        for ds in target_ds:
            target_plot.remove_dataset(ds)
    else:
        for ds in target_ds:
            if ds.title == dataset:
                target_plot.remove_dataset(ds)

# The preference system: 
def load_user_prefs(prefix = ''):
    """Load preferences, optionally prepending the value of
    prefix in the preference search.  This is typically used
    to load an alternative set of preferences"""
    # Run through our parameters, looking for the corresponding
    # preferences
    g = globals()
    p = g.keys()
    for name in p:
        if hasattr(g[name], 'value'):
            try:
                setattr(g[name], 'value', get_prof_value(prefix+name))
            except:
                print 'Failure setting %s to %s' % (name,str(get_prof_value(prefix+name)))
            print 'Set %s to %s' % (name, str(globals()[name].value))

def save_user_prefs(prefix=''):
    """Save user preferences, optionally prepending the value of
    prefix to the preferences. This prefix is typically used to
    save an alternative set of preferences.  Return lists of values
    as ASCII strings for logging purposes"""
    print 'In save user prefs'
    prof_names = []
    prof_vals = []
    # sneaky way to get all the preferences
    g = globals()
    p = g.keys()
    for name in p:
        if hasattr(g[name], 'value') and name[0] != '_':
            prof_val = getattr(g[name], 'value')
            set_prof_value(prefix+name,str(prof_val))
            print 'Set %s to %s' % (prefix+name,get_prof_value(prefix+name))
            prof_names.append(name)
            prof_vals.append(str(prof_val))
    return prof_names, prof_vals

""" Helper routines for run_script actions """
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

        return rs

def process_normalise_options():
    
    if norm_apply.value:
        
        # norm_ref is the source of information for normalisation
        # norm_tar is the value norm_ref should become,
        # by multiplication.  If 'auto', the maximum value of norm_ref
        # for the first dataset is used, otherwise any number may be entered.
        
        if len(str(norm_target.value))==0:
            norm_tar = -1
        else: 
            norm_tar = int(str(norm_target.value))
            
        # check if normalization reference is provided

        norm_ref = str(norm_reference.value)
        if len(norm_ref) == 0:
            norm_ref = None
            norm_tar = -1
            print 'WARNING: no reference for normalization was specified'
        else:     
            print 'utilized reference value for "' + norm_ref + '" is:', norm_tar
            
        # use provided reference value
        
        if norm_tar != -1:
            norm_tar = float(norm_tar)
            
    else:
        norm_ref = None
        norm_tar = None

    return norm_tar, norm_ref

def process_eff_options():
    if eff_apply.value:
        if not eff_map.value:
            eff = None
            print 'WARNING: no eff-map was specified'
        else:
            try:
                eff = Dataset(str(eff_map.value))
            except:
                raise ValueError, "Efficiency file %s not found" % str(eff_map.value)
    else:
        eff = None

    return eff

def get_detector_positions(ds):

    try:
        stth_value = sum(ds.stth)/len(ds.stth) # save for later
        all_stth = ds.stth[:] # also save for later
    except TypeError:
        stth_value = ds.stth
        all_stth= [stth_value]

    print "all_stth is %s" % repr(all_stth) 

    return all_stth

def get_euler_positions(ds):

    phi_locs = ("/entry1/sample/euler_phi", "/entry1/sample/phi")
    chi_locs = ("/entry1/sample/euler_chi", "/entry1/sample/chi")
    omega_locs = ("/entry1/sample/euler_omega", "/entry1/sample/omega")
    
    ret_vals = []
    for angle in omega_locs, chi_locs, phi_locs:
        result = None
        for locs in angle:
            try:
                result = ds[locs][:]
                break
            except:
                pass
        if result == None:
            print "Failed to find %s" % repr(angle)
            result = ["None"] * len(ds)
        print "Have a result for %s: %s" % (locs, result)
        ret_vals.append(result)

    return ret_vals

def get_frame_range(ds):
    """
    A single a:b restriction on frames to use
    """
    restrict_spec = str(output_restrict.value)
    if ':' in restrict_spec:
        current_frame_start,end_frame = map(int,restrict_spec.split(':'))
    else:
        end_frame = len(ds)
        current_frame_start = 0
    return current_frame_start, end_frame

def create_stem_template(ds, df, fn, frame_no):
    """Calculate the filename string by looking for special wildcards.
    %s = samplename, %t1,%t2, %vf for cryo and vacuum furnace temperatures.
    If frame_no is >= 0, a particular step will be chosen
    """
    import re
    temp_table = {"%t1":("/entry1/sample/tc1/sensor/sensorValueA","%.1fK"),
                  "%t2":("/entry1/sample/tc1/sensor/sensorValueB","%.1fK"),
                  "%vf":("/entry1/sample/tc1/sensor","%.0fC"),
                  "%phi":("/entry1/sample/phi","%04.0fd")
                  }
    stem_template = str(output_stem.value)
    stem_template = re.sub(r'[^\w+=()*^@~:{}\[\].%-]','_',stem_template)
    
    if '%s' in stem_template:
        samplename = ds.harvest_metadata("CIF")['_pd_spec_special_details']
        name_front = re.sub(r'[^\w+=()*^@~:{}\[\].%-]','_',samplename)
        stem_template = stem_template.replace('%s',name_front)

    for wildcard in temp_table.keys():
        if wildcard in stem_template:
            loc, fmt = temp_table[wildcard]

            try:
                temperature = df[fn][loc]
            except AttributeError:
                print "Unable to determine temperature for %s" % wildcard[1:]
                continue
            else:
                print `temperature`
                if frame_no >= 0 and hasattr(temperature,"__len__") and len(temperature) > 0:
                    temperature = temperature[frame_no]
                elif hasattr(temperature,"__len__"):
                    temperature = sum(temperature)/len(temperature)
            stem_template = stem_template.replace(wildcard, fmt % temperature)
 
    if stem_template != "": stem_template = "_"+stem_template
    print 'Filename stem is now ' + stem_template
    return stem_template

def process_straighten(cs, stth, bottom, top):
    from Reduction import straightening

    radius = float(cs.harvest_metadata("CIF")["_pd_instr_dist_spec/detc"])

    start_angles = stth
    print "Total length %d" % len(cs.axes[-1])
    if len(cs.axes[-1]) == cs.shape[-1] + 1:
        print "First bin is %f to %f" % (cs.axes[-1][0], cs.axes[-1][1])
        wires = getCenters(cs.axes[-1])
    else:
        print "First wire is at %f" % cs.axes[-1][0]
        wires = cs.axes[-1]
    print "First wire at offset %f, stth angles %s" % (wires[0], `start_angles`)
    vert_size = len(cs.axes[-2]) - 1
    vert_pos = getCenters(cs.axes[-2]) - cs.axes[-2][vert_size/2]
    vert_pos.title = "Vertical offset"
    wires.title = cs.axes[-1].title
    new_ds, new_contribs = straightening.correctGeometryjv(cs, radius, start_angles, wires, vert_pos,
                                                         bottom, top)
    # Add metadata record

    info_string = """Geometry was corrected by dividing pixel intensity and variance between ideal true two-theta bins based on deviation from ideal bin centre."""
    new_ds.add_metadata('_pd_proc_info_data_reduction', info_string, append=True)
    return new_ds, new_contribs
    
def process_vertical_sum(cs, stth_values, segment, contribs=None):
    """
    The detector is divided into three regions: top, middle, bottom. The value of
    segment determines which of these is summed.
    """
    from Reduction import reduction

    # Work out region
    
    regions = {"top":(86,127), "middle":(42,85), "bottom":(0,41)}
    bottom, top = regions[segment]

    # fix the axes

    print "Stth values: " + `stth_values`
    
    cs.set_axes([stth_values,cs.axes[1],cs.axes[2]],anames=["Texture step",
                                                         "Vertical Pixel",
                                                         "Two theta"],
                aunits=["None","mm","Degrees"])

    print 'cs axes: ' + cs.axes[0].title + ' ' + cs.axes[1].title + ' ' + cs.axes[2].title
    reduction.boundaries_to_mdpts(cs)
    print 'cs axes: ' + cs.axes[0].title + cs.axes[1].title
    Plot1.set_dataset(cs)
    print `cs.__dict__['ms']`
    gs = reduction.getVerticalIntegrated(cs, axis=1, bottom = bottom, top= top)
    return gs
    
''' Script Actions '''

# This function is called when pushing the Run button in the control UI.
def __run_script__(fns):
    
    from Reduction import reduction, AddCifMetadata
    from os.path import basename
    from os.path import join
    from Formats import output
    
    df.datasets.clear()
    
    # save user preferences

    prof_names,prof_values = save_user_prefs()

    elapsed = time.clock()

    # check input
    
    if (fns is None or len(fns) == 0) :
        print 'no input datasets'
        return

    # Get processing parameters
    
    norm_tar, norm_ref = process_normalise_options()

    # The error dialog only works at this level it seems, so anything that
    # could raise an error we group together here
    
    try:
        eff = process_eff_options()
    except ValueError as e:
        open_error(str(e))
        return

    # iterate through input datasets
    # note that the normalisation target (an arbitrary number) is set by
    # the first dataset unless it has already been specified.
    
    for fn in fns:

        Plot1.clear()   #otherwise we get an error
        
        # load dataset

        ds = df[fn]
                       
        # extract and store basic metadata

        ds = reduction.AddCifMetadata.extract_metadata(ds)
        reduction.AddCifMetadata.store_reduction_preferences(ds, prof_names, prof_values)

        # Get detector positions

        all_stth = get_detector_positions(ds)
        all_omega, all_chi, all_phi = get_euler_positions(ds)
 
        # Prepare dataset
        
        if ds.ndim > 3:
            rs = ds.get_reduced()
        else:
            rs = ds

        rs = rs * 1.0  #convert to float
        rs.copy_cif_metadata(ds)
        
        # Do normalisation

        if norm_ref:
            norm_tar = reduction.applyNormalization(rs, reference=norm_table[norm_ref], target=norm_tar)

        # check if efficiency correction is required
        
        assert rs.dtype == Array([1.2,1.3]).dtype

        if eff:
            ds = reduction.getEfficiencyCorrected(rs, eff)
        else:
            ds = rs

        # Keep axis information
        
        ds.axes = rs.axes
        
        # Calculate filename string

        stem_template = create_stem_template(ds, df, fn, -1)
        
        # restrict output set of frames

        current_frame_start, end_frame = get_frame_range(ds)

        # Now select these frames

        cs = ds[current_frame_start:end_frame]
        
        # We have a simple approach: each frame is a separate setting of
        # chi and/or phi, so there is no need to group
            
        cs.copy_cif_metadata(ds)
        stth_values = all_stth[current_frame_start:end_frame]
        all_omega = all_omega[current_frame_start:end_frame]
        all_chi = all_chi[current_frame_start:end_frame]
        all_phi = all_phi[current_frame_start:end_frame]

        # Extract temperature if requested

        stem_template = create_stem_template(ds, df, fn, -1)
        contribs = None
            
        # Always straighten
            
        cs, contribs = process_straighten(cs, stth_values, 1, 126)
            
        print 'Finished straightening'

        # Apply 2th offset
        print 'offset: ' + `stth_values[0]`
        cs.axes[-1] += stth_values[0]

        # Vertical summation

        data_parts = {}

        for region in ("bottom", "middle", "top"):

            gs = process_vertical_sum(cs, stth_values, region, contribs = contribs)

            data_parts[region] = gs
            print "Cs.shape %s, Gs shape: %s" % (cs.shape, gs.shape)
            if region == "middle":
                send_to_plot(gs,Plot1,add=False,title="Integrated data",quantity="Counts")
        for frameno in range(current_frame_start,end_frame):

            # get angles
            this_chi = all_chi[frameno]
            this_phi = all_phi[frameno]
            this_om = all_omega[frameno]
            
            #output.write_cif_data(gs,filename_base)

            this_frame = (data_parts["bottom"][frameno], data_parts["middle"][frameno],
                          data_parts["top"][frameno])
            for i,r in enumerate(("bottom", "middle", "top"),):
                this_frame[i].copy_cif_metadata(data_parts[r])
            eta = (-5.25, 0.0, 5.25)
            # Output datasets
            filename_base = join(str(out_folder.value), basename(str(fn))[:-7] + stem_template + "_chi" + `round(this_chi)` + "_phi" + `round(this_phi)`)

            print "Writing %s" % filename_base
            output.write_esg_data(this_chi, this_phi, this_om, this_frame, eta, filename_base)
            
            
''' Utility functions for plots '''
def send_to_plot(dataset,plot,add=False,title="",add_timestamp=True,quantity=""):
    """This routine appends a timestamp to the dataset title
    in order to keep uniqueness of the title for later 
    identification purposes. It also maintains plot
    consistency in terms of displaying d-spacing."""
    from datetime import datetime
    from Reduction import reduction
    if add_timestamp:
        timestamp = datetime.now().strftime("%H:%M:%S")
        dataset.title = dataset.title + timestamp
    # Check d-spacing status
    if plot.ds:
        if plot.ds[0].axes[0].name == 'd-spacing':
            reduction.convert_to_dspacing(dataset)
        elif plot.ds[0].axes[0].name == 'Two theta':
            reduction.convert_to_twotheta(dataset)
    if add:
        plot.add_dataset(dataset)
    else:
        plot.set_dataset(dataset)
    if title:
        plot.title = title
    #Vertical axis
    plot.set_y_label(quantity)
    #Update any widgets that keep a track of the plots
    plh_plot_changed()

# dispose
def __dispose__():
    global Plot1,Plot2,Plot3
    Plot1.clear()
    Plot2.clear()
    Plot3.clear()

''' Quick-Fix '''

def run_action(act):
    act.set_running_status()
    try:
        exec(act.command)
        act.set_done_status()
    except:
        act.set_error_status()
        traceback.print_exc(file = sys.stdout)
        raise Exception, 'Error in running <' + act.text + '>'
    
''' Execute this each time it is loaded to reload user preferences '''

load_user_prefs()

