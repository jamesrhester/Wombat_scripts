# Script control setup area
__script__.title     = 'WOM Reduction'
__script__.version   = '1.1'

import sys
# For direct access to the selected filenames
__datasource__ = __register__.getDataSourceViewer()

''' User Interface '''

# Output Folder
out_folder = Par('file')
out_folder.dtype = 'folder'
out_folder.title = 'Directory'
output_xyd = Par('bool','True')
output_xyd.title = 'XYD'
output_cif = Par('bool','True')
output_cif.title = 'pdCIF'
output_fxye = Par('bool','False')
output_fxye.title = 'GSAS FXYE'
output_topas = Par('bool','False')
output_topas.title = 'Topas'
output_stem = Par('string','reduced')
output_stem.title = 'Append to filename'
grouping_options = {"Frame":"run_number","TC1 setpoint":"/sample/tc1/sensor/setpoint1","None":None}
output_grouping = Par('string','None',options=grouping_options.keys())
output_grouping.title = 'Split frames'
output_restrict = Par('string','All')
output_restrict.title = 'These frames only (n:m)'
Group('Output File').add(output_xyd,output_cif,output_fxye,output_topas,output_stem,out_folder,
                           output_grouping,output_restrict)

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

# Background Correction
bkg_apply = Par('bool', 'False')
bkg_apply.title = 'Apply'
bkg_map   = Par('file', '')
bkg_map.ext = '*.hdf'
bkg_map.title = 'Background file'
bkg_show  = Act('bkg_show_proc()', 'Show') 
Group('Background Correction').add(bkg_apply, bkg_map, bkg_show)

# Efficiency Correction
eff_apply = Par('bool', 'True')
eff_apply.title = 'Apply'
eff_map   = Par('file', '')
eff_map.ext = '*.*'
eff_map.title = 'Efficiency File'
eff_show  = Act('eff_show_proc()', 'Show') 
Group('Efficiency Correction').add(eff_apply, eff_map, eff_show)

# Vertical Integration
vig_lower_boundary = Par('int', '0')
vig_lower_boundary.title = 'Lower limit'
vig_upper_boundary = Par('int', '127')
vig_upper_boundary.title = 'Upper limit'
vig_apply_rescale  = Par('bool', 'True')
vig_apply_rescale.title = 'Rescale'
vig_rescale_target = Par('float', '10000.0')
vig_rescale_target.title = 'Rescale to:'
Group('Vertical Integration').add(vig_lower_boundary, vig_upper_boundary, vig_apply_rescale, vig_rescale_target)

# Recalculate gain
regain_apply = Par('bool','False')
regain_apply.title = 'Apply'
regain_iterno = Par('int','5')
regain_iterno.title = 'Iterations'
regain_store = Par('bool','False')
regain_store.title = 'Store gain result'
regain_store_filename = Par('file')
regain_store_filename.title = 'Store in file:'
regain_load = Par('bool','False')
regain_load.title = 'Load gain from file'
regain_load_filename = Par('file')
regain_load_filename.title = 'Gain file'
#regain_dump_tubes = Par('bool','False')
#regain_dump_tubes.title = 'Dump values by tube'
regain_sum = Par('bool','False')
regain_sum.title = 'Sum before refinement'
Group('Recalculate Gain').add(regain_apply,regain_iterno,regain_store,regain_store_filename,
                              regain_load,regain_load_filename,regain_sum)

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

# Plot settings
ps_plotname = Par('string','Plot 2',options=['Plot 2','Plot 3'])
ps_plotname.title = 'Plot Name'
ps_dspacing = Par('bool',False,command='dspacing_change()')
ps_dspacing.title = 'd spacing'
Group('Plot settings').add(ps_plotname,ps_dspacing)

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

''' New global plots'''
if 'Plot4' not in globals() or 'Disposed' in str(Plot4.pv):
    Plot4 = Plot(title='Chi-squared behaviour')
    Plot4.close = noclose
if 'Plot5' not in globals() or 'Disposed' in str(Plot5.pv):
    Plot5 = Plot(title='Final gain values')
    Plot5.close = noclose

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
            
# show Background Correction Map
def bkg_show_proc():
    global Plot1
    show_helper(bkg_map.value, Plot1, "Background Map: ")

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

def dspacing_change():
    """Toggle the display of d spacing on the horizontal axis"""
    global Plot2,Plot3
    from Reduction import reduction
    plot_table = {'Plot 2':Plot2, 'Plot 3':Plot3}
    target_plot = plot_table[str(ps_plotname.value)]
    if target_plot.ds is None:
        return
    # Preliminary check we are not displaying something
    # irrelevant, e.g. monitor counts
    for ds in target_plot.ds:
        if ds.axes[0].name not in ['Two theta','d-spacing']:
            return
    change_dss = copy(target_plot.ds)
    # Check to see what change is required
    need_d_spacing = ps_dspacing.value
    # target_plot.clear() causes problems; use 'remove' instead
    # need to set the xlabel by hand due to gplot bug
    if need_d_spacing: target_plot.x_label = 'd-spacing (Angstroms)'
    elif not need_d_spacing: target_plot.x_label = 'Two theta (Degrees)'
    target_plot.y_label = 'Intensity'
    for ds in change_dss:
        current_axis = ds.axes[0].name
        print '%s has axis %s' % (ds.title,current_axis)
        if need_d_spacing:    
            result = reduction.convert_to_dspacing(ds)
        elif not need_d_spacing:
            result = reduction.convert_to_twotheta(ds)
        if result == 'Changed':
            target_plot.remove_dataset(ds)
            target_plot.add_dataset(ds)

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
        if eval('isinstance('+ name + ',Par)'):
            try:
               setattr(g[name],'value',get_prof_value(prefix+name))
            except:
                print 'Failure setting %s to %s' % (name,str(get_prof_value(prefix+name)))
            print 'Set %s to %s' % (name,str(eval(name+'.value')))

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
        if eval('isinstance('+ name + ',Par)'):
            prof_val = getattr(g[name], 'value')
            set_prof_value(prefix+name,str(prof_val))
            print 'Set %s to %s' % (prefix+name,get_prof_value(prefix+name))
            prof_names.append(name)
            prof_vals.append(str(prof_val))
    return prof_names,prof_vals

""" Helper routines for run_script actions """

def process_normalise_options():
    
    normalise_output = False
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
            location = norm_table[norm_ref]     
            print 'utilized reference value for "' + norm_ref + '" is:', norm_tar
            
        # use provided reference value
        
        if norm_tar != -1:
            norm_tar = float(norm_tar)
            normalise_output = True
            
    else:
        norm_ref = None
        norm_tar = None

    return norm_tar, norm_ref, normalise_output, location

def process_bkg_options():
    if bkg_apply.value:
        if not bkg_map.value:
            bkg = None
            print 'WARNING: no bkg-map was specified'
        else:
            try:
                bkg = Dataset(str(bkg_map.value))
            except:
                raise ValueError, "Background file %s not found" % str(bkg_map.value)
            
            # to avoid complaints in routines that expect it
            reduction.AddCifMetadata.add_metadata_methods(bkg)
    else:
        bkg = None

    return bkg

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

def process_rescale_options():

    if vig_apply_rescale.value:
        vig_normalisation = float(vig_rescale_target.value)
    else:
        vig_normalisation = -1
    group_val = grouping_options[str(output_grouping.value)]

    return vig_normalisation

def process_regain_options():
    regain_data = []
    pre_ignore = 0
    if regain_load.value:
        if not regain_load_filename.value:
            raise ValueError, "You have requested loading of gain correction from a file but no file has been specified"
        rlf = str(regain_load_filename.value)
        regain_data, pre_ignore = reduction.load_regain_values(rlf)
        print "Loaded gain data from %s, first/last %d wires ignored" % (rlf,pre_ignore)

    return regain_data, pre_ignore

''' Script Actions '''

# This function is called when pushing the Run button in the control UI.
def __run_script__(fns):
    
    from Reduction import reduction, AddCifMetadata
    from os.path import basename
    from os.path import join
    from Formats import output
    import re
    
    df.datasets.clear()
    
    # save user preferences

    prof_names,prof_values = save_user_prefs()

    elapsed = time.clock()

    # check input
    
    if (fns is None or len(fns) == 0) :
        print 'no input datasets'
        return

    # set the title for Plot2
    # Plot2.title = 'Plot 2'

    norm_tar, norm_ref, normalise_output, location = process_normalise_options()
    vig_normalisation = process_rescale_options()
    group_val = grouping_options[str(output_grouping.value)]

    # The error dialog only works at this level it seems, so anything that
    # could raise an error we group together here
    
    try:
        bkg = process_bkg_options()
        eff = process_eff_options()
        regain_data, pre_ignore = process_regain_options()
    except ValueError as e:
        open_error(str(e))
        return
    
    # iterate through input datasets
    # note that the normalisation target (an arbitrary number) is set by
    # the first dataset unless it has already been specified.
    
    for fn in fns:
        # load dataset
        ds = df[fn]
        # extract basic metadata
        ds = reduction.AddCifMetadata.extract_metadata(ds)
        reduction.AddCifMetadata.store_reduction_preferences(ds,prof_names,prof_values)
        try:
             stth_value = sum(ds.stth)/len(ds.stth) # save for later
             all_stth = ds.stth[:] # also save for later
        except TypeError:
            stth_value = ds.stth
            all_stth= [stth_value]
        if ds.ndim > 3:
            rs = ds.get_reduced()
        else:
            rs = ds
        print "all_stth is %s" % repr(all_stth) 
        rs = rs * 1.0  #convert to float
        rs.copy_cif_metadata(ds)
        # check if normalized is required 
        if norm_ref:
            norm_tar = reduction.applyNormalization(rs, reference=norm_table[norm_ref], target=norm_tar)
        if bkg:
            rs = reduction.getBackgroundCorrected(rs, bkg, norm_table[norm_ref], norm_tar)
        # check if efficiency correction is required
        assert rs.dtype == Array([1.2,1.3]).dtype
        if eff:
            ds = reduction.getEfficiencyCorrected(rs, eff)
        else:
            ds = rs
        # Calculate inserted string: %s for sample name, %t for temperature
        stem_template = str(output_stem.value)
        stem_template = re.sub(r'[^\w+=()*^@~:{}\[\].%-]','_',stem_template)
        if '%s' in stem_template:
             samplename = ds.harvest_metadata("CIF")['_pd_spec_special_details']
             name_front = re.sub(r'[^\w+=()*^@~:{}\[\].%-]','_',samplename)
             stem_template = stem_template.replace('%s',name_front)
        if '%t1' in stem_template and group_val is None:
             # get tc1
             temperature = df[fn]["/entry1/sample/tc1/sensor/sensorValueA"]
             print `temperature`
             try:
                 avetemp = sum(temperature)/len(temperature)
             except TypeError:
                 avetemp = temperature
             stem_template = stem_template.replace('%t1',"%.0fK" % avetemp)
        if '%t2' in stem_template and group_val is None:
             # get tc1
             temperature = df[fn]["/entry1/sample/tc1/sensor/sensorValueB"]
             print `temperature`
             try:
                 avetemp = sum(temperature)/len(temperature)
             except TypeError:
                 avetemp = temperature
             stem_template = stem_template.replace('%t2',"%.1fK" % avetemp)
        if '%vf' in stem_template and group_val is None:
             # get tc1
             temperature = df[fn]["/entry1/sample/tc1/sensor"]
             print `temperature`
             try:
                 avetemp = sum(temperature)/len(temperature)
             except TypeError:
                 avetemp = temperature
             stem_template = stem_template.replace('%vf',"%.0fC" % avetemp)
        if stem_template != "": stem_template = "_"+stem_template
        print 'Filename stem is now ' + stem_template
        # restrict output set of frames
        restrict_spec = str(output_restrict.value)
        if ':' in restrict_spec:
            first,last = map(int,restrict_spec.split(':'))
            start_frames = last
            current_frame_start = first
            frameno = first
        else:
            start_frames = len(ds)
            current_frame_start = 0
            frameno = 0
        # perform grouping of sequential input frames 
        # we accumulate the equivalent total monitor 
        # counts for requested normalisation later  
        while frameno <= start_frames:
            if regain_apply.value:   #take them all
                frameno = start_frames
                target_val = ""
            elif group_val is None:
                target_val = ""
                final_frame = start_frames-1
                frameno = start_frames
            else:
                stth_value = all_stth[current_frame_start]
                target_val = ds[group_val][current_frame_start]
                try:
                    if df[fn][group_val][frameno] == target_val:
                        frameno += 1
                        continue
                except:   #Assume an exception is due to too large frameno
                    print 'Exiting frame loop due to error'
            # frameno is the first frame with the wrong values
            cs = ds.get_section([current_frame_start,0,0],[frameno-current_frame_start,ds.shape[1],ds.shape[2]])
            cs.copy_cif_metadata(ds)
            # Store temperature if requested
            if "%t1" in stem_template:
                try:
                    temperature = df[fn]["/entry1/sample/tc1/sensor/sensorValueA"][current_frame_start]
                except AttributeError:
                    print 'Unable to determine temperature at step %d' % current_frame_start
                else:
                    stem_template = stem_template.replace('%t1',"%.0fK" % temperature)
                    print 'Filename is now ' + stem_template    
            elif "%vf" in stem_template:
                try:
                    temperature = df[fn]["/entry1/sample/tc1/sensor"][current_frame_start]
                except (AttributeError,TypeError):
                    print 'Unable to determine temperature at step %d' % current_frame_start
                else:
                    stem_template = stem_template.replace('%vf',"%.0fC" % temperature)
                    print 'Filename is now ' + stem_template    

            # check if we are recalculating gain 
            if regain_apply.value:
                # fix the axes
                cs.set_axes([all_stth,cs.axes[1],cs.axes[2]],anames=["Azimuthal angle",
                                                                    "Vertical Pixel",
                                                                    "Two theta"],
                            aunits=["Degrees","mm","Degrees"])
                bottom = int(vig_lower_boundary.value)
                top = int(vig_upper_boundary.value)
                dumpfile = None
                #if regain_dump_tubes.value:
                #    dumpfile = filename_base+".tubes"
                gs,gain,esds,chisquared,no_overlaps,ignored = reduction.do_overlap(cs,regain_iterno.value,bottom=bottom,top=top,
                                                                          drop_frames="",drop_wires="0:6",use_gains=regain_data,dumpfile=dumpfile,
                                                                                   do_sum=regain_sum.value, fix_ignore=pre_ignore)
                if gs is not None:
                    print 'Have new gains at %f' % (time.clock() - elapsed)
                    fg = Dataset(gain)
                    fg.var = esds**2
                    # set horizontal axis (ideal values)
                    Plot1.set_dataset(reduction.getStepSummed(cs))
                    Plot4.set_dataset(Dataset(chisquared))   #chisquared history
                    Plot5.set_dataset(fg)   #final gain plot
                    # now save the file if requested
                    if regain_store.value and not regain_load.value:
                        gain_comment = "Gains refined from file %s" % fn
                        reduction.store_regain_values(str(regain_store_filename.value),gain,gain_comment,
                                                      ignored=ignored)
                else:
                    open_error("Cannot do gain recalculation as the scan ranges do not overlap.")
                    return
            if not regain_apply.value:  #already done
                print 'Summing frames from %d to %d, shape %s, start 2th %f' % (current_frame_start,frameno-1,cs.shape,stth_value)
                if target_val != "":
                    print 'Corresponding to a target value of ' + `target_val`
                    # sum the input frames
                print 'cs axes: ' + cs.axes[0].title + ' ' + cs.axes[1].title + ' ' + cs.axes[2].title
                # es = cs.intg(axis=0)
                es = reduction.getSummed(cs,applyStth=stth_value)  # does axis correction as well
                es.copy_cif_metadata(cs)
                print 'es axes: ' + `es.axes[0].title` + es.axes[1].title
                Plot1.set_dataset(es)

                gs = reduction.getVerticalIntegrated(es, axis=0, normalization=vig_normalisation,
                                                     bottom = int(vig_lower_boundary.value),
                                                     top=int(vig_upper_boundary.value))
            if target_val != "":
                gs.title = gs.title + "_" + str(target_val)
            try:
                send_to_plot(gs,Plot2,add=True,title="Integrated data",quantity="Counts")
            except IndexError:  #catch error from GPlot ??
                send_to_plot(gs,Plot2,add=False,title="Integrated data",quantity="Counts")
            # Output datasets
            try:
                val_for_output = int(target_val)
                format_for_output = "%03d"
            except:
                val_for_output = str(target_val)
                format_for_output = "%s"
            if target_val != "": format_for_output = "_" + format_for_output
            filename_base = join(str(out_folder.value),basename(str(fn))[:-7] + stem_template +(format_for_output % val_for_output))
            if output_cif.value:
                output.write_cif_data(gs,filename_base)
            if output_xyd.value:
                output.write_xyd_data(gs,filename_base,comment="#",extension="xyd")
            if output_fxye.value:
                output.write_fxye_data(gs,filename_base)
            if output_topas.value:
                output.write_xyd_data(gs,filename_base,comment="'",extension="xye")
            #loop to next group of datasets
            current_frame_start = frameno
            frameno += 1
            
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

