import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from kstar import plot

def kstar_plot(activities, fpr, colormap, background_color = 'white', samples = None, x_label_dict = None, sig_restrict = False, kinases_to_drop = None, binary_sig = True, sort_kinases = True, sort_samples = True, binary_evidence = None, evidence_type = 'dots', context = None, phosphotype = None, fontsize = 10, figsize = (4,10)):
    #log transform activities
    log_results = -np.log10(activities)

    if samples is not None:
        log_results = log_results[samples]
        fpr = fpr[samples]

    if binary_evidence is not None:
        include_evidence = True
    else:
        include_evidence = False
    nrows = 1 + sort_samples*1 + include_evidence*1
    ncols = 1 + sort_kinases*1

    #get width/height ratios
    if sort_samples and include_evidence:
        height_ratios = [0.05,1,0.05]
    elif include_evidence:
        height_ratios = [1,0.05]
    elif sort_samples:
        height_ratios = [0.05,1]
    else:
        height_ratios = [1]

    #get width ratios
    if sort_kinases:
        width_ratios = [0.1,1]
    else:
        width_ratios = [1]

    fig, axes = plt.subplots(figsize = figsize,
                             nrows = nrows,
                             ncols = ncols,
                             height_ratios = height_ratios,
                             width_ratios = width_ratios,
                             sharex = 'col', sharey = 'row')
    fig.subplots_adjust(wspace = 0, hspace = 0)
    #instantiate the dotplot object
    dots = plot.DotPlot(log_results,
                fpr,
                legend_title = '-log10(p-value)',
                colormap = colormap,
                facecolor=background_color,
                binary_sig = binary_sig,
                x_label_dict=x_label_dict)

    if sig_restrict:
        dots.drop_kinases_with_no_significance()

    if kinases_to_drop:
        dots.drop_kinases(kinases_to_drop)

    #cluster the kinases with hierarchical clustering if desired
    if sort_kinases:
        ax = axes[1,0] if sort_samples else axes[0]
        dots.cluster(orientation = 'left', ax = ax, method='ward')
        ax.set_xticks([])

    if sort_samples:
        ax = axes[0,1] if sort_kinases else axes[0]
        dots.cluster(orientation = 'top', ax = ax, method='ward')
        axes[0,0].axis('off')
        ax.set_yticks([])

    if binary_evidence is not None:
        ax = axes[2,1] if sort_samples and sort_kinases else axes[2] if sort_samples else axes[1]
        ax.set_yticks([])
        dots.evidence_count(ax = ax, binary_evidence = binary_evidence, plot_type = evidence_type,
                    include_recommendations = True, phospho_type = phosphotype)
        ax.tick_params(labelsize = fontsize)
        ax.tick_params(axis = 'x', rotation = 90)
        if sort_kinases:
            axes[nrows-1,0].axis('off')



    #generate the dotplot
    ax = axes[1,1] if sort_samples and sort_kinases else axes[1] if sort_samples or sort_kinases else axes
    dots.dotplot(ax = ax)
    ax.tick_params(axis = 'y', labelsize = fontsize)
    ax.tick_params(axis = 'x', labelsize = fontsize, rotation = 90)
    st.pyplot(fig)
    return fig


st.title("Plotting KSTAR Results")

st.header("Load Data")


activities_fileupload = st.file_uploader("Upload KSTAR activities file:", type = 'tsv', key = 'activities', help = 'Find the file in the KSTAR results file labeled "..._mann_whitney_activities.tsv"')
fpr_fileupload = st.file_uploader("Upload KSTAR false positive rate file:", type = 'tsv', key = 'fpr', help = 'Find the file in the KSTAR results file labeled "..._mann_whitney_fpr.tsv"')

if activities_fileupload is not None and fpr_fileupload is not None:
    activities = pd.read_csv(activities_fileupload, sep = '\t', index_col = 0)
    fpr = pd.read_csv(fpr_fileupload, sep = '\t', index_col = 0)
    st.success('Activities and fpr file uploaded')


st.header('Plot')

#plot_type = st.selectbox('Select Plot Type:', ['Dotplot', 'Barplot'], key = 'plot_type')

t1, t2, t3, t4 = st.tabs(['Figure Parameters', 'Filter Kinases and Samples', 'Sort Kinases/Samples', 'Add Additional Context'])
with t1:
    #establish figure size
    st.markdown('*Figure Size*')
    col1, col2 = st.columns(2)
    fwidth = col1.number_input('Width (in)', min_value=1, max_value=20, value=4, key = 'width')
    fheight = col2.number_input('Height (in)', min_value=1, max_value=20, value=10, key = 'height')
    fontsize = st.number_input('Fontsize', min_value = 6, max_value = 18, value = 10, key = 'fontsize')

    
    st.markdown('*Color Scheme*')
    #determine how activity is shown, binary or gradient
    binary_sig = st.radio('How to show significant activity:', ['binary', 'gradient'], horizontal = True, key = 'binary')
    if binary_sig == 'binary':
        binary_sig = True
    else:
        binary_sig = False
    #establish color scheme
    col1, col2, col3 = st.columns(3)
    background_color = col1.color_picker('Background Color', value = '#FFFFFF', key = 'background_color')
    activity_color = col2.color_picker('Color to Represent Activity', value = '#FF3300', key = 'activity_color')
    noactivity_color = col3.color_picker('Color to Represent Lack of Activity', value = '#6b838f', key = 'lackactivity_color')


with t2:
    if 'activities' in locals():
        st.subheader('Filter Kinases')
        sig_restrict = st.checkbox('Restrict to kinases with significance in at least one sample', value = False, key = 'sig_restrict')
        #ask if user wants to drop any kinases
        manual_kinase_edit = st.selectbox('Manually edit kinases in plot:', ['No', 'Remove Specific Kinases', 'Select Specific Kinases'])
        if manual_kinase_edit == 'Remove Specific Kinases':
            kinases_to_drop = st.multiselect('Manually pick kinases to remove from analysis', activities.index, key = 'kinases_to_drop')
        elif manual_kinase_edit == 'Select Specific Kinases':
            kinases_to_keep = st.multiselect('Manually pick kinases to include in analysis', activities.index, key = 'kinases_to_keep')
            kinases_to_drop = [kin for kin in activities.index if kin not in kinases_to_keep]
        else:
            kinases_to_drop = None

        st.subheader('Filter Samples')
        samples_to_plot = st.multiselect('Manually select samples to plot:', options = activities.columns, default = list(activities.columns), key = 'samples_to_plot')
        
        #rename 
        sample_rename = st.checkbox('Change the xlabel for each sample',value = False, key = 'sample_rename')
        if sample_rename:
            x_label_dict = {}
            for samp in samples_to_plot:
                x_label_dict[samp] = st.text_input(label = samp, value = samp, key = samp + '_label')
        else:
            x_label_dict = None


with t3:
    sort_kinases = st.radio('Sort Kinases:', ['No Sorting', 'By Activity', 'By Hierarchical Clustering'], key = 'sort_kinases', horizontal = True)
    if sort_kinases == 'By Activity':
        sample_to_sort_by = st.selectbox('Sort Kinases by Activity in Sample:', activities.columns, key = 'sample_to_sort_by')
        activities = activities.sort_values(by = sample_to_sort_by, ascending = False)
        sort_kinases = False
    elif sort_kinases == 'By Hierarchical Clustering':
        sort_kinases = True
    else:
        sort_kinases = False

    sort_samples = st.radio('Sort Samples:', ['No Sorting', 'By Activity', 'By Hierarchical Clustering'], key = 'sort_samples', horizontal = True)
    if sort_samples == 'By Activity':
        kinase_to_sort_by = st.selectbox('Sort Samples by Activity for Kinase:', activities.index, key = 'kinase_to_sort_by')
        activities = activities.sort_values(by = kinase_to_sort_by, axis = 1, ascending = False)
        sort_samples = False
    elif sort_samples == 'By Hierarchical Clustering':
        sort_samples = True
    else:
        sort_samples = False

with t4:
    st.subheader('Add evidence size used for each sample')
    include_evidence = st.checkbox('Include Evidence Size for Each Sample', key = 'include_evidence')
    if include_evidence:
        phosphotype = st.radio("Type of kinase being plotted:", ['Y', 'ST'], horizontal = True)
        evidence_type = st.radio("How to display evidence size:", ['dots','bars'], horizontal = True)
        evidence_upload = st.file_uploader('Upload binary evidence file:', type = 'tsv', key = 'binary_evidence', help = 'This should be the file labeled "..._evidence_binary.tsv"')
        if evidence_upload is not None:
            binary_evidence = pd.read_csv(evidence_upload, sep = '\t', index_col = 0)
            st.success('Binary evidence uploaded')

    st.subheader('Add additional context for the sample, such as treatment or patient type')
    add_context = st.checkbox('Add Additional Context', key = 'add_context')
    if add_context:
        context = None


st.subheader('Construct DotPlot')
plot_now = st.checkbox('Show Plot', value = False, key = 'show_plot')
if plot_now:
    if 'activities' in locals() and 'fpr' in locals():
        if include_evidence and 'binary_evidence' in locals() and add_context and 'context' in locals():
            fig = kstar_plot(activities, fpr, colormap = {0: noactivity_color, 1: activity_color}, background_color = background_color, figsize = (fwidth, fheight), sig_restrict = sig_restrict, kinases_to_drop = kinases_to_drop, binary_sig = binary_sig, sort_kinases = sort_kinases, sort_samples = sort_samples, binary_evidence=binary_evidence, evidence_type = evidence_type, phosphotype = phosphotype, fontsize = fontsize, context = context, samples = samples_to_plot, x_label_dict=x_label_dict)
        elif include_evidence and 'binary_evidence' not in locals() and add_context and 'context' not in locals():
            st.error('Please upload binary evidence file and context file')
        elif include_evidence and 'binary_evidence' in locals() and add_context and 'context' not in locals():
            fig = st.error('Please upload context file')
        elif include_evidence and 'binary_evidence' in locals():
            fig = kstar_plot(activities, fpr, colormap = {0: noactivity_color, 1: activity_color}, background_color = background_color, figsize = (fwidth, fheight), sig_restrict = sig_restrict, kinases_to_drop = kinases_to_drop, binary_sig = binary_sig, sort_kinases = sort_kinases, sort_samples = sort_samples, binary_evidence=binary_evidence, evidence_type = evidence_type, phosphotype=phosphotype, fontsize = fontsize, samples = samples_to_plot, x_label_dict = x_label_dict)
        elif add_context and 'context' in locals():
            fig = kstar_plot(activities, fpr, colormap = {0: noactivity_color, 1: activity_color}, background_color = background_color, figsize = (fwidth, fheight), sig_restrict = sig_restrict, kinases_to_drop = kinases_to_drop, binary_sig = binary_sig, sort_kinases = sort_kinases, sort_samples = sort_samples, context = context, fontsize = fontsize, samples = samples_to_plot, x_label_dict = x_label_dict)
        elif include_evidence and 'binary_evidence' not in locals():
            st.error('Please upload binary evidence file')
        else:
            fig = kstar_plot(activities, fpr, colormap = {0: noactivity_color, 1: activity_color}, background_color = background_color, figsize = (fwidth, fheight), sig_restrict = sig_restrict, kinases_to_drop = kinases_to_drop, binary_sig = binary_sig, sort_kinases = sort_kinases, sort_samples = sort_samples, fontsize = fontsize, samples = samples_to_plot, x_label_dict = x_label_dict)
    else:
        st.error('Please upload both files')

        
    file_name = st.text_input('Figure name:', value = 'KSTAR_dotplot')
    extension = st.radio('Type of file to save as:', ['.png','.jpg', '.pdf', '.svg'], horizontal = True)
    itype = extension.split('.')[1]
    file_name = file_name + extension
    plt.savefig(file_name,  bbox_inches = 'tight')
    with open(file_name, 'rb') as img:
        st.download_button(label = 'Download Figure', data = img, file_name = file_name, mime='image/'+itype)
    
    #file_name = st.text_input('Figure name:', value = 'KSTAR_dotplot')
    #extension = st.radio('Type of file to save as:', ['.png','.jpg', '.pdf', '.svg'], horizontal = True)
    #st.write(file_name + extension)
    #st.write(extension.strip('.'))
    #plt.savefig(file_name + extension, format = extension.strip('.'), bbox_inches = 'tight')
    #with open(file_name+extension, 'rb') as img:
    #    st.download_button(label = 'Download Figure', data = img, file_name = file_name, mime='image/png')
