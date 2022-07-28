## Libraries
import click ## Pretty Commandline
import glob ## Directory Parsing
import os ## To check files, directories
from  skimage import io, color ## to make interpretable plot of segmentation
from matplotlib import pyplot as plt ## to show plots not needed in final
import numpy as np ## Image processing
import pandas as pd ## Read text files
from tqdm import tqdm ## progress bar
from PIL import Image ## IDK why we are using this package honestly
from natsort import natsorted ## IDK why use this sort function
import pandas as pd
import sys
import multiprocessing as mp 
from multiprocessing import Pool
from pkg_resources import get_distribution
import pkg_resources

## Function checks directory structure, quits if it is not there, If it is true then return all jpg images in /CellComposite
def check_smi_dirs(path):
    
    ## 7/20/21 
    morph = os.path.isdir(path + "/RawMorphologyImages")
    comp = os.path.isdir(path + "/CellComposite")
    over = os.path.isdir(path + "/CellOverlay")

    check = morph & comp & over 
    
    if not check: 
        click.echo("Nanostring Directory Structure (/CellComposite & /CellOverlay & /RawMorphologyImages) not detected quitting....")
        exit()
    else: 
                          
        ## NOW MODIFIED TO RETURN NOTHING ALL FUNCTIONS PROCESS ON THEIR OWN                    
        #comp_paths = glob.glob("{}/CellComposite/*.jpg".format(path))
        #comp_paths = sorted(comp_paths) ## Sort to give a logical structure. 
        #return(comp_paths)
        return() ## Return nothing if successful
    
    
    
    
## Accepts one jpg input and returns numpy array ready for segmentation 
## Remember R people, python is 0 indexed so 2 == 3rd position
def prepare_composite(one,dapi_idx = 2):
    
    one_fov = io.imread(one) ## read in input jpeg with all its beautiful channels
    
    ## V1: use basic 2 channels... 
    ## Extract only Blue or Green
    ## I have noticed weird performance when using the opposite indexing. 
    #fov_GB = one_fov[:, :, (1, 2)] ## slice out only the second and third indexes (Green & Blue)
    #fov_GB = one_fov[:, :, (2, 1)] ## Flip indexes to see if that ruins segmentation? Expects nuclear DAPI dim1 then Membrane

    ## V2: Accept a dapi channel index and then shove everything else into one channel and switch
    ## I have convinced myself I will average the not dapi channels to merge. https://e2eml.school/convert_rgb_to_grayscale.html
    dapi_channel = one_fov[:,:,(dapi_idx)]
    squanch_idx = list(range(one_fov.shape[2])) ## should be able to accept a tiff file in theory and squish it together. 
    membrane_channel = np.mean( one_fov[:,:,(squanch_idx)], axis=2).astype(np.uint8) ## makes a float but force into int
    ## Stack them along the last axes.
    fov_GB = np.dstack((dapi_channel,membrane_channel))
    
    
    ## Add BATCH dim on the first axes to match mesmer input.
    fov_GB = np.expand_dims(fov_GB,0)
    
    return(fov_GB)
  
  


## Function Accepts an_fov (F001), wd (path to smi dir), mesm_seg (numpy array of M x N x2) 
## JG 7/26/22: Updated to overlay mesmer segmentation in yellow in one image
def vis_fov(an_fov,wd,title=""):
    print("Visualizing: {}".format(an_fov))
    
    ## Read in default segment image
    overlay_path = "{}/CellOverlay/CellOverlay_{}.jpg".format(wd,an_fov)

    default_im = io.imread(overlay_path)
    
    ## Composite path (Note I could use the path parsed before but I just want to pass the two above. 
    #comp_path = "{}/CellComposite/CellComposite_{}.jpg".format(wd,an_fov)
    #comp_fov = io.imread(comp_path) ## read in input jpeg with all its beautiful channels

    ## Make only Green Blue Image for segmentation result plotting :) 
    #img_GB = comp_fov.copy() ## Normally make a copy to not mess but we don't care
    #comp_fov[:, :, (0,1)] = 0 ## Fill RED & GREEN channel with 0's to visualize only DAPI

    #img_BW = color.rgb2gray(comp_fov)
    
    
    ## Generalized to read in the output from our specific mesmer pipeline. 
    #mesm_segs_path = "{}/MesmerSegmentation/{}_nuclear_membrane_mesmer.npy".format(wd,an_fov)
    mesm_segs_path = "{}/MesmerSegmentation/{}_nuclear_membrane_mesmer.npz".format(wd,an_fov)

    
    mesm_segs = np.load(mesm_segs_path)
    #('wholecell', 'nuclear')
    #print(mesm_segs['wholecell'].shape())
    #print(mesm_segs['wholecell'].shape)
    
    ## Split Segmentation
    #cell_seg = mesm_segs[:, :, (0)] ## Whole Cell
    #nuc_seg = mesm_segs[:, :, (1)] ## Nuclear
    cell_seg = mesm_segs['wholecell']
    nuc_seg = mesm_segs['nuclear']

    
    plot_out = "{}/MesmerSegmentation/{}_Mesmer_Comparision.jpg".format(wd,an_fov)
    
    ## Start Plot and Save
    fig, axes = plt.subplots(1, 1, figsize=(8, 6), sharey=True)
    
    ## Plot Default from nanostring
    #axes[0].imshow(default_im)
    #axes[0].set_title("Nanostring vs. Mesmer")
    #axes[0].contour(cell_seg, [0.5], linewidths=.2, colors='goldenrod')
    
    axes.imshow(default_im)
    axes.set_title("Nanostring vs. Mesmer {}:{}".format(title,an_fov))
    axes.contour(cell_seg, [0.5], linewidths=.4, colors='goldenrod')
    axes.contour(nuc_seg, [0.5], linewidths=.4, colors='indigo')
    axes.axis('off')
    
    
    ## Plot WHOLE CELL
    #axes[1].imshow(img_BW, cmap='gray')
    #axes[1].contour(cell_seg, [0.5], linewidths=.3, colors='teal')
    #axes[1].set_title("Mesmer Segmentation")
    

    #for a in axes:
        #a.axis('off')

    plt.tight_layout()

    plt.savefig(plot_out,dpi=300)
    
    plt.close()
    
    return()



## infiles are the direct paths to 
def vis_wrapper(basedir, title, parallel = True):
    
    
    
    ## REMOVE THE FINAL '/' BC basename returns empty string :( 
    basedir = basedir.rstrip('/') ## this does nothing if it is missing
    sample_name = os.path.basename(basedir)
    
    pos = pd.read_csv(os.path.join(basedir, sample_name+'_fov_positions_file.csv'), index_col=0)
    
    #fovs = ['F{0:03d}'.format(i) for i in pos.index]
    
    if parallel is True: 
        pool_arg_list = [('F{0:03d}'.format(i), basedir,sample_name) for i in pos.index]
        
        with Pool(processes=mp.cpu_count()) as pool:
                 pool.starmap(vis_fov,pool_arg_list )
        
    else:
        [vis_fov('F{0:03d}'.format(i), basedir,title = sample_name) for i in pos.index]


### Peter Du's TIFF processing pipeline
### Must try to generalize as a plug-n-play function for future datasets... 

## Peter's workflow is as follows: 1) process_tif (this calls composite_zstack) 2) run_mesmer on "_nuclear_membrane.npy" as input
def composite_zstack(filelist, n_channels=5):
    stacks = []
    for f in filelist:
        stack = []
        im = Image.open(f)
        for i in range(n_channels):
            im.seek(i)
            stack.append(np.array(im))
        stacks.append(np.stack(stack, axis=-1))
    composite = np.stack(stacks, axis=-1)
    return composite
    
    
def composite_zstack_base(f,n_channels):
    stack = []
    im = Image.open(f)
    for i in range(n_channels):
        im.seek(i)
        stack.append(np.array(im))
   
    stack = np.stack(stack, axis=-1) ## Do the stack in parallel
    return(stack)
    
    
def composite_zstack_parallel(filelist, n_channels=5):
    
    pool_arg_list = [(K, n_channels) for K in filelist]
        
    with Pool(processes=mp.cpu_count()) as pool:
             stacks = pool.starmap(composite_zstack_base,pool_arg_list )

            
    composite = np.stack(stacks, axis=-1)
    return composite



## JG 7/20/22: Modifying function to spit out the nuclear_membrane.npy paths to make easier
def process_tif(basedir, outdir):
    
    ## REMOVE THE FINAL '/' BC basename returns empty string :( 
    basedir = basedir.rstrip('/') ## this does nothing if it is missing
    sample_name = os.path.basename(basedir)

    tifdir = os.path.join(basedir, 'RawMorphologyImages')
    pos = pd.read_csv(os.path.join(basedir, sample_name+'_fov_positions_file.csv'), index_col=0)
    
    
    out_mem = []
    for i in tqdm(pos.index, desc='fov'):
        fov = 'F{0:03d}'.format(i)
        outpath = os.path.join(outdir, fov+'_nuclear_membrane.npy')
        
        ## if file exists then skip the computationally expensive part
        if os.path.exists(outpath):
            print(f'FOV previously processed {fov}')
            ## Add path to output
            out_mem.append(outpath)
            
        else:      
            image_files = [os.path.join(tifdir, f) for f in os.listdir(tifdir) if f'_{fov}_' in f]
            if len(image_files) == 0:
                print(f'No files found for {fov}')
                continue
            #composite = composite_zstack(image_files) ## One thread approach
            composite = composite_zstack_parallel(image_files) ## parallelized read in function
            max_stack = composite.max(axis=-1)

            # (nuclear, membrane)
            nuc_mem = max_stack[:, :, [4, 0]] # 4 == DAPI, 0 == membrane        
            np.save(outpath, nuc_mem)

            ## Add path to output
            out_mem.append(outpath)
        
    return out_mem
        

## https://github.com/vanvalenlab/deepcell-tf/blob/master/notebooks/applications/Mesmer-Application.ipynb        
def run_mesmer(infiles, outdir, mes_app,n_batches=2, image_mpp=0.18, compart = "both",maxima_threshold = 0.8,interior_threshold = 0.5):
    # infiles
    stack = []
    for f in natsorted(infiles):
        stack.append(np.load(f))
    data = np.stack(stack)
    
    preds = []
    batches = np.digitize(range(len(stack)), np.linspace(0, len(stack), n_batches+1))
    
    
    # batch process
    for i in tqdm(range(1, n_batches+1), desc='batch'):
        idx = np.where(batches == i)[0]

        pred = mes_app.predict(data[idx, ...], image_mpp=image_mpp,compartment=compart,postprocess_kwargs_whole_cell={'maxima_threshold': maxima_threshold, 'interior_threshold': interior_threshold },postprocess_kwargs_nuclear={'maxima_threshold': maxima_threshold, 'interior_threshold': interior_threshold })
        
        # save results
        for ii, j in enumerate(idx):
            infile = os.path.basename(infiles[j])
            #outpath = os.path.join(outdir, os.path.splitext(infile)[0] + '_mesmer.npy')
            #np.save(outpath, pred[ii, ...])
            
            outpath = os.path.join(outdir, os.path.splitext(infile)[0] + '_mesmer.npz')

            #save_path = "{}/mesmer_{}.npz".format(seg_dir, fov)
            np.savez_compressed(
                    outpath, wholecell= pred[ii, :, :, (0)], nuclear=pred[ii,:, :, (1)]
                )

    
   
# wd (path to smi dir)
# seg_files (iter of paths to segmentation, .csv or .npz)
# seg_type ('wholecell' or 'nuclear' for seg.npz)
# return_dict (bool; return results as dictionary instead of writing to disk) 
def assign_molecules(wd, seg_files, seg_type='wholecell', return_dict=False):
    assert seg_type in ('wholecell', 'nuclear'), click.echo('seg_type must be wholecell or nuclear')
    
    sample_name = os.path.basename(wd)
    
    # load molecules
    tx_file = os.path.join(wd, f'{sample_name}_tx_file.csv')
    tx = pd.read_csv(tx_file)
    
    # load segmentations
    counts_dict = {}
    for f in seg_files:
        if not os.path.exists(f):
            click.echo(f'{f} not found, skipping')
            continue
            
        fov = os.path.splitext(f)[0].split('_')[-1]
        fov_int = int(fov[1:]) ## remove the leading F char
        local = tx[tx['fov'] == fov_int].copy()
        
        if f.endswith('.npz'):
            seg = np.load(f)[seg_type]
        elif f.endswith('.csv'):
            seg = np.loadtxt(f)
        else:
            click.echo('Segmentation file extension not recognized')
            continue
    
        # make counts matrix
        xcoord, ycoord = local[['x_local_px', 'y_local_px']].round().astype(int).values.T
        local['mesmer_cell_ID'] = seg[ycoord, xcoord]
        counts = local.groupby(['mesmer_cell_ID', 'target']).size().reset_index() ## count molecules
        counts = counts.pivot(index='mesmer_cell_ID', columns='target', values=0) ## pivot into matrix
        counts = counts.replace(np.nan, 0).astype(int).reset_index()
        if return_dict:
            counts_dict[os.path.basename(f)] = counts
        else:
            counts.to_csv(os.path.join(wd, f'{sample_name}_{fov}_mesmer_counts.csv'), index=False)
    
    if return_dict:
        return counts_dict
    
    
    
    
    
