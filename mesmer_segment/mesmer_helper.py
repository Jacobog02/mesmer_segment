## Libraries
import click ## Pretty Commandline
import glob ## Directory Parsing
import os ## To check files, directories
from mesmer_helper import * ## import all my functions to do the things
from skimage import io, color ## to make interpretable plot of segmentation
from matplotlib import pyplot as plt ## to show plots not needed in final
import numpy as np ## Image processing

## Function checks directory structure, quits if it is not there, If it is true then return all jpg images in /CellComposite
def check_smi_dirs(path):
    
    comp = os.path.isdir(path + "/CellComposite")
    over = os.path.isdir(path + "/CellOverlay")

    check = comp & over
    
    if not check: 
        click.echo("Nanostring Directory Structure (/CellComposite & /CellOverlay) not detected quitting....")
        exit()
    else: 
        comp_paths = glob.glob("{}/CellComposite/*.jpg".format(path))
        comp_paths = sorted(comp_paths) ## Sort to give a logical structure. 
        
        return(comp_paths)
    
    
    
    
## Accepts one jpg input and returns numpy array ready for segmentation 
def prepare_composite(one):
    
    one_fov = io.imread(one) ## read in input jpeg with all its beautiful channels
    
    ## Extract only Blue or Green
    ## I have noticed weird performance when using the opposite indexing. 
    #fov_GB = one_fov[:, :, (1, 2)] ## slice out only the second and third indexes (Green & Blue)
    fov_GB = one_fov[:, :, (2, 1)] ## Flip indexes to see if that ruins segmentation? Expects nuclear DAPI dim1 then Membrane

    ## Add BATCH dim on the first axes to match mesmer input.
    fov_GB = np.expand_dims(fov_GB,0)
    
    return(fov_GB)


## Function Accepts an_fov (F001), wd (path to smi dir), mesm_seg (numpy array of M x N x2) 
def vis_fov(an_fov,wd,mesm_segs):
    
    ## Read in default segment image
    overlay_path = "{}/CellOverlay/CellOverlay_{}.jpg".format(wd,an_fov)

    default_im = io.imread(overlay_path)
    
    ## Composite path (Note I could use the path parsed before but I just want to pass the two above. 
    comp_path = "{}/CellOverlay/CellOverlay_{}.jpg".format(wd,an_fov)
    comp_fov = io.imread(comp_path) ## read in input jpeg with all its beautiful channels

    ## Make only Green Blue Image for segmentation result plotting :) 
    #img_GB = comp_fov.copy() ## Normally make a copy to not mess but we don't care
    comp_fov[:, :, (0,1)] = 0 ## Fill RED & GREEN channel with 0's to visualize only DAPI

    img_BW = color.rgb2gray(comp_fov)
    
    
    ## Split Segmentation
    cell_seg = mesm_segs[:, :, (0)] ## Whole Cell
    nuc_seg = mesm_segs[:, :, (1)] ## Neuclear
    
    
    plot_out = "{}/MesmerSegmentation/Mesmer_comparision_{}.jpg".format(wd,an_fov)
    
    ## Start Plot and Save
    fig, axes = plt.subplots(1, 3, figsize=(12, 8), sharey=True)
    
    ## Plot Default from nanostring
    axes[0].imshow(comp_fov)
    axes[0].set_title("Nanostring Default")

    ## Plot WHOLE CELL
    axes[1].imshow(img_BW, cmap='gray')
    axes[1].contour(cell_seg, [0.5], linewidths=.8, colors='teal')
    axes[1].set_title("Nuclear Segmentation")
    
    ## Plot 
    axes[2].imshow(img_BW, cmap='gray')
    axes[2].contour(nuc_seg, [0.5], linewidths=.8, colors='teal')
    axes[2].set_title("Whole Cell Segmentation")


    for a in axes:
        a.axis('off')

    plt.tight_layout()

    plt.savefig(plot_out,dpi=300)
    
    return()

    
    
    
    
    
    
    