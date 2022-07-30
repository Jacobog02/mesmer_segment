## Jacob Gutierrez | Satpathy Lab | jacobog@stanford.edu | 5/5/22
## This script is a callable function which is input with a directory path to nanostring SMI results. Extracts the CellComposite Images for MESMER segmentation. MESMER segmentation masks are written to file for downstream use. Additional plotting functions to visualize the mask output on the DAPI channel compared directly to the Nanostring CellSegmentation images.

## Libraries

from cmath import nan
import click  ## Pretty Commandline
import os

from mesmer_segment.mesmer_helper import *  ## import all my functions to do the things
from mesmer_segment.process_masks import convert_masks

from skimage import io, color  ## to make interpretable plot of segmentation
from matplotlib import pyplot as plt  ## to show plots not needed in final
from deepcell.applications import Mesmer  ## The good stuff.
import numpy as np  ## Image processing


@click.command()
@click.option(
    "--input",
    "-i",
    required=True,
    help="Path to SMI output directory (not the CellComposite Directory!) for resegmentation with MESMER",
)
@click.option(
    "--visual",
    "-v",
    is_flag=True,
    default=False,
    required=False,
    help="Logical Flag to generate plots of new segmentation compared to default SMI segmentation",
)
# @click.option(
#     "--save-npz",
#     is_flag=True,
#     default=False,
#     required=False,
#     help="Logical Flag to determine way to save segmentation masks",
# )
@click.option(
    "--build-outs",
    is_flag=True,
    default=False,
    required=False,
    help="Logical Flag to whether to convert segmentation masks into CosMx SMI-like output files",
)
@click.option(
    "--skip-segmentation",
    is_flag=True,
    default=False,
    required=False,
    help="Logical Flag to whether to run segmentation (if its already been run)",
)
def cli(input, visual, build_outs, skip_segmentation):
    """
    mesmer_segment: Process SMI CellComposite images for segmentation.\n
    Jacob Gutierrez
    """
    #segment(input, visual, save_npz, build_outs, skip_segmentation)
    segment(input, visual, build_outs, skip_segmentation) ## For TIFF input 


## JG 7/27/22: Forcing input as TIF images. 
def segment(input, visual, build_outs, skip_segmentation):

    ## set working dir as input.
    wd = input

    ## Check if directory valid and return input cellcomposite paths for analysis
    check_smi_dirs(input) ## stops if the structure is incorrect

    
    ## Create output path
    seg_dir = "{}/MesmerSegmentation/".format(wd)
    os.makedirs(seg_dir, exist_ok=True)
    
    
    
    ## Process TIF input
    ## NOTE: the way this function works requires the base directory to be named the exact sample name following a common naming convention.
    input_np = process_tif(wd, seg_dir)

    if not skip_segmentation:
          ## We have files now so Spark up Mesmer
          app = Mesmer()
    
          ## Run Mesmer 
          ## Updated to write both nuclear and membrane 
          whole_param = {"maxima_threshold": .3,"interior_threshold":.6}
          nuc_param = {"maxima_threshold" : .3,"interior_threshold":.8}
          run_mesmer(input_np, seg_dir,mes_app = app,image_mpp=0.18,whole_param = whole_param, nuc_param = nuc_param) ## update to npz compressed & for both nuclear + whole cell 
          click.echo("All FOV's segmented DONE")

## Logical for Visualization
    if visual:
       ## Evaluate over list and directly write results
       vis_wrapper(wd,title = os.path.basename(wd.rstrip('/')),parallel = True)

       #vis_fov(fov, wd, seg_fov)
       click.echo("VISUALIZATION WRITTEN")    
    
     
    ## Convert masks into Nanostring CosMx SMI-like output files
    if build_outs:
        click.echo("Converting FOV masks to Nanostring CosMx SMI-like output files")
        ## Create output path
        output_dir = "{}/mesmer_nanostring_outs/".format(wd)
        os.makedirs(output_dir, exist_ok=True)
        convert_masks(
            nanostring_outs_path=input,
            mask_directory=seg_dir,
            output_prefix="mesmer",
            save_path=output_dir,
            fovs=None,
            expand_pixels=0,
            verbose=True,
            num_processes=16,
        )


def segment_OLD(input, visual, save_npz, build_outs, skip_segmentation):
    ## set working dir as input.
    wd = input

    ## Check if directory valid and return input cellcomposite paths for analysis
    all_in_paths = check_smi_dirs(input)

    ## Extract FOV names for easy writing and reading.
    all_fovs = [
        os.path.splitext(os.path.basename(n))[0].split("_")[1] for n in all_in_paths
    ]

    ## Create output path
    seg_dir = "{}/MesmerSegmentation/".format(wd)
    os.makedirs(seg_dir, exist_ok=True)

    if not skip_segmentation:
        ## We have files now so Spark up Mesmer
        app = Mesmer()

        ## Loop over each FOV and create segmentation
        ## In theory could parallelize across images here depending on MESMER size which I am not pythonic enough to answer
        ## i is a string...
        for i, fov in enumerate(all_fovs):
            click.echo("Processing: {}".format(fov))
            fov_GB = prepare_composite(all_in_paths[int(i)])

            ## this gives an N x M x 2 matrix. Index 0 == Whole-Cell | Index 1 == Nueclear
            seg_fov = app.predict(
                fov_GB, compartment="both"
            )  ## Do the thing! Takes ~80-120 seconds per image

            ## Squash the extra batch dim
            seg_fov = np.squeeze(seg_fov)

            ## Save whole cell and nuclear masks
            if save_npz:
                save_path = "{}/mesmer_{}.npz".format(seg_dir, fov)
                np.savez_compressed(
                    save_path, wholecell=seg_fov[:, :, (0)], nuclear=seg_fov[:, :, (1)]
                )
                # load the masks with:
                # data = np.load(save_path)
                # nuclear_mask = data["nuclear"]
                # wholecell_mask = data["wholecell"]
            else:
                ## Save Whole Cell Mask
                whole_cell_path = "{}/mesmer_whole_cell_{}.csv".format(seg_dir, fov)
                np.savetxt(whole_cell_path, seg_fov[:, :, (0)], delimiter=",")

                ## Save Nuclear Mask
                nuc_path = "{}/mesmer_nucleus_{}.csv".format(seg_dir, fov)
                np.savetxt(nuc_path, seg_fov[:, :, (1)], delimiter=",")

            ## Logical for Visualization
            if visual:
                vis_fov(fov, wd, seg_fov)
                click.echo("VISUALIZATION WRITTEN")

        click.echo("All FOV's segmented DONE")

    ## Convert masks into Nanostring CosMx SMI-like output files
    if build_outs and save_npz:
        click.echo("Converting FOV masks to Nanostring CosMx SMI-like output files")
        ## Create output path
        output_dir = "{}/mesmer_nanostring_outs/".format(wd)
        os.makedirs(output_dir, exist_ok=True)
        convert_masks(
            nanostring_outs_path=input,
            mask_directory=seg_dir,
            output_prefix="mesmer",
            save_path=output_dir,
            fovs=None,
            expand_pixels=12,
            verbose=True,
            num_processes=16,
        )
