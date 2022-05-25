"""Functions to compute Nanostring CosMX SMI-like output files from segmentation masks which can be processed by downstream tools such as Seurat."""

import numpy as np
import pandas as pd
from skimage.segmentation import expand_labels, find_boundaries

import random
import os
import time
import multiprocessing


### Constants
TX_FILE_COLNAMES = []
POLYGONS_FILE_COLNAMES = ["fov", "cellID", "x_local_px", "y_local_px", "x_global_px", "y_global_px"]
EXPRMAT_FILE_COLNAMES = []
METADATA_FILE_COLNAMES = [
        "fov",
        "cell_ID",
        "CenterX_local_px",
        "CenterY_local_px",
        "CenterX_global_px",
        "CenterY_global_px"
    ]
METADATA_X_GLOBAL = "CenterX_global_px"
METADATA_Y_GLOBAL = "CenterY_global_px"
POLYGONS_X_GLOBAL = "x_global_px"
POLYGONS_Y_GLOBAL = "y_global_px"
POLYGONS_X_LOCAL = "x_local_px"
POLYGONS_Y_LOCAL = "y_local_px"
FOV_SIZE_X = 3638
FOV = "fov"
CELL_ID = "cell_ID"


def convert_masks(
    nanostring_outs_path,
    mask_directory,
    output_prefix,
    save_path,
    fovs=None,
    expand_pixels=12,
    verbose=True,
    num_processes=16
):
    """Rebuild the gene expression matrix, generate segmentation vertices, and reassign each transcripts CellComp attribute.

    Args:
        nanostring_outs_path (str): Path to the default Nanostring CosMx SMI outputs.
        mask_directory (str): Path to directory containing compressed .npz segmentation masks for each FOV.
        output_prefix (str): Prefix for each newly generated output file.
        save_path (str): Path to directory to save new output files
        expand_pixels (int): Number of pixels to expand segmentations by. Defaults to 12.
        verbose (bool): Defaults to True.
    """
    if verbose:
        print("Starting timer...")
        start = time.time()

    for file in os.listdir(nanostring_outs_path):
        if file.endswith("fov_positions_file.csv"):
            fov_positions_file = os.path.join(nanostring_outs_path, file)
        elif file.endswith("tx_file.csv"):
            tx_file = os.path.join(nanostring_outs_path, file)

    # check paths
    files = [fov_positions_file, tx_file, mask_directory, save_path]
    # dirs = [mesmer_masks, save_path]
    for file in files:
        if not os.path.exists(file):
            msg = f"{file} does not exist"
            raise Exception(msg)  

    expr_path = os.path.join(save_path, output_prefix + "_expreMat_file.csv")
    tx_path = os.path.join(save_path, output_prefix + "_tx_file.csv")
    polygons_path = os.path.join(save_path, output_prefix + "_polygons_file.csv")
    metadata_path = os.path.join(save_path, output_prefix + "_metadata_file.csv")

    # find compressed masks
    masks = [
        x for x in os.listdir(mask_directory) if x.endswith("npz")
    ]
    mask_files = dict()
    for mask in masks:
        # get FOV number from filename
        mask_number = int(mask[-6:-4])
        if fovs:
            if mask_number in fovs:
                mask_files[mask_number] = os.path.join(mask_directory, mask)
        else:
            mask_files[mask_number] = os.path.join(mask_directory, mask)
    if verbose:
        print(mask_files)

    # construct outputs in parallel by FOV
    results = build_parallel(
        tx_file=tx_file,
        fov_positions_file=fov_positions_file,
        mask_files=mask_files,
        expand_pixels=expand_pixels,
        verbose=verbose,
        num_processes=num_processes
    )

    # save dfs to csv
    results["transcripts"].to_csv(tx_path, index=False)
    results["expression_matrix"].to_csv(expr_path, index=False)
    results["segmentation_vertices"].to_csv(polygons_path, index=False)
    results["metadata"].to_csv(metadata_path, index=False)

    if verbose:
        print(results["transcripts"].head())
        print(results["expression_matrix"].head())
        print(set(results["expression_matrix"][FOV]))
        print(results["expression_matrix"].shape)
        print(results["transcripts"].shape)
        print(f"The code took {time.time() - start} to run.")


def build_parallel(
    tx_file,
    mask_files,
    fov_positions_file,
    expand_pixels=12,
    verbose=True,
    num_processes = 16
):
    """Build an expression matrix from a segmentation mask and generate an updated transcripts file.

    Args:
        tx_file (str): The path to the default Nanostring transcripts file. This file contains the position of each transcript detected.
        mask_files (dict): Dictionary where keys represent the name or FOV number for a particular mask and the value is a path to the mask.
        expand_pixels (int): Number of pixels to expand segmentations by. Defaults to 12.
        verbose (bool): Defaults to True.

    Returns:
        The output is a dictionary containing the updated transcripts dataframe and a requantified expression matrix

            {
                'transcripts': transcripts_df,
                'expression_matrix': expression_matrix
            }
    """
    txs = pd.read_csv(tx_file)
    fov_shifts = get_fov_shifts(fov_positions_file=fov_positions_file)

    # map each transcript to an index
    tx_map = dict()
    for i, t in enumerate(list(set(txs.target))):
        tx_map[t] = i
        
    # get set of unique genes
    genes = set(txs.target)

    # initialize returned dataframes
    new_txs = pd.DataFrame(columns = txs.columns)
    new_expr = pd.DataFrame(columns = genes)

    segmentation_vertices = pd.DataFrame(columns=POLYGONS_FILE_COLNAMES)
    metadata_df = pd.DataFrame(columns=METADATA_FILE_COLNAMES)

    pool = multiprocessing.Pool(num_processes)
    start_time = time.perf_counter()
    processes = [
        pool.apply_async(
            construct_fov_files, args=(txs, fov, mask_files[fov], fov_shifts, )
        ) for fov in mask_files.keys()
    ]
    results = [p.get() for p in processes]
    finish_time = time.perf_counter()
    print(f"Program finished in {finish_time-start_time} seconds")
    print("Concatenating results...")
    for result in results:
        # append the updated transcripts and expression for current FOV
        new_txs = pd.concat([new_txs, result["transcripts"]])
        new_expr = pd.concat([new_expr, result["expression_matrix"]])

        # get polygon vertices
        segmentation_vertices = pd.concat([segmentation_vertices, result["segmentation_vertices"]])
        metadata_df = pd.concat([metadata_df, result["metadata"]])

    return {
        "transcripts": new_txs,
        "expression_matrix": new_expr,
        "segmentation_vertices": segmentation_vertices,
        "metadata": metadata_df
    }


def construct_fov_files(
    txs,
    fov,
    mask_file,
    fov_shifts,
    expand_pixels=12,
    verbose=True
):
    # map each transcript to an index
    tx_map = dict()
    for i, t in enumerate(list(set(txs.target))):
        tx_map[t] = i
        
    # get set of unique genes
    genes = set(txs.target)

    print(f"Building the expression matrix for FOV {fov}.")
    if verbose:
        start = time.time()
    txs_subset = txs[txs.fov == fov].copy()

    # check file extension and read masks into np.ndarrays
    if mask_file.endswith("npz"):
        mask = np.load(mask_file)["wholecell"]
    elif mask_file.endswith("csv"):
        mask = np.genfromtxt(mask_file, delimiter=",")
    else:
        raise Exception("Unknown file extension")

    # expand segmentations
    mask = expand_labels(mask, distance=expand_pixels)

    # create empty expression matrix
    n_cells = len(np.unique(mask)) - 1
    exp_mtx = np.zeros((n_cells, len(genes)), dtype=int)

    # flip y coordinates
    max_y = max(txs_subset.y_local_px)
    txs_subset["tmp_y_local_px"] = abs(txs_subset.y_local_px - max_y)

    cellcomp, cell_id = [], []
    for i, row in txs_subset.iterrows():
        if verbose and i % 100000 == 0:
            print(f"\t\t{i} transcripts processed in FOV {fov}")
        x_pos = int(row.x_local_px)
        y_pos = int(row.tmp_y_local_px)

        # Check if position is non-zero, indicating transcript is within a cell.
        v = mask[y_pos, x_pos]
        if v > 0:
            exp_mtx[int(v) - 1, tx_map[row.target]] += 1
            cellcomp.append("Cytoplasm")
        else:
            cellcomp.append("0")
        cell_id.append(v)

    txs_subset.drop(["tmp_y_local_px"], inplace=True, axis=1)

    # replace CellComp column with updated CellComp reflecting mesmer segmentations
    txs_subset["CellComp"] = cellcomp
    txs_subset[CELL_ID] = cell_id
    txs_subset[CELL_ID] =  txs_subset[CELL_ID].astype(int)

    # convert expression array to dataframe
    cell_ids = [int(i+1) for i in range(n_cells)]
    matrix = pd.DataFrame(exp_mtx, columns=genes)
    
    # add cell_ID and fov columns to the matrix
    matrix[CELL_ID] = cell_ids
    matrix[CELL_ID] = matrix[CELL_ID].astype(int)
    matrix[FOV] = fov
    matrix[FOV] = matrix[FOV].astype(int)

    # get polygon vertices
    vertice_files = get_polygon_vertices(fov, mask)

    # shift global coordinates
    shift = fov_shifts[fov]
    x_shift = shift["x"]
    y_shift = shift["y"]

    vertice_files["metadata_df"][METADATA_X_GLOBAL] += x_shift
    vertice_files["metadata_df"][METADATA_Y_GLOBAL] += y_shift
    vertice_files["segmentations_df"][POLYGONS_X_GLOBAL] += x_shift
    vertice_files["segmentations_df"][POLYGONS_X_GLOBAL] += y_shift

    if verbose:
        print(f"Completed FOV {fov} in {round(time.time() - start)} seconds")

    return {
        "expression_matrix": matrix,
        "transcripts": txs_subset,
        "segmentation_vertices": vertice_files["segmentations_df"],
        "metadata": vertice_files["metadata_df"]
    }


############# Utility functions ###################


def get_fov_shifts(fov_positions_file):
    with open(fov_positions_file) as f:
        lines = f.readlines()
    shifts = dict()
    for line in lines[1:]:
        line = line.strip()
        vals = line.split(",")
        shifts[int(vals[0])] = {
            "x": float(vals[1]),
            "y": float(vals[2])
        }
    return shifts


def get_metadata_df(fov, boundary_dict):
    metadata_df = pd.DataFrame(columns=METADATA_FILE_COLNAMES)
    for k in boundary_dict.keys():
        v = boundary_dict[k]
        x, y = np.average(v, axis=0)
        metadata_df.loc[len(metadata_df.index)] = [int(fov), k, y, -x + FOV_SIZE_X, y, -x + FOV_SIZE_X]
    metadata_df[FOV]=metadata_df[FOV].astype(int)
    metadata_df[CELL_ID]=metadata_df[CELL_ID].astype(int)
    return metadata_df


def get_polygon_vertices(fov, expanded, num_vertices=30):
    """
    Loop over each pixel in boundaries. Boundaries has the same dimensions as expanded
    For every non-zero entry in boundaries, grab the corresponding pixel in expanded.
    Add the expanded pixel to a list describing the edges for each segmentation
    """

    # find the boundaries of each segmentation
    boundaries = find_boundaries(expanded, mode='inner')

    # create empty lists for each cell to store segmentations vertices
    u = np.unique(expanded)
    boundary_dict = {}
    for v in u[1:]:
        boundary_dict[v] = []
        
    assert expanded.shape == boundaries.shape
    for r in range(expanded.shape[0]):
        for c in range(expanded.shape[1]):
            v = boundaries[r, c]
            if v > 0:
                cell = expanded[r, c]
                boundary_dict[cell].append((r, c))

    # get metadata dataframe before downsampling the numer of vertices (since the downsampling is random its possible it could occasionally lead to poorly chosen cell center positions)
    metadata_df = get_metadata_df(fov=fov, boundary_dict=boundary_dict)

    # randomly sample vertices if cell has more than num_vertices vertices describing it
    for v in boundary_dict.keys():
        d = boundary_dict[v]
        if len(d) < 4: # monitor if there are any segmentations with an oddly low number of vertices
            print(v)
        elif len(d) > num_vertices:
            boundary_dict[v] = random.sample(d, num_vertices)

    # TODO: there's gotta be a much faster less hacky way to properly order segmentation vertices and it should be much faster
    segmentations_df = pd.DataFrame(columns=POLYGONS_FILE_COLNAMES)
    for k in boundary_dict.keys():
        if k % 500 == 0:
            print(f"\t\tOrdering vertices for cell {k} in fov {fov}")
        series_list = []
        for x, y in boundary_dict[k]:
            # this is where we need to sort the vertices
            # x and y are intentionally flipped here
            series_list.append(pd.Series([fov, k, y, -x + 3638, y, -x + 3638], index=POLYGONS_FILE_COLNAMES))
            
        cur_seg_df = pd.DataFrame(series_list)
        
        pts = cur_seg_df[[POLYGONS_X_LOCAL, POLYGONS_Y_LOCAL]]
        neighbors = np.zeros((pts.shape[0], pts.shape[0]))
        for idx1, row1 in pts.iterrows():
            for idx2, row2 in pts.iterrows():
                if idx1 == idx2:
                    neighbors[idx1, idx2] = 100000
                else:
                    d = np.sqrt((row1[POLYGONS_X_LOCAL] - row2[POLYGONS_X_LOCAL])**2 + (row1[POLYGONS_Y_LOCAL] - row2[POLYGONS_Y_LOCAL])**2)
                    neighbors[idx1, idx2] = d

        order = [0]
        mask = [100000] + [1] * (pts.shape[0] - 1)
        current = 0
        for _ in range(pts.shape[0] - 1):
            val, idx = min((val, idx) for (idx, val) in enumerate(neighbors[current, :] * mask))
            mask[idx] *= 10000
            order.append(idx)
            current = idx

        cur_seg_df = cur_seg_df.reindex(order)
        segmentations_df = pd.concat([segmentations_df, cur_seg_df])
    segmentations_df.columns = POLYGONS_FILE_COLNAMES

    return {
       "segmentations_df": segmentations_df,
       "metadata_df": metadata_df
    }


if __name__ == "__main__":
    # test code
    convert_masks(
        "/brahms/hartmana/spatial_vignette_data/nanostring/lung5_rep1",
        "/brahms/hartmana/mesmer/lung5_rep1/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/MesmerSegmentation",
        "test_outs",
        "/brahms/hartmana/mesmer/testing/test_output_directory",
        fovs=[1,3],
    )
