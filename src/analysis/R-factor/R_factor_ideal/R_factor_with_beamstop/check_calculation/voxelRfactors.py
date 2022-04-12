#!/usr/bin/env python

import numpy as np
from tqdm import tqdm
import scipy.io
import click
import h5py
from pathlib import Path


def add_beamstop(array, radius=8):
    center = int((array.shape[0] - 1) / 2)
    i_pixel = np.indices(array.shape)
    r_pixel = np.sqrt(
        (i_pixel[0] - center) ** 2
        + (i_pixel[1] - center) ** 2
        + (i_pixel[2] - center) ** 2
    )
    array[r_pixel <= radius] = 0


def extractGroup(fn):
    with h5py.File(fn, "r") as h5:
        for grp in h5:
            if "output" in grp:
                return grp


def get_mat_array(fn, grp="output_ori_corr"):
    with h5py.File(fn, "r") as h5:
        array = h5[grp][()]
    return array


def crop_array(array, size=550):
    center = int((array.shape[0] - 1) / 2)
    offset = int((size-1) // 2)
    array_new = array[
        center - offset : center + offset,
        center - offset : center + offset,
        center - offset : center + offset,
    ]
    return array_new


def get_Ds(q_map):
    qs = np.linspace(q_map.min(), q_map.max(), q_map.shape[0])
    return 1 / qs  # Angstrom


def R_r_vol(r_val, img, img_ref, r_map):
    """Get R factor of a specific r for an ndarray.

    Args:
        r_val: r value in pixel unit.
        img: The array to be inspected.
        img_ref: The reference array
        r_map: The map of r

    Returns:
        The R_factor of a specific r_val
    """

    N_real = img
    N_ideal = img_ref

    q_range = r_map <= r_val
    # print(r_val, np.count_nonzero(q_range))

    N_real_sub = N_real[q_range]
    N_real_sub[N_real_sub < 0] = 0
    N_ideal_sub = N_ideal[q_range]
    N_ideal_sub[N_ideal_sub < 0] = 0

    N_d = np.sum(np.sqrt(N_real_sub))
    N_ideal_d = np.sum(np.sqrt(N_ideal_sub))
    print(r_val, N_d, N_ideal_d)

    # import pdb; pdb.set_trace()
    eqs = abs(np.sqrt(N_real_sub) / N_d - np.sqrt(N_ideal_sub) / N_ideal_d)

    return np.sum(eqs)


def R_d_vol(d, img, img_ref, q_map):
    """Get R factor of a specific d for a reciprocal volume.

    Args:
        d: resolution (Angstrom).
        img: q list of the pattern to be analyzed. (Angstrom^-1).
        img_ref: The pattern to be analyzed.
        q_map: The map of |q| in the unit of 1/Angstrom.

    Returns:
        The R_factor of a specific d
    """

    N_real = img
    N_ideal = img_ref

    q_range = q_map <= 1 / d

    N_real_sub = N_real[q_range]
    N_real_sub[N_real_sub < 0] = 0
    N_ideal_sub = N_ideal[q_range]
    N_ideal_sub[N_ideal_sub < 0] = 0

    N_d = np.sum(np.sqrt(N_real_sub))
    N_ideal_d = np.sum(np.sqrt(N_ideal_sub))

    eqs = abs(np.sqrt(N_real_sub) / N_d - np.sqrt(N_ideal_sub) / N_ideal_d)

    return np.sum(eqs)


def get_r_map(img):
    r_idx = np.indices(img.shape)
    center = np.array(img.shape) // 2
    r_map2 = 0
    for i, indices in enumerate(r_idx):
        r_map2 += (indices - center[i]) ** 2
    return np.sqrt(r_map2)


def R_factor_vol(
    img,
    img_ref,
):
    """Get residual factor for a reciprocal space volume.

    Args:
        img: Diffraction pattern array
        img_ref: Referece diffraction pattern array

    Returns:
        The residual factor of each pixel
    """

    assert img.shape == img_ref.shape
    print("Data shape:", img.shape)
    r_map = get_r_map(img)
    r_range = range(int(r_map.max()))
    R_factors = []
    for r_val in r_range:
        R_factors.append(R_r_vol(r_val, img, img_ref, r_map))

    return R_factors


def get_q_map(crop_size):
    # Get q map
    ijk_map = np.indices((crop_size, crop_size, crop_size))
    map_pixel_unit = ijk_map + 0.5
    voxel_dis_map = np.linalg.norm(map_pixel_unit + 0.5 - crop_size / 2, axis=0)
    q_voxel_size = 0.000307438  # in 1/Angstrom
    q_map = voxel_dis_map * q_voxel_size
    return q_map


def get_ref(fn):
    """Get reference volume"""
    scaled_respace = scipy.io.loadmat(fn)
    volume_ref = scaled_respace["output_ori_corr_clip"]
    return volume_ref


def analyze(fn, volume_ref, crop_size):

    sample = get_mat_array(fn, extractGroup(fn))
    sample = crop_array(sample, crop_size)
    sample = sample.transpose()
    add_beamstop(sample, 10)

    # Run for the analysis for a run
    print(f"Analyzing {fn}")
    # R_factors = R_factor_vol(volume, volume_ref, voxel_dis_map)
    return R_factor_vol(sample, volume_ref)


@click.command()
@click.argument("input", type=str)
@click.argument("output", type=str)
def main(input, output):
    print(f"input: {input}")
    print(f"output: {output}")
    scaled_respace_path = "/beegfs/desy/user/kimchan/SPI_Dragonfly_AGIPD_noise/original_3d_dt/GPU_RV500_new_OC_result.mat"
    volume_ref = get_ref(scaled_respace_path)
    crop_size = volume_ref.shape[0]
    # crop_size = 200
    print("crop_size:", crop_size)
    # assert crop_size == 444
    volume_ref = crop_array(volume_ref, crop_size)
    add_beamstop(volume_ref, 10)
    r_map = get_r_map(volume_ref)
    r_range = range(int(r_map.max()))

    with h5py.File(output, "w") as h5:
        R_factors = analyze(input, volume_ref, crop_size)
        h5["R_factors"] = R_factors
        h5["resolutions"] = r_range


if __name__ == "__main__":
    main()
