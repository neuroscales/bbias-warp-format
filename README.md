# connects-warp-format

* Status of this document: alpha
* Editor: Yael Balbastre <y.balbastre at ucl.ac.uk>
* Version: 1.0a

## Abstract

A scalable format for large displacement or coordinates fields ("warps").

## Table of Content

0. [References](#0-references)
1. [Introduction](#1-introduction)
2. [Format Specification](#2-format-specification)
3. [Main differences with NIfTI and/or OME-NGFF](#3-main-differences-with-nifti-andor-ome-ngff)
4. [Conversion tables](#4-conversion-tables)
5. [Reference implementations](#5-reference-implementations)

## 0. References

* [__Zarr__](https://zarr.readthedocs.io) is a format for the storage of
  chunked, compressed, N-dimensional arrays inspired by HDF5, h5py and bcolz.
* [__OME-NGFF__](https://ngff.openmicroscopy.org) (Next Generation File
  Format) is a format based on zarr for the storage of biomedical imaging data.
  It's version 0.6 specifies chains of geometric transformations, including coordinates
  and displacement fields.

## 1. Introduction

The CONNECTS consortium generates and analyses peta-bytes-scale images of the human brain.
Some of the analysis steps involve estimating and applying geometric transformations to these images.
Many transformation formats exist in the neuroimaging world, and none seems to have been accepted as 
a standard by the community. Furthermore, no format seems to scale to TB- ro PB-scale datasets.

The Big Brain Imaging Analysis and Standardization Working Group (BBIAS), within the CONNECTS consortium
is tasked with drafting and implementing a scalable format for large non linear transformations ("warps"),
to be adopted by the consortium and the wider community.

Because it also targets large images and does specify a transformation format, this draft aligns as best as possible
with the OME-NGFF v0.6 specification. However, it does highlight limitations and propose extensions to this format.

## 2. Existing warp formats

### NIfTI-based warps

The [NIfTI-1 format](https://github.com/NIFTI-Imaging/nifti_clib/blob/master/niftilib/nifti1.h) does list two intent
codes related to coordinates and displacements fields:

```C
#define NIFTI_INTENT_DISPVECT  1006   /* specifically for displacements */
#define NIFTI_INTENT_VECTOR    1007   /* for any other type of vector */
```

