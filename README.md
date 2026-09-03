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

[ITK](https://itk.org/) uses these intent codes when storing either displacements (`1006`) or coordinates (`1007`) fields, 
as do all software built in top of ITK ([Slicer](https://www.slicer.org/), [ANTs](https://github.com/antsx/ants), etc).
ITK stores the coordinates or displacements in the **LPS** model space (x=right-to-left, y=anterior-to-posterior, z=inferior-to-superior), 
with the NIfTI's affine used to encode a mapping from the vector field voxel grid to the LPS model space.

A RAS coordinate is transformed by an ITK transform via:

```text
# coordinates field
[ x ]    [ -1  0  0 ]                      ( [ Axx Axy Axz ]   [ x ]   [ Tx ] )
[ y ] =  [  0 -1  0 ] x LPSCoordinatesField( [ Ayx Ayy Ayz ] x [ y ] + [ Ty ] )
[ z ]    [  0  0 +1 ]                      ( [ Azx Azy Azz ]   [ z ]   [ Tz ] )
  v            v                                    v            v       v
 RAS        LPS2RAS                            RAS2VOX(Lin)     RAS    RAS2VOX(Tr)

# displacement field
[ x ]   [ x ]   [ -1  0  0 ]                       ( [ Axx Axy Axz ]   [ x ]   [ Tx ] )
[ y ] = [ y ] + [  0 -1  0 ] x LPSDisplacementField( [ Ayx Ayy Ayz ] x [ y ] + [ Ty ] )
[ z ]   [ z ]   [  0  0 +1 ]                       ( [ Azx Azy Azz ]   [ z ]   [ Tz ] )
  v       v           v                                     v            v       v
 RAS     RAS       LPS2RAS                             RAS2VOX(Lin)     RAS    RAS2VOX(Tr)
```

[SPM](https://www.fil.ion.ucl.ac.uk/spm/) also uses NIfTI files to store its warps, although it does not set the
corresponding vector codes, and in contrast with ITK, it saves its coordinates in terms of the **RAS** model space
(x=left-to-right, y=posterior-to-anterior, z=inferior-to-superior).

A RAS coordinate is transformed by an SPM transform via:

```text
[ x ]                      ( [ Axx Axy Axz ]   [ x ]   [ Tx ] )
[ y ] = RASCoordinatesField( [ Ayx Ayy Ayz ] x [ y ] + [ Ty ] )
[ z ]                      ( [ Azx Azy Azz ]   [ z ]   [ Tz ] )
  v                                 v            v       v
 RAS                           RAS2VOX(Lin)     RAS    RAS2VOX(Tr)
```

[FSL](https://fsl.fmrib.ox.ac.uk/) possesses its own NIfTI intent codes:

```C
#define NIFTI_INTENT_FSL_FNIRT_DISPLACEMENT_FIELD       2006
#define NIFTI_INTENT_FSL_CUBIC_SPLINE_COEFFICIENTS      2007
#define NIFTI_INTENT_FSL_DCT_COEFFICIENTS               2008
#define NIFTI_INTENT_FSL_QUADRATIC_SPLINE_COEFFICIENTS  2009
```

FSL always stores displacements (not absolute coordinates), and these displacements are in terms of
"scaled voxels", not in terms of an abstract model space.

A RAS coordinate is transformed by an FSL transform via:

```text
[ x ]   [ Axx Axy Axz ]   [ Sx  0  0 ] ^(-1)                          ( [ Axx Axy Axz ]   [ x ]   [ Tx ] )   [ Tx ]
[ y ] = [ Ayx Ayy Ayz ] x [  0 Sy  0 ]       x ScaledDisplacementField( [ Ayx Ayy Ayz ] x [ y ] + [ Ty ] ) + [ Ty ]
[ z ]   [ Azx Azy Azz ]   [  0  0 Sz ]                                ( [ Azx Azy Azz ]   [ z ]   [ Tz ] )   [ Tz ]
  v           v                 v                                              v            v       v          v
 RAS      VOX2RAS(Lin)    Inverse Scale                                   RAS2VOX(Lin)     RAS  RAS2VOX(Tr)  VOX2RAS(Tr)
```

### OME-NGFF v0.6

The v0.6 (currently v0.6rc0) of OME-NGFF introduces coordinates spaces and coordinates transformations. It defines
a large range of transformation types, two of them being qualified as "warps": coordinates fields and displacement fields.

Because OME-NGFF does not commit to a specific organ, it stays clear of any mandatory model space (e.g. RAS vs LPS vs voxels). 
It is therefore the user's  (or data creator) job to ensure the consistency of its spaces. 
In general, the common space of two transformations that compose are always exactly identical. 
This is slightly more complicated for displacement and coordinates field, because the input space is 
a voxel grid and, especially in the case of displacement fields, can lead to inconsistencies. 
Cleaning up these inconsistencies has  been one of the main working points during the later stages of the 0.6 development.

When applied in a chain of transformations, i.e., within the `"coordinatesTransformations"` of an OME multiscales dataset, it has the form:

```yaml
# coordinates
{
  "type": "coordinates",
  "path": "s3://bucket/path/to/coords.ome.zarr",
  "interpolation": "bspline-cubic"
},
# displacements
{
  "type": "displacements",
  "path": "s3://bucket/path/to/disp.ome.zarr",
  "interpolation": "bspline-cubic"
}
```

The coordinates or displacements pointed to by `"path"` are full-fledged multiscales datasets, with metadata:

```yaml
{
  "ome": {
    "version": "0.6rc0",
    "name": "displacements",
    "multiscales": [
      {
        "coordinateSystems": [
          {
            "name": "physical",
            "axes": [
              {"name": "c", "type": "displacement", "discrete": true},
              {"name":"z", "type": "space", "unit": "micrometer"},
              {"name":"y", "type": "space", "unit": "micrometer"},
              {"name":"x", "type": "space", "unit": "micrometer"}
            ]
          },
          {
            "name": "ras+",
            "axes": [
              {"name": "c", "type": "displacement", "discrete": true},
              {"name":"z", "type": "space", "unit": "micrometer"},
              {"name":"y", "type": "space", "unit": "micrometer"},
              {"name":"x", "type": "space", "unit": "micrometer"}
            ]
          }
        ],
        "coordinateTransformations": [
          {
            "name": "phys2ras",
            "input": "physical",
            "output": "ras+",
            "type": "affine",
            "matrix": [ ... ]
          }
        ],
        "datasets": [
          {
            "path": "s0",
            "coordinateTransformations": [
              {
                "type": "scale",
                "scale": [2.0, 2.0],
                "input" : {"path": "s0"},
                "output" : {"name": "physical"}
              }
            ]
          }
        ]
      }
    ]
  }
}
```
