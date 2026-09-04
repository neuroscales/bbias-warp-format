# connects-warp-format

* Status of this document: alpha
* Editor: Yael Balbastre <y.balbastre at ucl.ac.uk>
* Version: 1.0a

## Abstract

A scalable format for large displacement or coordinates fields ("warps").

## Table of Content

0. [References](#0-references)
1. [Introduction](#1-introduction)
2. [Existing Warp Formats](#2-existing-warp-formats)
3. [Format Specification](#3-format-specification)
4. [Limitations](#4-limitations)

## 0. References

* [__Zarr__](https://zarr.readthedocs.io) is a format for the storage of
  chunked, compressed, N-dimensional arrays inspired by HDF5, h5py and bcolz.
* [__OME-NGFF__](https://ngff.openmicroscopy.org) (Next Generation File
  Format) is a format based on zarr for the storage of biomedical imaging data.
  It's version 0.6 specifies chains of geometric transformations, including coordinates
  and displacement fields.
* [__ITK__](https://itk.org/) is a C++ library for biomedical imaging processing, which
  implements its own set of transformation formats and algorithms. It is used by popular
  biomedical software, such as ANTs or 3DSlicer, making its transformation formats ubiquitous
  in the field.
* [__FSL__](https://fsl.fmrib.ox.ac.uk/) is a package for the processing and analysis of
  neuroimaging data. It defines its own set of affine and nonlinear transformation formats.
  It is heavily used by the CMC project within CONNECTS.

## 1. Introduction

The NIH-funded BRAIN CONNECTS consortium generates and analyses peta-bytes-scale images of the human brain.
Some of the analysis steps involve estimating and applying geometric transformations to these images.
Many transformation formats exist in the neuroimaging world, and none seems to have been accepted as 
a standard by the community. Furthermore, no format seems to scale to TB- or PB-scale datasets.

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
[ x ]   [ x ]   [ Axx Axy Axz ]   [ Sx  0  0 ] ^(-1)                          ( [ Axx Axy Axz ]   [ x ]   [ Tx ] )   [ Tx ]
[ y ] = [ y ] + [ Ayx Ayy Ayz ] x [  0 Sy  0 ]       x ScaledDisplacementField( [ Ayx Ayy Ayz ] x [ y ] + [ Ty ] ) + [ Ty ]
[ z ]   [ z ]   [ Azx Azy Azz ]   [  0  0 Sz ]                                ( [ Azx Azy Azz ]   [ z ]   [ Tz ] )   [ Tz ]
  v       v            v                v                                              v            v       v          v
 RAS     RAS      VOX2RAS(Lin)    Inverse Scale                                   RAS2VOX(Lin)     RAS  RAS2VOX(Tr)  VOX2RAS(Tr)
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

The only interpolation methods currently supported by the specification are 
`"nearest"`, `"linear"`, `"bspline-cubic"`, and `"windowed sinc"`.

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
                "scale": [1.0, sz, sy, sx],
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

Multiple (and a least one) `coordinateSystems` can be defined in the same dataset.
The **first** system in the list is the "intrinsic" system, which is advised to be a scaled voxel space
(eventually translated -- but only to account for shifts across resolution levels). The **last** system in the list
is the default model space (e.g. used for visualisation by neuroglancer). Together with all `coordinateTransformations`
(at the within-dataset level, and the across-datasets level), they define a transformation from the 
field's voxel space to the model space. This transformation **must be invertible**.

Note the the fact that the model space is named (and interpreted as) "ras+" is a choice we've made here, not a requirement 
from the specification. Assuming that the the model space is indeed RAS, and that the implied "voxel2tas" transformation is affine, 
an input RAS coordinate is transformed by an OME transform via:

```text
# displacement

[ z ]   [ z ]                    ( [ Azz Azy Azx ]   [ z ]   [ Tz ] )
[ y ] = [ y ] + DisplacementField( [ Ayz Ayy Ayx ] x [ y ] + [ Ty ] )
[ x ]   [ x ]                    ( [ Axz Axy Axx ]   [ x ]   [ Tx ] )
  v       v                               v            v       v       
 RAS     RAS                         RAS2VOX(Lin)     RAS  RAS2VOX(Tr)

# coordinates

[ z ]                   ( [ Azz Azy Azx ]   [ z ]   [ Tz ] )
[ y ] = CoordinatesField( [ Ayz Ayy Ayx ] x [ y ] + [ Ty ] )
[ x ]                   ( [ Axz Axy Axx ]   [ x ]   [ Tx ] )
  v                              v            v       v       
 RAS                        RAS2VOX(Lin)     RAS  RAS2VOX(Tr)
```

> [!WARNING]
> Array dimensions in an OME-NGFF dataset **must** be ordered ([t], [c], [z], y, x).
> Axes within a coordinates system must also have this order. In other words, dimensions
> are transposed compared to conventions used in NIfTI and other neuroimaging tools.
> This means that the components of the field are also ordered as {z, y, x} (the first channel
> of the displacement field contains the displacement along the z axis).
>
> We could either decide that (R, A, S) == (z, y, x) so that the input and output coordinates
> in the equation above remain "RAS", or keep (R, A, S) == (x, y, z) and learn to reverse
> our RAS vectors when computing the transformations. I lean towards the second solution.

> [!NOTE]
> Although we have supposed here that the "voxel to model" transformation is affine, nothing mandates
> it in the OME-NGFF spec. The list of coordinates transformation in the displacement field's metadata could very
> well contain a (invertible) displacement field as well. The more general transformation equation becomes
> 
> ```
> output_model_coord = CoordinatesField( inv(voxel_to_model)(input_model_coord) )
> ```

## 3. Format Specification

This document proposes to leverage OME-NGFF displacements and coordinates fields for the CONNECTS common warp format, with
the additional requirement that the output model space is (R, A, S) == (x, y, z). It **recommends** that the OME-NGFF 
[RFC 4](https://ngff.openmicroscopy.org/rfc/4/index.html) be used to explictely encode this convention. It **does not 
recommend** that a specific unit be used (e.g. "mmRAS"), meaning that downstream users **must** check and adapt units
when composing transformations that act on RAS coordinates spaces with different units.

Consequently, CONNECTS displacement and coordinates fields will be multiscales OME-NGFF datasets with metadata of the form:

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
              {
                "name": "z", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "inferior-to-superior"}
              },
              {
                "name": "y", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "posterior-to-anterior"}
              },
              {
                "name": "x", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "left-to-right"}
              }
            ]
          }
        ],
        "coordinateTransformations": [
          {
            "name": "phys2ras",
            "input": "physical",
            "output": "ras+",
            "type": "byDimension",
            "transformations": [
              {
                "inputAxes": [1, 2, 3],
                "outputAxes": [1, 2, 3],
                "type": "affine",
                "matrix": [
                  [Azz, Azy, Azx, Tz],
                  [Ayz, Ayy, Ayx, Ty],
                  [Axz, Axy, Axx, Tx]
                ]
              }
            ]
          }
        ],
        "datasets": [
          {
            "path": "s0",
            "coordinateTransformations": [
              {
                "type": "scale",
                "scale": [1.0, sz, sy, sx],
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

This is fully compatible with the OME-NGFF v0.6 specification, with the constraint that 
the model space is RAS+, and the voxel-to-RAS transformation is (at most) an affine transformation.

### Multi-resolution extension

Displacement and coordinates transformations as defined in the v0.6 specification assume that 
the fields themselves are rather coarse, since a single resolution level is mandated and used.
However, because they are full-fledged multiscales datasets, they easily lend themselves to multi-resolution
extensions. Multi-resolution displacement fields may be needed when the algorithm that estimates them refines them 
at increasingly higher resolution. In such cases, the resoluting displacement fields can reach a resolution on par
with that of the images that have been collected. Applying such large fields is computationnaly demanding; however, 
during visualization, only contant-size chunks of the deformed image ever need to be calculated, albeit at different 
resolution levels. The availability of different resolution levels of the same displacement field allows the visualization
software to pick the transformation scale most suited to its current viewport, thereby greatly minimizing its memory
footprint.

Such a multi-resolution displacement fields would have OME metadata of the form


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
              {
                "name": "z", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "inferior-to-superior"}
              },
              {
                "name": "y", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "posterior-to-anterior"}
              },
              {
                "name": "x", "type": "space", "unit": "micrometer", 
                "orientation": {"type": "anatomical", "value": "left-to-right"}
              }
            ]
          }
        ],
        "coordinateTransformations": [
          {
            "name": "phys2ras",
            "input": "physical",
            "output": "ras+",
            "type": "byDimension",
            "transformations": [
              {
                "inputAxes": [1, 2, 3],
                "outputAxes": [1, 2, 3],
                "type": "affine",
                "matrix": [
                  [Azz, Azy, Azx, Tz],
                  [Ayz, Ayy, Ayx, Ty],
                  [Axz, Axy, Axx, Tx]
                ]
              }
            ]
          }
        ],
        "datasets": [
          {
            "path": "s0",
            "coordinateTransformations": [
              {
                "input" : {"path": "s0"},
                "output" : {"name": "physical"}
                "type": "scale",
                "scale": [1.0, sz, sy, sx],
              }
            ]
          },
          {
            "path": "s1",
            "coordinateTransformations": [
              {
                "input" : {"path": "s1"},
                "output" : {"name": "physical"}
                "type": "sequence",
                "transformations": [
                  {
                    "type": "scale",
                    "scale": [1.0, sz*2, sy*2, sx*2],
                  },
                  {
                    "type": "translation",
                    "scale": [0.0, sz/2, sy/2, sx/2],
                  },
                ]
              }
            ]
          },
          ...

        ]
      }
    ]
  }
}
```

Note that because displacements and coordinates are defined in terms of model space (RAS) coordinates, 
downsampled levels can be computed using simple moving-average or gaussian-pyramid downsampling. The
downsampled displacements do not need to be re-scaled.

## 4. Limitations

- Boundary conditions when interpolating a warp are not specified in the specification. This matters particularly for cubic splines.
- Some of the interpolation orders used by FSL are not available in OME-NGFF (quadratic spline, DCT), which means that fields
  that use these orders cannot be converted losslessly to OME-NGFF.
- Furthermore, displacement [coordinates] fields must be saved in terms of "model space" displacements [coordinates], which again differs
  from FSL's "voxel space" displacements [coordinates]. This will also lead to numerical differences when interpolating the field.
