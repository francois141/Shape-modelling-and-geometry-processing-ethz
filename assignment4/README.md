# Assignment 4

Name: François Costa

Legi-Nr: 19-931-989

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Mandatory Tasks

1) Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off, Octo_cut2.obj)

2) Several examples of the distortion visualizations.


## Reports

## The code was compiled with cpp17!

This line needs to be added in the cmake file.

```cmake
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
```

### (mandatory-1) parameterization and checkerboard texture models
#### cathead
| Method            | checkerboard textured models          | Parameterization                                     |
| :--------------:  | ------------------------------------- |------------------------------------------------------|
| Uniform (fixed)   |<img align="center" src="./images/cathead_1.png" >| <img align="center"  src="./images/cat_flat_1.png" > |
| Cotangent (fixed) |<img align="center" src="./images/cathead_2.png" >| <img align="center"  src="./images/cat_flat_2.png" > |
| LSCM (fixed)      |<img align="center" src="./images/cathead_3.png" >| <img align="center"  src="./images/cat_flat_3.png" > |
| ARAP (fixed)      |<img align="center" src="./images/cathead_4.png" >| <img align="center"  src="./images/cat_flat_4.png" > |
| LSCM (free)       |<img align="center" src="./images/cathead_3_free.png" >| <img align="center"  src="./images/cat_flat_5.png" > |
| ARAP (free)       |<img align="center" src="./images/cathead_4_free.png" >| <img align="center"  src="./images/cat_flat_6.png" > |

#### hemisphere
| Method            | checkerboard textured models          | Parameterization                                        |
| :--------------:  | ------------------------------------- |---------------------------------------------------------|
| Uniform (fixed)   |<img align="center" src="./images/hemisphere_1.png" >| <img align="center"  src="./images/sphere_flat_1.png" > |
| Cotangent (fixed) |<img align="center" src="./images/hemisphere_2.png" >| <img align="center"  src="./images/sphere_flat_2.png" > |
| LSCM (fixed)      |<img align="center" src="./images/hemisphere_3.png" >| <img align="center"  src="./images/sphere_flat_3.png" > |
| ARAP (fixed)      |<img align="center" src="./images/hemisphere_4.png" >| <img align="center"  src="./images/sphere_flat_4.png" > |
| LSCM (free)       |<img align="center" src="./images/hemisphere_3_free.png" >| <img align="center"  src="./images/sphere_flat_5.png" > |
| ARAP (free)       |<img align="center" src="./images/hemisphere_4_free.png" >| <img align="center"  src="./images/sphere_flat_6.png" > |


#### hemisphere_non_convex_boundary
| Method            | checkerboard textured models          | Parameterization                                                   |
| :--------------:  | ------------------------------------- |--------------------------------------------------------------------|
| Uniform (fixed)   |<img align="center" src="./images/hemisphere_non_convex_1.png" >| <img align="center"  src="./images/sphere_non_convex_flat_1.png" > |
| Cotangent (fixed) |<img align="center" src="./images/hemisphere_non_convex_2.png" >| <img align="center"  src="./images/sphere_non_convex_flat_2.png" > |
| LSCM (fixed)      |<img align="center" src="./images/hemisphere_non_convex_3.png" >| <img align="center"  src="./images/sphere_non_convex_flat_3.png" > |
| ARAP (fixed)      |<img align="center" src="./images/hemisphere_non_convex_4.png" >| <img align="center"  src="./images/sphere_non_convex_flat_4.png" > |
| LSCM (free)       |<img align="center" src="./images/hemisphere_non_convex_3_free.png" >| <img align="center"  src="./images/sphere_non_convex_flat_5.png" > |
| ARAP (free)       |<img align="center" src="./images/hemisphere_non_convex_4_free.png" >| <img align="center"  src="./images/sphere_non_convex_flat_6.png" > |

#### Octo_cut2
| Method            | checkerboard textured models          | Parameterization                                      |
| :--------------:  | ------------------------------------- |-------------------------------------------------------|
| Uniform (fixed)   |<img align="center" src="./images/octocut_1.png" >| <img align="center"  src="./images/octo_flat_1.png" > |
| Cotangent (fixed) |<img align="center" src="./images/octocut_2.png" >| <img align="center"  src="./images/octo_flat_2.png"  > |
| LSCM (fixed)      |<img align="center" src="./images/octocut_3.png" >| <img align="center"  src="./images/octo_flat_3.png"  > |
| ARAP (fixed)      |<img align="center" src="./images/octocut_4.png" >| <img align="center"  src="./images/octo_flat_4.png"  > |
| LSCM (free)       |<img align="center" src="./images/octocut_3_free.png" >| <img align="center"  src="./images/octo_flat_5.png" > |
| ARAP (free)       |<img align="center" src="./images/octocut_4_free.png" >| <img align="center"  src="./images/octo_flat_6.png" > |



### (mandatory-2) distortion visualization
#### cathead
|     mtd \ metric     | Conformal (angle)                                                     | Authalic (area)                                                       | Isometric  (length)                                                   |
|:--------------------:|-----------------------------------------------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------------------------|
| LSCM (free) - View 1 | <img align="center" src="./images/cat_dist_0_view_1.png" >            | <img align="center" src="./images/cat_dist_2_view_1.png" >            | <img align="center" src="./images/cat_dist_1_view_1.png" >            |
| LSCM (free) - View 2 | <img align="center" src="./images/cat_dist_0_view_2.png" >            | <img align="center" src="./images/cat_dist_2_view_2.png" >            | <img align="center" src="./images/cat_dist_1_view_2.png" >            |
| ARAP (free) - View 1 | <img align="center" src="./images/cat_arap_distortion_0_view_1.png" > | <img align="center" src="./images/cat_arap_distortion_2_view_1.png" > | <img align="center" src="./images/cat_arap_distortion_1_view_1.png" > |
| ARAP (free) - View 2 | <img align="center" src="./images/cat_arap_distortion_0_view_2.png" > | <img align="center" src="./images/cat_arap_distortion_2_view_2.png" > | <img align="center" src="./images/cat_arap_distortion_1_view_2.png" > |


#### hemisphere

| mtd \ metric      | Conformal (angle) |    Authalic (area)  |  Isometric  (length)    |
| :--------------:  | ----------------- | ------------------- | ----------------------- |
| LSCM (free)       | <img align="center" src="./images/hemisphere_dist_0.png" >     | <img align="center" src="./images/hemisphere_dist_2.png" >     | <img align="center" src="./images/hemisphere_dist_1.png" >     | 
| ARAP (free)       | <img align="center" src="./images/hemisphere_arap_dist_0.png" >| <img align="center" src="./images/hemisphere_arap_dist_2.png" >|<img align="center" src="./images/hemisphere_arap_dist_1.png" >|


