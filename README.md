# Evaluate Lie group functions safely in Fortran
This project contains several modules that allow easy and safe evaluation of Lie group functions related to the following Lie groups of unit quaternions ๐ยณ and its semidirect product ๐ยณโโยณ. In order to do so, there is a module that evaluates singular functions, minimizing round-off errors. Additionally, some functionality of quaternions are included.

## Contents

### Module `cross_functions`
Contains just three functions related to cross products:

 * `cross(v,w)`: Calculates the cross product ๐ฃร๐ค.
 * `skw(v)`: Calculates the skew-symmetric matrix skw(๐ฃ) such that skw(๐ฃ)โ๐ค=๐ฃร๐ค.
 * `dynadic_product(x,y)`: Calculates the dyadic product ๐ฅโ๐ฆแต of column vectors ๐ฅ, ๐ฆ of arbitrary (but equal) length.

### Module `quaternion_functions`
Contains functions related to quaternions in general. Quaternions are represented as 4D vectors.

 * `qp(p,q)`: Calculates the quaternion product ๐*๐.
 * `apply_quat(p,v)`: Applies the rotation represented by a (unit) quaternion ๐ to a 3D vector ๐ฃ.
 * `apply_conj_quat(p,v)`: Applies the inverse rotation represented by a (unit) quaternion ๐ to a 3D vector ๐ฃ. The name refers to the fact that the inverse of a unit quaternion is its conjugate.
 * `conj_quat(p)`: Calculate the conjugate of a quaternion.

### Module `s3_functions`
Contains functions related to the Lie group of unit quaternions ๐ยณ:

 * `lp_s3(p1,p2)`: Calculates the Lie product ๐โ*๐โ of two unit quaternions. This is just the quaternion product.
 * `inv_s3(p)`: Calculates the inverse of a unit quaternion. This is just the conjugate.
 * `lie_bracket_s3(v1,v2)`: Calculates the result of ฬ๐ฃโโ๐ฃโ with the hat operator. The name comes from the fact that the hat operator is defined by the Lie bracket of the Lie algebra ๐ฐยณ and the tilde operator chosen for ๐ยณ. Also, the mapping (๐ฃโ,๐ฃโ)โฆฬ๐ฃโโ๐ฃโ fulfills the definition of a Lie bracket on โยณ.
 * `expt_s3(v)`: Calculates the exponential of the tilde operator expt(๐ฃ)=exp(ฬ๐ฃ).
 * `tan_op_s3(v)`: Calculates the tangent operator ๐(๐ฃ).
 * `tan_tr_mult_s3(v,w)`: Calculates the result of ๐แต(๐ฃ)โ๐ค.
 * `tan_op_inv_s3(v)`: Calculates the inverse matrix ๐โปยน(๐ฃ)=(๐(๐ฃ))โปยน of the tangent operator.
 * `tan_tr_inv_s3(v)`: Calculates ๐โปแต(๐ฃ)=((๐(๐ฃ))แต)โปยน, the inverse matrix of the transposed tangent operator.
 * `tan_op_inv_mult_s3(v,w)`: Calculates ๐โปยน(๐ฃ)โ๐ค.
 * `tan_tr_inv_mult_s3(v,w)`: Calculates ๐โปแต(๐ฃ)โ๐ค.
 * `logt_s3(p)`: Calculates the logarithm logt(๐) with results in โยณ of a quaternion ๐ sufficiently close to the identity quaternion. Note that it holds tilde(logt(๐))=log(๐).
 * `d_tan_tr_inv_s3(Om, A)`: Calculates the Jacobi matrix of ๐โปแต(ฮฉ)โ๐ด with respect to the vector ฮฉโโยณ, where ๐ดโโยณ.

It is recommended to use the functions with `_mult` in the function name whereever possible for increased speed.

### Module `s3sdr3_functions`
Contains functions related to the semidirect product Lie group ๐ยณโโยณ:

 * `lp_s3sdr3(q1,q2)`: Calculates the Lie product ๐โ*๐โ for ๐โ,๐โโ๐ยณโโยณ.
 * `inv_s3sdr3(q)`: Calculates the inverse of an element ๐โ๐ยณโโยณ.
 * `lie_bracket_s3(v1,v2)`: Calculates the result of ฬ๐ฃโโ๐ฃโ with the hat operator and ๐ฃโ,๐ฃโโโโถ. The name is explained in the module `s3_functions`.
 * `hat_tr_s3sdr3(v)`: Calculates the 6ร6-matrix ฬ๐ฃแต.
 * `hat_tr_mult_s3sdr3(v,w)`: Calculates the result of ฬ๐ฃแตโ๐ค for ๐ฃ,๐คโโโถ.
 * `expt_s3sdr3(v)`: Calculates the exponential of the tilde operator expt(๐ฃ)=exp(ฬ๐ฃ).
 * `tan_op_s3sdr3(v)`: Calculates the tangent operator ๐(๐ฃ).
 * `tan_op_inv_s3sdr3(v)`: Calculates the inverse matrix ๐โปยน(๐ฃ)=(๐(๐ฃ))โปยน of the tangent operator.
 * `tan_tr_inv_s3sdr3(v)`: Calculates ๐โปแต(๐ฃ)=((๐(๐ฃ))แต)โปยน, the inverse matrix of the transposed tangent operator.
 * `tan_tr_mult_s3sdr3(v,w)`: Calculates the result of ๐แต(๐ฃ)โ๐ค.
 * `tan_op_inv_mult_s3sdr3(v,w)`: Calculates ๐โปยน(๐ฃ)โ๐ค.
 * `tan_tr_inv_mult_s3sdr3(v,w)`: Calculates ๐โปแต(๐ฃ)โ๐ค.
 * `logt_s3sdr3(q)`: Calculates the logarithm logt(๐) with results in โโถ of a ๐โ๐ยณโโยณ sufficiently close to the identity. Note that it holds tilde(logt(๐))=log(๐).
 * `d_tan_tr_inv_s3sdr3(v, w)`: Calculates the Jacobi matrix of ๐โปแต(ฮฉ)โ๐ด with respect to the vector ฮฉโโยณ, where ๐ดโโยณ.

### Module `singular_functions`
This is a module that implements functions that are used in the other modules of this project, which have a singularity. In all of these functions, we assume that the input variable ๐ฅ is nonnegative, although there will be no error if ๐ฅ is actually negative. Usually, ๐ฅ will be the norm of some vector. All function names try to describe the function it evaluates by stringing together all letters and numbers by underscores.
If ๐ฅ is near a singular value, the function will be evaluated by a Taylor approximation. All thresholds and Taylor polynomial degrees are determined to minimize the errors as well as the loss-of-significance errors, which could occasionally become very large for ๐ฅ very close to a singular value and no Taylor approximation would be used.
Note that there is a file `singular_functions_1e-6.F90`, which implements the same module but with higher precision. This will, however, result in more computing time.

### Testing
The file `test.F90` contains a bunch of tests. At this point in time, all tests should pass, except for the test of `d_tan_tr_inv_s3sdr3` which is slightly above the error threshold.

## Usage
Call `make` in order to build the project. Object and module files will be created in the directory `obj/`. Call `make test` in order to run the test.

## Validation and finding error bounds
In order to validate the implementations there is the Mathematica file `misc/Lie_group_functions.nb`. It implements most of the Lie group functions.
Finding the error bounds that are used in the module `singular_functions` is not easy. Here, the Mathematica file `misc/singular_coefficients.nb` was used. It defines a function `singAnal[f,OptionsPattern[]]`, where a function `f` can be passed. In the options pattern, the critical point, the threshold, the order to the Taylor polynomial (among other things) can be specified. Calling `singAnal` on a function will give an overview of the function near the critical value. For this, the function is evaluated as given with machine precision and with extremely high precision. Futhermore, the Taylor approximation given and is evaluated and the three are compared. The goal is to adjust the threshold and the order of the Taylor polynomial in such a way that neither the Taylor approximation is too imprecise inside the threshold nor the as-is-evaluation outside the threshold has loss of significance that is too big.

## Related projects
Integrators:

 * [The Lie group generalized-ฮฑ method `gena`](https://github.com/StHante/gena)
 * [The Lie group BDF method `BLieDF`](https://github.com/StHante/BLieDF)
 * [The Lie group RATTLE method `RATTLie`](https://github.com/StHante/RATTLie)
 * [The Lie group SHAKE method `SHAKELie`](https://github.com/StHante/SHAKELie)
 * [The nonholonomic RATTLie method `RATTLie_nonhol`](https://github.com/StHante/RATTLie_nonhol)

Test problems:

 * [The heavy top example `heavy_top`](https://github.com/StHante/heavy_top)
 * [The constrained Cosserat beam model `crmS3R3`](https://github.com/StHante/crmS3R3)
 * [The rolling disk example `rolling_disk`](https://github.com/StHante/rolling_disk)

Miscellaneous:

 * [Implementation of Lie group functions `liegroup`](https://github.com/StHante/liegroup)
 * [Expand a config file with different configurations to several files `expandconfig`](https://github.com/StHante/expandconfig)
 * [Read lua files in Matlab and Octave `readLua`](https://github.com/StHante/readLua-for-Matlab-and-Octave)

Third party projects:

 * [Reading lua files in Fortran `aotus`](https://geb.sts.nt.uni-siegen.de/doxy/aotus/)
 * [GFortran](https://gcc.gnu.org/fortran/)
 * [GNU Parallel](https://www.gnu.org/software/parallel/)
