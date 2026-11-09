# Structure for Unifying Multiple Modeling Alternatives: SUMMA

[![Build Status](https://travis-ci.org/NCAR/summa.svg?branch=develop)](https://travis-ci.org/NCAR/summa)
[![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/NCAR/SUMMA/master/COPYING)
[![Documentation Status](https://readthedocs.org/projects/summa/badge/?version=latest)](http://summa.readthedocs.org/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.800772.svg)](https://doi.org/10.5281/zenodo.800772)

SUMMA (Clark et al., [2015a](#clark_2015a);[b](#clark_2015b);[c](#clark_2015c); [2021](#clark_2021)) is a hydrologic modeling framework that can be used for the systematic analysis of alternative model conceptualizations with respect to flux parameterizations, spatial configurations, and numerical solution techniques. It can be used to configure a wide range of hydrological model alternatives and we anticipate that systematic model analysis will help researchers and practitioners understand reasons for inter-model differences in model behavior. When applied across a large sample of catchments, SUMMA may provide insights in the dominance of different physical processes and regional variability in the suitability of different modeling approaches. An important application of SUMMA is selecting specific physics options to reproduce the behavior of existing models – these applications of "**model mimicry**" can be used to define reference (benchmark) cases in structured model comparison experiments, and can help diagnose weaknesses of individual models in different hydroclimatic regimes.

SUMMA is built on a common set of conservation equations and a common numerical solver, which together constitute the  “**structural core**” of the model. Different modeling approaches can then be implemented within the structural core, enabling a controlled and systematic analysis of alternative modeling options, and providing insight for future model development.

The important modeling features are:

 1. The formulation of the conservation model equations is cleanly separated from their numerical solution;

 2. Different model representations of physical processes (in particular, different flux parameterizations) can be used within a common set of conservation equations; and

 3. The physical processes can be organized in different spatial configurations, including model elements of different shape and connectivity (e.g., nested multi-scale grids and HRUs).


## Documentation
SUMMA documentation is available [online](http://summa.readthedocs.io/) and remains a work in progress. Additional SUMMA information including publications, test data sets, and sample applications can be found on the [SUMMA web site](http://www.ral.ucar.edu/projects/summa) at NCAR.

## Nested Newton Development Features
In this development branch, Newton iteration features are controlled using hard-coded parameters. In future work, these parameters will be migrated to SUMMA's model decision and parameter input files. For consistency with other SUMMA branches, the source code is set up to use the default standard Newton solver. However, nested Newton (and hybrid Newton) features may be controlled by setting the following parameters in `build/source/engine/varSubstep.f90`.

`nested_Newton_flag` (Line 219): `.false.` for SUMMA's default standard solver (all parameters below are ignored), `.true.` for nested Newton or hybrid Newton iterations

`nested_Newton % dynamic` (Line 739): `.false.` for nested Newton, `.true.` for hybrid Newton

`nested_Newton % kmax_classical` (Line 717, $\mathcal{K}_{\text{SN}}-1$): max number of standard Newton iterations minus one for hybrid Newton

`nested_Newton % kmax` (Line 749, $\mathcal{K}-1$): max number of outer iterations minus one for nested or hybrid Newton

`nested_Newton % lmax` (Line 749, $\mathcal{L}-1$): max number of inner iterations minus one for nested or hybrid Newton

`nested_Newton % order_min` (Line 740, $o_{\text{min}}$): computational convergence order threshold for hybrid Newton

Additionally, within `build/source/nested_Newton_library/Newton_solvers.f90` we have the following parameter.

`k_check` (Line 66, $k^{*}$): number of standard Newton iterations minus one prior to computational convergence order check for hybrid Newton

## Credits
SUMMA's initial implementation is described in two papers published in [Water Resources Research](http://onlinelibrary.wiley.com/journal/10.1002/(ISSN)1944-7973). If you use SUMMA, please credit these two publications.

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, 2015a: A unified approach for process-based hydrologic modeling: Part 1. Modeling concept. _Water Resources Research_, [doi:10.1002/2015WR017198](http://dx.doi.org/10.1002/2015WR017198).<a id="clark_2015a"></a>

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015b: A unified approach for process-based hydrologic modeling: Part 2. Model implementation and case studies. _Water Resources Research_, [doi:10.1002/2015WR017200](http://dx.doi.org/10.1002/2015WR017200).<a id="clark_2015b"></a>
 
 * Clark, M. P., Zolfaghari, R., Green, K. R., Trim, S., Knoben, W. J. M., Bennett, A., Nijssen, B., Ireson, A., Spiteri, R. J., 2021: The Numerical Implementation of Land Models: Problem Formulation and Laugh Tests. _Journal of Hydrometeorology_, [doi:10.1175/JHM-D-20-0175.1](http://dx.doi.org/10.1175/JHM-D-20-0175.1).<a id="clark_2021"></a>

In addition, an NCAR technical note describes the SUMMA implementation in detail:

 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, L. D. Brekke, J. R. Arnold, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015c: The structure for unifying multiple modeling alternatives (SUMMA), Version 1.0: Technical Description. _NCAR Technical Note NCAR/TN-514+STR_, 50 pp., [doi:10.5065/D6WQ01TD](http://dx.doi.org/10.5065/D6WQ01TD).<a id="clark_2015c"></a>


## License
SUMMA is distributed under the GNU Public License Version 3. For details see the file `COPYING` in the SUMMA root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
