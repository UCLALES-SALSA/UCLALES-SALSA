# UCLALES-SALSA

This description covers UCLALES-SALSA branch **IceDevelOrg** (https://github.com/UCLALES-SALSA/UCLALES-SALSA/tree/IceDevelOrg).

## Description
UCLALES-SALSA (Tonttila et al., 2017; Ahola et al., 2020) is a Large-Eddy Simulator (LES) coupled with a detailed sectional aerosol microphysics module SALSA (Kokkola et al., 2008) extended for cloud, rain, and ice microphysics. This version of SALSA includes semi-volatile organics (Prank et al., 2022) described with the Volatility Basis Set (VBS) approach. SALSA and all modifications and components related to the coupled UCLALES-SALSA have been designed and developed by the Finnish Meteorological Institute (see the MIT license in this directory).

The core model UCLALES is the work of Stevens et al. (Stevens et al., 1999, 2005; Stevens and Seifert, 2008). The original two-moment warm cloud microphysics by Seifert and Beheng (2001) and the mixed-phase and ice clouds by Seifert et al. (2006, 2008, 2012, 2014) are optionally included. UCLALES code is originally from https://github.com/uclales/uclales (last access: Mar 12, 2014, for the core model and May 24, 2018, for the ice microphysics). UCLALES code is copyrighted by Bjorn Stevens and protected by the General Public License, version 3, see gpl-3.0.txt in this directory

#### References
- Ahola, J., Korhonen, H., Tonttila, J., Romakkaniemi, S., Kokkola, H., and Raatikainen, T.: Modelling mixed-phase clouds with the large-eddy model UCLALES-SALSA, Atmos. Chem. Phys., 20, 11639-11654, https://doi.org/10.5194/acp-20-11639-2020, 2020.
- Kokkola, H., Korhonen, H., Lehtinen, K. E. J., Makkonen, R., Asmi, A., Jarvenoja, S., Anttila, T., Partanen, A.-I., Kulmala, M., Jarvinen, H., Laaksonen, A., and Kerminen, V.-M.: SALSA - a Sectional Aerosol module for Large Scale Applications, Atmos. Chem. Phys., 8, 2469-2483, https://doi.org/10.5194/acp-8-2469-2008, 2008.
- Prank, M., Tonttila, J., Ahola, J., Kokkola, H., Kuhn, T., Romakkaniemi, S., and Raatikainen, T.: Impacts of marine organic emissions on low-level stratiform clouds - a large eddy simulator study, Atmos. Chem. Phys., 22, 10971-10992, https://doi.org/10.5194/acp-22-10971-2022, 2022.
- Seifert, A., and Beheng, K. D.: A double-moment parameterization for simulating autoconversion, accretion and selfcollection, Atmos. Res., 59-60, 265-281, https://doi.org/10.1016/S0169-8095(01)00126-0, 2001.
- Seifert, A., and Beheng, K. D.: A two-moment cloud microphysics parameterization for mixed-phase clouds. Part 1: Model description, Meteorol. Atmos. Phys., 92, 45-66, https://doi.org/10.1007/s00703-005-0112-4, 2006.
- Seifert, A.: On the parameterization of evaporation of raindrops as simulated by a one-dimensional rainshaft model, J. Atmos. Sci., 65, 3608-3619, https://doi.org/10.1175/2008JAS2586.1, 2008.
- Seifert, A., Kohler, C., and Beheng, K. D.: Aerosol-cloud-precipitation effects over Germany as simulated by a convective-scale numerical weather prediction model, Atmos. Chem. Phys., 12, 709-725, https://doi.org/10.5194/acp-12-709-2012, 2012.
- Seifert, A., Blahak, U., and Buhr, R.: On the analytic approximation of bulk collision rates of non-spherical hydrometeors, Geosci. Model Dev., 7, 463-478, https://doi.org/10.5194/gmd-7-463-2014, 2014.
- Stevens, B., Moeng, C.-H., and Sullivan, P. P.: Large-Eddy Simulations of Radiatively Driven Convection: Sensitivities to the Representation of Small Scales, J. Atmos. Sci., 56, 3963-3984, https://doi.org/10.1175/1520-0469(1999)056<3963:LESORD>2.0.CO;2, 1999.
- Stevens, B., Moeng, C.-H., Ackerman, A. S., Bretherton, C. S., Chlond, A., de Roode, S., Edwards, J., Golaz, J.-C., Jiang, H., Khairoutdinov, M., Kirkpatrick, M. P., Lewellen, D. C., Lock, A., Muller, F., Stevens, D. E., Whelan, E., and Zhu, P.: Evaluation of Large-Eddy Simulations via Observations of Nocturnal Marine Stratocumulus, Mon.Weather Rev., 133, 1443-1462, https://doi.org/10.1175/MWR2930.1, 2005.
- Stevens, B. and Seifert, A.: Understanding macrophysical outcomes of microphysical choices in simulations of shallow cumulus convection, J. Meteorol. Soc. Jpn. Ser. II, 86A, 143-162, https://doi.org/10.2151/jmsj.86A.143, 2008.
- Tonttila, J., Maalick, Z., Raatikainen, T., Kokkola, H., Kuhn, T., and Romakkaniemi, S.: UCLALES-SALSA v1.0: a large-eddy model with interactive sectional microphysics for aerosol, clouds and precipitation, Geosci. Model Dev., 10, 169-188, https://doi.org/10.5194/gmd-10-169-2017, 2017.

### Compilation
A sample makefile *Makefile_Puhti_Intel* is valid for CSC's Puhti supercomputer when using Intel Fortran compiler environment. Dependencies (netCDF, HDF5, mpi) needs to be modified according to each computing environment. The compilation with flag *seq* or *mpi* produces either single processor sequential (./bin/les.seq) or parallel (./bin/les.mpi) executable.

### Run the model
The bin directory contains a sample NAMELIST (*NAMELIST_sample*) and sounding (*sound_in_sample*) for a test run. Just remove "_sample" from the file names. For a simple sequential run on a command prompt just execute **./les.seq**, but otherwise running the model depends on the computing environment.

## UCLALES-SALSA publications
- Tonttila, J., Maalick, Z., Raatikainen, T., Kokkola, H., Kuhn, T., and Romakkaniemi, S.: UCLALES-SALSA v1.0: a large-eddy model with interactive sectional microphysics for aerosol, clouds and precipitation, Geosci. Model Dev., 10, 169-188, https://doi.org/10.5194/gmd-10-169-2017, 2017.
- Stevens, R. G., Loewe, K., Dearden, C., Dimitrelos, A., Possner, A., Eirund, G. K., Raatikainen, T., Hill, A. A., Shipway, B. J., Wilkinson, J., Romakkaniemi, S., Tonttila, J., Laaksonen, A., Korhonen, H., Connolly, P., Lohmann, U., Hoose, C., Ekman, A. M. L., Carslaw, K. S., and Field, P. R.: A model intercomparison of CCN-limited tenuous clouds in the high Arctic, Atmos. Chem. Phys., 18, 11041-11071, https://doi.org/10.5194/acp-18-11041-2018, 2018.
- Ahola, J., Korhonen, H., Tonttila, J., Romakkaniemi, S., Kokkola, H., and Raatikainen, T.: Modelling mixed-phase clouds with the large-eddy model UCLALES-SALSA, Atmos. Chem. Phys., 20, 11639-11654, https://doi.org/10.5194/acp-20-11639-2020, 2020.
- Prank, M., Tonttila, J., Ahola, J., Kokkola, H., Kuhn, T., Romakkaniemi, S., and Raatikainen, T.: Impacts of marine organic emissions on low-level stratiform clouds - a large eddy simulator study, Atmos. Chem. Phys., 22, 10971-10992, https://doi.org/10.5194/acp-22-10971-2022, 2022.
- Ahola, J., Raatikainen, T., Alper, M. E., Keskinen, J.-P., Kokkola, H., Kukkurainen, A., Lipponen, A., Liu, J., Nordling, K., Partanen, A.-I., Romakkaniemi, S., Raisanen, P., Tonttila, J., and Korhonen, H.: Technical note: Parameterising cloud base updraft velocity of marine stratocumuli, Atmos. Chem. Phys., 22, 4523-4537, https://doi.org/10.5194/acp-22-4523-2022, 2022.
- Raatikainen, T., Prank, M., Ahola, J., Kokkola, H., Tonttila, J., and Romakkaniemi, S.: The effect of marine ice-nucleating particles on mixed-phase clouds, Atmos. Chem. Phys., 22, 3763-3778, https://doi.org/10.5194/acp-22-3763-2022, 2022.
