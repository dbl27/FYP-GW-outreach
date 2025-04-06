FYP - Dancing With the Neutron Stars

Diego Baird-Ludlow & Alec Shackleton

Gravitational wave outreach project meant to help illustrate the creation process of gravitaional waves and GW strain amplitude.

-> lib:
    Contains fromACC2.py file, which details the creation of the MassPoint data class and the subsequent calculation of the time-varying mass quadrupole moment, from predefined input trajectories.

-> results:
    data_animated.ipynb creates animated plots of trajectories and GW strain data for a given .mp4 recording. The csv output files in /blender/data are created externally with blender, with help of CSVscript.txt in /blender.

-> validating_results:
    The accuracy of the time-varying QM calculations in fromACC2.py are verified with GWOSC data and trjectories from SXS simulation database. compare_GWOSC.ipynb plots trajectories from SXS data, calculates and plots strain with functions from fromACC2.py, and plots observed strain data from the GWOSC database for comparison. /GWOSC_data contains H1 and L1 strain data from merger event GW150914, with H1 and L1 representing data taken from the Hanford and Livingstone observatories respectively (https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/). /SXS_data contains simulated trajectories for black hole merger BBH:0305 (https://data.black-holes.org/waveforms/catalog.html).

-> CMU:
    Unfinished- contains motion tracking data used for testing purposes from Carnegie Mellon University Graphics Lab Motion Capture Database. Includes cartwheel and Lindy Hop trajectories.