#!/usr/bin/env python

import os,shutil

# Helpers
from SimEx.Utilities.Units import meter, electronvolt, joule, radian

# PMI
from SimEx.Calculators.XMDYNDemoPhotonMatterInteractor import XMDYNDemoPhotonMatterInteractor

# Simple Beam Parameters
from SimEx.Parameters.PhotonBeamParameters import PhotonBeamParameters

# Diffraction
from SimEx.Parameters.DetectorGeometry import DetectorGeometry, DetectorPanel
from SimEx.Parameters.SingFELPhotonDiffractorParameters import SingFELPhotonDiffractorParameters
from SimEx.Calculators.SingFELPhotonDiffractor import SingFELPhotonDiffractor

def cleanUp(out_path):
    dirs_to_remove=[out_path]
    files_to_remove=[out_path+'.h5']

    for d in dirs_to_remove:
        if os.path.isdir(d):
            shutil.rmtree(d)
    for f in files_to_remove:
        if os.path.isfile(f):
            os.remove(f)
# Detector setup
p0 = DetectorPanel(ranges={
    'fast_scan_min': 0,
    'fast_scan_max': 109,
    'slow_scan_min': 0,
    'slow_scan_max': 109
},
                   pixel_size=1000e-6 * meter,
                   photon_response=1.0,
                   distance_from_interaction_plane=0.25 * meter,
                   corners={
                       'x': 54.5,
                       'y': 54.5
                   },
                   fast_scan_xyz='1.0x',
                   slow_scan_xyz='1.0y')

detector_geometry = DetectorGeometry(panels=[p0])

beam = PhotonBeamParameters(
    photon_energy=6.0e3 * electronvolt,
    beam_diameter_fwhm=226.0e-9 * meter, #sqrt(250*160/pi)*2
    pulse_energy=1.52e-3 * joule, #4.0 * 0.38 mJ
)

diffraction_parameters = SingFELPhotonDiffractorParameters(
                                               uniform_rotation=False,
                                               calculate_Compton = False,
                                               number_of_diffraction_patterns=6240,
                                               detector_geometry=detector_geometry,
                                               beam_parameters = beam,
                                               forced_mpi_command='mpirun -np 40'
                                               )

out_path = "diffr"


cleanUp(out_path)

diffractor = SingFELPhotonDiffractor(parameters=diffraction_parameters,
                                     input_path='./4V7V.pdb',
                                     output_path=out_path)


from timeit import default_timer as timer
start = timer()

diffractor.backengine()

end = timer()
print(end - start,'s') # Time in seconds


# Link the .h5 output generated from mpirun together
diffractor.saveH5()



