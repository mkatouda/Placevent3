#!/bin/bash

#pwd
guv=complex_r3d-kh_guv.O.0.dx
thresh_p=5.0
watOs_pdb=watOs_placevent5.0.pdb

placevent ${guv} 55.5 ${thresh_p} 2.4 > ${watOs_pdb}
