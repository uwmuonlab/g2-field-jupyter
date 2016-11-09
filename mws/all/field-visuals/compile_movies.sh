#!/bin/bash

ffmpeg -f concat -i calibration_dipole.txt -c copy vid/calibration_dipole.mpg
ffmpeg -f concat -i calibration_multipoles.txt -c copy vid/calibration_multipoles.mpg
ffmpeg -f concat -i pole_moves_dipole.txt -c copy vid/pole_moves_dipole.mpg
ffmpeg -f concat -i pole_moves_multipoles.txt -c copy vid/pole_moves_multipoles.mpg
ffmpeg -f concat -i shim_moves_dipole.txt -c copy vid/shim_moves_dipole.mpg
ffmpeg -f concat -i shim_moves_multipoles.txt -c copy vid/shim_moves_multipoles.mpg
ffmpeg -f concat -i pole_laminations_dipole.txt -c copy vid/pole_laminations_dipole.mpg
ffmpeg -f concat -i pole_laminations_multipoles.txt -c copy vid/pole_laminations_multipoles.mpg
