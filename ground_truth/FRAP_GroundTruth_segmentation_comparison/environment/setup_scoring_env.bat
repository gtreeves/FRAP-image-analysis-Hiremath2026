@echo off
rem Creates the conda environment "frapscore" used to score the segmentations
rem against the ground truth (code\scoring). Versions are the ones the paper's
rem numbers were produced with. Progress is written to setup_scoring_log.txt.

set LOG=%~dp0setup_scoring_log.txt
echo START %date% %time% > "%LOG%"
call "%USERPROFILE%\anaconda3\Scripts\activate.bat" >> "%LOG%" 2>&1
call conda create -y -n frapscore -c conda-forge --override-channels python=3.13 >> "%LOG%" 2>&1
call conda activate frapscore >> "%LOG%" 2>&1
python -m pip install --upgrade pip >> "%LOG%" 2>&1
python -m pip install numpy==2.5.3 scipy==1.18.1 pandas==3.0.5 matplotlib==3.11.2 tifffile==2026.9.15 openpyxl==3.1.5 pillow==12.3.0 h5py==3.16.0 >> "%LOG%" 2>&1
python -c "from importlib.metadata import version as v; print(*[p + ' ' + v(p) for p in ['numpy','scipy','pandas','matplotlib','tifffile','openpyxl','pillow','h5py']], sep=' | ')" >> "%LOG%" 2>&1
echo DONE %date% %time% >> "%LOG%"
