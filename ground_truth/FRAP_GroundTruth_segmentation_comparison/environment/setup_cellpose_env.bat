@echo off
rem Creates the conda environment used for the paper:
rem   cellpose3  - Cellpose 3.1.1.3 (nuclei model), Python 3.11,
rem                PyTorch 2.14.1 built for CUDA 12.6
rem Without an NVIDIA GPU, replace the cu126 index URL with
rem https://download.pytorch.org/whl/cpu (Cellpose then runs on the CPU, more slowly).
rem Progress is written to setup_log.txt next to this file.

set LOG=%~dp0setup_log.txt
echo START %date% %time% > "%LOG%"
call "%USERPROFILE%\anaconda3\Scripts\activate.bat" >> "%LOG%" 2>&1

call conda create -y -n cellpose3 -c conda-forge --override-channels python=3.11 >> "%LOG%" 2>&1
call conda activate cellpose3 >> "%LOG%" 2>&1
python -m pip install torch==2.14.1 --index-url https://download.pytorch.org/whl/cu126 >> "%LOG%" 2>&1
python -m pip install cellpose==3.1.1.3 tifffile >> "%LOG%" 2>&1

python -c "import torch; from importlib.metadata import version; print('cellpose3: cellpose', version('cellpose'), '| torch', torch.__version__, '| GPU', torch.cuda.is_available())" >> "%LOG%" 2>&1
echo DONE %date% %time% >> "%LOG%"
