@echo off
rem Creates the conda environment "stardist" used for the paper:
rem   Python 3.11, TensorFlow 2.15.1 (CPU), StarDist 0.9.2, CSBDeep 0.8.2, NumPy 1.26.4
rem TensorFlow is kept at 2.15 because newer versions switch to Keras 3, which
rem StarDist does not support. Progress is written to setup_log.txt next to this file.

set LOG=%~dp0setup_log.txt
echo START %date% %time% > "%LOG%"
call "%USERPROFILE%\anaconda3\Scripts\activate.bat" >> "%LOG%" 2>&1
call conda create -y -n stardist -c conda-forge --override-channels python=3.11 >> "%LOG%" 2>&1
call conda activate stardist >> "%LOG%" 2>&1
python -m pip install --upgrade pip >> "%LOG%" 2>&1
python -m pip install tensorflow==2.15.1 numpy==1.26.4 stardist==0.9.2 csbdeep==0.8.2 tifffile >> "%LOG%" 2>&1
python -c "from importlib.metadata import version as v; print('stardist', v('stardist'), '| csbdeep', v('csbdeep'), '| tensorflow', v('tensorflow'), '| numpy', v('numpy'))" >> "%LOG%" 2>&1
echo DONE %date% %time% >> "%LOG%"
