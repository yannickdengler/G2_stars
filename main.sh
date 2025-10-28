#!/usr/bin/bash

cp src/EoS_light_I.py src/EoS.py
python3 src/run.py 500 0.001 light_EoS_I
python3 src/run.py 1000 0.001 light_EoS_I
python3 src/run.py 2000 0.003 light_EoS_I
python3 src/run.py 4000 0.003 light_EoS_I
cp src/EoS_light_II.py src/EoS.py
python3 src/run.py 500 0.001 light_EoS_II
python3 src/run.py 1000 0.001 light_EoS_II
python3 src/run.py 2000 0.003 light_EoS_II
python3 src/run.py 4000 0.003 light_EoS_II
cp src/EoS_light_III.py src/EoS.py
python3 src/run.py 500 0.001 light_EoS_III
python3 src/run.py 1000 0.001 light_EoS_III
python3 src/run.py 2000 0.003 light_EoS_III
python3 src/run.py 4000 0.003 light_EoS_III
cp src/EoS_heavy_I.py src/EoS.py
python3 src/run.py 500 0.001 heavy_EoS_I
python3 src/run.py 1000 0.001 heavy_EoS_I
python3 src/run.py 2000 0.003 heavy_EoS_I
python3 src/run.py 4000 0.003 heavy_EoS_I
cp src/EoS_heavy_II.py src/EoS.py
python3 src/run.py 500 0.001 heavy_EoS_II
python3 src/run.py 1000 0.001 heavy_EoS_II
python3 src/run.py 2000 0.003 heavy_EoS_II
python3 src/run.py 4000 0.003 heavy_EoS_II
cp src/EoS_heavy_III.py src/EoS.py
python3 src/run.py 500 0.001 heavy_EoS_III
python3 src/run.py 1000 0.001 heavy_EoS_III
python3 src/run.py 2000 0.003 heavy_EoS_III
python3 src/run.py 4000 0.003 heavy_EoS_III