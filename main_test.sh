#!/usr/bin/bash

cp src/EoS_light_I.py src/EoS.py
python3 src/run.py 500 0.01 light_EoS_I
python3 src/run.py 1000 0.01 light_EoS_I
python3 src/run.py 2000 0.03 light_EoS_I
python3 src/run.py 4000 0.03 light_EoS_I
cp src/EoS_light_II.py src/EoS.py
python3 src/run.py 500 0.01 light_EoS_II
python3 src/run.py 1000 0.01 light_EoS_II
python3 src/run.py 2000 0.03 light_EoS_II
python3 src/run.py 4000 0.03 light_EoS_II
cp src/EoS_light_III.py src/EoS.py
python3 src/run.py 500 0.01 light_EoS_III
python3 src/run.py 1000 0.01 light_EoS_III
python3 src/run.py 2000 0.03 light_EoS_III
python3 src/run.py 4000 0.03 light_EoS_III
cp src/EoS_heavy_I.py src/EoS.py
python3 src/run.py 500 0.01 heavy_EoS_I
python3 src/run.py 1000 0.01 heavy_EoS_I
python3 src/run.py 2000 0.03 heavy_EoS_I
python3 src/run.py 4000 0.03 heavy_EoS_I
cp src/EoS_heavy_II.py src/EoS.py
python3 src/run.py 500 0.01 heavy_EoS_II
python3 src/run.py 1000 0.01 heavy_EoS_II
python3 src/run.py 2000 0.03 heavy_EoS_II
python3 src/run.py 4000 0.03 heavy_EoS_II
cp src/EoS_heavy_III.py src/EoS.py
python3 src/run.py 500 0.01 heavy_EoS_III
python3 src/run.py 1000 0.01 heavy_EoS_III
python3 src/run.py 2000 0.03 heavy_EoS_III
python3 src/run.py 4000 0.03 heavy_EoS_III