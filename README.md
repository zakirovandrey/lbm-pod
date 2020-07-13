
To run.

Prepare: 

```
> export LD_RUN_PATH=$LD_RUN_PATH:./                                            
> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./
```

Run:

```
> make clean; ARCH=zu vglrun ./shear_wave_test.py --test > log.log
```

Shear wave diagnostics:

```
> for n in log*; do echo "#".$n; grep Max $n | python3 shear_wave_diag.py ; done
```

Output Max Velocity graph:

```
> grep Max log.log|awk '{print $10}'>u.dat     
```

