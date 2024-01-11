### runing ###
three steps:
```
python run_scaler.py <dir_name>
python run_fitter.py <dir_name>
python run_plotter.py <dir_name>
```
where `<dir_name>` is a directory with a file called `summed_egamma.root`

then for BAT:
```
python process_for_BAT.py <dir_name>
```

### to get root and python2 and scipy
```
$ conda install -n base conda-forge::mamba
$ conda create --name <env_name>
$ mamba install root=6.20.0=py27h97dbdcd_0 scipy
```
