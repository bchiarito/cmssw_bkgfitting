#! /bin/bash
if [[ $# -ne 1 ]]; then
    echo "Must supply input rootfile" >&2
    exit 2
fi
rm -rf loose_fit_hists/
rm -rf tight_templates/
rm plots.pdf
## <<< run scaler here
python fitter.py $1 --testBin "noniso_sym barrel 100" --createLooseFits
python fitter.py $1 --testBin "noniso_sym barrel 100" --saveTightTemplates --savePoly --ftest "3 2"
python plotter.py $1 --testBin "noniso_sym barrel 100"
xdg-open plots.pdf
rm -rf loose_fit_hists/
rm -rf tight_templates/
