#! /bin/bash
if [[ $# -ne 1 ]]; then
    echo "Must supply input rootfile" >&2
    exit 2
fi
rm -rf scaled_tight_hists/
rm -rf loose_fit_hists/
rm -rf tight_templates/
rm plots.pdf 2> /dev/null
python tightscaler.py $1 --testBin "noniso_sym barrel 100"
python fitter.py $1 --testBin "noniso_sym barrel 100" --createLooseFits
python fitter.py $1 --testBin "noniso_sym barrel 100" --saveTightTemplates --savePoly --ftest "3 2" --useScaledTight
python plotter.py $1 --testBin "noniso_sym barrel 100" --name "test"
xdg-open test.pdf
rm plots.pdf
rm -rf scaled_tight_hists/
rm -rf loose_fit_hists/
rm -rf tight_templates/
