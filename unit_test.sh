#! /bin/bash
if [[ $# -ne 1 ]]; then
    echo "Must supply input directory" >&2
    exit 2
fi
rm -rf $1/scaled_tight_hists/
rm -rf $1/loose_fit_hists/
rm -rf $1/tight_templates/
python run_scaler.py $1 --testBin "noniso_sym barrel 100"
python run_fitter.py $1 --testBin "noniso_sym barrel 100" --createLooseFits
python run_fitter.py $1 --testBin "noniso_sym barrel 100" --saveTightTemplates --savePoly --ftest "3 2" --useScaledTight
python run_plotter.py $1 --testBin "noniso_sym barrel 100" --name "test"
rm $1/plots.pdf
rm -rf $1/scaled_tight_hists/
rm -rf $1/loose_fit_hists/
rm -rf $1/tight_templates/
xdg-open $1/test.pdf
