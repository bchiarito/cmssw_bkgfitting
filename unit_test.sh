#! /bin/bash
if [[ $# -ne 1 ]]; then
    echo "Must supply input directory" >&2
    exit 2
fi
rm -rf $1/scaled_tight_hists/
rm -rf $1/scaled_phislice_tight_hists/
rm -rf $1/loose_fit_hists/
rm -rf $1/tight_templates/
rm $1/test.pdf
python run_scaler.py $1 --testBin "noniso_sym barrel 100"
python run_fitter.py $1 --testBin "noniso_sym barrel 100" --ftest "3 1" --useScaledTight
python run_fitter.py $1 --testBin "noniso_sym barrel 100" --ftest "3 1"
python run_plotter.py $1 --testBin "noniso_sym barrel 100" --name "test" --useScaledTight
rm $1/plots.pdf
xdg-open $1/test.pdf
