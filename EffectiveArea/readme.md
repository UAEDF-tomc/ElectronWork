# Notes on computing the effective areas for cut-based electron ID tuning.

=== March 2018 ===

   Simply run:

```
./compileAndRun.sh computeEffectiveAreaWithIsoCutoffs
```

The effective areas are based on 90% contours, like for 
the 2016 electron ID. One can run this script in "verbose"
mode and on a smaller event sample by changing the appropriate
boolean flag in the script.

To validate the effective area constants, use the script
validateCorrectedIsolationV2.C. In the beginning, edit the
line with the effective area constants, cut and paste there
the array of the slopes printed out by the previous script,
and run this validate... script twice, with the arguments true (for barrel)
and false (for the endcap). The efficiency for the isolation with
the effective area correction on the resulting plot should be nearly flat.

One can also compute effective areas that are mean-based with the
script computeEffectiveAreaV2.C, and validate them with the above validate...
script.


### Input tuples
The ntuples that are fed into this script (defined in the beginning
of the script as const TString) come from running the SimpleElectronNtupler,
available from the code below:

```
cmsrel CMSSW_10_6_2
cd CMSSW_10_6_2/src/
cmsenv
git clone https://github.com/UAEDF-tomc/EgammaWork.git
cd EgammaWork/
git checkout ntupler_106X
scram build -j 10
```

The cmsRun executables are in EgammaWork/ElectronNtupler/test:
runElectrons_AOD.py
runElectrons.py
(they only differ by a boolean)

Submit them using crabSubmit.py



# Quick notes on the scripts in this directory:

- computeEffectiveArea.C
- validateCorrectedIsolation.C
     These scripts calibrate the effective area based on the mean of the
     isolation and then draw the efficiency as a function of Nvtx.
     No generator level weights are used.

- computeEffectiveAreaV2.C
- validateCorrectedIsolationV2.C
     These scripts calibrate the effective area based on the mean of the
     isolation and then draw the efficiency as a function of Nvtx.
     Generator level weights are used in these scripts.

- validateCorrectedIsolationFlat.C
     Here we look at the central barrel region only and use different effective
     area constants.
     No generator weights.

- computeIsolationCutoffs.C
     Computes the cutoff line in the iso(rho) space for a given
     XX% efficiency.

- drawIsoWithCutoffs.C
     The script plots 2D uncorrected (pho + neu.had.) isolation vs rho
     and draws the XX% contours and the mean. The constants for contours 
     and the mean come from another script.

- computeEffectiveAreaWithIsoCutoffs.C
     The script computes the effective areas based on the XX-% efficiency contours.

- validateCorrectedIsolationAbs.C
     A variation of the validate script, draws the eff(Nvtx)
     not for relative but for the absolute total neutral isolation.

- drawChIsoConvProgression.C
- drawChIsoConvProgressionTest.C
     These scripts take charged isolation and convolute it with 
     a model pile-up isolation. The first works with the Poisson model
     and the second works with the Gaussian model.

- drawChIsoConvolution.C
     This script draws plots demonstrating the convolution of the
     charged isolation with two models for a particular rho slice:
     a Poisson and a Gaussian.

- drawIsolationROC.C
    Draws a single variable ROC: one for the default effective area
    correction and another for the alternative case (such as XX% cutoff-based).

- drawSigBgIso.C
    Draw (neu.had. + pho. isolation) vs (rho) with and without the effective
    area correction, with and without a cut line superimposed.

- computeEffectiveArea_v2.C
    A test of computing effective areas with via the mean, but
    after applying some upper cut-off to the values that go into the mean.

- testCutoffCalculation.C
    A script that computes error on the cutoff using the efficiency curve method
    and the toy MC method, for comparisons.

- hoeEnergyPrepare2D.C
- hoeEnergyFindDependence.C
- hoeRhoPrepare2D.C
- hoeRhoFindDependence.C
- hoeROC.C
- drawROCs.C
- drawFlatPtEff.C
    This set of scripts is for studying H/E cut dependence on energy
    and pileup.


