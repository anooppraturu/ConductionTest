(*Import initial radial temperature profiles and radial temperature \
profiles after the central temperature has increased by ~17% *)

iT64 = Drop[Import["iT64.dat", "Table"], 2][[All, 3]];
iT128 = Drop[Import["iT128.dat", "Table"], 2][[All, 3]];
iT256 = Drop[Import["iT256.dat", "Table"], 2][[All, 3]];
iT512 = Drop[Import["iT512.dat", "Table"], 2][[All, 3]];

T64 = Drop[Import["T64.dat", "Table"], 2][[All, 3]];
T128 = Drop[Import["T128.dat", "Table"], 2][[All, 3]];
T256 = Drop[Import["T256.dat", "Table"], 2][[All, 3]];
T512 = Drop[Import["T512.dat", "Table"], 2][[All, 3]];


(*Plot Showing the initial condition compared to the profile after it \
has been affcted by effective isotropic conduction. Notice that for \
all of our resolutions the shape of the profile is the same once the \
central temperature has changed by a certain amount. The only thing \
which changes is how long it took to get there!*)
ListPlot[{iT512, 
  T512, T128, T256, T64}, DataRange -> {0, 0.5}, ImageSize -> 512, 
 AxesLabel -> {"r", "T"}, 
 LabelStyle -> Directive[Black, Bold, Medium], Joined -> True]

(*Initialize the analytic solution to isotropic thermal conduction in \
2-D polar using Bessel's Functions*)
\[Sigma] = 0.15;
f0[x_] := Exp[-((x - 0.25)^2/(2 \[Sigma]^2))];
a[n_] := With[{\[Alpha] = BesselJZero[0, n] // N},
   2/BesselJ[1, \[Alpha]]^2 NIntegrate[
     x f0[x] BesselJ[0, \[Alpha] x], {x, 0, 1}]];
f[\[Kappa]_, t_] := Table[
    With[{\[Alpha] = N[BesselJZero[0, n]]},
     a[n] BesselJ[0, x \[Alpha]] Exp[-\[Kappa] \[Alpha]^2 t]], {n, 1, 
     20}] // Total;

(*Plot of the actual profile compared with analytic solutions within \
a  range of well fitting conductivities*)
ListPlot[
 Prepend[Table[
   Table[Evaluate[f[n, 0.02*5] + 1], {x, 0, 0.5, 0.01}], {n, 0.01, 
    0.03, 0.01}], T64], DataRange -> {0, 0.5}, ImageSize -> 512, 
 AxesLabel -> {"r", "T"}, 
 LabelStyle -> Directive[Black, Bold, Medium], Joined -> True]

(*Clearly, the analytic solution is not a solution to the equations \
Athena solves which gives rise to the deviations from perfectly \
anisotropic conduction. However The eigenmodes of isotropic \
conduction still give a fairly good fit, and we can exploit this to \
define a effective isotropic conductivity derived from the conduction \
timescale it takes for the central temperature to rise by ~ 17%*)

(*Kappas calculated based on conduction times*)

Kappa = {0.04/2, 0.04/7, 0.04/29, 0.04/100};
(*The initial perturbation was a gaussian with \[Sigma]=0.15, so here \
we give the number of cells per perturbation size (2\[Sigma]). \
(Domain size in physical units is 1)*)
res = 0.3*{64, 128, 256, 512};
(*fitting parameters for conductivities on log-log plot to obtain \
power law*)

bfit = Covariance[res // Log, Kappa // Log]/Variance[res // Log];
afit = Mean[Kappa // Log] - bfit*Mean[res // Log];
kapfit = Table[
   afit + bfit*i, {i, res[[1]] // Log, res[[4]] // Log, 0.05}];
fitdom = Table[i, {i, res[[1]] // Log, res[[4]] // Log, 0.05}];

(*Plots on log-log and normal scale of conductivity versus number of \
cells per perturbation.*)
(*Call cells per perturbation C*)
\
ListPlot[{Transpose[{res // Log, Kappa // Log}], 
  Transpose[{fitdom, kapfit}]}, PlotRange -> All, ImageSize -> 512, 
 AxesLabel -> {"Log[C]", "Log[\[Kappa]]"}, 
 LabelStyle -> Directive[Black, Bold, Medium]]
ListPlot[{Transpose[{res, Kappa}], 
  Transpose[{Exp[fitdom], Exp[kapfit]}]}, PlotRange -> All, 
 ImageSize -> 512, AxesLabel -> {"C", "\[Kappa]"}, 
 LabelStyle -> Directive[Black, Bold, Medium]]
