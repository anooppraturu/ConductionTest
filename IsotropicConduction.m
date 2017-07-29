(* IsotropicConduction.m: makes plots for the convergence analysis in
   the appendix of Anoop's paper. *)

(* This script both documents our work and enables you to reproduce
   it; run it using Mathematica. *)



(* ========================================================================== *)
(* preliminaries *)

(* start by making a directory to store our plots *)
CreateDirectory["plots"];

(* use a uniform padding so plots line up *)
pad = {{60,10}, {40,10}};



(* ========================================================================== *)
(* analytic solution to the heat equation... *)
(* ... this can be compared to the numerical results to see if there's
   an effective perpendicular conductivity *)

(* shape of the perturbation at time t=0 *)
sigma = 0.15;

f0[x_] := Exp[-(x-0.25)^2 / (2*sigma^2)];

(* nth Fourier-Bessel coefficient of the initial condition *)
a[n_] := With[{alpha = BesselJZero[0, n] // N},
              (2/BesselJ[1, alpha]^2) * NIntegrate[
                  x f0[x] BesselJ[0, alpha x], {x, 0, 1}]];

(* Fourier-Bessel expansion of the solution at time t *)
f[kappa_, t_] := 
    Table[
        With[{alpha = N[BesselJZero[0, n]]},
             a[n] BesselJ[0, x alpha] Exp[-kappa alpha^2 t]], 
        {n, 1, 20}] // Total;



(* ========================================================================== *)
(* temperature profiles from numerical simulations. *)
(* profiles taken at the time when the central temperature has
   increased by ~17% *)

(* import a file and drop the header *)
myImport[file_] := 
    Block[{tmp = Import[file, "Table"]},
          (* drop header *)
          tmp = Select[tmp, VectorQ[#, NumberQ] &];

          (* return {r,T} pairs *)
          tmp[[All, {2,3}]] ];

ics = Map[myImport, 
          {"iT64.dat",
           "iT128.dat",
           "iT256.dat",
           "iT512.dat"}];

vals = Map[myImport, 
           {"T64.dat",
            "T128.dat",
            "T256.dat",
            "T512.dat"}];



(* ========================================================================== *)
(* plot comparing the initial temperature perturbation to the profile
   after it has diffused out due to numerical errors.  notice that for
   all of our resolutions the shape of the profile is the same once
   the central temperature has changed by a certain amount.  the only
   thing which changes is how long it took to get there! *)

p0 = ListPlot[Last[ics],
              PlotStyle    -> {{Thick, Dashed, Black}},
              PlotRange    -> {{0.0, 0.5}, {1.0, 2.1}},
              ImageSize    -> 400,
              Frame        -> True,
              FrameLabel   -> {"r", "T"},
              ImagePadding -> pad,
              BaseStyle    -> {"FontSize" -> 14, "FontFamily" -> "Helvetica"},
              Joined       -> True];

p = ListPlot[vals, 
             Joined    -> True,
             PlotStyle -> Thick];

Export["plots/temperatures.pdf", Show[{p0,p}], "PDF"];



(* ========================================================================== *)
(* Plot of the actual profile compared with analytic solutions within
   a range of well fitting conductivities *)

kappas = {0.01, 0.02, 0.03};

sols = Map[f[#, 0.02*5] + 1.0 &, 
           kappas];

p0 = ListPlot[Last[vals],
              PlotStyle    -> {{Thick, Dashed, Black}},
              PlotRange    -> {{0.0, 0.5}, {1.0, 2.1}},
              ImageSize    -> 400,
              Frame        -> True,
              FrameLabel   -> {"r", "T"},
              ImagePadding -> pad,
              BaseStyle    -> {"FontSize" -> 14, "FontFamily" -> "Helvetica"},
              Joined       -> True];

p = Plot[sols, {x, 0, 0.5},
         PlotStyle -> Thick];

Export["plots/compare-analytic-conduction.pdf", 
       Show[{p0, p}], 
       "PDF"];


(* the fit isn't exact, which makes sense since Athena is not solving
   a heat equation in the radial direction -- in fact, the temperature
   profile should not evolve at all.  we're seeing the effect of
   numeric errors, which are under no obligation to look like the heat
   equation.  the fit looks reasonable, however, so we can interpret
   numeric errors as an effective perpendicular conductivity. *)



(* ========================================================================== *)
(* try defining kappa from the conduction timescale *)

kappas = {0.04/2, 0.04/7, 0.04/29, 0.04/100};

res = (2*sigma/0.5) * {64, 128, 256, 512};


(* fit data to a line in log-space *)
data = Transpose[{Log[res], Log[kappas]}];

nlm = NonlinearModelFit[data,
                        a x + b,
                        {a, b}, x];

p = LogLogPlot[Exp[nlm[Log[x]]], {x, 8, 512},
               Epilog       -> {PointSize[Large], Point[data]},
               FrameLabel   -> {"\[CapitalDelta]x \[Del]T", 
                                "\!\(\*SubscriptBox[\"\[Kappa]\", \"perp\"]\)/\!\(\*SubscriptBox[\" \[Kappa]\", \"para\"]\)"},
               ImageSize    -> 400,
               Frame        -> True,
               FrameLabel   -> {"r", "T"},
               ImagePadding -> pad,
               BaseStyle    -> {"FontSize" -> 14, "FontFamily" -> "Helvetica"}];

Export["plots/convergence.pdf", p, "PDF"];
