(* ========================================================================== *)
(* analytic solution to the heat equation... *)
(* ... this can be compared to the numerical results to see if there's
   an effective perpendicular conductivity *)

(* shape of the perturbation at time t=0 *)
sigma = 0.15;

f0[x_] := Exp[-(x-0.5)^2 / (2*sigma^2)];

(* nth Fourier-Bessel coefficient of the initial condition *)
a[n_] := With[{alpha = BesselJZero[0, n] // N},
              (2/BesselJ[1, alpha]^2) * NIntegrate[x f0[x] BesselJ[0, alpha x], {x, 0, 1}]];

(* Fourier-Bessel expansion of the solution at time t *)
f[kappa_, t_] :=
    Table[
        With[{alpha = N[BesselJZero[0, n]]},
             a[n] BesselJ[0, x alpha] Exp[-kappa alpha^2 t]],
        {n, 1, 20}] // Total;



(* ========================================================================== *)
(* print out analytic solutions to the heat equation at *)

tcond = sigma^2/1.0;

x  = Range[0, 100]/100.0;

y1 = Map[1.0 + 0.1 f0[#] &, x];

y2 = 1 + 0.1 f[1, 0.5 tcond];
y3 = 1 + 0.1 f[1, 1.0 tcond];

myfmt[num_] := ToString[PaddedForm[num, {10, 6}]];

data = Transpose[{x, y1, y2, y3}];
data = Map[myfmt, data, {2}];

header = {
"# analytic.dat",
"#",
"# analytic solutions to the heat equation for our initial condition",
"# assuming an effective perpendicular conductivity",
"#",
"# [1] = r, [2] = T(t=0), [3] = T(t=0.5 tcond), [4] = T(t=1.0 tcond)",
"#",
"# tcond is defined as (sigma^2)/kappa",
"#"
}

data = Join[header, data];

Export["analytic.dat",
       Append[data, {}]];



(* ========================================================================== *)
(* print out central temperatures at different times *)

ClearAll[x];

centraltemp[t_] := 1 + 0.1 f[1, t tcond] /. {x -> 0};

ts    = {0.0, 0.05, 0.1, 0.5, 1.0};
temps = Map[centraltemp, ts];

data = Transpose[{ts, temps}];
data = Map[myfmt, data, {2}];

header = {
"# central-temp.dat",
"#",
"# central temperature as a function of time from the analytic",
"# solution to the head equation for our initial condition",
"# assuming an effective perpendicular conductivity",
"#",
"# [1] = t/tcond, [2] = T(r=0)",
"#",
"# tcond is defined as (sigma^2)/kappa",
"#"
}

data = Join[header, data];

Export["central-temp.dat",
       Append[data, {}]];

Exit[0];
