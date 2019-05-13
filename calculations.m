ClearAll["Global`*"]
(* Solve With Trial Wave Functions *)

\[Lambda] = 1; (*set to 1 for incident electron -1 for incident hole*) 
rCoeffs = {roe, roh, \[Omega]oe, \[Omega]oh};
tCoeffs = {toe, toh, \[Omega]oe, \[Omega]oh};
consts = {c1, c1, c2, c2};
kvalues = {ke, kh, ke, kh};

If[\[Lambda] == 1,
 wfsD1 = {Exp[I*k*x] + roe*Exp[-I*k*x], toe*Exp[I*k*x]};
 wfsB1 = {roh*Exp[I*k*x], toh*Exp[-I*k*x]};
 wfsD2 = {\[Omega]oe*Exp[-I*k*x], \[Omega]oe*Exp[I*k*x]};
 wfsB2 = {\[Omega]oh*Exp[I*k*x], \[Omega]oh*Exp[-I*k*x]};,
 If[\[Lambda] == -1,
  wfsD1 = {roe*Exp[-I*k*x], toe*Exp[I*k*x]};
  wfsB1 = {Exp[-I*k*x] + roh*Exp[I*k*x], toh*Exp[-I*k*x]};
  wfsD2 = {\[Omega]oe*Exp[-I*k*x], \[Omega]oe*Exp[I*k*x]};
  wfsB2 = {\[Omega]oh*Exp[I*k*x], \[Omega]oh*Exp[-I*k*x]};,
  {Print["error: \[Lambda] must be \[PlusMinus]1"], Abort[]},
  {Print["error: \[Lambda] must be \[PlusMinus]1"], Abort[]}],
 {Print["error: \[Lambda] must be \[PlusMinus]1"], Abort[]}]

listOfWfs = {wfsD1, wfsB1, wfsD2, wfsB2};

diffEqFunction[\[Phi]_] := D[\[Phi], {x, 2}] + k^2*\[Phi] + c1*DiracDelta[x]

testFunction[\[Phi]Lis_] := Table[
  If[FullSimplify[
    diffEqFunction[\[Phi]Lis[[i]]] == 0,
    If[i == 1, x > 0, x < 0]],
   Nothing,
   Print["error: wfs do not solf diffeq"],
   Print["error: wfs do not solf diffeq"]],
  {i, 1, Length[\[Phi]Lis]}]
If[
 ParallelMap[testFunction, listOfWfs] == Table[{Nothing},
   {Length[listOfWfs]}],
 Null,
 {Print[ParallelMap[testFunction, listOfWfs]], 
  Abort[]}] (*check that wfs solve diffeq*)

zValuePaper = enM^2 - (en + I * (t1^2 / ve + t1^2/vh)) * (en + I * (t2^2 / ve + t2^2 / vh)); (*form of Z given in Paper*)
c1Paper = \[Lambda] * t1 * (en + I * (t2^2 / ve + t2^2 / vh)) / zValuePaper;   (*form of c1 given in Paper*)
c2Paper = \[Lambda]*(I * enM * t1) / zValuePaper;   (*form of c2 given in Paper*)

(* Match wfs for first BC *)
boundaryConds1 = Flatten[Table[
   Solve[FullSimplify[listOfWfs[[i, 1]] == listOfWfs[[i, 2]], x == 0], tCoeffs[[i]]],
   {i, Length[listOfWfs]}], 1]

(* Second BC derivatives have discontinuity due to bound state *)
(* \!\(TraditionalForm\` \*SubscriptBox[\(\[PartialD]\), \(x\)] \*SubscriptBox[\(\[Phi]\), \(+\)] \*SubscriptBox[\(|\), \(x = \[Epsilon]\)]\ -\) \!\(TraditionalForm\`\(\(\*SubscriptBox[\(\[PartialD]\), \(x\)]\*SubscriptBox[\(\[Phi]\), \(-\)]\) \*SubscriptBox[\(|\), \(x = \(-\[Epsilon]\)\)]\) + \*SubscriptBox[\(c\), \(1\)] \*SubscriptBox[\(t\), \(1\)]\) = 0 *)

probsInit = Flatten[Table[
   Solve[Limit[(D[listOfWfs[[i, 2]], x] /. x -> \[Epsilon]) - (D[listOfWfs[[i, 1]], x] /. x -> \[Epsilon]), \[Epsilon] -> 0] + consts[[i]] == 0 /. boundaryConds1[[i]] /. k -> kvalues[[i]], rCoeffs[[i]]],
   {i, 1, Length[listOfWfs]}]]

tunProbs = Table[probsInit[[i, 2]], {i, 1, 4}];

(* Solve for c1 and c2 by combining other two equations *)
paperConsts = FullSimplify[Solve[
      {t1*Limit[listOfWfs[[2, 2]] - listOfWfs[[1, 2]], x -> 0] - I*enM*const2 - en*const1 == 0,
       t2*Limit[listOfWfs[[4, 2]] - listOfWfs[[3, 2]], x -> 0] + I*enM*const1 - en*const2 == 0},
      {const1, const2}] /. boundaryConds1[[1]] /. boundaryConds1[[2]]];

constsPaperForm = FullSimplify[ Solve[{Evaluate[paperConsts[[1, 1, 2]] /. probsInit] == \[Beta]*c1/ t1, Evaluate[paperConsts[[1, 2, 2]] /. probsInit] == \[Beta]* c2/t2}, {c1, c2}] /. {ke -> ve/(2*\[Beta]), kh -> vh/(2*\[Beta])}];

If[
 And[Reduce[FullSimplify[\[Beta]/t1*constsPaperForm[[1, 1, 2]] == c1Paper]], FullSimplify[\[Beta]/t2*constsPaperForm[[1, 2, 2]] == c2Paper]],
 Print["Math is Correct"]]

{roe -> (I c1)/(2 ke), roh -> -((I c1)/(2 kh)), \[Omega]oe -> (I c2)/(2 ke), \[Omega]oh -> -((I c2)/(2 kh))};
{roe -> (I c1)/(2 ke), roh -> -((I c1)/(2 kh)), \[Omega]oe -> (I c2)/(2 ke), \[Omega]oh -> -((I c2)/(2 kh))};
