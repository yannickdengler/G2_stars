#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsubsection:: *)
(*Polytropic interpolation for G2 EoS*)


(* ::Input:: *)
(*dheavy={{0,0,0},{0.18,0.00110637,0.000123486},{0.2,0.00151811,0.000111776},{0.22,0.00208304,0.000079017},{0.24,0.00309918,0.000144038},{0.26,0.00429562,0.000260806},{0.28,0.00387751,0.000195464},{0.3,0.00786396,0.000404354},{0.32,0.0137287,0.000650502},{0.34,0.0197046,0.000598608},{0.36,0.0243828,0.000372009},{0.4,0.0271982,0.000182171},{0.42,0.0264943,0.000181594},{0.44,0.0274103,0.000173365},{0.48,0.0273374,0.000204433},{0.5,0.0295112,0.000382481},{0.5,0.0308065,0.000260144},{0.52,0.0352139,0.000335109},{0.54,0.0352143,0.000894907},{0.56,0.0374233,0.000881497},{0.6,0.0537042,0.000975687},{0.7,0.106987,0.00046004},{0.8,0.202169,0.000862469},{0.9,0.359465,0.00151906},{1.,0.693259,0.0012576},{1.1,1.5767,0.00192981},{1.2,3.22497,0.00198903},{1.3,5.51716,0.00143838},{1.4,8.44217,0.00112608},{1.5,11.5558,0.00174252},{1.6,13.6381,0.00275352}};*)
(*(**,{1.7,13.9944,0.000736476},{1.8,14.0011,0.000458305},{1.9,14.0019,0.000562317},{2.,14.0011,0.000318698})) cs>1**)*)


(* ::Input:: *)
(*hbarc = 197.32698;(**MeVfm**)*)
(*alat  = 0.357 ;(**fm**)*)
(*amn = 1.7;(**a**)*)
(*amn/alat*hbarc*)


(* ::Input:: *)
(*\[Mu]= dheavy[[All,1]]*3/amn*)
(*n=dheavy[[All,2]]/(3*amn^3) *)
(*\[Mu][[20]]*)
(*\[Mu][[21]]*)


(* ::Input:: *)
(*(**Length[latdat]*)
(*nsat = latdat[[38,2]];*)
(*\[Mu]sat = latdat[[38,1]];**)*)
(*\[Mu]arr=\[Mu][[21;;]]*)
(*narr = n[[21;;]]*)


(* ::Input:: *)
(*\[Mu]PP[n_, c_, \[CapitalGamma]_, K_]:=c+\[CapitalGamma]*K*n^(\[CapitalGamma]-1)/(\[CapitalGamma]-1)*)
(*PPP[n_, \[CapitalGamma]_, K_]:=K*n^\[CapitalGamma]*)
(*epsPP[n_,c_,\[CapitalGamma]_,K_]:=c*n+K*n^\[CapitalGamma]/(\[CapitalGamma]-1)*)
(*csPP[n_,c_,\[CapitalGamma]_,K_]:=(\[CapitalGamma]*PPP[n, \[CapitalGamma], K])/(PPP[n, \[CapitalGamma], K]+epsPP[n,c, \[CapitalGamma], K])*)


(* ::Input:: *)
(*(**next\[CapitalGamma][\[Mu]0_,\[Mu]1_,K0_,n0_,n1_,\[CapitalGamma]0_]:=NSolve[0==(\[CapitalGamma] K0 n0^(-\[CapitalGamma]+\[CapitalGamma]0) (-n0^(-1+\[CapitalGamma])+n1^(-1+\[CapitalGamma])))/(1-\[CapitalGamma])-\[Mu]0+\[Mu]1,\[CapitalGamma]]**)*)
(*next\[CapitalGamma][\[Mu]0_,\[Mu]1_,K0_,n0_,n1_,\[CapitalGamma]0_]:=\[CapitalGamma]/.FindRoot[(\[CapitalGamma] K0 n0^(-\[CapitalGamma]+\[CapitalGamma]0) (-n0^(-1+\[CapitalGamma])+n1^(-1+\[CapitalGamma])))/(1-\[CapitalGamma])-\[Mu]0+\[Mu]1,{\[CapitalGamma],1.5},PrecisionGoal->6]*)
(*nextK[n0_,\[CapitalGamma]1_,P0_]:=n0^-\[CapitalGamma]1 P0*)
(*nextc[\[Mu]1_,n1_,K1_,\[CapitalGamma]1_]:=\[Mu]1-(K1 n1^(-1+\[CapitalGamma]1) \[CapitalGamma]1)/(\[CapitalGamma]1-1)*)
(*nextP[n1_,\[CapitalGamma]1_,K1_]:=K1 n1^\[CapitalGamma]1*)


(* ::Input:: *)
(*\[CapitalGamma]00=5/3//N;*)
(*c00=1;*)
(*\[CapitalGamma]arr = Table[0,{i,Length[narr]}];*)
(*carr = Table[0,{i,Length[narr]}];*)
(*Karr = Table[0,{i,Length[narr]}];*)
(*Parr = Table[0,{i,Length[narr]}];*)
(*epsarr = Table[0,{i,Length[narr]}];*)
(*\[CapitalGamma]arr[[1]]=\[CapitalGamma]00;*)
(*carr[[1]]=c00;*)
(*Karr[[1]]=(\[Mu]arr[[1]]-carr[[1]])*(\[CapitalGamma]arr[[1]]-1)*narr[[1]]^(1-\[CapitalGamma]arr[[1]])/\[CapitalGamma]arr[[1]];*)
(*Parr[[1]]=Karr[[1]]*narr[[1]]^\[CapitalGamma]arr[[1]];*)
(*epsarr[[1]]=epsPP[narr[[1]],carr[[1]],\[CapitalGamma]arr[[1]],Karr[[1]]];*)


(* ::Input:: *)
(*For[i=2,i<Length[narr]+1,i++,*)
(*\[CapitalGamma]arr[[i]]=next\[CapitalGamma][\[Mu]arr[[i-1]],\[Mu]arr[[i]],Karr[[i-1]],narr[[i-1]],narr[[i]],\[CapitalGamma]arr[[i-1]]];*)
(*Karr[[i]]=nextK[narr[[i-1]],\[CapitalGamma]arr[[i]],Parr[[i-1]]];*)
(*carr[[i]]=nextc[\[Mu]arr[[i]],narr[[i]],Karr[[i]],\[CapitalGamma]arr[[i]]];*)
(*Parr[[i]]=nextP[narr[[i]],\[CapitalGamma]arr[[i]],Karr[[i]]];*)
(*epsarr[[i]]=epsPP[narr[[i]],carr[[i]],\[CapitalGamma]arr[[i]],Karr[[i]]];*)
(*]*)


(* ::Input:: *)
(*\[Mu]arr*)
(*narr;*)
(*Parr;*)
(*carr*)
(*Karr*)
(*\[CapitalGamma]arr*)
(*Parr*)
(*Export[NotebookDirectory[]<>"output/cKgP_DM_heavy.csv",Table[{carr[[i]],Karr[[i]],\[CapitalGamma]arr[[i]],Parr[[i]]},{i,1,Length[Parr]}], "CSV"]*)


(* ::Input:: *)
(*Plotlist={};*)
(*For[i=2,i<Length[narr]+1,i++,*)
(*(**AppendTo[Plotlist,Plot[\[Mu]PP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Blue]];*)
(*AppendTo[Plotlist,ListPlot[Transpose@{narr,\[Mu]arr}]];**)*)
(*AppendTo[Plotlist,Plot[csPP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Purple]];*)
(*AppendTo[Plotlist,Plot[1,{n,narr[[i-1]],narr[[i]]},PlotStyle->Black]];*)
(*AppendTo[Plotlist,Plot[PPP[n, \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Red]];*)
(*AppendTo[Plotlist,Plot[epsPP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Green]];*)
(*AppendTo[Plotlist,ListPlot[Transpose@{narr,Parr},PlotStyle->Pink]];*)
(*]*)
(*Show[Plotlist,PlotRange->Automatic,PlotLabel->Automatic]*)


(* ::Input:: *)
(*convmN4=101425.86;*)


(* ::Input:: *)
(*ListLogLogPlot[Table[{epsarr[[i]]*convmN4,Parr[[i]]*convmN4},{i,1,Length[Parr]}]]*)


(* ::Input:: *)
(*epsarrout=epsarr[[1;;9]]*convmN4 (** 1GeV^4 in MeV/fm^3**)*)
(*Parrout=Parr[[1;;9]]*convmN4 (** 1GeV^4 in MeV/fm^3**)*)
(*Export[NotebookDirectory[]<>"output/heavy_points.dat",{epsarrout,Parrout}]*)



