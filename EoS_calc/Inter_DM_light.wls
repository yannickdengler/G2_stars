#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsubsection:: *)
(*Polytropic interpolation for G2 EoS*)


(* ::Input:: *)
(*(**latdat=Import["/Users/yannickdengler/Documents/Uni/PhD_Graz/Code/G2_stars/EoS/data/EoS_G2_prime"]**)*)
(*latdat=Import[NotebookDirectory[]<>"input/EoS_G2_prime"]*)


(* ::Input:: *)
(*Length[latdat]*)
(*nsat = latdat[[38,2]]*)
(*\[Mu]sat = latdat[[38,1]]*)
(*\[Mu]gt1=latdat[[30;;37]];*)
(*Length[%]*)


(* ::Input:: *)
(*\[Mu]PP[n_, c_, \[CapitalGamma]_, K_]:=c+\[CapitalGamma]*K*n^(\[CapitalGamma]-1)/(\[CapitalGamma]-1)*)
(*PPP[n_, \[CapitalGamma]_, K_]:=K*n^\[CapitalGamma]*)
(*epsPP[n_,c_,\[CapitalGamma]_,K_]:=c*n+K*n^\[CapitalGamma]/(\[CapitalGamma]-1)*)
(*csPP[n_,c_,\[CapitalGamma]_,K_]:=(\[CapitalGamma]*PPP[n, \[CapitalGamma], K])/(PPP[n, \[CapitalGamma], K]+epsPP[n,c, \[CapitalGamma], K])*)


(* ::Input:: *)
(*nFD[\[Mu]_, a_,b_]:=nsat/(Exp[a-b*\[Mu]]+1)*)
(*\[Mu]FD[n_,a_,b_]:=(a-Log[nsat/n-1])/b*)


(* ::Input:: *)
(*n1=latdat[[35,2]];*)
(*n2=latdat[[37,2]];*)
(*\[Mu]1=latdat[[35,1]];*)
(*\[Mu]2=latdat[[37,1]];*)
(*a=(\[Mu]1 Log[(n1 (n2-nsat))/(n2 (n1-nsat))])/(\[Mu]1-\[Mu]2)+Log[-1+nsat/n1];*)
(*b=Log[(n1 (n2-nsat))/(n2 (n1-nsat))]/(\[Mu]1-\[Mu]2);*)
(*\[Mu]arrFD = Table[x,{x,1.5,3.5,0.5}];*)
(*narrFD = Table[nFD[\[Mu]arrFD[[i]],a,b],{i,Length[\[Mu]arrFD]}];*)


(* ::Input:: *)
(*Plotlist={};*)
(*AppendTo[Plotlist,Plot[nFD[\[Mu],a,b],{\[Mu],\[Mu]1,\[Mu]sat}]];*)
(*AppendTo[Plotlist,ListPlot[Transpose@{\[Mu]gt1[[All,1]],\[Mu]gt1[[All,2]]}]];*)
(*AppendTo[Plotlist,ListPlot[{{\[Mu]sat,nsat}},PlotStyle->Orange]];*)
(*AppendTo[Plotlist,ListPlot[Transpose@{\[Mu]arrFD,narrFD},PlotStyle->Red]];*)
(*Show[Plotlist, PlotRange->All]*)


(* ::Input:: *)
(*(**next\[CapitalGamma][\[Mu]0_,\[Mu]1_,K0_,n0_,n1_,\[CapitalGamma]0_]:=NSolve[0==(\[CapitalGamma] K0 n0^(-\[CapitalGamma]+\[CapitalGamma]0) (-n0^(-1+\[CapitalGamma])+n1^(-1+\[CapitalGamma])))/(1-\[CapitalGamma])-\[Mu]0+\[Mu]1,\[CapitalGamma]]**)*)
(*next\[CapitalGamma][\[Mu]0_,\[Mu]1_,K0_,n0_,n1_,\[CapitalGamma]0_]:=\[CapitalGamma]/.FindRoot[(\[CapitalGamma] K0 n0^(-\[CapitalGamma]+\[CapitalGamma]0) (-n0^(-1+\[CapitalGamma])+n1^(-1+\[CapitalGamma])))/(1-\[CapitalGamma])-\[Mu]0+\[Mu]1,{\[CapitalGamma],1.5},PrecisionGoal->6]*)
(*nextK[n0_,\[CapitalGamma]1_,P0_]:=n0^-\[CapitalGamma]1 P0*)
(*nextc[\[Mu]1_,n1_,K1_,\[CapitalGamma]1_]:=\[Mu]1-(K1 n1^(-1+\[CapitalGamma]1) \[CapitalGamma]1)/(\[CapitalGamma]1-1)*)
(*nextP[n1_,\[CapitalGamma]1_,K1_]:=K1 n1^\[CapitalGamma]1*)


(* ::Input:: *)
(*narr = Join[\[Mu]gt1[[All,2]],narrFD];*)
(*\[Mu]arr = Join[\[Mu]gt1[[All,1]],\[Mu]arrFD];*)
(*Length[\[Mu]arr]*)
(*Length[\[Mu]gt1[[All,1]]]*)


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
(*Export[NotebookDirectory[]<>"output/cKgP_DM_light.csv",Table[{carr[[i]],Karr[[i]],\[CapitalGamma]arr[[i]],Parr[[i]]},{i,1,Length[Parr]}], "CSV"]*)


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
(*Export[NotebookDirectory[]<>"output/light_points.dat",{epsarrout,Parrout}]*)



