#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*\[Mu]PP[n_, c_, \[CapitalGamma]_, K_]:=c+\[CapitalGamma]*K*n^(\[CapitalGamma]-1)/(\[CapitalGamma]-1)*)
(*PPP[n_, \[CapitalGamma]_, K_]:=K*n^\[CapitalGamma]*)
(*epsPP[n_,c_,\[CapitalGamma]_,K_]:=c*n+K*n^\[CapitalGamma]/(\[CapitalGamma]-1)*)
(*csPP[n_,c_,\[CapitalGamma]_,K_]:=(\[CapitalGamma]*PPP[n, \[CapitalGamma], K])/(PPP[n, \[CapitalGamma], K]+epsPP[n,c, \[CapitalGamma], K])*)


(* ::Input:: *)
(*next\[CapitalGamma][\[Mu]0_,\[Mu]1_,K0_,n0_,n1_,\[CapitalGamma]0_]:=\[CapitalGamma]/.FindRoot[(\[CapitalGamma] K0 n0^(-\[CapitalGamma]+\[CapitalGamma]0) (-n0^(-1+\[CapitalGamma])+n1^(-1+\[CapitalGamma])))/(1-\[CapitalGamma])-\[Mu]0+\[Mu]1,{\[CapitalGamma],1.5},PrecisionGoal->6]*)
(*nextK[n0_,\[CapitalGamma]1_,P0_]:=n0^-\[CapitalGamma]1 P0*)
(*nextc[\[Mu]1_,n1_,K1_,\[CapitalGamma]1_]:=\[Mu]1-(K1 n1^(-1+\[CapitalGamma]1) \[CapitalGamma]1)/(\[CapitalGamma]1-1)*)
(*nextP[n1_,\[CapitalGamma]1_,K1_]:=K1 n1^\[CapitalGamma]1*)


(* ::Subsubsection:: *)
(*Polytropic interpolation for Hebeler soft plus 1*)


(* ::Input:: *)
(*n0=0.16;(**fm^-3**)*)
(*neutronmass=939.5654205;(**MeV**)*)
(*hbarc=197.32698;(**MeV fm**)*)


(* ::Input:: *)
(*Nucldat=Transpose[Import[NotebookDirectory[]<>"input/EoS_Hebeler.csv","CSV"]];*)
(*Length[Nucldat];*)
(*EoSI=Transpose[Import[NotebookDirectory[]<>"input/EoS_I.csv","CSV"]];*)
(*nEoSI=EoSI[[1]]*n0/(neutronmass/hbarc)^3;*)
(*PEoSI=EoSI[[2]]*hbarc^3/neutronmass^4;*)
(*\[Epsilon]EoSI=EoSI[[3]]*hbarc^3/neutronmass^4;*)
(*\[Mu]EoSI=EoSI[[4]]/(neutronmass/1000);*)


(* ::Input:: *)
(*novern0=Nucldat[[1]][[;;6]];*)
(*narr=novern0*n0/(neutronmass/hbarc)^3;*)
(*PNucl=Nucldat[[3]][[;;6]]*hbarc^3/neutronmass^4;*)
(*\[Epsilon]Nucl=Nucldat[[4]][[;;6]]*hbarc^3/neutronmass^4;*)
(*\[Mu]arr=Table[(\[Epsilon]Nucl[[i]]+PNucl[[i]])/narr[[i]],{i,Length[narr]}];*)
(*(** EoS I **)*)
(*AppendTo[narr,1.1*n0/(neutronmass/hbarc)^3];*)
(*AppendTo[narr,5.9*n0/(neutronmass/hbarc)^3];*)
(*AppendTo[narr,7.5*n0/(neutronmass/hbarc)^3];*)
(*AppendTo[\[Mu]arr,0.9657/(neutronmass/1000)];*)
(*AppendTo[\[Mu]arr,1.653/(neutronmass/1000)];*)
(*AppendTo[\[Mu]arr,1.772/(neutronmass/1000)];*)


(* ::Input:: *)
(*\[CapitalGamma]00=5/3//N;*)
(*c00=1;*)
(*Clear[\[CapitalGamma]arr,carr,Karr,Parr]*)
(*\[CapitalGamma]arr = Table[0,{i,Length[narr]+1}];*)
(*carr = Table[0,{i,Length[narr]+1}];*)
(*Karr = Table[0,{i,Length[narr]+1}];*)
(*Parr = Table[0,{i,Length[narr]+1}];*)
(*epsarr = Table[0,{i,Length[narr]+1}];*)
(*\[CapitalGamma]arr[[1]]=\[CapitalGamma]00;*)
(*carr[[1]]=c00;*)
(*Karr[[1]]=(\[Mu]arr[[1]]-carr[[1]])*(\[CapitalGamma]arr[[1]]-1)*narr[[1]]^(1-\[CapitalGamma]arr[[1]])/\[CapitalGamma]arr[[1]];*)
(*Parr[[1]]=Karr[[1]]*narr[[1]]^\[CapitalGamma]arr[[1]];*)
(*epsarr[[1]]=epsPP[narr[[1]],carr[[1]],\[CapitalGamma]arr[[1]],Karr[[1]]];*)
(*For[i=2,i<=Length[narr],i++,*)
(*\[CapitalGamma]arr[[i]]=next\[CapitalGamma][\[Mu]arr[[i-1]],\[Mu]arr[[i]],Karr[[i-1]],narr[[i-1]],narr[[i]],\[CapitalGamma]arr[[i-1]]];*)
(*Karr[[i]]=nextK[narr[[i-1]],\[CapitalGamma]arr[[i]],Parr[[i-1]]];*)
(*carr[[i]]=nextc[\[Mu]arr[[i]],narr[[i]],Karr[[i]],\[CapitalGamma]arr[[i]]];*)
(*Parr[[i]]=nextP[narr[[i]],\[CapitalGamma]arr[[i]],Karr[[i]]];*)
(*epsarr[[i]]=epsPP[narr[[i]],carr[[i]],\[CapitalGamma]arr[[i]],Karr[[i]]];*)
(*Print[\[CapitalGamma]arr[[i]], " ",Karr[[i]], " ",carr[[i]], " ",Parr[[i]]]*)
(*]*)
(*nadd =0.035;*)
(*\[CapitalGamma]arr[[Length[\[CapitalGamma]arr]]]=\[CapitalGamma]arr[[Length[\[CapitalGamma]arr]-1]];*)
(*Karr[[Length[Karr]]]=Karr[[Length[Karr]-1]];*)
(*carr[[Length[carr]]]=carr[[Length[carr]-1]];*)
(*Parr[[Length[Parr]]]=PPP[nadd, \[CapitalGamma]arr[[Length[\[CapitalGamma]arr]]], Karr[[Length[Karr]]]];*)
(*AppendTo[narr,nadd];*)
(*AppendTo[\[Mu]arr,\[Mu]PP[nadd, carr[[Length[carr]]], \[CapitalGamma]arr[[Length[carr]]], Karr[[Length[carr]]]]];*)
(*Print[\[CapitalGamma]arr[[Length[\[CapitalGamma]arr]]], " ",Karr[[Length[Karr]]], " ",carr[[Length[carr]]], " ",Parr[[Length[Parr]]]];*)


(* ::Input:: *)
(*\[Mu]arr*neutronmass*)
(*narr*(neutronmass/hbarc)^3*)
(*Parr;*)
(*carr;*)
(*Karr;*)
(*\[CapitalGamma]arr;*)
(*Parr;*)
(*Export[NotebookDirectory[]<>"output/cKgP_EoSI.csv",Table[{carr[[i]],Karr[[i]],\[CapitalGamma]arr[[i]],Parr[[i]]},{i,1,Length[Parr]}], "CSV"]*)


(* ::Input:: *)
(*Plotlist={};*)
(*For[i=2,i<Length[narr]+1,i++,*)
(*(**AppendTo[Plotlist,Plot[\[Mu]PP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Blue]];*)
(*(**AppendTo[Plotlist,ListPlot[Transpose@{narr,\[Mu]arr}]];**)*)
(*AppendTo[Plotlist,ListPlot[Transpose@{nEoSI,\[Mu]EoSI}]];**)*)
(*AppendTo[Plotlist,Plot[PPP[n, \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Red]];*)
(*AppendTo[Plotlist,Plot[epsPP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Green]];*)
(*(**AppendTo[Plotlist,Plot[1,{n,narr[[i-1]],narr[[i]]},PlotStyle->Black]];*)
(*AppendTo[Plotlist,Plot[csPP[n, carr[[i]], \[CapitalGamma]arr[[i]], Karr[[i]]],{n,narr[[i-1]],narr[[i]]},PlotStyle->Purple]];**)]*)
(*AppendTo[Plotlist,ListPlot[Transpose@{nEoSI,PEoSI}]];*)
(*AppendTo[Plotlist,ListPlot[Transpose@{nEoSI,\[Epsilon]EoSI}]];*)
(*(**Show[Plotlist,PlotRange->{{0,0.004},{0,0.0001}},PlotLabel->Automatic]**)*)
(*Show[Plotlist,PlotRange->All,PlotLabel->Automatic]*)


(* ::Input:: *)
(*convmN4=101425.86*)


(* ::Input:: *)
(*ListLogLogPlot[Table[{epsarr[[i]]*convmN4,Parr[[i]]*convmN4},{i,1,Length[Parr]}]]*)


(* ::Input:: *)
(*epsarrout=epsarr[[1;;9]]*convmN4 (** 1GeV^4 in MeV/fm^3**)*)
(*Parrout=Parr[[1;;9]]*convmN4 (** 1GeV^4 in MeV/fm^3**)*)
(*Export[NotebookDirectory[]<>"output/EoS_I_points.dat",{epsarrout,Parrout}]*)



