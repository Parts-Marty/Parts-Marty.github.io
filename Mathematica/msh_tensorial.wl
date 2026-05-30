(* ::Package:: *)

(* :Context: mshtensorial` *)

(* :Title: Tensor formulation *)
(* :Sub Titile:Using Mathematica for tensor analysis*)
(* :Author:Lic. Renan Cabrera *)
(* :Improved by Shinabe Munehiro (Nickneme:Parts Marty)*)

(* :Summary: This is package for Tensor calculus *)
(* :Package Version: 2.6       2025/07/02*)

(*
This package was originally developed by author Lic. Renan Cabrera in collaboration with
 "Academia Nacional de Ciencias de Bolivia" and improved by author Parts Marty.
  It may be freely used and distributed for non-commercial purposes, but you must credit
  the author in each use.   
This software package and its accompanying documentation are provided without guarantee
of support or maintenace.Mathematica V.13.3 can be acceptable.
 If you want to help me to improve this package or you have any suggestion or you want
   the last version please contact me at "shinabe.munehiro@hotmail.co.jp"
*)

(*  Mathematica is a registered trademark of Wolfram Research, Inc. *)


BeginPackage["mshtensorial`"]

(*&&&&&&&&&&&&&&&&&&&& Matrix Function &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
MatTen::usage ="MatTen[A_,ind_] represents a Matrix with label A and indexes ind_"
ExpandMat::usage = "ExpandMat[ MatTen[A,indn], indt] ,  Expand a Matrix over with indt and indn" 

(*&&&&&&&&&&&&&&&&&&&& Tensor Function &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
Tensor::usage = "Tensor[A,sub,sup] represents a tensor with label A and indexes _sub and  ^sup"
x::usage ="Tensor[x,{Void},{p}]"
(*Overscript[x, _]::usage ="Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{p}]"*)
OrderTensor::usage ="OrderTensor[w,g], metric Tensor[g,{j,i},{t,s}] to Tensor[g,{i,j},{s,t}]"
SetTensorRule::usage = "SetTensorRule[Tensor,Values], Defines a set of rules for a expanded 
 Tensor with particular Values"
ExpandIndex::usage = "ExpandIndex[Tensor , ind ], Expand a tensor over free indexes ind" 
ExpandDoubleIndex::usgae = "ExpandDoubleIndex[Tensor,ind], Expand double indexes like a Sum 
over indexes ind and contract"
TensorSimplify::usage = "MetricSimplify [w,g], simplify Metric Tensor g "
DummyVariableSimplify::usage = 
             "DummyVariableSimplify[ T  ], Simplify dummy variables in Tensor expresions T"
             
(*&&&&&&&&&&&&&&&&&&&& Kronecker,SetSystem of Metric,LeviCivita &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
Kronecker::usage = "Kronecker[i,j] ,is the ordinary Kroneckers Delta,
                   Kronecker[ i,j... , k,l,m..... ] is the generalized Kroneckers Delta"
IndexDerivativeSimplify::usage = "IndexDerivativeSimplify[ Tensor ],Convert to Kronecker"
KroneckerSimplify::usage = "KroneckerSimplify[ Tensor ], simplify Kronecker"                   
MetricSimplify::usage = "MetricSimplify [w,g], simplify Metric Tensor g "
LeviCivitaRule::usage = " LeviCivitaRule "
SetSystem::usage = "SetSystem[x,g,GM] set up the geometrical space. x is the label for the
coordinates, g is the label for the Metric Tensor and GM is the matrix containing the 
particular values of the Metric Tensor"

(*&&&&&&&&&&&&&&&&&&&&  a lot of label &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
Void::usage =""
\[Epsilon] ::usage = " LeviCivita label"
\[Delta] ::usage = " Delta Kronecker label & Dummy Coeff & AbsoluteDerivative"
\[CapitalGamma] ::usage = " Christoffel label"
\[Rho] ::usage = " Dummy Coeff"
\[Sigma] ::usage = " Dummy Coeff"
\[Xi] ::usage = " Dummy Ricci Coeff"
NDim::usage = "ExpandIndex dimensin"

(*&&&&&&&&&&&&&&&&&&&&  Derivative &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
ToCovariantIndex::usage = " ToCovariantIndex[w,ind,newind] Brings Tensor w to 
   Covariant Index ind with new label newind"
ToContravariantIndex::usage = " ToCovariantIndex[w,ind,newind] Brings Tensor w to 
   Contravariantvariant Index ind with new label newind"
   
IndexDerivative::usage = " IndexDerivative[ T , Tensor[x,{Void},{i}] ] , 
Represents the partial derivative of T respect the coordinate i "

IndexDerivativeUp::usage = " IndexDerivative[ T , Tensor[x,{Void},{i}], m , g , x], 
Represents the partial derivative of T respect the coordinate i "

ExpandIndexDerivativeUp::usage ="ExpandIndexDerivativeUp[w],
Expand IndexDerivativeUp to Metric Tensor*IndexDerivative "

CovariantDerivative::usage = "CovariantDerivative[Tensor[T,{r},{Void}],x,s,g],
Represents the Covariant Derivative respect the index s"
 
AbsoluteDerivative::usage = "AbsoluteDerivative[Tensor[T,{Void},{i}],r,t,g,x], 
AbsoluteDerivative of Tensor[T,{Void},{i}] over with the index r,parameter t,coordinates x and MetricTensor g"

ExpandABSDer::usage ="Expand AbsoluteDerivative to CovariantDerivative*IndexDerivative "(*by msh *)

EvaluateIndexDerivative::usage = "EvaluateIndexDerivative[ T ], Evaluate IndexDerivatives in T"

OrderIndexDeri::usage="IndexDeri[w,x],IndexDerivative[ T , x[i],x[j]] to IndexDerivative[ T , x[j],x[i]]"


(*&&&&&&&&&&&&&&&&&&&&  Christoffel symbol &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)
Christoffel1::usage = " Cristoffel1[g,x][{i,k},{j}] symbols of class 1 with metric Tensor g 
and coordinates x"
Christoffel2::usage = " Christoffel2[g,x][{i,k},j}] symbols of class 2 with metric Tensor g 
and coordinates x"
ExpandChristoffel1::usage ="Expand Cristoffel1 to metric Tensor g and coordinates x"
ExpandChristoffel2::usage ="Expand Cristoffel2 to metric Tensor g and coordinates x"
OrderChristoffel::usage ="{b,a} of Cristoffel1,2 to {a,b}"
         
(*&&&&&&&&&&&&&&&&&&&&  Riemman-Chritoffel curvature Tensor &&&&&&&&&&&&&&&&&&&&&&&&&*)
RieChrisCur::usage = "RieChrisCur[{j},{k,l},{i},g,x] Riemman-Chritoffel curvature Tensor
					  RieChrisCur[{i,j},{k,l},{Void},g,x] Riemman-Chritoffel curvature Tensor
                      RieChrisCur[{i},{j},{Void},g,x] Ricci tensor of the first kind
                      RieChrisCur[{Void},{j},{i},g,x] Ricci tensor of the second kind
                      RieChrisCur[{Void},{Void},{i,j},g,x] Ricci tensor 
                      RieChrisCur[{Void},{Void},{Void},g,x] Ricci curvature"  (*by msh*)
ExpandRCC::usage = "Expand RieChrisCur Tensor"   (*by msh*)
ExpandCODer::usage = "Expand CovariantDerivative"   (*by msh*)
Expand2ndCODer::usage = "Expand 2nd CovariantDerivative"   (*by msh*)  
RiccFormula::usage ="RiccFormula[{k,l},T,g,x]]" (*by msh*) 
ExpandRiccFor::usage ="Expand Ricc Formula" (*by msh*) 
ExpandRicci::usage ="Expand to Ricci tensor" (*by msh*) 
EinsteinTensor::usage = "EinsteinTensor"   
(*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*)


(*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*)

Begin["`Private`"]

(*.............................................*)
(*...............Format Display................*)
(*.............................................*)
(*mshfont="Cambria Math";*)
mshfont= "Times New Roman";
(*mshfont="Yu Gothic UI";*)
mshsize=16;
Void=" "

Format[ Tensor[A_] ]=Style[A,Black,Italic,mshsize,FontFamily->mshfont];

Format[ Tensor[A_,sub_,sup_] ]:=   
Style[Subsuperscript[A, SequenceForm@@sub , SequenceForm@@sup ],Black,Italic,mshsize,FontFamily->mshfont];

Unprotect[Power];
Format[ Power[A_,sup_] ]:=Style[Superscript[ SequenceForm["(",A ,")"],sup],Black,Plain,14,FontFamily->mshfont];

(*Format[ Power[A_,sup_] ]:=If[sup==1/2,Style[SqrtBox["A"],Black,Plain,16,FontFamily->mshfont],
Style[Superscript[ SequenceForm["(",A ,")"],sup],Black,Plain,16,FontFamily->mshfont]];*)

(*Format[ Power[A_,sup_] ]:=Superscript[ SequenceForm["(",A ,")"],sup];*)
Protect[Power];

(*.............................................*)
(*.................SetSystem...................*)
(*.............................................*)
(*x[p_]:=Tensor[x,{Void},{p}];*)
x[p_]:=Tensor[x,{Void},{p}];
\!\(\*OverscriptBox[\(x\), \(_\)]\)[p_]:=Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{p}];
NDim=3; (*by msh*)
 (*SetSystem\:306b\:3088\:308a\:8a08\:91cf\:30c6\:30f3\:30bd\:30ebg\:3092\:8a2d\:5b9a\:3059\:308bSetSystem\:306b\:3088\:308a\:8a08\:91cf\:30c6\:30f3\:30bd\:30ebg\:3092\:8a2d\:5b9a\:3057\:3001TensorSimplify\:304c\:4f7f\:3048\:308b\:3002*)
SetSystem[var_,g_,mat_]:= Module[{i,j,w},
   NDim = Length[mat];
      
   gExpanded  = ExpandIndex[ Tensor[g,{i,j},{Void,Void}]  , {i,j} ]; 
   gInvExpanded = ExpandIndex[ Tensor[g,{Void,Void},{i,j}]  , {i,j} ]; 
    
   Evaluate[ gExpanded ]=    mat ;
   Evaluate[  gInvExpanded  ] = Inverse[mat] ;    



(*IDXX[ IndexDerivative[F_,Tensor[var,{Void},{i_}] ],
IndexDerivative[G_,Tensor[var,{Void},{j_}]]] :=  Kronecker[i,j];*)



TensorMetricRuleSimplifyg = { 
(*     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,i,j]/;MemberAllQ[w,i]&& FreeAllQ[w,IndexDerivative] ,
     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,j,i]/;MemberAllQ[w,j]&& FreeAllQ[w,IndexDerivative]  ,*)
(*\:964d\:683c by msh*)     
     Tensor[w_,sub_,sup_]Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[Tensor[w,sub,sup],i,j]/;MemberAllQ[sup,i],
     Tensor[w_,sub_,sup_]Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[Tensor[w,sub,sup],j,i]/;MemberAllQ[sup,j],
(*\:6607\:683c by msh*) 
	Tensor[w_,sub_,sup_] Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[Tensor[w,sub,sup],i,j]/;MemberAllQ[sub,i] ,
    Tensor[w_,sub_,sup_] Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[Tensor[w,sub,sup],j,i]/;MemberAllQ[sub,j],
      
	Tensor[g,{a_,i_},{Void,Void}]Tensor[g,{Void,Void},{a_,j_}] :> Kronecker[i,j],
	Tensor[g,{i_,a_},{Void,Void}]Tensor[g,{Void,Void},{j_,a_}] :> Kronecker[i,j],
	Tensor[g,{a_,i_},{Void,Void}]Tensor[g,{Void,Void},{j_,a_}] :> Kronecker[i,j],
	(*Kronecker[i_,i_]:>1,*)
	
    Tensor[g,{i_,Void},{Void,j_}] :> Kronecker[i,j],
    Tensor[g,{Void,i_},{j_,Void}] :> Kronecker[i,j],
    
   (IndexDerivative[Tensor[g,{i_,p_},{Void,Void}],x[s_]]-IndexDerivative[Tensor[g,{s_,i_},{Void,Void}],x[p_]])
   Tensor[g,{Void,Void},{s_,p_}]:>0
        
(*    IndexDerivative[Tensor[g,{i_,p_},{Void,Void}],x[s_]]Tensor[g,{Void,Void},{s_,p_}]
   -IndexDerivative[Tensor[g,{s_,i_},{Void,Void}],x[p_]]Tensor[g,{Void,Void},{s_,p_}]:>0      *)  
   
                       
        } ;

(*CovariantDerivative[ Tensor[g,{i_,j_},{Void,Void}]  ,  __ ] = 0;
CovariantDerivative[ Tensor[g,{Void,Void},{i_,j_}] ,  __ ] = 0; by msh *)

ToContravariantIndex[ w_  ,ind_ ,toind_] := 
   w Tensor[g,{Void,Void},{ind,toind}];

            ]	


		
(*.............................................*)
(*...............ExpandIndex...................*)
(*.............................................*)

ExpandIndex[ T_ , ind_] :=  
  Module[{ni,dum},
		
  dum = {#,Sequence@@{1,NDim}} &/@  ind ;  
  (*"&":(3+#)&[x] 3+x   "/@":Map "@@":Apply \:5f0f expr \:306e\:982d\:90e8\:3092 f \:3067\:7f6e\:63db\:3059\:308b *)
  (*dum={{n,1,3},{m,1,3},{k,1,3}}*)
   (*Evaluate[Sequence@@dum]={n,1,3},{m,1,3},{k,1,3}*)
  s1=Table[ T, Evaluate[Sequence@@dum] ]
  (*MatrixForm[s1]*)
]/;VectorQ[ind];

(*............................................*)	  
(*............ExpandDoubleIndex...............*)	
(*............................................*)

ExpandDoubleIndex[w_,ind_]:=ExpandDoubleIndexINT[ Expand[w] , ind ];
(*ExpandDoubleIndexINT[Power[w_,n_],ind_List]:=Power[Fold[ExpandDoubleIndexINT[#1,#2]&,w,ind],n];*)
ExpandDoubleIndexINT[Power[w_,n_],ind_List]:=Power[ExpandDoubleIndexINT[w,ind],n];
ExpandDoubleIndexINT[w_,ind_List]:= Fold[ExpandDoubleIndexINT[#1,#2]& , w , ind];
ExpandDoubleIndexINT[w_Plus,ind_]:=  Map[ExpandDoubleIndexINT[#,ind]& , w ];
ExpandDoubleIndexINT[w_,ind_] := w/;FreeAllQ[w,ind];

ExpandDoubleIndexINT[w_,ind_]:=Module[{q},
       Plus@@Table[ w /.ind->q, {q,1,NDim} ]  
				]; 
(*...........................................*)
(*...............TensorSimplify..............*)
(*....KroneckerRuleSimplify..................*)
(*....TensorMetricRuleSimplifyg.............*)

TensorSimplify[w_]:=Module[{Regla},
	Regla = {Sequence@@KroneckerRuleSimplify,
                 Sequence@@TensorMetricRuleSimplifyg
                 };
        w //.Regla
		]

(*...........................................*)
(*....IndexDerivativeSimplify..? by msh........*)
(*...........................................*)
IndexDerivativeSimplify[w_]:=Module[{Regla},
	Regla = {IndexDerivative[Tensor[x_,{Void},{j_}],Tensor[x_,{Void},{i_}]]->Kronecker[i,j],
	IndexDerivative[Tensor[x_,{j_},{Void}],Tensor[x_,{i_},{Void}]]->Kronecker[i,j],
	IndexDerivative[ F_ ,v_] IndexDerivative[ v_ , u_ ]-> IndexDerivative[F, u],
	IndexDerivative[ F_ , n_, m_ ] IndexDerivative[ n_ , u_ ]-> IndexDerivative[F,u,m],
	IndexDerivative[ F_ , n_, m_ ] IndexDerivative[ m_ , u_ ]-> IndexDerivative[F,n,u],
	IndexDerivative[Kronecker[j_,i_],u_]->0,	
	IndexDerivative[ IndexDerivative[F_,Tensor[x_,{Void},{i_}] ],IndexDerivative[F_,Tensor[x_,{Void},{j_}]]]->Kronecker[j,i]
	(*IndexDerivative[F_,Tensor[x_,{Void},{i_}] ]IndexDerivative[Tensor[x_,{Void},{j_}],F_]->Kronecker[i,j]*)
	(*\:2460\:3068\:2461\:3067\:5b9f\:73fe*)	
			};
       w //.Regla

		](*by msh*)


(*...........................................*)
(*....MetricSimplify[w_,g_].................*)
(*....\:8a08\:91cf\:30c6\:30f3\:30bd\:30eb\:306b\:3088\:308b\:964d\:968e\:3001\:6607\:968e..\:5bfe\:79f0\:30c6\:30f3\:30bd\:30eb......*)
(*...........................................*)
	
MetricSimplify[w_,g_]:=Module[{Regla},
TensorMetricRuleSimplify = { 
(*     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,i,j]/;MemberAllQ[w,i]&& FreeAllQ[w,IndexDerivative] ,
     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,j,i]/;MemberAllQ[w,j]&& FreeAllQ[w,IndexDerivative]  ,*)
(*\:964d\:683c by msh*)     
     Tensor[y_,sub_,sup_]Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[Tensor[y,sub,sup],i,j]/;MemberAllQ[sup,i],
     Tensor[y_,sub_,sup_]Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[Tensor[y,sub,sup],j,i]/;MemberAllQ[sup,j],
(*\:6607\:683c by msh*) 
	Tensor[y_,sub_,sup_] Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[Tensor[y,sub,sup],i,j]/;MemberAllQ[sub,i] ,
    Tensor[y_,sub_,sup_] Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[Tensor[y,sub,sup],j,i]/;MemberAllQ[sub,j],
      
	Tensor[g,{a_,i_},{Void,Void}]Tensor[g,{Void,Void},{a_,j_}] :> Kronecker[i,j],
	Tensor[g,{i_,a_},{Void,Void}]Tensor[g,{Void,Void},{j_,a_}] :> Kronecker[i,j],
	Tensor[g,{a_,i_},{Void,Void}]Tensor[g,{Void,Void},{j_,a_}] :> Kronecker[i,j],
	(*Kronecker[i_,i_]:>1,*)
	
    Tensor[g,{i_,Void},{Void,j_}] :> Kronecker[i,j],
    Tensor[g,{Void,i_},{j_,Void}] :> Kronecker[i,j],
     
(*    IndexDerivative[Tensor[g,{i_,p_},{Void,Void}],x[s_]]Tensor[g,{Void,Void},{s_,p_}]
   -IndexDerivative[Tensor[g,{s_,i_},{Void,Void}],x[p_]]Tensor[g,{Void,Void},{s_,p_}]:>0,*)
       (IndexDerivative[Tensor[g,{i_,p_},{Void,Void}],x[s_]]
   -IndexDerivative[Tensor[g,{s_,i_},{Void,Void}],x[p_]])Tensor[g,{Void,Void},{s_,p_}]:>0,
          (IndexDerivative[Tensor[g,{i_,p_},{Void,Void}],x[s_]]
   -IndexDerivative[Tensor[g,{s_,i_},{Void,Void}],x[p_]])Tensor[g,{Void,Void},{p_,s_}]:>0,
          (IndexDerivative[Tensor[g,{p_,i_},{Void,Void}],x[s_]]
   -IndexDerivative[Tensor[g,{i_,s_},{Void,Void}],x[p_]])Tensor[g,{Void,Void},{s_,p_}]:>0,
          (IndexDerivative[Tensor[g,{p_,i_},{Void,Void}],x[s_]]
   -IndexDerivative[Tensor[g,{i_,s_},{Void,Void}],x[p_]])Tensor[g,{Void,Void},{p_,s_}]:>0,
    IndexDerivative[Tensor[g,{\[Mu]_,\[Nu]_},{Void,Void}],x[j_]]*IndexDerivative[Tensor[g,{Void,Void},{\[Mu]_,\[Nu]_}],x[i_]]
   -IndexDerivative[Tensor[g,{\[Mu]_,\[Nu]_},{Void,Void}],x[i_]]*IndexDerivative[Tensor[g,{Void,Void},{\[Mu]_,\[Nu]_}],x[j_]]:>0,
   
   RieChrisCur[{m_},{k_,l_},{k_},g,x_]:>RieChrisCur[{m},{l},{Void},g,x],
   RieChrisCur[{m_},{l_,k_},{k_},g,x_]:>-RieChrisCur[{m},{l},{Void},g,x],

   RieChrisCur[{m_},{l_},{Void},g,x_]Tensor[g,{Void,Void},{m_,l_}] :>RieChrisCur[{Void},{Void},{Void},g,x],
   RieChrisCur[{j_},{k_,l_},{n_},g,x_]Tensor[g,{i_,n_},{Void,Void}]:>RieChrisCur[{i,j},{k,l},{Void},g,x],
   
   RieChrisCur[{m_},{n_},{Void},g,x_]Tensor[g,{Void,Void},{i_,m_}]Tensor[g,{Void,Void},{j_,n_}]:>
   RieChrisCur[{Void},{Void},{i,j},g,x],
   
   RieChrisCur[{k_},{j_},{Void},g,x_] Tensor[g,{Void,Void},{k_,i_}]:>RieChrisCur[{Void},{j},{i},g,x],
   RieChrisCur[{k_},{j_},{Void},g,x_]Tensor[g,{Void,Void},{i_,k_}]:>RieChrisCur[{Void},{j},{i},g,x],
   
   RieChrisCur[{j_},{k_,l_},{m_},g,x_]Tensor[g,{i_,m_},{Void,Void}]:>RieChrisCur[{i,j},{k,l},{Void},g,x],
   RieChrisCur[{j_},{l_,k_},{m_},g,x_]Tensor[g,{i_,m_},{Void,Void}]:>RieChrisCur[{i,j},{k,l},{Void},g,x],
   
   RieChrisCur[{h_,m_},{i_,k_},{Void},g,x_]Tensor[g,{Void,Void},{m_,k_}]:>RieChrisCur[{h},{i},{Void},g,x],
   RieChrisCur[{h_,m_},{i_,k_},{Void},g,x_]Tensor[g,{Void,Void},{k_,m_}]:>RieChrisCur[{h},{i},{Void},g,x],
    RieChrisCur[{h_,m_},{k_,i_},{Void},g,x_]Tensor[g,{Void,Void},{m_,k_}]:>-RieChrisCur[{h},{i},{Void},g,x],
   RieChrisCur[{h_,m_},{k_,i_},{Void},g,x_]Tensor[g,{Void,Void},{k_,m_}]:>-RieChrisCur[{h},{i},{Void},g,x],
   
   RieChrisCur[{m_,h_},{i_,k_},{Void},g,x_]Tensor[g,{Void,Void},{m_,k_}]:>-RieChrisCur[{h},{i},{Void},g,x],
   RieChrisCur[{m_,h_},{i_,k_},{Void},g,x_]Tensor[g,{Void,Void},{k_,m_}]:>-RieChrisCur[{h},{i},{Void},g,x],
    RieChrisCur[{m_,h_},{k_,i_},{Void},g,x_]Tensor[g,{Void,Void},{m_,k_}]:>RieChrisCur[{h},{i},{Void},g,x],
   RieChrisCur[{m_,h_},{k_,i_},{Void},g,x_]Tensor[g,{Void,Void},{k_,m_}]:>RieChrisCur[{h},{i},{Void},g,x],
   
   
   RieChrisCur[{Void},{i_},{l_},g,x_]Tensor[g,{Void,Void},{i_,j_}]:>RieChrisCur[{Void},{Void},{l,j},g,x],
   RieChrisCur[{Void},{i_},{l_},g,x_]Tensor[g,{Void,Void},{j_,i_}]:>RieChrisCur[{Void},{Void},{l,j},g,x]
                           
        };
	Regla = {Sequence@@TensorMetricRuleSimplify};
        w //.Regla
		](*by msh*)        



(*............................................*)
(*...........KroneckerSimplify............*)
(*............................................*)
KroneckerSimplify[w_]:=Module[{Regla},
	Regla = {Sequence@@KroneckerRuleSimplify};
        w //.Regla
		](*by msh*)

Format[ Kronecker[x_,y_] ]:=  
 Tensor[\[Delta],{x},{y}]; 

Kronecker[i_?NumericQ, i_?NumericQ] :=  1;
Kronecker[i_?NumericQ, j_?NumericQ] := 0;
(*Kronecker[i_, i_] := 1;*)
(*Kronecker[i_,j_] := 1/; i-j==0; i j \:304c\:540c\:3058\:30b7\:30f3\:30dc\:30eb\:306e\:5834\:5408\:3000\:ff11\:3068\:3059\:308b\:3000by msh*)
(*Kronecker[i_,j_] := 0/;i-j!=0; i i+1 \:304c\:540c\:3058\:30b7\:30f3\:30dc\:30eb\:306e\:5834\:5408\:30000\:3068\:3059\:308b\:3000by msh*)		
KroneckerRuleSimplify = 
{
(*T_ Kronecker[i_, i_] :> T,*)
T_ Kronecker[i_,j_]:> ( T/. i->j)/; TensorRank[Position[   T , i]]>=2,
T_ Kronecker[i_,j_]:> ( T/. j->i)/;TensorRank[Position[   T , j]]>=2
		};
(*..................*)

VectorNumericQ[w_List]:= And@@Map[ NumericQ , w ];

Tensor[ \[Epsilon] ,x_ ,{Void...}]  := 
      Signature[x]/; VectorNumericQ[x];

Tensor[ \[Epsilon] ,{Void...},x_]  := 
      Signature[x]/; VectorNumericQ[x];

Kronecker[x_List,y_List] := KroneckerDelta[x,y]
(*Tensor[ \[Epsilon],x, Table[Void,{NDim}] ]*
Tensor[\[Epsilon],Table[Void,{NDim}],y]/;VectorNumericQ[x]&&VectorNumericQ[y];
*)
(*...........................................*)
(*...............SetTensorRule...............*)
(*..........{Subsuperscript[A, Void, 1]\[Rule]X,Subsuperscript[A, Void, 2]\[Rule]Y,Subsuperscript[A, Void, 3]\[Rule]Z}..........*)

(*SetTensorRule[x_,y_]:= MapThread[Rule,{x,y}, TensorRank[y]  ]//Flatten*)	        
SetTensorRule[x_,y_]:= MapThread[Rule,{x,y}, TensorRank[x]  ]//Flatten



(*............................................*)
(*..........ToCovariantIndex..................*)
(*.........\:5171\:5909\:30c6\:30f3\:30bd\:30eb\:306b\:5909\:66f4.....................*)

ToCovariantIndex[w_,ind_,toind_]:= w /;FreeAllQ[w,ind];

ToCovariantIndex[ Tensor[w_,sub_,sup_],ind_,toind_ ]:=Tensor[w,sub,sup]/;FreeQ[sup,ind];
ToCovariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Module[{pos,SUB,SUP},
 pos =		Position[sup,ind];
 SUB = ReplacePart[sub,toind,pos];
 SUP = ReplacePart[sup,Void,pos];
		 Tensor[w,SUB,SUP]
		];

ToCovariantIndex[ T_Times  ,ind__ ]:=
  Map[ToCovariantIndex[ #  ,ind ]&, T ];

ToCovariantIndex[ w_Plus  ,ind__ ]:= 
    Map[ToCovariantIndex[ #  ,ind ]&,  Expand[w] ];

ToCovariantIndex[ w_  ,ind_List ,toind_List]:=
Fold[
ToCovariantIndex[ #1  ,Part[#2,1],Part[#2,2] ]&,w, Transpose[{ind ,toind}]];


(* ::Code::Initialization::"Tags"-><|"UppercasePattern" -> <|Enabled -> False|>|>:: *)
(*............................................*)
(*..........ToContravariantIndex..................*)
(*.........\:53cd\:5909\:30c6\:30f3\:30bd\:30eb\:306b\:5909\:66f4.....................*)
ToContravariantIndex[w_,ind_,toind_]:= w /;FreeAllQ[w,ind];

ToContravariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Tensor[w,sub,sup]/;FreeQ[sub,ind];

ToContravariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Module[{pos,SUB,SUP},
 pos =		Position[sub,ind];
 SUB = ReplacePart[sub,Void,pos];
 SUP = ReplacePart[sup,toind,pos];
		 Tensor[w,SUB,SUP]
		]/; MemberQ[sub,ind];

ToContravariantIndex[ T_Times  ,ind__ ]:=
  Map[ToContravariantIndex[ #  ,ind ]&, T ];

ToContravariantIndex[ w_Plus  ,ind__ ]:= 
    Map[ToContravariantIndex[ #  ,ind ]&,  Expand[w] ];

ToContravariantIndex[ w_  ,ind_List ,toind_List]:=
Fold[
ToContravariantIndex[ #1  ,Part[#2,1],Part[#2,2] ]&,w, Transpose[{ind ,toind}]];


(*............................................*)	
(*...............IndexDerivativeUp.............*)
(*............................................*)
mshsize2=20;
(*Format[ IndexDerivativeUp[ w_ , ind_]  ]:=*)   
(*Style[Subsuperscript["\[PartialD]", SequenceForm@@w , ind],Black,Italic,mshsize2,FontFamily->mshfont];*)
Format[ IndexDerivativeUp[ w_ , ind_,g_,x_]]:=   
(*Style[Superscript["\[PartialD]",Row[{ind}]]w,Black,Italic,mshsize2,FontFamily->mshfont];*)
Style[SequenceForm[Superscript["\[PartialD]",Row[{ind},","]],w],Black,Italic,mshsize2,FontFamily->mshfont];
ExpandIndexDerivativeUp[w_]:=  w /. IndexDerivativeUp -> IDXUP;
IDXUP[ w_, m_,g_,x_]:=Tensor[g,{Void,Void},{m,\[Rho]}]IndexDerivative[w,x[\[Rho]]]


(*............................................*)	
(*...............IndexDerivative.............*)
(*............................................*)
(*Attributes[IndexDerivative]={HoldFirst};(*by msh *)
Attributes[IDXX]={HoldFirst};(*by msh *)*)
(*Format[ IndexDerivative[ w_ , ind__] ] :=
Style[SequenceForm[Subscript["\[PartialD]",Row[{ind},","]],w],Black,Italic,mshsize2,FontFamily->mshfont];*)
(*Format[ IndexDerivative[ w_ , ind__] ] := Style[Times[Row[{"\[PartialD]",w}],Power[Row[{"\[PartialD]",Row[{ind},"\[PartialD]"]}],-1]]
,Black,Italic,18,FontFamily->"Times New Roman"];*)
(*Format[ IndexDerivative[ w_ , ind__] ] = HoldForm[ D[ w , ind ] ];*)


Format[ IndexDerivative[ w_ , ind__] ] := Style[SequenceForm["",HoldForm[ D[ w , ind ] ],""],Black,Italic,mshsize2,FontFamily->mshfont];
IndexDerivative[IndexDerivative[F_,u__], v__ ] := IndexDerivative[ F,v,u];
IndexDerivative[z_,z_] := 1;
IndexDerivative[(w_)^n_,z_] := n w^(n-1) IndexDerivative[w,z];
IndexDerivative[Norm[w_],t_] := w . IndexDerivative[w,t]/Norm[w];
(*IndexDerivative[(w_)^-1,z_] := -1/w^2 IndexDerivative[w,z];*)
IndexDerivative[  w_Plus , z_] := IndexDerivative[ # , z]& /@ w ;(*by msh *)
IndexDerivative[ u_*v_ , z_] :=  IndexDerivative[u,z]v + u IndexDerivative[v,z];(*by msh *)
IndexDerivative[ Dot[u_,v_] , z_] :=  v . IndexDerivative[u,z] + u . IndexDerivative[v,z];(*by msh *)
IndexDerivative[a_List,b_]:=Map[IndexDerivative[#,b]&,a];(*by msh *)
IndexDerivative[ F_[u_] , v_ ] := Derivative[1][F][u] IndexDerivative[u,v] /; (u=!=v && F=!=Tensor); 
IndexDerivative[ Exp[u_] , v_ ] := Exp[u] IndexDerivative[u,v] /; (u=!=v); 
IndexDerivative[ Exp[u_] , u_ ] := Exp[u] ; 
(*EvaluateIndexDerivative[w__]:=  w /. IndexDerivative -> D;*)
EvaluateIndexDerivative[w_]:=  w /. IndexDerivative -> IDXX /. IDXX -> D;
IndexDerivative[ x_?NumericQ ,w_] = 0;
IDXX[u_?NumericQ,w_]=0;



(* IDXX[u_,u_]=1;

IDXX[ F_[u_] , v_ ] := 
   IDXX[F[u],u] IDXX[u,v] /; (u=!=v && F=!=Tensor); 
   	
IDXX[ u_*v_ , z_] =  IDXX[u,z]v + u IDXX[v,z];

IDXX[ u_ + v_ , z_] =  IDXX[u,z] + IDXX[v,z];

IDXX[u_ , x_, y__] := IDXX[ IDXX[u,x] , y ];deleate by msh *)
(*		
IDXX[Power[u_,m_] ,v_ ] = m Power[u,m-1]IDXX[u,v]+ Power[u,m]Log[u]IDXX[m,v];

IDXX[ ArcCos[u_] , v_ ] = -1/Sqrt[1-u^2] IDXX[u,v];

IDXX[ ArcSin[u_] , v_ ] = 1/Sqrt[1-u^2] IDXX[u,v];

IDXX[ Sin[u_] , v_ ] = Cos[u] IDXX[u,v];

IDXX[ Cos[u_] , v_ ] = -Sin[u] IDXX[u,v];

IDXX[ Log[u_] , v_ ] = 1/u IDXX[u,v]; 
IDXX[ u_ , x_ ]:=  0/; Position[u,Tensor]=={};    deleate by msh *)


(*............................................*)
(*.............Covariant Derivative...........*)
(*...........................modified by msh..*)

Attributes[CovariantDerivative]={HoldFirst};
(*Attributes[CDCov]={HoldFirst};*)

Format[ CovariantDerivative[T_, x_ , ind_ , g_]  ]:= 
(*Subscript[  T  ,SequenceForm[ ",",SequenceForm@@ind ] ];*)
  SequenceForm[ Tensor["\[Del]",{SequenceForm@@ind},Void],T];(*by msh*)
  
  Format[ CovariantDerivative[T_, x_ , {k_,l_} , g_]  ]:= 
   Style[SequenceForm[ Tensor["\[Del]",{SequenceForm@@k},{Void}],"(",Tensor["\[Del]",{SequenceForm@@l},{Void}],T,")"],Black,Italic,18,FontFamily->mshfont];(*by msh*)
(*Format[ 
CovariantDerivative[  CovariantDerivative[ T_, x_ , ind1_ , g_] , x_ , ind2_ , g_] ]:=
(*Subscript[  T  ,SequenceForm[ ",",SequenceForm@@{ind1,ind2} ] ];*)  
 (* SequenceForm[ Tensor["\[Del]",{SequenceForm@@{ind1,ind2}},{Void,Void}],T];*)(*by msh*)
    SequenceForm[ Tensor["\[Del]",{SequenceForm@@ind2},{Void}],"(",Tensor["\[Del]",{SequenceForm@@ind1},{Void}],T,")"];(*by msh*)

Format[ 
CovariantDerivative[  
   CovariantDerivative[ 
        CovariantDerivative[ T_, x_ , ind0_ , g_], x_ , ind1_ , g_] , x_ , ind2_ , g_] ]:=
 (*Subscript[  T  ,SequenceForm[ ",",SequenceForm@@{ind0,ind1,ind2} ] ];*) 
  SequenceForm[ Tensor["\[Del]",{SequenceForm@@{ind0,ind1,ind2}},Void],T]; (*by msh*)*)
(*.....................*)
(*(*Fold:\:95a2\:6570 f \:3092\:9023\:7d9a\:7684\:306b\:9069\:7528\:3057\:3066 x \:304a\:3088\:3073\:30ea\:30b9\:30c8\:8981\:7d20\:3092\:30b7\:30fc\:30c9\:3059\:308b*)
CovariantDerivative[  w_ ,  x_ , j_List, g_ ] :=
     Fold[ CovariantDerivative[ #1 , x , #2 ,g ]& , w , j ];*)
     
CovariantDerivative[  CovariantDerivative[ T_, x_ , ind1_ , g_] , x_ , ind2_ , g_]:=                
CovariantDerivative[ T, x , {ind1,ind2} , g];   
(*w_Plus: a_+b_ -> CovariantDerivative[a,x,j,g]+CovariantDerivative[b,x,j,g]*)
CovariantDerivative[  w_Plus , x_ ,j_ ,g_ ]:= 
     CovariantDerivative[ # , x,j,g]& /@ w;

CovariantDerivative[  w_*u_ ,  x_ , j_ , g_ ] := 
     CovariantDerivative[  w ,  x ,j , g   ]* u + 
    CovariantDerivative[  u ,  x ,j , g   ]* w ;

(*CovariantDerivative[  w_?NumericQ , Tensor[x_,{Void},{j_}], g_ ] = 0; *)
CovariantDerivative[ w_?NumericQ , x_,k_, g_ ] = 0; 
(*CovariantDerivative[-1,x_,k_,g_]=0;*)
CovariantDerivative[  Kronecker[i_,j_] ,  __ ] = 0;

CovariantDerivative[ Tensor[g_,{i_,j_},{Void,Void}] ,_,_,g_ ] = 0;(*by msh *)
CovariantDerivative[ Tensor[g_,{Void,Void},{i_,j_}] ,_,_,g_ ] = 0;
CovariantDerivative[ Tensor[g_,{i_,j_},{" "," "}]   ,_,_,g_ ] = 0;
CovariantDerivative[ Tensor[g_,{" "," "},{i_,j_}]   ,_,_,g_ ] = 0;
CovariantDerivative[ Tensor[T_,{},{}] , x_,i_,g_] := 
         IndexDerivative[Tensor[T,{},{}] , Tensor[x,{Void},{i}] ] ;
                                                                        
(*......................*)
(*......................*)
(*CDCov[IndexDerivative[Tensor[T_,sub_,sup_],Tensor[x_,{_},{l_}]],x_,k_,g_]:=
IndexDerivative[Tensor[T,sub,sup],x[l],x[k]]+
Christoffel2[g,x][{k,\[Sigma]},sup](IndexDerivative[Tensor[T,Void,\[Sigma]],x[l]]);
*)
(*CDCov[Christoffel2[g_,x_][{l_,m_},{n_}],x_,k_,g_]:=
(IndexDerivative[Christoffel2[g,x][{l,m},{n}],x[k]])*)   


ExpandCODer[w_]:= w /. CovariantDerivative->CDCov; 
CDCov[ Tensor[w_,{Void},{Void}] , x_ , m_ , g_]:=IndexDerivative[Tensor[w,{Void},{Void}],x[m]];

CDCov[ IndexDerivative[Tensor[w_,{Void},{Void}],Tensor[x_,{Void},{n_}]] , x_ , m_ , g_]:=
IndexDerivative[Tensor[w,{Void},{Void}],x[m],x[n]]-
IndexDerivative[Tensor[w,{Void},{Void}],x[\[Rho]]] Christoffel2[g,x][{m,n},{\[Rho]}];

CDCov[ IndexDerivative[Tensor[w_,{Void},{Void}],t_] , x_ , m_ , g_]:=
IndexDerivative[Tensor[w,{Void},{Void}],x[m],t]/;FreeTensorQ[t];

CDCov[ IndexDerivative[Tensor[w_,{Void},{i_}],Tensor[x_,{Void},{n_}]] , x_ , j_ , g_]:=
IndexDerivative[Tensor[w,{Void},{i}],x[j],Tensor[x,{Void},{n}]]-
IndexDerivative[Tensor[w,{Void},{i}],Tensor[x,{Void},{\[Rho]}]]Christoffel2[g,x][{j,n},{\[Rho]}]+
IndexDerivative[Tensor[w,{Void},{\[Rho]}],Tensor[x,{Void},{n}]]Christoffel2[g,x][{j,\[Rho]},{i}];

CDCov[ IndexDerivative[Tensor[w_,{Void},{i_}],t_] , x_ , j_ , g_]:=IndexDerivative[Tensor[w,{Void},{i}],x[j],t]+
IndexDerivative[Tensor[w,{Void},{\[Rho]}],t]Christoffel2[g,x][{j,\[Rho]},{i}]/;FreeTensorQ[t];
(*/;FreeTensorQ[t_]*)
(*\:6ce8\:610f\:3000\:639b\:3051\:7b97\:306e\:6642\:540c\:3058\:30c0\:30df\:30fc\:5909\:6570\[Rho]\:3068\:306a\:308b\:306e\:3067\:5909\:66f4\:3059\:308b\:5fc5\:8981\:304c\:3042\:308a*)
CDCov[ Tensor[T_,sub_,sup_] , x_ , i_ , g_]:=
Module[{DP,xx,ind},
	ind = Transpose[{ sub , sup , Range[Length[sub]](*,Table[mshtensorial`dk[p],{p,1,Length[sub]}]*) }];
	(*ind={{a,Void,1},      sub={a,Void,c}  sup={Void,b,Void} 
	       {Void,b,2},
	       {c,Void,3}}   *)
	
	xx = Tensor[x,{Void},{i}];
	(*\[PartialD]x[i]+CDCov:\[CapitalGamma]m km *)
   IndexDerivative[Tensor[T,sub,sup] , Tensor[x,{Void},{i}]  ] +
   (*+CDCov:\[CapitalGamma]m km *)
   Plus@@Map[sCDCov[ Tensor[T,sub,sup] ,#,xx,g]& ,ind]
		]/;ListQ[i]==False;
		(*/;AllNumberQ[{sub,sup,i}];*)(*\:6570\:5024\:3067\:3042\:308b\:306a\:3089\:5b9f\:884c*)
(************)
sCDCov[ Tensor[T_,subT_,supT_] , {sub_,Void,n_}, Tensor[x_,_,{i_}] ,g_]:=
      -Tensor[T,ReplacePart[subT,\[Rho],n],supT] * Christoffel2[g,x][{i,sub},{\[Rho]}]; 
       (*ReplacePart[subT,k,n]  {sub}\:306en\:756a\:76ee\:3092k\:306b\:7f6e\:304d\:63db\:3048\:308b *)    
sCDCov[ Tensor[T_,subT_,supT_] , {Void,sup_,n_},Tensor[x_,_,{i_}],g_]:=
      Tensor[T,subT, ReplacePart[supT,\[Rho],n] ]*Christoffel2[g,x][{i,\[Rho]},{sup}] ;
     


Expand2ndCODer[w_]:= w /. CovariantDerivative->Covariant2ndDerivative;

Covariant2ndDerivative[ Tensor[T_,sub_,sup_] ,  x_ , {k_,l_}, g_ ] :=
Module[{DP,xx,ind},
	ind = Transpose[{ sub , sup , Range[Length[sub]]}];(*,Table[mshtensorial`dk[p],{p,1,Length[sub]}]*) 
	(*ind={{a,Void,1},      sub={a,Void,c}  sup={Void,b,Void} 
	       {Void,b,2},
	       {c,Void,3}}   *)
	
	xx = Tensor[x,{Void,Void},{k,l}];
	(*\[PartialD]x[i]+CDCov:\[CapitalGamma]m km *)
IndexDerivative[CovariantDerivative[Tensor[T,sub,sup],x, l, g],Tensor[x,{Void},{k}]]-
Christoffel2[g,x][{k,l},{\[Sigma]}]CovariantDerivative[Tensor[T,sub,sup],x, \[Sigma], g]+
   (*+CDCov:\[CapitalGamma]m km *)
   Plus@@Map[SndCDCov[ Tensor[T,sub,sup] ,#,xx,g]& ,ind]
		];
          
SndCDCov[ Tensor[T_,subT_,supT_] , {sub_,Void,n_}, Tensor[x_,_,{k_,l_}] ,g_]:=
      - CovariantDerivative[Tensor[T,ReplacePart[subT,\[Sigma],n],supT],x,l,g ]* Christoffel2[g,x][{k,sub},{\[Sigma]}]; 
       (*ReplacePart[subT,k,n]  {sub}\:306en\:756a\:76ee\:3092k\:306b\:7f6e\:304d\:63db\:3048\:308b *)    
SndCDCov[ Tensor[T_,subT_,supT_] , {Void,sup_,n_},Tensor[x_,_,{k_,l_}],g_]:=
      CovariantDerivative[Tensor[T,subT, ReplacePart[supT,\[Sigma],n]],x,l,g]* Christoffel2[g,x][{k,\[Sigma]},{sup}] ;
                                                                
(*Covariant2ndDerivative[ Tensor[T_,{Void},{n_}] ,  x_ , {k_,l_}, g_ ] :=
IndexDerivative[CovariantDerivative[Tensor[T,{Void},{n}],x, l, g],Tensor[x,{Void},{k}]]-
Christoffel2[g,x][{k,l},{\[Sigma]}]CovariantDerivative[Tensor[T,{Void},{n}],x, \[Sigma], g]+   
Christoffel2[g,x][{k,\[Sigma]},{n}]CovariantDerivative[Tensor[T,{Void},{\[Sigma]}],x, l, g];    

Covariant2ndDerivative[ Tensor[T_,{n_},{Void}] ,  x_ , {k_,l_}, g_ ] :=
IndexDerivative[CovariantDerivative[Tensor[T,{n},{Void}],x, l, g],Tensor[x,{Void},{k}]]-
Christoffel2[g,x][{k,l},{\[Sigma]}]CovariantDerivative[Tensor[T,{n},{Void}],x, \[Sigma], g]-   
Christoffel2[g,x][{k,n},{\[Sigma]}]CovariantDerivative[Tensor[T,{\[Sigma]},{Void}],x, l, g];      *)         
      


(*............................................*)
(*.............Absolut Derivative.............*)
(*...........................modified by msh..*)

Attributes[AbsoluteDerivative]={HoldFirst};
Format[AbsoluteDerivative[ w_ ,r_ ,s_,g_,x_]]:= 
(*  Times[ SequenceForm[ \[Delta] ,w ] , Power[ SequenceForm[\[Delta],s],-1]];*)
Style[Times[ SequenceForm[ \[Delta] ,w ] , Power[ SequenceForm[\[Delta],s],-1]],Black,Plain,18,FontFamily->mshfont];

AbsoluteDerivative[u_*w_,r_ ,s_,g_,x_]= 
  u AbsoluteDerivative[w,r,s,g,x]+w AbsoluteDerivative[u,r,s,g,x];
AbsoluteDerivative[  w_Plus , r_ ,s_,g_,x_] := AbsoluteDerivative[ # ,r,s,x,g]& /@ w ;

ExpandABSDer[w_]:= w /. AbsoluteDerivative -> ABCov;
ABCov[ w_ ,r_,s_,g_,x_]:= Module[{aa},
	CovariantDerivative[w,x,r,g]IndexDerivative[Tensor[x,{Void},{r}],s]
	];


(*............................................*)
(*................Dummy Variables.............*)
(*............................................*)

DoubleVoidRule = 
  Tensor[P_,{subP__},{supP1___,Void,supP2___}] -> Tensor[P,{subP},{supP1,fXX[Void],supP2}] ;

DoubleIndexSimplifyRule01=	
a_.*Tensor[A_,{subA__},{supA1___,i_,supA2___}]*Tensor[B_,{subB1___,i_,subB2___},{supB___}] +
b_.*Tensor[A_,{subAA__},{supAA1___,j_,supAA2___}]*Tensor[B_,{subBB1___,j_,subBB2___},{supBB___}] :> 
a*Tensor[A,{subA},{supA1,i,supA2}]  Tensor[B,{subB1,i,subB2},{supB}]+
b*Tensor[A,{subAA},{supAA1,fXX[i],supAA2} ]*Tensor[B,{subBB1,i,subBB2},{supBB}] ;

(*............*)
DoubleIndexSimplifyRule02 =
a_.*Tensor[ A_,{subA1___,i_,subA2___},{supA1___,i_,supA2___}] + 
b_.*Tensor[ A_,{subAA1___,k_,subAA2___},{supAA1___,k_,supAA2___}]	:>
a*Tensor[ A,{subA1,i,subA2},{supA1,i,supA2}]+ 
b*Tensor[ A,{subAA1,fXX[i],subAA2},{supAA1,i,supAA2}]/;Length[{subA2}]===Length[{subAA2}];

(*..............*)
DoubleIndexSimplifyRule03=	
a_. Tensor[A_,{subA__},{supA1___,i_,supA2___}]*Tensor[B_,{subB1___,i_,subB2___},{supB___}]+
b_. Tensor[A_,{subAA__},{supAA1___,j_,supAA2___}]*Tensor[B_,{subBB1___,j_,subBB2___},{supBB___}]* 
Tensor[P_,{subP__},{supP1___,i_,supP2___}]*Tensor[Q_,{subQ1___,i_,subQ2___},{supQ___}] :> 
a*Tensor[A,{subA},{supA1,i,supA2}]  Tensor[B,{subB1,i,subB2},{supB}] +
b*Tensor[A,{subAA},{supAA1,fXX[i],supAA2} ]*Tensor[B,{subBB1,i,subBB2},{supBB}]*Tensor[P,{subP},{supP1,j,supP2}]*  
Tensor[Q,{subQ1,j,subQ2},{supQ}] ;

(*............*)
DoubleIndexSimplifyRule04 =
a_.*Tensor[ A_,{subA1___,i_,subA2___},{supA1___,i_,supA2___}]+ 
b_.* Tensor[ A_,{subAA1___,k_,subAA2___},{supAA1___,k_,supAA2___}]*Tensor[P_,{subP__},{supP1___,i_,supP2___}]*  
Tensor[Q_,{subQ1___,i_,subQ2___},{supQ___}]  :> a*Tensor[ A,{subA1,i,subA2},{supA1,i,supA2}]+ 
b*Tensor[ A,{subAA1,fXX[i],subAA2},{supAA1,i,supAA2}]* Tensor[P,{subP},{supP1,k ,supP2}]*  
Tensor[Q,{subQ1,k,subQ2},{supQ}]/;Length[{subA2}]===Length[{subAA2}];	
DummyVariableSimplify[w_]:=Module[{},
Expand[w]//.DoubleVoidRule//.DoubleIndexSimplifyRule04//.
DoubleIndexSimplifyRule02//.DoubleIndexSimplifyRule03//.DoubleIndexSimplifyRule01/.fXX->Times];


(*............................................*)
(*.............Christofel....................*)
(*............................................*)

Format[Christoffel1[g_,x_][{i_,j_},{k_}]]=
(*  Style[SequenceForm[ Subscript[\[CapitalGamma],k,SequenceForm[i,j]],"[",g,"]"],Black,Italic,mshsize,FontFamily->mshfont];*)
    Style[SequenceForm[ Subscript[\[CapitalGamma],k,SequenceForm[i," ",j]]],Black,Plain,mshsize,FontFamily->mshfont];
    
Format[\!\(\*OverscriptBox[\(Christoffel1\), \(_\)]\)[g_,x_][{i_,j_},{k_}]]=
    Style[SequenceForm[ Subscript[\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(_\)]\),k,SequenceForm[i," ",j]]],Black,Plain,mshsize,FontFamily->mshfont];  

Format[\!\(\*OverscriptBox[\(Christoffel2\), \(_\)]\)[g_,x_][{i_,j_},{k_}]]=
  Style[Subsuperscript[\!\(\*OverscriptBox[\(\[CapitalGamma]\), \(_\)]\), SequenceForm[i," ",j], SequenceForm[" ",k]],Black,Plain,mshsize,FontFamily->mshfont];   
 (*       SequenceForm[ Tensor[Overscript[\[CapitalGamma], _],{i," ",j},{Void,k}]];*) 

Format[Christoffel2[g_,x_][{i_,j_},{k_}]]=
   Style[Subsuperscript[\[CapitalGamma], SequenceForm[i," ",j], SequenceForm[" ",k]],Black,Plain,mshsize,FontFamily->mshfont];   
       (* SequenceForm[ Tensor[\[CapitalGamma],{i,j},{Void,k}]];*)
                 
                        
Christoffel1[g_,x_][{i_,j_},{k_}]:= Module[{p},
(*x[p_] = Tensor[x,{Void},{p}];*)
		
1/2( IndexDerivative[  Tensor[g,{i,k},{Void,Void}]  , Tensor[x,{Void},{j}]] +
     IndexDerivative[  Tensor[g,{j,k},{Void,Void}] ,Tensor[x,{Void},{i}]]-
     IndexDerivative[  Tensor[g,{i,j},{Void,Void}] , Tensor[x,{Void},{k}]] )		
		]/;IntegerQ[i]&&IntegerQ[j]&&IntegerQ[k];
(*..............*)

Christoffel2[g_,x_][{i_,j_},{k_}]:= Module[{aa},
ExpandDoubleIndex[
	Tensor[g,{Void,Void},{k,aa}] Christoffel1[g,x][{i,j},{aa}]      
       ,{aa}]
		]/;IntegerQ[i]&&IntegerQ[j]&&IntegerQ[k];
		
ExpandChristoffel2[w_]:=  w /.{Christoffel2 -> ExpCri2,\!\(\*OverscriptBox[\(Christoffel2\), \(_\)]\) -> ExpCri2};(*by msh*)
	
ExpCri2[g_,x_][{i_,j_},{k_}]:= Module[{aa},
	Tensor[g,{Void,Void},{k,\[Sigma]}] Christoffel1[g,x][{i,j},{\[Sigma]}]      
       ];(*\:6ce8\:610f\:3000\:639b\:3051\:7b97\:306e\:6642\:540c\:3058\:30c0\:30df\:30fc\:5909\:6570\[Sigma]\:3068\:306a\:308b\:306e\:3067\:5909\:66f4\:3059\:308b\:5fc5\:8981\:304c\:3042\:308a*)
       
ExpandChristoffel1[w_]:=  w /.{Christoffel1 -> ExpCri1,\!\(\*OverscriptBox[\(Christoffel1\), \(_\)]\) -> ExpCri1};(*by msh*)	
	
ExpCri1[g_,x_][{i_,j_},{k_}]:=Module[{p},
(*x[p_] = Tensor[x,{Void},{p}];*)
		
1/2( IndexDerivative[  Tensor[g,{i,k},{Void,Void}],Tensor[x,{Void},{j}]] +
     IndexDerivative[  Tensor[g,{j,k},{Void,Void}],Tensor[x,{Void},{i}]]-
     IndexDerivative[  Tensor[g,{i,j},{Void,Void}],Tensor[x,{Void},{k}]] )		
		];
		
OrderChristoffel[w_]:=Module[{Regla},
	Regla = {Sequence@@OrdCri};
        w /.Regla
		];
		
OrdCri={Christoffel2[g_,x_][{j_,k_},{i_}] :> Christoffel2[g,x][{k,j},{i}]/;Order[j,k]==-1,
Christoffel1[g_,x_][{j_,k_},{i_}] :> Christoffel1[g,x][{k,j},{i}]/;Order[j,k]==-1}		
(*OrdCri2={Christoffel2[g_,x_][{j_,k_},{i_}] :> Christoffel2[g,x][{k,j},{i}]/;Order[j,k]==-1};(*by msh*)
OrdCri2={Christoffel1[g_,x_][{j_,k_},{i_}] :> Christoffel1[g,x][{k,j},{i}]/;Order[j,k]==-1};(*by msh*)*)

OrderTensor[w_,T_]:=Module[{Regla},
	OrdTen={Tensor[T,{i_,j_},{s_,t_}] :> Tensor[T,{j,i},{s,t}]/;Order[i,j]==-1,
			Tensor[T,{i_,j_},{s_,t_}] :> Tensor[T,{j,i},{t,s}]/;Order[s,t]==-1};
	Regla = {Sequence@@OrdTen};
    w //.Regla
		];

OrderIndexDeri[w_,x_]:=Module[{Regla},
	OrdIndDer={
	IndexDerivative[T_,Tensor[x,{Void},{i_}],Tensor[x,{Void},{j_}]]
	 :> IndexDerivative[T,Tensor[x,{Void},{j}],Tensor[x,{Void},{i}]]/;Order[i,j]==-1,
	IndexDerivative[T_,i_,j_]
	 :>  IndexDerivative[T,j,i]/;Order[i,j]==-1, 
	IndexDerivative[T_,Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{i_}],Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{j_}]]
	 :> IndexDerivative[T,Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{j}],Tensor[\!\(\*OverscriptBox[\(x\), \(_\)]\),{Void},{i}]]/;Order[i,j]==-1	 
	 };			
	Regla = {Sequence@@OrdIndDer};
    w //.Regla
		];
						
(*OrdTen={Tensor[T,{i_,j_},{s_,t_}] :> Tensor[T,{j,i},{s,t}]/;Order[i,j]==-1,
	Tensor[T,{i_,j_},{s_,t_}] :> Tensor[T,{j,i},{t,s}]/;Order[s,t]==-1
};*)(*by msh*)			
																								
(*............................................*)
(*.............Member Functions...............*)
(*............................................*)

FreeTensorQ[w_]:= Position[w,Tensor]=={};
FreeAllQ[w_,i_]:= Position[w,i] == {};
FreeAllQ[w_,i_List]:=  And@@Map[FreeAllQ[w,#]&,i] ;
MemberAllQ[w_,in_]:= Position[w,in]!= {};
MemberAllQ[w_,in_List]:= And@@Map[ MemberAllQ[w,#]& , in ];
AllNumberQ[w_]:= NumericQ[  Apply[Plus,Flatten[w]/.Void->0]    ];(*\:6570\:5024\:3067\:3042\:308b\:304b\:5224\:65ad*)	


(*..........................................*)
(*......Reimman-Chritoffel curvature........*)
(*.......(curvature tensor).......by msh....*)

(*Attributes[RieChrisCur]={HoldFirst};*)
(*Attributes[RCCur]={HoldFirst};*)
(*Format[RieChrisCur[ind_,{j_,k_},{a_},g_,x_]]:=
 SequenceForm[ Tensor[mshtensorial`RCC,{SequenceForm@@{Row[ind],",",j,k}},{Void,a}]];*)
(* Format[RieChrisCur[ind1_,ind2_,{a_},g_,x_]]:=
 SequenceForm[ Tensor[mshtensorial`RCC,{SequenceForm@@{Row[ind1],",",Row[ind2]}},{Void,a}]];*)
 
  Format[RieChrisCur[ind1_,ind2_,ind3_,g_,x_]]:=Style[
 (*SequenceForm[ Tensor[mshtensorial`RCC,{SequenceForm@@{Row[ind1],",",Row[ind2]}},{Void,Row[ind3]}]],Black,Italic,18,FontFamily->mshfont];*)
  SequenceForm[ Tensor["RCC",{SequenceForm@@{Row[ind1],",",Row[ind2]}},{Void,Row[ind3]}]],Black,Italic,18,FontFamily->mshfont];
 

(*RieChrisCur[{i_},{j_},{Void},g_,x_]:=RieChrisCur[{i},{\[Sigma],j},{\[Sigma]},g,x]
RieChrisCur[{Void},{j_},{i_},g_,x_]:=Tensor[g,{Void},{\[Rho],i}]RieChrisCur[{\[Rho]},{\[Sigma],j},{\[Sigma]},g,x]*)
ExpandRicci[w_]:= w /. RieChrisCur -> RicTen;(*by msh*)

RicTen[{Void},{Void},{Void},g_,x_]:=Tensor[g,{Void,Void},{\[Xi],\[Delta]}]RieChrisCur[{\[Xi]},{\[Sigma],\[Delta]},{\[Sigma]},g,x];
RicTen[{Void},{Void},{i_,j_},g_,x_]:=Tensor[g,{Void,Void},{i,\[Rho]}]Tensor[g,{Void,Void},{j,\[Sigma]}]RieChrisCur[{\[Rho]},{\[Xi],\[Sigma]},{\[Xi]},g,x];
RicTen[{i_},{j_},{Void},g_,x_]:=RieChrisCur[{i},{\[Sigma],j},{\[Sigma]},g,x];
RicTen[{Void},{j_},{i_},g_,x_]:=Tensor[g,{Void,Void},{\[Xi],i}]RieChrisCur[{\[Xi]},{\[Sigma],j},{\[Sigma]},g,x];
RicTen[{i_,j_},{k_,l_},{Void},g_,x_]:=RieChrisCur[{j},{k,l},{\[Sigma]},g,x]Tensor[g,{i,\[Sigma]},{Void,Void}];

ExpandRCC[w_]:=  w /. RieChrisCur -> RCCur;(*by msh*)
RCCur[{i_},{j_,k_},{a_},g_,x_] := 
     (*ExpandDoubleIndex[*)
IndexDerivative[ Christoffel2[g,x][{k,i},{a}] , Tensor[x,{Void},{j}] ]-
IndexDerivative[ Christoffel2[g,x][{j,i},{a}] , Tensor[x,{Void},{k}] ]+
Christoffel2[g,x][{k,i},{\[Rho]}] Christoffel2[g,x][{j,\[Rho]},{a}] -
Christoffel2[g,x][{j,i},{\[Rho]}] Christoffel2[g,x][{k,\[Rho]},{a}];

RCCur[{i_,j_},{k_,l_},{Void},g_,x_] := 
     (*ExpandDoubleIndex[*)
IndexDerivative[ Christoffel1[g,x][{j,l},{i}] , Tensor[x,{Void},{k}] ]-
IndexDerivative[ Christoffel1[g,x][{j,k},{i}] , Tensor[x,{Void},{l}] ]+
Christoffel1[g,x][{i,l},{\[Rho]}] Christoffel2[g,x][{j,k},{\[Rho]}] -
Christoffel1[g,x][{i,k},{\[Rho]}] Christoffel2[g,x][{j,l},{\[Rho]}];
      (*,{\[Beta]}];*)



(*..........................................*)
(*......Ricci Formula.......................*)
(*................................by msh....*)
(*Attributes[RiccFormula]={HoldFirst};*)
(*Attributes[RiFor]={HoldFirst};*)
Format[RiccFormula[{k_,l_},T_,g_,x_]]:=
 Style[SequenceForm["[", Tensor["\[Del]",k,Void], " ,", Tensor["\[Del]",l,Void], "]", T ],Black,Italic,18,FontFamily->mshfont];(*by msh*)
 ExpandRiccFor[w_]:=  w /.RiccFormula -> RiccFor;(*by msh*)

RiccFor[{k_,l_},Tensor[T_,sub_,sup_],g_,x_] := 
 Module[{ind},
	ind = Transpose[{ sub , sup , Range[Length[sub]] }];
	Plus@@Map[ RCF[{k,l},Tensor[T,sub,sup] ,#,g,x]& ,ind]
	];
	
RCF[{k_,l_}, Tensor[T_,subT_,supT_],{sub_,Void,n_},g_,x_]:=
	-RieChrisCur[{sub},{k,l},{\[Sigma]},g,x]*Tensor[T,ReplacePart[subT,\[Sigma],n],supT]; 
RCF[{k_,l_}, Tensor[T_,subT_,supT_],{Void,sup_,n_},g_,x_]:=
	RieChrisCur[{\[Sigma]},{k,l},{sup},g,x]*Tensor[T,subT,ReplacePart[supT,\[Sigma],n]];	


(*..........................................*)
(*......Matrix function....................*)
(*................................by msh....*)
Format[ MatTen[A_,indn_] ]:=
Style[Subscript[Row[{"[",A,"]"}],Row[indn,","]],Black,Italic,16]
ExpandMat[ MatTen[A_,indn_], indt_] :=  
  Module[{ni,dum},		
  dum =  Sequence@{#[[1]],1,#[[2]]}&/@ Evaluate[{indt,indn}\[Transpose]] ;  
  (*"/@":Map "@@":Apply \:5f0f expr \:306e\:982d\:90e8\:3092 f \:3067\:7f6e\:63db\:3059\:308b 
f/@{a,b,c,d}--->{f[a],f[b],f[c],f[d]}
f@@{a,b,c,d}--->f[a,b,c,d]
dum={{n,1,3},{m,1,3},{k,1,3}}*)
  ss=Table[ A, Evaluate[Sequence@@dum] ]
  (*MatrixForm[s1]*)
]/;VectorQ[indt]&&VectorQ[indn];
ExpandMat[  w_Plus , indt_ ]:= ExpandMat[  # , indt ]& /@ w
ExpandMat[  w_*u_ ,  indt_ ] := ExpandMat[  w,  indt ]ExpandMat[  u ,  indt ]


(*..........................................*)
(*..............Protect.....................*)
(*..........................................*)
End[]
EndPackage[]
