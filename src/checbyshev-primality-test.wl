polmul[f_, g_, r_, n_] := Mod[f.NestList[RotateRight, g, r - 1], n]

matmul[a_, b_, r_, n_] :=  Mod[
 {{polmul[a[[1, 1]], b[[1, 1]], r, n] + polmul[a[[1, 2]], b[[2, 1]], r, n], 
   polmul[a[[1, 1]], b[[1, 2]], r, n] + polmul[a[[1, 2]], b[[2, 2]], r, n]}, 
  {polmul[a[[2, 1]], b[[1, 1]], r, n] + polmul[a[[2, 2]], b[[2, 1]], r, n], 
   polmul[a[[2, 1]], b[[1, 2]], r, n] + polmul[a[[2, 2]], b[[2, 2]], r, n]}}, n]

matsq[a_, r_, n_] := matmul[a, a, r, n]

matpow[a_, k_, r_, n_] := If[k == 1, a,
 If[EvenQ[k],
  matpow[matsq[a, r, n], k/2, r, n], 
  matmul[a, matpow[matsq[a, r, n], (k - 1)/2, r, n], r, n]
 ]
]

xmat[r_, n_] :=
 {{PadRight[{0, 2}, r], PadRight[{n - 1}, r]},
  {PadRight[{1}, r], ConstantArray[0, r]}}

isprime[n_] := With[{r = smallestr[n]}, 
 If[r == 0, n == 2,
  With[{xp = matpow[xmat[r, n], n - 1, r, n]},
   Mod[RotateRight[xp[[1, 1]]] + xp[[1, 2]], n]
    === PadRight[Append[ConstantArray[0, Mod[n, r]], 1], r]
  ]
 ]
]
