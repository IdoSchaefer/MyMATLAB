[result, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 7, 1e-5)
[result, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-5)
rresult = Ftest1((0:10:200).', 0:0.01:0.05, 9, 21)
result_1 = Ftest1((0:10:200).', 0:0.01:0.05, 9, 21)
(result-result1)./result1
abs(result-result_1)./abs(result_1)
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13)
abs(result2-result_1)./abs(result_1)
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
abs(result2-result_1)./abs(result_1)
result_1 = Ftest1((0:10:200).', 0:0.01:0.05, 9, 21)
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
is_big
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
abs(result_minus1(is_big)
abs(result_minus1(is_big))
abs(result_minus1)
eps./abs(result_minus1)
eps./abs(result_minus1(is_big))>tol
(eps./abs(result_minus1(is_big)))>tol
tol
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
eps./abs(result_minus1(is_big))>tol
eps./abs(result_minus1)>tol
abs(result_minus1(is_big))
eps./abs(result_minus1(is_big))
eps./abs(result_minus1(is_big))>tol
is_big(ans)
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
epsM(is_big)./abs(result_minus1(is_big)
epsM(is_big)./abs(result_minus1(is_big))
epsM(is_big)./abs(result_minus1(is_big))>tol
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
(result2-result1)./result1
result_1 = Ftest1((0:10:200).', 0:0.01:0.05, 9, 21)
(result2-result_1)./result_1
[result, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-5)
[result, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13)
[result, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-5)
[result2, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, eps)
(result-result2)./result2
abs((result-result2)./result2)
[result3, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13)
abs((result3-result)./result)
[result3, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13)
abs(result_minus1(is_big)
abs(result_minus1(is_big))
tic;[result3, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13) toc
tic;[result3, polydeg] = Ftest3((0:10:200).', 0:0.01:0.05, 9, 1e-13); toc
tic,result_1 = Ftest1((0:10:200).', 0:0.01:0.05, 9, 21); toc
%-- 07/07/2015 14:57 --%
[result3, polydeg] = Ftest3((0:1:20).', 0:0.01:0.05, 9, 1e-13)
(0:10:200).'*(0:0.01:0.05)
[yt, all_polydeg] = Ftaylor((0:0.01:10).', 1, 9);
yd = Fdirect((0:0.01:10).', 1, 9);
relE = abs(yt-yd)./abs(yt)
[(0:0.01:20).', relE(1:21), all_polydeg(1:21)]
[(0:0.01:2).', relE(1:21), all_polydeg(1:21)]
[(0:0.01:2).', relE(1:201), all_polydeg(1:201)]
[result3, polydeg] = Ftest3((0:0.01:2).', 1, 9, 1e-5)
abs(result_minus1(is_big))
result(1)
result(2)
polyi
factorial(8)./minus_ixt
abs(result_minus1(is_big))
factorial(9)./minus_ixt
[result3, polydeg] = Ftest3((0:0.01:2).', 1, 9, 1e-13)
[yt, all_polydeg] = Ftaylor((0:0.01:10).', 1, 9);
term
[yt, all_polydeg] = Ftaylor((0:0.01:0.2).', 1, 9);
term
[yt, all_polydeg] = Ftaylor((0:0.01:0.2).', 1, 9);
term*minus_izt^9/factorial(9)
term.*minus_izt.^9/factorial(9)
abs(term).*abs(minus_izt).^9/factorial(9)
[yt, all_polydeg] = Ftaylor((0:0.01:2).', 1, 9);
term.*minus_izt.^9/factorial(9)
ans/eps
eps./(term.*minus_izt.^9/factorial(9))
[yt, all_polydeg] = Ftaylor((0:0.01:2).', 1, 9);
yd = Fdirect((0:0.01:2).', 1, 9);
relE = abs(yt-yd)./abs(yt)
term = -1i*(0:0.01:2).'/10
Eest = eps./(term.*minus_izt.^9/factorial(9));
Eest = eps./abs(term.*(-1i*(0:0.01:2).').^9/factorial(9));
[Erel, Eest]
[relE, Eest]
yt8 = Ftaylor((0:0.01:2).', 1, 8);
yt8-1
relE9=eps./abs(yt8 - 1)
relE9*factorial(9)./((0:0.01:2).').^9
Eest2 = relE9*factorial(9)./((0:0.01:2).').^9;
[relE, Eest, Eest2]
Eest./Eest2
Eest2 = relE9*factorial(9)./((0:0.01:2).').^9;
E = abs(yt-yd)
relE
abs(yt)
save Ferror