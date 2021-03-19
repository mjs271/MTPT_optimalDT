function v2 = calcV2(v1, Lv, tDef)

invV = @(v1, v2, Lv, x) 1 ./ (v1 + ((v2 - v1) ./ Lv) * x);

TT = @(v2) integral(@(x) invV(v1, v2, Lv, x), 0, Lv) - tDef;

v2 = fsolve(TT, 0.0);

end