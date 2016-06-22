function [result] = calJS(P,Q)
%%
% º∆À„D_JS(P||Q)
R = (P+Q)/2;
result = calKL(P,R)/2 + calKL(Q,R)/2;
end
