function [result] = calKL(A, B)
%%
% º∆À„D_KL(A|B)
% topic A, topic Bs
M = length(A);
result = 0;
for i=1:M
   result = result +  A(i) * log(A(i)/B(i));
end

end