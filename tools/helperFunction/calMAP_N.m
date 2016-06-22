function [result]=calMAP_N(x_ori, x_recom, n)
%%
% ¼ÆËãMAP@NµÄº¯Êý
% value = calMAP_N(x_ori, x_recom,n)
result = 0;
sum = 0;
for i=1:n
    temX = x_recom(i);
    oriIndex = find(x_ori==temX);
    if (oriIndex<=n)
        sum = sum + 1;
        result = result + sum/i; 
    end
end
result = result/n;
end