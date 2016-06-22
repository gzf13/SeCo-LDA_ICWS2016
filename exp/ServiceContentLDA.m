%%
%Content-based LDA for Service


%%
%读入数据及参数初始化

load('../Description.mat');


%%
%读取每一个单词的信息
fid = fopen('wordcount_service.txt','wt');
for i=1:1:sNum
    %Service_i
    temIndex = find(SW(i,:));
    temL = length(temIndex);
    for j=1:1:temL
        %输出到文件 
        fprintf(fid,'%d %d %d\n',i,temIndex(j),SW(i,temIndex(j)));
    end
end
fclose(fid);