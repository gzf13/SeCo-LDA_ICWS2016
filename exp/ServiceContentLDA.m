%%
%Content-based LDA for Service


%%
%�������ݼ�������ʼ��

load('../Description.mat');


%%
%��ȡÿһ�����ʵ���Ϣ
fid = fopen('wordcount_service.txt','wt');
for i=1:1:sNum
    %Service_i
    temIndex = find(SW(i,:));
    temL = length(temIndex);
    for j=1:1:temL
        %������ļ� 
        fprintf(fid,'%d %d %d\n',i,temIndex(j),SW(i,temIndex(j)));
    end
end
fclose(fid);