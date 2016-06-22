%%
%experiment1
%edit by gzf 1209
%STEP2: Draw service-service relation
%��mashup-service���ù�ϵ�ó�service-service֮������ù�ϵ

%%
%�������
load('../Record.mat');

%%
[sNum, mNum] = size(useRecord);
ssRelation = zeros(sNum, sNum);

%%
for i=1:1:sNum
    %��ÿ��service����ͳ�ƣ�ͳ�����������ͬһ��mashup�е�service������������������
    relatedMashup = find(useRecord(i,:));           %�ҳ���ص�mashup
    relatedMashupNum = length(relatedMashup);
    if (relatedMashupNum~=0)
       for j=1:1:relatedMashupNum                   %����ÿ��mashup���ҳ����а���������service
           relatedService = find(useRecord(:,relatedMashup(j)));
           relatedServiceNum = length(relatedService);
           if (relatedServiceNum>1)                 %>1��˵������������serviceҲ�����mashup����
               for k=1:1:relatedServiceNum
                  if (relatedService(k)~=i)
                     ssRelation(i,relatedService(k)) = ssRelation(i,relatedService(k)) + 1;
                  end
               end
           end
       end
    end
    i
end

%%

flag_ssRelation = sum(ssRelation,2);    %��ssRelation�������
leftServiceSet = find(flag_ssRelation>0);   %�������ù�ϵ��service���
leftServiceNum = length(leftServiceSet);    %�������ù�ϵ��service����

ssRelationLeft = ssRelation(leftServiceSet,leftServiceSet);

%%
%��������ϵ������ĵ���,���ssRelationLeft�е�����
fid = fopen('step2_wordcount_ssRelationLeft.txt','wt');
for i=1:1:leftServiceNum
    %Service_i
    temIndex = find(ssRelationLeft(i,:));
    temL = length(temIndex);
    for j=1:1:temL
        %������ļ�
        fprintf(fid,'%d %d %d\n',i,temIndex(j),ssRelationLeft(i,temIndex(j)));
    end
    i
end
fclose(fid);

%%
% �������������ļ�
save STEP2DATA ssRelation ssRelationLeft leftServiceSet leftServiceNum

