%%
%experiment1
%edit by gzf 1209
%STEP2: Draw service-service relation
%由mashup-service调用关系得出service-service之间的引用关系

%%
%读入参数
load('../Record.mat');

%%
[sNum, mNum] = size(useRecord);
ssRelation = zeros(sNum, sNum);

%%
for i=1:1:sNum
    %对每个service进行统计，统计与其出现在同一个mashup中的service个数（不包括它本身）
    relatedMashup = find(useRecord(i,:));           %找出相关的mashup
    relatedMashupNum = length(relatedMashup);
    if (relatedMashupNum~=0)
       for j=1:1:relatedMashupNum                   %对于每个mashup，找出其中包含的其他service
           relatedService = find(useRecord(:,relatedMashup(j)));
           relatedServiceNum = length(relatedService);
           if (relatedServiceNum>1)                 %>1才说明还有其他的service也被这个mashup调用
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

flag_ssRelation = sum(ssRelation,2);    %将ssRelation按行求和
leftServiceSet = find(flag_ssRelation>0);   %存在引用关系的service编号
leftServiceNum = length(leftServiceSet);    %存在引用关系的service数量

ssRelationLeft = ssRelation(leftServiceSet,leftServiceSet);

%%
%将关联关系输出到文档中,输出ssRelationLeft中的内容
fid = fopen('step2_wordcount_ssRelationLeft.txt','wt');
for i=1:1:leftServiceNum
    %Service_i
    temIndex = find(ssRelationLeft(i,:));
    temL = length(temIndex);
    for j=1:1:temL
        %输出到文件
        fprintf(fid,'%d %d %d\n',i,temIndex(j),ssRelationLeft(i,temIndex(j)));
    end
    i
end
fclose(fid);

%%
% 保存结果到数据文件
save STEP2DATA ssRelation ssRelationLeft leftServiceSet leftServiceNum

