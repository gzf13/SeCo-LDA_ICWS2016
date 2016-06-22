function WriteTopicKeywords(T2_W, WO, topicNum, outKeywordNum, FILENAME)
%%function WriteTopicKeywords
%%
%��������-����Ƶ������T2_W���������ļ�FILENAME�����ÿ�������ǰoutKeywordNum���ع��Ĺؼ���
%WOΪ�ֵ�

%%
fid = fopen( FILENAME , 'W' );

for i=1:topicNum
   fprintf( fid , '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'); 
   fprintf( fid , 'Topic_%d Keywords:\n',i);
   [~, indexOrder] = sort(-T2_W(i,:));
   for j=1:outKeywordNum
       %���������ǰ�����ÿһ���ؼ���
       fprintf( fid , '%10s' , char(WO(indexOrder(j))));
       if (mod(j,10)==0)
           fprintf( fid , '\n');
       end
   end
   fprintf( fid , '\n');
end

fclose( fid ); 
end
