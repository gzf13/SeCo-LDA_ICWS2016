function WriteTopicKeywords(T2_W, WO, topicNum, outKeywordNum, FILENAME)
%%function WriteTopicKeywords
%%
%根据主题-单词频数矩阵T2_W，依次在文件FILENAME中输出每个主题的前outKeywordNum个重构的关键词
%WO为字典

%%
fid = fopen( FILENAME , 'W' );

for i=1:topicNum
   fprintf( fid , '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n'); 
   fprintf( fid , 'Topic_%d Keywords:\n',i);
   [~, indexOrder] = sort(-T2_W(i,:));
   for j=1:outKeywordNum
       %依次输出当前主题的每一个关键词
       fprintf( fid , '%10s' , char(WO(indexOrder(j))));
       if (mod(j,10)==0)
           fprintf( fid , '\n');
       end
   end
   fprintf( fid , '\n');
end

fclose( fid ); 
end
