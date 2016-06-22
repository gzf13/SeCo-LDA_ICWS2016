%%
%experiment1
%edit by gzf 1209
%STEP1:Content-based LDA for Service


%%
%读入数据及参数初始化

load('../Description.mat');

%%
%读取每一个单词的信息并保存到文档中
fid = fopen('step1_wordcount_service.txt','wt');
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

%%
% Content-based LDA
% 参数准备
[WS1, DS1] = importworddoccounts('step1_wordcount_service.txt');    %利用importworddoccounts读入所需参数
WO_S = WD';                                                   %转成1*W的单词字典矩阵

% Set the number of topics
T=35; 

% Set the hyperparameters
BETA=0.01;
ALPHA=50/T;

% The number of iterations
N = 400; 

% The random seed
SEED = 3;

% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 2;

% This function might need a few minutes to finish
tic
%WP DP Z
%DIC: WO_S
[ W_T_S, S_T, Z_S ] = GibbsSamplerLDA( WS1 , DS1 , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
%预备工作一：整理D1相关数据
%DIC: WO_S
%[WP,DP,Z]
%[ W_T_S, S_T, Z_S ] = GibbsSamplerLDA( WS1 , DS1 , T , N , ALPHA , BETA , SEED , OUTPUT );
%D1：service-topic-word

%每个主题的概率
W_S = size( W_T_S , 1 );
T_S = size( W_T_S , 2 );
sumW_T_S = sum( W_T_S , 1 ) + BETA*W_S;
probtopic = sumW_T_S / sum( sumW_T_S ); 

%每个主题下单词的概率
%前100个
K_S = 100;                              %显示前100个单词
Sorted_P1_T_W = zeros( K_S , T_S );     %D1：topic-word概率分布的具体值大小（从大到小）
Index_P1_T_W = zeros( K_S , T_S );      %D1：topic-word概率分布的word标签
%所有单词
for t=1:T_S
   [ temp1 , temp2 ] = sort( -W_T_S( : , t ) );
   Sorted_P1_T_W( : , t )  = ( full( -temp1( 1:K_S )) + BETA ) ./ ( repmat( sumW_T_S( t ) , K_S , 1 ));
   Index_P1_T_W( : , t )   = temp2( 1:K_S );
end

P1_T_W = zeros( length(WO_S) , T_S );             %D1，每个主题中word的分布情况，D1_topic_word值的大小（未排序）,按列看
for t=1:T_S
    temp1 = W_T_S(:,t);
    P1_T_W( : , t )  = ( temp1 + BETA ) / sumW_T_S( t );
end

P1_S_T = zeros(sNum, T_S);                %D1，每个文档(service)的主题分布情况，D1_service_topic值得大小，按行看
for i=1:sNum
    P1_S_T(i,:) = (S_T(i,:) + ALPHA)/(sum(S_T(i,:)) + ALPHA*T_S);
    i
end

%%
% 保存输出结果到文件
save STEP1DATA Sorted_P1_T_W Index_P1_T_W P1_T_W P1_S_T WO_S sNum probtopic



