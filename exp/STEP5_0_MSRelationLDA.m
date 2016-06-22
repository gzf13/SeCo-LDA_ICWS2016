%%
%experiment3
%edit by gzf 1220
%STEP5_0 MSRelation LDA (直接用m-s调用关系做LDA)

%%
% Content-based LDA
% 参数准备

load('../Record.mat');

%%
% 除去引用service次数为0和1的mashup
msum = find(sum(useRecord)>1);
msMNum_new = length(msum);
msRelation = useRecord(:,msum);

% 再除去没有在mashup中出现过的service
ssum = find(sum(msRelation,2)>0);
msSNum_new = length(ssum);
msRelation = msRelation(ssum,:);

%%
% 输出文件以供LDA函数使用
fid = fopen('step5_0_wordcount_msRelation.txt','wt');
for i=1:1:msMNum_new
    %Mashup i
    temIndex = find(msRelation(:,i));
    temL = length(temIndex);
    for j=1:1:temL
        %输出到文件
        fprintf(fid,'%d %d %d\n',i,temIndex(j),1);
    end
    i
end
fclose(fid);

%%

[WS, DS] = importworddoccounts('step5_0_wordcount_msRelation.txt');    %利用importworddoccounts读入所需参数

%使用清洗过的数据

%WO为service列表
%选取存在引用关系的service作为词典

WO = service(ssum,1);
WO_MSR = WO';

% Set the number of topics
T=35; 

% Set the hyperparameters
BETA=0.01;
ALPHA=50/T;

% The number of iterations
N = 500; 

% The random seed
SEED = 3;

% What output to show (0=no output; 1=iterations; 2=all output)
OUTPUT = 2;

% This function might need a few minutes to finish
tic
%[WP,DP,Z]
[ W_T_MSR,S_T_MSR,Z_MSR ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
toc

%%
%预备工作二：整理D2相关数据
%DIC: WO_R
%[WP,DP,Z]
%[ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
%D2:service-topic-service

%每个主题的概率
W_MSR = size( W_T_MSR , 1 );                    %单词的数量，D2中为service的数量
T_MSR = size( W_T_MSR , 2 );                    %主题的数量，50
sumW_T_MSR = sum( W_T_MSR , 1 ) + BETA*W_MSR;
probtopic_MSR = sumW_T_MSR / sum( sumW_T_MSR );   %D2，每个主题的概率分布情况

%每个主题下单词的概率
K_R = 100;                                     %显示前50个单词(此处的“单词”为service)
Sorted_P3_T_S = zeros( K_R , T_MSR );           %D2，每个主题中service的分布情况，D2_topic_service值的大小(已排序)
Index_P3_T_S = zeros( K_R , T_MSR );            %D2，每个主题中service的下标分布，D2_topic_service具体service的下标(已排序)

for t=1:T_MSR
   [ temp1 , temp2 ] = sort( -W_T_MSR( : , t ) );
   Sorted_P3_T_S( : , t )  = ( full( -temp1( 1:K_R )) + BETA ) ./ ( repmat( sumW_T_MSR( t ) , K_R , 1 ));
   Index_P3_T_S( : , t )   = temp2( 1:K_R );
end

sNum_R = size(W_T_MSR,1);
P3_T_S = zeros( sNum_R , T_MSR );             %D2，每个主题中service的分布情况，D3_topic_service值的大小（未排序）
for t=1:T_MSR
    temp1 = W_T_MSR(:,t);
    P3_T_S( : , t )  = ( temp1 + BETA ) / sumW_T_MSR( t );
end

P3_M_T = zeros(size(S_T_MSR,1), T_MSR);                %D2，每个文档(service)的主题分布情况，D3_mashup_topic值得大小
for i=1:size(S_T_MSR,1)
    for t=1:T_MSR
        P3_M_T(i,t) = (S_T_MSR(i,t) + ALPHA)/(sum(S_T_MSR(i,:)) + ALPHA*T_MSR);
    end
end
%%
% 保存结果到文件
save STEP5_0_DATA Sorted_P3_T_S Index_P3_T_S P3_M_T P3_T_S WO_MSR W_T_MSR S_T_MSR Z_MSR K_R probtopic_MSR msum ssum msMNum_new msSNum_new msRelation T_MSR
%%

