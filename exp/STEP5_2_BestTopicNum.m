%%
% 12.23
% 最优主题数目的选择

%%
% service-service


load('../Record.mat');
load('STEP2DATA.mat');

%%
% 待测试的主题数量
T_N = 5:5:70;
LValue = zeros(1,length(T_N));
minD_JS = zeros(1,length(T_N));     % 平均最小主题相似度
perplexity = zeros(1,length(T_N));  % 困惑度


for t_i = 1:length(T_N)
   disp('**************************************************************')
   disp(['当前主题数目:',num2str(T_N(t_i))]);
   disp('1. 开始LDA模型计算...');
   % LDA模型
   [WS, DS] = importworddoccounts('step2_wordcount_ssRelationLeft.txt');    %利用importworddoccounts读入所需参数

    WO = service(leftServiceSet,1);
    WO_R = WO';

    % Set the number of topics
    % 第t_i组实验
    T = T_N(t_i);

    % Set the hyperparameters
    BETA=0.01;
    ALPHA=1;

    % The number of iterations
    N = 400; 
    % The random seed
    SEED = 3;
    % What output to show (0=no output; 1=iterations; 2=all output)
    OUTPUT = 0;
    % This function might need a few minutes to finish
    tic
    %[WP,DP,Z]
    [ W_T_R,S_T_R,Z_R ] = GibbsSamplerLDA( WS , DS , T , N , ALPHA , BETA , SEED , OUTPUT );
    toc
    
    % 计算P_S_T && P_T_S
    disp('2. 开始相关概率分布计算...');
    
    %每个主题的概率
    W_R = size( W_T_R , 1 );                    %单词的数量，D2中为service的数量
    T_R = size( W_T_R , 2 );                    %主题的数量
    %T_R
    sumW_T_R = sum( W_T_R , 1 ) + BETA*W_R;
    probtopic_R = sumW_T_R / sum( sumW_T_R );   %D2，每个主题的概率分布情况

    %P_T_S
    sNum_R = size(S_T_R,1);
    P2_T_S = zeros( sNum_R , T_R );             %D2，每个主题中service的分布情况，D2_topic_service值的大小（未排序）
    for t=1:T_R
        temp1 = W_T_R(:,t);
        P2_T_S( : , t )  = ( temp1 + BETA ) ./ ( repmat( sumW_T_R( t ) , sNum_R , 1 ));
    end
    %P_S_T
    P2_S_T = zeros(sNum_R, T_R);                %D2，每个文档(service)的主题分布情况，D2_service_topic值得大小
    for i=1:sNum_R
        P2_S_T(i,:) = (S_T_R(i,:) + ALPHA)/(sum(S_T_R(i,:)) + ALPHA*T_R);
    end

    % 计算似然函数值
    disp('3. 开始似然函数值计算...')
    temLValue = 0;
    for i=1:leftServiceNum
       % each doucment(service here)
       word_index = find(ssRelationLeft(i,:));  %找出di包含的所有单词
       for w_i = 1:length(word_index)
           n_wi = ssRelationLeft(i,word_index(w_i));
           temLValue = temLValue + n_wi * log(sum(P2_T_S(word_index(w_i),:).*P2_S_T(i,:)));     %sigma[z]P(w|z)P(z|d)
       end
    end
    LValue(t_i) = temLValue;
    disp(['似然函数值为:',num2str(temLValue)]);
    
    % 计算平均最小主题相似度
    disp('4. 开始平均最小主题相似度计算...');
    D_JS = zeros(T,T);
    temAverJS = 0;
    for a = 1:T
        for b = 1: T
            if (a==b) 
                D_JS(a,b) = 9999;
            else
               D_JS(a,b) =  calJS(P2_T_S(:,a),P2_T_S(:,b));
            end
        end
    end
    for i = 1:T
        temAverJS = temAverJS + min(D_JS(i,:));
    end
    temAverJS = temAverJS/T;
    minD_JS(t_i) = temAverJS;
    
    % 计算困惑度
    disp('5. 开始困惑度计算...');
    sigma_log_pwd = 0;
    for i=1:leftServiceNum
        % each document(service here)
        word_index = find(ssRelationLeft(i,:));  %找出di包含的所有单词
        for w_i = 1:length(word_index)
           sigma_log_pwd = sigma_log_pwd + log(sum(P2_T_S(word_index(w_i),:).*P2_S_T(i,:)));     %sigma[z]P(w|z)P(z|d)
       end
    end
    sigma_Nd = sum(sum(ssRelationLeft));
    perplexity(t_i) = exp(-sigma_log_pwd/sigma_Nd);
end

save STEP5_0_scitation minD_JS perplexity;

% % 绘制曲线
% plot(T_N,LValue,'-b*');
% [~,index] = max(LValue);
% disp('******************************************************')
% disp(['最佳主题数为：',num2str(T_N(index))]);