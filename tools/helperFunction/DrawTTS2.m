function DrawTTS2(TTS,average_TTS,i,i2,totalT,sigma)
%%
% 绘制第i个主题的时间特性曲线

%%
temT = TTS(i,:);
[~, x0] = find(temT);
y0 = temT(x0);

%%
% 绘制原始数据点

% 1. 直接绘出

% 2. 根据权重大小设置不同的makersize
% red,cubic
% for j=1:length(x0)
%     plot(x0(j),y0(j),'rs','markersize',20*y0(j)/max(y0));
%     hold on;
% end

% 3. 插值绘制
rTTS = zeros(1,totalT);
for j=1:totalT
    tem = 0;
    for k=1:totalT
        % norm N(0,sigma)
        tem = tem + normpdf(k-j,0,sigma)/normpdf(0,0,sigma)*temT(k);
    end
    rTTS(j) = tem;
    j
end
rTTS = rTTS/sum(rTTS)*120;
plot(rTTS,'b');                                             %绘制插值平滑后的时间特性曲线
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*2],'b');  %绘制均值
axis([0 3800 0 max(temT)*1.2]);                             %设置坐标轴范围

%%
% 绘制第二条曲线
temT = TTS(i2,:);
[~, x02] = find(temT);
y02 = temT(x02);
%%
% 绘制原始数据点

% 1. 直接绘出

% 2. 根据权重大小设置不同的makersize
% red,cubic
% for j=1:length(x0)
%     plot(x0(j),y0(j),'rs','markersize',20*y0(j)/max(y0));
%     hold on;
% end

% 3. 插值绘制
rTTS2 = zeros(1,totalT);
for j=1:totalT
    tem = 0;
    for k=1:totalT
        % norm N(0,sigma)
        tem = tem + normpdf(k-j,0,sigma)/normpdf(0,0,sigma)*temT(k);
    end
    rTTS2(j) = tem;
    j
end

rTTS2 = rTTS2/sum(rTTS2)*120;

%% 绘图
figure();
plot(x0,y0,'ks','MarkerFaceColor','b');
hold on;

f1=plot(rTTS,'b','LineWidth',1.6);                                             %绘制插值平滑后的时间特性曲线
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*2],'b');  %绘制均值
axis([0 3800 0 max(temT)*1.2]);                             %设置坐标轴范围

plot(x02,y02,'ks','MarkerFaceColor','r');


f2=plot(rTTS2,'r','LineWidth',1.6);                                             %绘制插值平滑后的时间特性曲线
plot([average_TTS(i2) average_TTS(i2)],[0,max(temT)*2],'r');  %绘制均值

%设置横坐标显示内容
set(gca,'XTick',120:365:3600);
set(gca,'XTickLabel',{'2006','2007','2008','2009','2010','2011','2012','2013','2014','2015'});

str1 = strcat('Topic',num2str(i));
str2 = strcat('Topic',num2str(i2));
legend([f1,f2],str1,str2);

xlabel('Time');
ylabel('P(z|t)');

end
