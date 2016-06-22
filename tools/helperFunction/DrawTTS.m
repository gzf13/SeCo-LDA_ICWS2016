function DrawTTS(TTS,average_TTS,i,totalT,sigma)
%%
% 绘制第i个主题的时间特性曲线

%%
temT = TTS(i,:);
[~, x0] = find(temT);
y0 = temT(x0);

%%
% 绘制原始数据点

% 1. 直接绘出
figure();
plot(x0,y0,'rs');
hold on;
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
title(i);
rTTS = rTTS/max(rTTS)*max(temT);
plot(rTTS,'b');                                             %绘制插值平滑后的时间特性曲线
plot([average_TTS(i) average_TTS(i)],[0,max(temT)*2],'g');  %绘制均值
axis([0 3800 0 max(temT)*1.2]);                             %设置坐标轴范围

end
