legendText = [];
distance = [45 60 80 100 200 300 400 500 1000];
distance = [100];
xsAP = round(sqrt(distance.^2-(32-1.5).^2));
scenarios = [3 4 10 11 12 14];
scenarios = [3];
%% Строим графики
for scenario=scenarios
    for xAP=xsAP
        if ((scenario==11)||(scenario==12))&&(xAP==xsAP(1))
            continue
        end;
        load(strcat('PG_Scenario',int2str(scenario),'_QAM1024_x',...
            int2str(xAP),'.mat'))
        Npath = size(PG{2,1},3);
        color = getRainbow(Npath);

        % Покажем поближе несколько (Пока только первый луч):
        legendText = {};
        for i=1:Npath
            legendText{end+1} = strcat('МС1, луч №',int2str(i));
             if (i==1)%||(mod(i,10)==0)
                 PathGain = PG{2,1}(1,1,i,:); % Интересующая МС 1
                 f1 = figure;
                 scatter(real(PathGain), imag(PathGain), 5,color(i,:), 'filled');
                 grid on
                 xlabel('In-Phase')
                 ylabel('Qudrature')
                 legend(strcat('МС1, луч №',int2str(i)))
                 title(strcat('Коэффициенты пути для МС1, Сценарий: ',int2str(scenario), ', x = ',int2str(xAP),' м'))
                 saveas(gcf,strcat('PathGains_path',int2str(i),'_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
                 
             end;
             close all
        end;
        % Построим общий
        figure
        for i=1:Npath
            PathGain = PG{2,1}(1,1,i,:);
            scatter(real(PathGain), imag(PathGain), 40,color(i,:), 'filled');
            hold on
        end;
        grid on
        title(strcat('Коэффициенты пути для МС1, Сценарий: ',int2str(scenario), ', x = ',int2str(xAP),' м'))
        xlabel('In-Phase')
        ylabel('Qudrature')
        legend(legendText)
        saveas(gcf,strcat('AllPathGains_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllPathGains_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
        
        % Амплитуда и фаза 1 луча:
        figure
        absPG = abs(squeeze(PG{2,1}(1,1,1,:)));
        anglePG = angle(squeeze(PG{2,1}(1,1,1,:)));
        vecLength = numOFDM * (Ncp+Nfft);
        [hAx,hLine1,hLine2] = plotyy(1:vecLength,absPG,1:vecLength,anglePG);
        title(strcat('Амплитуда и фаза для Path Gains МС 1, луча №1, Сценарий: ',...
            int2str(scenario), ', x = ',int2str(xAP),' м'))
        xlabel('Отсчет')
        ylabel(hAx(1),'abs') % left y-axis
        ylabel(hAx(2),'angle') % right y-axis
        hLine1.LineWidth = 2;
        hLine2.LineWidth = 2;
        grid on
        saveas(gcf,strcat('AbsAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AbsAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
        
        % Амплитуда и фаза каждого луча:
        figure;
        absAllPG = abs(squeeze(PG{2,1}(1,1,:,1)));
        angleAllPG = angle(squeeze(PG{2,1}(1,1,:,1)));
        f = stem(absAllPG, 'filled');
        title(strcat('Амплитудный коэффициент пути МС 1, Сценарий: ',...
            int2str(scenario), ', x = ',int2str(xAP),' м'))
        f.MarkerSize = 6;
        f.LineWidth = 1.5;
        f.Color = [0.8 0.1 0.05];
        xlabel('Порядковый номер луча')
        ylabel('Abs')
        grid on
        saveas(gcf,strcat('AllAbs_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllAbs_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
    
        figure
        f = stem(angleAllPG, 'filled');
        title(strcat('Фазовый коэффициент пути МС 1, Сценарий: ',...
            int2str(scenario), ', x = ',int2str(xAP),' м'))
        f.MarkerSize = 6;
        f.LineWidth = 1.5;
        f.Color = [0 0.35 0.5];
        xlabel('Порядковый номер луча')
        ylabel('Angle')
        grid on
        saveas(gcf,strcat('AllAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
    end;
end;