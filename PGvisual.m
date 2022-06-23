legendText = [];
distance = [45 60 80 100 200 300 400 500 1000];
distance = [100];
xsAP = round(sqrt(distance.^2-(32-1.5).^2));
scenarios = [3 4 10 11 12 14];
scenarios = [3];
%% ������ �������
for scenario=scenarios
    for xAP=xsAP
        if ((scenario==11)||(scenario==12))&&(xAP==xsAP(1))
            continue
        end;
        load(strcat('PG_Scenario',int2str(scenario),'_QAM1024_x',...
            int2str(xAP),'.mat'))
        Npath = size(PG{2,1},3);
        color = getRainbow(Npath);

        % ������� ������� ��������� (���� ������ ������ ���):
        legendText = {};
        for i=1:Npath
            legendText{end+1} = strcat('��1, ��� �',int2str(i));
             if (i==1)%||(mod(i,10)==0)
                 PathGain = PG{2,1}(1,1,i,:); % ������������ �� 1
                 f1 = figure;
                 scatter(real(PathGain), imag(PathGain), 5,color(i,:), 'filled');
                 grid on
                 xlabel('In-Phase')
                 ylabel('Qudrature')
                 legend(strcat('��1, ��� �',int2str(i)))
                 title(strcat('������������ ���� ��� ��1, ��������: ',int2str(scenario), ', x = ',int2str(xAP),' �'))
                 saveas(gcf,strcat('PathGains_path',int2str(i),'_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
                 
             end;
             close all
        end;
        % �������� �����
        figure
        for i=1:Npath
            PathGain = PG{2,1}(1,1,i,:);
            scatter(real(PathGain), imag(PathGain), 40,color(i,:), 'filled');
            hold on
        end;
        grid on
        title(strcat('������������ ���� ��� ��1, ��������: ',int2str(scenario), ', x = ',int2str(xAP),' �'))
        xlabel('In-Phase')
        ylabel('Qudrature')
        legend(legendText)
        saveas(gcf,strcat('AllPathGains_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllPathGains_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
        
        % ��������� � ���� 1 ����:
        figure
        absPG = abs(squeeze(PG{2,1}(1,1,1,:)));
        anglePG = angle(squeeze(PG{2,1}(1,1,1,:)));
        vecLength = numOFDM * (Ncp+Nfft);
        [hAx,hLine1,hLine2] = plotyy(1:vecLength,absPG,1:vecLength,anglePG);
        title(strcat('��������� � ���� ��� Path Gains �� 1, ���� �1, ��������: ',...
            int2str(scenario), ', x = ',int2str(xAP),' �'))
        xlabel('������')
        ylabel(hAx(1),'abs') % left y-axis
        ylabel(hAx(2),'angle') % right y-axis
        hLine1.LineWidth = 2;
        hLine2.LineWidth = 2;
        grid on
        saveas(gcf,strcat('AbsAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AbsAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
        
        % ��������� � ���� ������� ����:
        figure;
        absAllPG = abs(squeeze(PG{2,1}(1,1,:,1)));
        angleAllPG = angle(squeeze(PG{2,1}(1,1,:,1)));
        f = stem(absAllPG, 'filled');
        title(strcat('����������� ����������� ���� �� 1, ��������: ',...
            int2str(scenario), ', x = ',int2str(xAP),' �'))
        f.MarkerSize = 6;
        f.LineWidth = 1.5;
        f.Color = [0.8 0.1 0.05];
        xlabel('���������� ����� ����')
        ylabel('Abs')
        grid on
        saveas(gcf,strcat('AllAbs_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllAbs_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
    
        figure
        f = stem(angleAllPG, 'filled');
        title(strcat('������� ����������� ���� �� 1, ��������: ',...
            int2str(scenario), ', x = ',int2str(xAP),' �'))
        f.MarkerSize = 6;
        f.LineWidth = 1.5;
        f.Color = [0 0.35 0.5];
        xlabel('���������� ����� ����')
        ylabel('Angle')
        grid on
        saveas(gcf,strcat('AllAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.png'))
        saveas(gcf,strcat('AllAngle_Scenario',int2str(scenario),'_QAM',int2str(M),'_h',int2str(xAP),'.fig'))
    end;
end;