clear; close all
%% ��������� �������������
numOFDM = 100; % ����� ������������ OFDM-��������
M = 1024; % ������������� ���������
minVal=0; % ����������� �������� ������������ ��������
maxVal=M-1; % ������������ �������� ������������ ��������
Nfft = 256; % ����������� FFT
NSamples = 245; % ����� ������������ ����������
Tofdm = 12.8*1e-6; % ������������ ������ OFDM-�������
Tcp = 1.6*1e-6; % ������������ ������������ ��������
fs = 1/(Tofdm/Nfft); % ������� �������������
Ncp = Tcp*fs; % ����� ��������, ��������� ��� ������������ ��������

%% ������������ ������ ���������� �����
clear Tcp; clear Tofdm
outMod = []; % ������������� ������� ������������������ � ���������� QAM-M
dataForChan = []; % ������������� ������� ������������������ � ���������� OFDM
for i=1:numOFDM
    % ��������� �������� �������������� ������������������:
    data = randi([minVal maxVal], NSamples, 1); 
    % OFDM-��������� ������������������:
    dataMod = qammod(data, M);
    outMod = [outMod dataMod(:).'];
    dataMod2 = zeros(Nfft,1);
    dataMod2(1:round(NSamples/2),:) = dataMod(round(NSamples/2):end,:);
    dataMod2(end-round(NSamples/2-2):end,:) = dataMod(1:round(NSamples/2-1),:);
    dataModIFFT = ifft(dataMod2,Nfft);
    % ����������� �������:
    dataModCP = vertcat(dataModIFFT(end-Ncp+1:end,:), dataModIFFT);
    dataForChan = [dataForChan dataModCP(:).'];
end;

%% ������������ ������ ������
% ������ ���������� ������� �� ������������:
numSTSVec = [1 2]; % ���� ����� �� ������ �� � ��� ������ ������ ��
numTx = sum(numSTSVec); % ����� ���������� �������
% ������ ������ ������������ ���������� �� �� �� ������ ��:
distance = [45 60 80 100 200 300 400 500 1000];
% ����������� ������ ���������� �� ��������� � ������ ����� �� � ��
% (1.5 � � 32 � ��������������):
xsAP = round(sqrt(distance.^2-(32-1.5).^2));
scenarios = [3 4 10 11 12 14]; % ��������, ������������ � ������
for scenario=scenarios
    switch scenario
        case 3
            title_end = 'B1';
        case 4
            title_end = 'B2';
        case 10
            title_end = 'C1';
        case 11
            title_end = 'C2';
        case 12
            title_end = 'C3';
        case 14
            title_end = 'D1';
    end;
    for xAP=xsAP
        % ��������� �������� �2 � �3 ��� ������� �� ���������� 45 � ��
        % ������� �������������� ������ ��������� ��� ������ WINNER II ��
        % ���� ����������:
        if ((scenario==11)||(scenario==12))&&(xAP==xsAP(1))
            continue
        end;
        % �������� UCA-3 (����������� �������� ������ 
        % ����������� �� 3 ��������� �� �������� 10 ��):
        AP  = winner2.AntennaArray('UCA', numTx, 0.1);
        % �������� ULA-1 (���� ����������) ��� �� 1:
        STA1 = winner2.AntennaArray('ULA', numSTSVec(1));
        % �������� ULA-2 (����������� �������� ������ 
        % ����������� �� ���� ��������� �� ������ 5 ��) ��� �� 2:
        STA2 = winner2.AntennaArray('ULA', numSTSVec(2), 0.05);
        
        % ������� Config Layout:
        layoutArray = [AP,STA1,STA2];
        % ������ ������ ��� ����������� �������� �� � layoutArray:
        STAIdx = [2 3];
        % ������ ������ �������� � �����, ������ ������� ��������
        % ������������ ���� ��:
        APIdx = {1};
        K = 2; % ����� ����� ����� ��� ����� �� (links)
        rndSeed = 12;
        cfgLayout = winner2.layoutparset(STAIdx,APIdx, ...
            K,layoutArray,[],rndSeed);
        % ������������� ��������� �� � ��
        cfgLayout.Stations(1).Pos(1:2) = [0, 0]; % ��������� ��
        % ���������� �� �� 1 �� �� �������� �� 45 � �� 1 ��
        % (�� 33 � �� 1 �� �� ���������):
        cfgLayout.Stations(2).Pos(1:2) = [xAP, 0]; % ��������� �� 1
        % ���������� �� �� 2 �� �� �������� 58 � (50 � �� ���������):
        cfgLayout.Stations(3).Pos(1:2) = [0, 50]; % ��������� �� 2
        % ������ �������� ����� ��� ������ ��:
        cfgLayout.ScenarioVector = [scenario scenario];
        
        % ������� Config Model
        cfgModel = winner2.wimparset;
        % ������ ���������� ��������� �������:
        cfgModel.NumTimeSamples = numOFDM*(Nfft+Ncp);
        % ������ ���������������� �������� � �������� ��� ����������
        % ���������:
        cfgModel.FixedPdpUsed = 'yes';
        % ������ ���������������� ���� ������ (AoD) � �������� (AoA) ���
        % ���������� ���������:
        cfgModel.FixedAnglesUsed = 'yes';
        % �������� ������ �� ���� ����� ������ ��������� �� ��� �����������
        % ��� ������ ����� ����� (��� ������ ��):
        cfgModel.IntraClusterDsUsed = 'yes';
        % ���������� ������� � ������� ������������:
        cfgModel.PolarisedArrays = 'yes';
        % ������� ������� ��������������� �� ���������������� ������������
        % LOS:
        cfgModel.UseManualPropCondition = 'no';
        % ��������� ������� ������� (��):
        cfgModel.CenterFrequency = 2.4e9;
        % ��������� ���������� ��������� ������� �� �������� ����� �����:
        cfgModel.SampleDensity = fs;
        % ������� ��������� � ������:
        cfgModel.ShadowingModelUsed = 'yes';
        % ������� ������ �� ���� ��������������� � ������:
        cfgModel.PathLossModelUsed = 'yes';
        cfgModel.RandomSeed = rndSeed;
        
        % ������� ��������� ������ ������:
        winChannel = comm.WINNER2Channel(cfgModel,cfgLayout);
        chanInfo = info(winChannel)
        % �������� ������ ������� ������:
        txSig = {[dataForChan.' dataForChan.' dataForChan.'],...
            [dataForChan.' dataForChan.' dataForChan.']}';
        % ���������� ������ ����� ����� WINNER II:
        [rxSig, PG] = winChannel(txSig);
        % rxSig - ���������� �������������������; PD - Path Gains
        
        %% ������������ ������ �������� ����� (������ ������������ �� 1)
        % ��������� ������ � ������� ������������: (numOFDM*(Ncp+Nfft))x(numTx)
        dataRx = [cell2mat(rxSig(2,1)) cell2mat(rxSig(1,1))];
        % ��������� �������� ������ � ������ (������ ������������ �� 1)
        Rx1 = dataRx(:,1).';
        % ��������� ������ � ������� ������������: (Ncp+Nfft)�(numOFDM)
        Rx2 = reshape(Rx1.',(Nfft+Ncp), numOFDM);
        % ������������� ������� ���������� �������������� ������:
        outRx = [];
        % �����������:
        for i=1:numOFDM
            Rx3 = fft(Rx2(33:end,i),Nfft);
            Rx4 = vertcat(Rx3(end-round(NSamples/2-2):end),...
                Rx3(1:round(NSamples/2)));
            outRx = [outRx Rx4(:).'];
        end;
        
        %% ������������ ������
        % ���������� ������������������� ��� ���������� ���:
        save(strcat('SAparams_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.mat'),'dataForChan','rxSig')
        save(strcat('PG_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.mat'),'PG')
        if (scenario==3)
            % ������������ ������������ �� � ��
            figure
            APPos  = cfgLayout.Stations(1).Pos;
            STA1Pos = cfgLayout.Stations(2).Pos;
            STA2Pos = cfgLayout.Stations(3).Pos;
            plot3(APPos(1),  APPos(2),  APPos(3),  'bo', ...
                STA1Pos(1), STA1Pos(2), STA1Pos(3), 'rs', ...
                STA2Pos(1), STA2Pos(2), STA2Pos(3), 'rd');
            grid on; xlim([0 1000]); ylim([0 500]); zlim([0 40]);
            xlabel('x, [�]');
            ylabel('y, [�]');
            zlabel('z, [�]');
            title(strcat('���������� ����� ��1 � ��: ',num2str(sqrt(xAP.^2+(32-1.5).^2)),' �'))
            legend(strcat('��(',num2str(APPos(1)),',',num2str(APPos(2)),',',...
                num2str(APPos(3)),')'),...
                strcat('��1(',num2str(STA1Pos(1)),',',num2str(STA1Pos(2)),',',...
                num2str(STA1Pos(3)),')'), ...
                strcat('��2(',num2str(STA2Pos(1)),',',num2str(STA2Pos(2)),',',...
                num2str(STA2Pos(3)),')'), ...
                'Location', 'northeast');
            saveas(gcf,strcat('Position_x',int2str(xAP),'.fig'))
        end;
        %% ������������ �������� ����� (������ ��� �� 1)
        % ���������� ����������� ���������
        figure
        f2 = scatter(real(outRx), imag(outRx), 5,'r', 'filled');
        f2.MarkerFaceAlpha = 0.1;
        grid on
        xlabel('In-Phase')
        ylabel('Qudrature')
        title(strcat('���������� ��������� ��������� ������� �� �� 1 ��� ��������: ', title_end))
        legend(strcat(int2str(M),'QAM. ���������� ����� ��1 � ��: ',num2str(sqrt(xAP.^2+(32-1.5).^2)),' �'))
        saveas(gcf,strcat('Constellation_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.png'))
        close all;
    end;
end;
% ������������ ������������ ��������� ������� ��:
figure
plot(cellfun(@(x) x(1),{AP.Element(:).Pos}),cellfun(@(x) x(2),{AP.Element(:).Pos}),...
    'o', 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
grid on
title('��������� ��������� ������� ��');
xlabel('x, [�]')
ylabel('y, [�]')
xlim([-0.08, 0.12])
saveas(gcf,strcat('AP_Elements.png'))
%% ������������ ���������� �����
% ���������� ����������� ���������
figure
fConst = scatter(real(outMod), imag(outMod), 30, 'k', 'filled');
grid on
xlabel('In-Phase')
ylabel('Qudrature')
title('���������� ��������� ��������� �������')
legend(strcat(int2str(M),'-QAM'))
clear f1;
hold on
ref = constellation(comm.RectangularQAMModulator(M));
f = scatter(real(ref), imag(ref), 8,'r', '+');
legend(strcat(int2str(M),'-QAM'), 'Reference Constellation')
xlim([-45, 45])
ylim([-45, 45])
saveas(gcf,strcat('Constellation_Tx.png'))
close all
