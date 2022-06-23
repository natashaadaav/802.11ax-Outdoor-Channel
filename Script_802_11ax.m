clear; close all
%% Параметры моделирования
numOFDM = 100; % Число передаваемых OFDM-символов
M = 1024; % Позиционность модуляции
minVal=0; % Минимальное значение генерируемых символов
maxVal=M-1; % Максимальное значение генерируемых символов
Nfft = 256; % Размерность FFT
NSamples = 245; % Число используемых поднесущих
Tofdm = 12.8*1e-6; % Длительность одного OFDM-символа
Tcp = 1.6*1e-6; % Длительность циклического префикса
fs = 1/(Tofdm/Nfft); % Частота дискретизации
Ncp = Tcp*fs; % Число символов, отводимых для циклического префикса

%% Имитационная модель передающей части
clear Tcp; clear Tofdm
outMod = []; % Инициализация входной последовательности с модуляцией QAM-M
dataForChan = []; % Инициализация входной последовательности с модуляцией OFDM
for i=1:numOFDM
    % Генерация исходной информационной последовательности:
    data = randi([minVal maxVal], NSamples, 1); 
    % OFDM-модуляция последовательности:
    dataMod = qammod(data, M);
    outMod = [outMod dataMod(:).'];
    dataMod2 = zeros(Nfft,1);
    dataMod2(1:round(NSamples/2),:) = dataMod(round(NSamples/2):end,:);
    dataMod2(end-round(NSamples/2-2):end,:) = dataMod(1:round(NSamples/2-1),:);
    dataModIFFT = ifft(dataMod2,Nfft);
    % Циклический префикс:
    dataModCP = vertcat(dataModIFFT(end-Ncp+1:end,:), dataModIFFT);
    dataForChan = [dataForChan dataModCP(:).'];
end;

%% Имитационная модель канала
% Задаем количество потоков на пользователя:
numSTSVec = [1 2]; % Один поток на первую МС и два потока вторую МС
numTx = sum(numSTSVec); % Общее количество потоков
% Задаем вектор используемых расстояний от ТД до второй МС:
distance = [45 60 80 100 200 300 400 500 1000];
% Расчитываем вектор расстояний на плоскости с учетом высот МС и ТД
% (1.5 м и 32 м соответственно):
xsAP = round(sqrt(distance.^2-(32-1.5).^2));
scenarios = [3 4 10 11 12 14]; % Сценарии, используемые в модели
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
        % Исключаем сценарии С2 и С3 для расчета на расстоянии 45 м по
        % причине неприменимости данных сценариев для модели WINNER II на
        % этом расстоянии:
        if ((scenario==11)||(scenario==12))&&(xAP==xsAP(1))
            continue
        end;
        % Создание UCA-3 (равномерный круговой массив 
        % излучателей из 3 элементов ТД радиусом 10 см):
        AP  = winner2.AntennaArray('UCA', numTx, 0.1);
        % Создание ULA-1 (Один излучатель) для МС 1:
        STA1 = winner2.AntennaArray('ULA', numSTSVec(1));
        % Создание ULA-2 (равномерный линейный массив 
        % излучателей из двух элементов МС длиной 5 см) для МС 2:
        STA2 = winner2.AntennaArray('ULA', numSTSVec(2), 0.05);
        
        % Создаем Config Layout:
        layoutArray = [AP,STA1,STA2];
        % Задаем вектор для обозначения индексов МС в layoutArray:
        STAIdx = [2 3];
        % Задаем массив столбцов и ячеек, каждый элемент которого
        % представляет одну БС:
        APIdx = {1};
        K = 2; % Число линий связи или число МС (links)
        rndSeed = 12;
        cfgLayout = winner2.layoutparset(STAIdx,APIdx, ...
            K,layoutArray,[],rndSeed);
        % Устанавливаем положения МС и ТД
        cfgLayout.Stations(1).Pos(1:2) = [0, 0]; % Положение ТД
        % Расстояние от МС 1 до ТД примерно от 45 м до 1 км
        % (от 33 м до 1 км по плоскости):
        cfgLayout.Stations(2).Pos(1:2) = [xAP, 0]; % Положение МС 1
        % Расстояние от МС 2 до ТД примерно 58 м (50 м по плоскости):
        cfgLayout.Stations(3).Pos(1:2) = [0, 50]; % Положение МС 2
        % Задаем сценарии связи для каждой МС:
        cfgLayout.ScenarioVector = [scenario scenario];
        
        % Создаем Config Model
        cfgModel = winner2.wimparset;
        % Задаем количество временных выборок:
        cfgModel.NumTimeSamples = numOFDM*(Nfft+Ncp);
        % Задаем предопределенные задержки и мощности для конкретных
        % сценариев:
        cfgModel.FixedPdpUsed = 'yes';
        % Задаем предопределенные углы вылета (AoD) и прибытия (AoA) для
        % конкретных сценариев:
        cfgModel.FixedAnglesUsed = 'yes';
        % Разделим каждый из двух самых мощных кластеров на три подкластера
        % для каждой линии связи (для каждой МС):
        cfgModel.IntraClusterDsUsed = 'yes';
        % Используем массивы с двойной поляризацией:
        cfgModel.PolarisedArrays = 'yes';
        % Зададим условия распространения из предопределенных вероятностей
        % LOS:
        cfgModel.UseManualPropCondition = 'no';
        % Установим несущую частоту (Гц):
        cfgModel.CenterFrequency = 2.4e9;
        % Установим количество временных выборок на половину длины волны:
        cfgModel.SampleDensity = fs;
        % Включим затенение в модель:
        cfgModel.ShadowingModelUsed = 'yes';
        % Включим потери на пути распространения в модель:
        cfgModel.PathLossModelUsed = 'yes';
        cfgModel.RandomSeed = rndSeed;
        
        % Создаем системный объект канала:
        winChannel = comm.WINNER2Channel(cfgModel,cfgLayout);
        chanInfo = info(winChannel)
        % Образуем массив потоков данных:
        txSig = {[dataForChan.' dataForChan.' dataForChan.'],...
            [dataForChan.' dataForChan.' dataForChan.']}';
        % Пропускаем сигнал через канал WINNER II:
        [rxSig, PG] = winChannel(txSig);
        % rxSig - полученная постледовательность; PD - Path Gains
        
        %% Имитационная модель приемной части (Только интересующей МС 1)
        % Помещение данных в матрицу размерностью: (numOFDM*(Ncp+Nfft))x(numTx)
        dataRx = [cell2mat(rxSig(2,1)) cell2mat(rxSig(1,1))];
        % Помещение принятых данных в вектор (только интересующей МС 1)
        Rx1 = dataRx(:,1).';
        % Помещение данных в матрицу размерностью: (Ncp+Nfft)х(numOFDM)
        Rx2 = reshape(Rx1.',(Nfft+Ncp), numOFDM);
        % Инициализация вектора полученных информационных данных:
        outRx = [];
        % Демодуляция:
        for i=1:numOFDM
            Rx3 = fft(Rx2(33:end,i),Nfft);
            Rx4 = vertcat(Rx3(end-round(NSamples/2-2):end),...
                Rx3(1:round(NSamples/2)));
            outRx = [outRx Rx4(:).'];
        end;
        
        %% Визуализация канала
        % Сохранение последовательностей для построения СПМ:
        save(strcat('SAparams_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.mat'),'dataForChan','rxSig')
        save(strcat('PG_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.mat'),'PG')
        if (scenario==3)
            % Визуализация расположения ТД и МС
            figure
            APPos  = cfgLayout.Stations(1).Pos;
            STA1Pos = cfgLayout.Stations(2).Pos;
            STA2Pos = cfgLayout.Stations(3).Pos;
            plot3(APPos(1),  APPos(2),  APPos(3),  'bo', ...
                STA1Pos(1), STA1Pos(2), STA1Pos(3), 'rs', ...
                STA2Pos(1), STA2Pos(2), STA2Pos(3), 'rd');
            grid on; xlim([0 1000]); ylim([0 500]); zlim([0 40]);
            xlabel('x, [м]');
            ylabel('y, [м]');
            zlabel('z, [м]');
            title(strcat('Расстояние между МС1 и ТД: ',num2str(sqrt(xAP.^2+(32-1.5).^2)),' м'))
            legend(strcat('ТД(',num2str(APPos(1)),',',num2str(APPos(2)),',',...
                num2str(APPos(3)),')'),...
                strcat('МС1(',num2str(STA1Pos(1)),',',num2str(STA1Pos(2)),',',...
                num2str(STA1Pos(3)),')'), ...
                strcat('МС2(',num2str(STA2Pos(1)),',',num2str(STA2Pos(2)),',',...
                num2str(STA2Pos(3)),')'), ...
                'Location', 'northeast');
            saveas(gcf,strcat('Position_x',int2str(xAP),'.fig'))
        end;
        %% Визуализация приемной части (только для МС 1)
        % Построение сигнального созвездия
        figure
        f2 = scatter(real(outRx), imag(outRx), 5,'r', 'filled');
        f2.MarkerFaceAlpha = 0.1;
        grid on
        xlabel('In-Phase')
        ylabel('Qudrature')
        title(strcat('Сигнальное созвездие принятого сигнала на МС 1 для сценария: ', title_end))
        legend(strcat(int2str(M),'QAM. Расстояние между МС1 и ТД: ',num2str(sqrt(xAP.^2+(32-1.5).^2)),' м'))
        saveas(gcf,strcat('Constellation_Scenario',int2str(scenario),'_QAM',int2str(M),'_x',int2str(xAP),'.png'))
        close all;
    end;
end;
% Визуализация расположения элементов антенны ТД:
figure
plot(cellfun(@(x) x(1),{AP.Element(:).Pos}),cellfun(@(x) x(2),{AP.Element(:).Pos}),...
    'o', 'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
grid on
title('Положение элементов антенны ТД');
xlabel('x, [м]')
ylabel('y, [м]')
xlim([-0.08, 0.12])
saveas(gcf,strcat('AP_Elements.png'))
%% Визуализация передающей части
% Построение сигнального созвездия
figure
fConst = scatter(real(outMod), imag(outMod), 30, 'k', 'filled');
grid on
xlabel('In-Phase')
ylabel('Qudrature')
title('Сигнальное созвездие исходного сигнала')
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
