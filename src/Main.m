clear; close all; clc;

Freq = 122.76e6; % Частота дискретизации
%% Формирование ПСП в NRZ
L1F_PT = zeros(1,511);
x = [1 1 1 1 1 1 1 1 1]; % Стартовая последовательность
for i = 1:511
   L1F_PT(i) = x(7);
   x_temporary = xor(x(5), x(9));
   x = circshift(x,[0 1]);
   x(1) = x_temporary;
end

L1F_PT(L1F_PT==1)=-1;
L1F_PT(L1F_PT==0)=1;

%% Баллистика НКА (все группировки ГЛОНАСС для 24 спутников - намеренно немного неправильно, но баллистика индивидуальна)
% Все координаты приведены в декартову систему координат
r_glonass = 19100000+6370000; % высота НС относительно поверхности Земли для ГЛОНАСС
incle_glonass = [64.8 184.8 304.8]*pi/180;
T_glonass = 11*3600+15*60+44;
t_modeling= 1/Freq;
omega_p_glonass = (1:8)*45*pi/180;
omega_set_glonass = [215.15 335.15 95.15]*pi/180;
NKA1 = 1; % Выберите космический аппарат от которого осуществляется передача сигнала наземному потребителю
k_p = rem(NKA1,8);
k_set = rem(NKA1,3);

phi = 60; % условная долгота потребителя
psi = 55; % условная широта потребителя
x_UE = 6370000*cos(phi*pi/180)*cos(psi*pi/180); % проекция потребителя на ось х в декартовой системе координат
y_UE = 6370000*sin(phi*pi/180)*cos(psi*pi/180); % проекция потребителя на ось y в декартовой системе координат
z_UE = 6370000*sin(psi*pi/180);                 % проекция потребителя на ось z в декартовой системе координат

L1F_PT = Sampling(1e-3,L1F_PT,Freq); % Дискретизируем входной сигнал
Signal_sampl = kron(ones(1,1e3),L1F_PT); % дополняем до 1 сек





% Сохраним в WV
wvfile = 'LF.wvd';
WV_file_id = fopen(wvfile,'a+'); % создаем файл для загрузки на внешнее файловое хранилище
for count = 1:T_glonass/4*Freq/1e3
     time_count = count*(1/Freq:1/Freq:1);
     % координаты спутника в декартовой системе координат
     x0_glonass = r_glonass*(cos(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*cos(omega_set_glonass(k_set)) - sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p)) * sin(omega_set_glonass(k_set))*cos(incle_glonass(k_set)));
     y0_glonass = r_glonass*(cos(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*sin(omega_set_glonass(k_set)) + sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p)) * cos(omega_set_glonass(k_set))*cos(incle_glonass(k_set)));
     z0_glonass = r_glonass*sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*sin(incle_glonass(k_set));
     % Взаимная дальность от наземного потребителя к КА
     R = sqrt((x0_glonass-x_UE).^2+(y0_glonass-y_UE).^2+(z0_glonass-z_UE).^2);
     % Задержка
     tau = R./physconst('LightSpeed');
     % вектор скорости для КА
     V_sat = [diff(x0_glonass-x_UE); diff(y0_glonass-y_UE); diff(z0_glonass-z_UE)];
     % относительный угол от плоскости вектора направленного к Земле к
     % вектору скорости в направлении движения КА
     teta = abs(acos((V_sat(1,:)+V_sat(2,:)+V_sat(1,:))./(sqrt(3)*sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2))));
     % проекция вектора скорости в декартову систему координат
     V_r_z = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*cos(teta);
     V_r_x = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*sin(teta);
     V_r_y = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*sin(teta);
     % доплеровское сдвиг частоты
     f_d = 2*pi*1602e6/physconst('LightSpeed')*sqrt(V_r_x.^2+V_r_y.^2+V_r_z.^2);
     
     % Сигнал к которому вводится доплеровский сдвиг и задержка
     Signal = Signal_sampl.*exp(2*pi*(f_d)*(time_count(1:end-1)+tau(1:end-1)));
     % Сохранение и представление отсчетов на устройстве векторного
     % генератора
     IQ_data_len = length(single(real(Signal)));
     IQ_data = single(zeros(1, 2*IQ_data_len));  
     IQ_data(1:2:end) = single(real(Signal));
     IQ_data(2:2:end) = single(imag(Signal));
     IQ_data = IQ_data /  max(abs(single(real(Signal)) + 1i*single(imag(Signal))));
     rms  = sqrt(mean(IQ_data(1:2:end).*IQ_data(1:2:end) + IQ_data(2:2:end).*IQ_data(2:2:end))) / 1.0;
     IQ_data = floor((IQ_data*32767+0.5)); 

     fwrite  (WV_file_id, IQ_data, 'int16');
end
fclose  (WV_file_id);



wvfile = 'LF.wvh'; 
WV_file_id = fopen(wvfile,'a+');
fprintf (WV_file_id, '%s','{TYPE:RAW16LE}' ); 
fprintf (WV_file_id, '%s','{COMPONENTS:IQ}' ); 
fprintf (WV_file_id, '%s',['{CLOCK:' num2str(Freq) '.000000}']); 
fprintf (WV_file_id, '%s','{RESOLUTION:16}');
fprintf (WV_file_id, '%s','{SAMPLES: ', num2str(count*IQ_data_len),'}');
fprintf (WV_file_id, '%s','{FREQUENCY:30.69000000.000000}');
fprintf (WV_file_id, '%s', ['{DATE: ' datestr(date,29) ';' datestr(clock,13) '}']);
fprintf (WV_file_id, '%s','{COPYRIGHT:2018 Rohde&Schwarz IQW}');
fclose  (WV_file_id);