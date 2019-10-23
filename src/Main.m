clear; close all; clc;

Freq = 122.76e6; % ������� �������������
%% ������������ ��� � NRZ
L1F_PT = zeros(1,511);
x = [1 1 1 1 1 1 1 1 1]; % ��������� ������������������
for i = 1:511
   L1F_PT(i) = x(7);
   x_temporary = xor(x(5), x(9));
   x = circshift(x,[0 1]);
   x(1) = x_temporary;
end

L1F_PT(L1F_PT==1)=-1;
L1F_PT(L1F_PT==0)=1;

%% ���������� ��� (��� ����������� ������� ��� 24 ��������� - ��������� ������� �����������, �� ���������� �������������)
% ��� ���������� ��������� � ��������� ������� ���������
r_glonass = 19100000+6370000; % ������ �� ������������ ����������� ����� ��� �������
incle_glonass = [64.8 184.8 304.8]*pi/180;
T_glonass = 11*3600+15*60+44;
t_modeling= 1/Freq;
omega_p_glonass = (1:8)*45*pi/180;
omega_set_glonass = [215.15 335.15 95.15]*pi/180;
NKA1 = 1; % �������� ����������� ������� �� �������� �������������� �������� ������� ��������� �����������
k_p = rem(NKA1,8);
k_set = rem(NKA1,3);

phi = 60; % �������� ������� �����������
psi = 55; % �������� ������ �����������
x_UE = 6370000*cos(phi*pi/180)*cos(psi*pi/180); % �������� ����������� �� ��� � � ���������� ������� ���������
y_UE = 6370000*sin(phi*pi/180)*cos(psi*pi/180); % �������� ����������� �� ��� y � ���������� ������� ���������
z_UE = 6370000*sin(psi*pi/180);                 % �������� ����������� �� ��� z � ���������� ������� ���������

L1F_PT = Sampling(1e-3,L1F_PT,Freq); % �������������� ������� ������
Signal_sampl = kron(ones(1,1e3),L1F_PT); % ��������� �� 1 ���





% �������� � WV
wvfile = 'LF.wvd';
WV_file_id = fopen(wvfile,'a+'); % ������� ���� ��� �������� �� ������� �������� ���������
for count = 1:T_glonass/4*Freq/1e3
     time_count = count*(1/Freq:1/Freq:1);
     % ���������� �������� � ���������� ������� ���������
     x0_glonass = r_glonass*(cos(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*cos(omega_set_glonass(k_set)) - sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p)) * sin(omega_set_glonass(k_set))*cos(incle_glonass(k_set)));
     y0_glonass = r_glonass*(cos(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*sin(omega_set_glonass(k_set)) + sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p)) * cos(omega_set_glonass(k_set))*cos(incle_glonass(k_set)));
     z0_glonass = r_glonass*sin(2*pi*1/T_glonass*time_count + omega_p_glonass(k_p))*sin(incle_glonass(k_set));
     % �������� ��������� �� ��������� ����������� � ��
     R = sqrt((x0_glonass-x_UE).^2+(y0_glonass-y_UE).^2+(z0_glonass-z_UE).^2);
     % ��������
     tau = R./physconst('LightSpeed');
     % ������ �������� ��� ��
     V_sat = [diff(x0_glonass-x_UE); diff(y0_glonass-y_UE); diff(z0_glonass-z_UE)];
     % ������������� ���� �� ��������� ������� ������������� � ����� �
     % ������� �������� � ����������� �������� ��
     teta = abs(acos((V_sat(1,:)+V_sat(2,:)+V_sat(1,:))./(sqrt(3)*sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2))));
     % �������� ������� �������� � ��������� ������� ���������
     V_r_z = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*cos(teta);
     V_r_x = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*sin(teta);
     V_r_y = sqrt(V_sat(1,:).^2+V_sat(2,:).^2+V_sat(3,:).^2).*sin(teta);
     % ������������ ����� �������
     f_d = 2*pi*1602e6/physconst('LightSpeed')*sqrt(V_r_x.^2+V_r_y.^2+V_r_z.^2);
     
     % ������ � �������� �������� ������������ ����� � ��������
     Signal = Signal_sampl.*exp(2*pi*(f_d)*(time_count(1:end-1)+tau(1:end-1)));
     % ���������� � ������������� �������� �� ���������� ����������
     % ����������
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