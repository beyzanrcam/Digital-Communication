clear 
close all
clc
% 1. soru İlk şık
A=2;
T=1;
Ts= 0.000001;
t=0: Ts : T-Ts;

s1 = @(t) A.*sin((2*pi*t)/T) .* (0<=t & t<=T/2);
s2 = @(t) -s1(t-T/2) .* (T/2<=t & t<=T);

figure;
subplot(2,1,1)
plot(t,s1(t))
xlabel('t(sn)')
ylabel('s1(t)')

subplot(2,1,2)
plot(t, s2(t))
xlabel('t(sn)')
ylabel('s2(t)')



Eb_calc = (A.^2)*T/2 ;
SNRdb_calc= 0:0.01:16;
Pb_calc= qfunc(sqrt(10.^(SNRdb_calc/10)));


figure;
semilogy(SNRdb_calc, Pb_calc);
xlabel('SNR(db) calculated');
ylabel('Pb calculated');
grid on;

% 1.soru 2. şık
Eb = 1 ;
SNRdb = 0:16;



a1= 1; % Eb=1 olduğundan bunlar hesaplanabilir
a2= -1;
gama= (a1+a2)/2;


bit_sayim=10^7;
bitler = randi([0,1], 1, bit_sayim);

bit= zeros(1, length(bitler));

bit(bitler==1) = a1;
bit(bitler==0) = a2;


Pb= zeros(1, length(SNRdb));
N0 = 10.^(-SNRdb./10);

for i= 1:length(SNRdb)
    
    n0 = sqrt(N0(i)).*randn(1,bit_sayim); %varyansı N0'ın karekökü olcak onla çarpılması lazım
    z = bit + n0;

    compare = z>=gama;
    bit_hata=xor(bitler,compare);
    Pb(i) = sum(bit_hata)/bit_sayim;
end


figure;
semilogy(SNRdb, Pb);
xlabel('SNR(db)');
ylabel('Pb');
grid on;

%%
% 2. soru ilk şık

gama_b= N0/2;

P1=0.75;
P2=0.25;

sigma_b=sqrt(N0);

SNRdb_calc_b= 0:16;
Pb_calc_b= (1-qfunc((gama_b-a1)./sigma_b)).*P2  + (qfunc((gama_b-a2)./sigma_b)).*P1  ;


figure;
semilogy(SNRdb_calc_b, Pb_calc_b);
xlabel('SNR(db) calculated b');
ylabel('Pb calculated b');
grid on;

%2. soru 2.kısım

bitler_b = randsrc(1,bit_sayim,[0 1; P1 P2]);

bit_b= zeros(1, length(bitler_b));

bit_b(bitler_b==1) = a1;
bit_b(bitler_b==0) = a2;

Pb= zeros(1, length(SNRdb_calc_b));
N0 = 10.^(-SNRdb_calc_b./10);

for i= 1:length(SNRdb_calc_b)
    
    n0 = sqrt(N0(i)).*randn(1,bit_sayim); %varyansı N0'ın karekökü olcak onla çarpılması lazım
    z_b = bit_b + n0;

    compare_b = z_b>=gama_b(i);
    bit_hata_b=xor(bitler_b,compare_b);
    Pb(i) = sum(bit_hata_b)/bit_sayim;
end


figure;
semilogy(SNRdb, Pb);
xlabel('SNR(db) b');
ylabel('Pb b');
grid on;


