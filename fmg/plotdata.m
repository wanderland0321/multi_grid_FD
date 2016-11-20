% plot 'data.txt'
clear;

A=importdata('data.txt');
x=A(:,1);
y1=A(:,2);
y2=A(:,3);

semilogy(x,y1,'r-*');
hold on;
semilogy(x,y2,'b-*');

xlabel('# of iterations');
ylabel('log(L2error)');
legend('L2error','L2residual');
title('FMG, Nx=Ny=64');
