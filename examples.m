% Some examples:

%% Plots of Turing patterns for different topologies:

figure(1);
makeGraph(10,'1Dlattice','model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoLattice','epsc')

figure(2);
makeGraph(10,'cycle','model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoCycle','epsc')

figure(3);
makeGraph(10,'complete','model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoComplete','epsc')

figure(4);
makeGraph(10,'regular','circulant', [0,1,0,0,0,1,0,0,0,1],'model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoCirculant','epsc')

figure(5);
makeGraph(25,'cycle','model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoCycle25','epsc')

figure(6);
makeGraph(25,'star','model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoStar25','epsc')

figure(7);
makeGraph(10,'random','erdos',0.4,'model',1,'parameters',[1,0.8,0.2],'d',[1,0.01]);
saveas(gcf,'homoRandom','epsc')

%% Convergence plots:

rng(20);
N1 = 10;
N2 = 10;
A1 = diag(ones(N1-1,1),1); A1 = A1+A1'; 
A2 = diag(ones(N2-1,1),1); A2 = A2+A2'; 

C0 = rand(N1,N2)<0.2;
if(C0 == zeros(N1,N2))
    C0(1,1) = 1;
end

fu1 = -1;
fu2 = -0.1;

ende = 5;
error = zeros(1,ende);
epsilon = zeros(1,ende);
for i = 1:ende
    e = 10^(-i);
    epsilon(i) = e;
[evalFull, evecFull, evalAsym0, evecAsym0, evalAsym1, lambdaAsym] = ComputeGraphAsymptotics(A1, A2, C0, e, fu1, fu2);
error(i) = norm(evalFull-lambdaAsym,Inf); 
end

p = polyfit(log(epsilon),log(error),1);
z = polyval(p,log(epsilon));
loglog(epsilon,exp(z),'r-','Linewidth',1.2);
hold on 
loglog(epsilon,error,'b.','Markersize',13);
set(gca,'Fontsize',13)
xlabel('$\varepsilon$','interpreter', 'latex', 'FontSize',20)
ylabel('$||\lambda - \lambda_a||_\infty$','interpreter', 'latex','FontSize',18)
h = legend(strcat('Fit $m=$ ', num2str(p(1))),'Error','location','southeast','FontSize',14);
set(h,'interpreter','latex')
hold off
saveas(gcf,'ConvergenceScalar','epsc')
saveas(gcf,'ConvergenceScalar','fig')

%% Plots of a Turing parameter space:

rng(100);
N1 = 10;
N2 = 10;
A1 = makeGraph(N1,'regular', 'circulant', [0,1,1,1,0,0,0,1,1,1]);
A2 = makeGraph(N2,'1Dlattice');

alpha1 = 1;
alpha2 = 1;
beta2 = 0.6;
zeta2 = 0.4;

e = 0.1;
d = [1,0.01,1,0.1];

betas = 2.5:-0.01:0.5;
zetas = 0:0.01:0.3;

C0 = rand(N1,N2)<0.2;

res = zeros(length(betas),length(zetas));
for i = 1:length(betas)
    beta1 = betas(i);
    for j = 1:length(zetas)
        zeta1 = zetas(j);
        J1 = [-(beta1+zeta1)^2/alpha1^2, -2*beta1*alpha1/(beta1+zeta1);(beta1+zeta1)^2/alpha1^2, 2*beta1*alpha1/(beta1+zeta1)-alpha1];
        J2 = [-(beta2+zeta2)^2/alpha2^2, -2*beta2*alpha2/(beta2+zeta2);(beta2+zeta2)^2/alpha2^2, 2*beta2*alpha2/(beta2+zeta2)-alpha2];
        evalJ = [eig(J1);eig(J2)];
        
        if(sum(real(evalJ)>0)>0)
            res(i,j) = 1;
        else
            evalFull = ComputeGraphAsymptoticsSystem(A1, A2, C0, e, J1, J2, d);
            if(sum(real(evalFull)>0)>0)
                res(i,j) = 2;
            end
        end
    end
end

imagesc( zetas,betas, res)
set(gca,'Fontsize',22)
xlabel('$\zeta$','interpreter', 'latex', 'FontSize',30);
ylabel('$\beta$','interpreter', 'latex', 'FontSize',30);
set(gca,'YDir','normal')
caxis([0, 2]);
saveas(gcf,'space_Circ6_Lat','epsc')
