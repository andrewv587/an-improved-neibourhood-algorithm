%%---------------------------------------------------------------
% NA反演程序
% 作者：黄卫航
% 日期：2015/12/31
% 工作单位：中国科学院地质与地球物理研究所
% 邮箱:huangweihang14@mails.ucas.ac.cn
% 参考文献 Sambridge M. Geophysical inversion with a neighbourhood 
% algorithm—I. Searching a parameter space[J]. 
% Geophysical Journal International, 1999, 138(2): 479-494.
%%---------------------------------------------------------------

%%---------------------------------------------------------------
%参数设置
%%---------------------------------------------------------------
clc;clear;%close all;
n_itr=5;     %迭代次数
ns=50;         %ns参数
nr=10;         %nr参数 可设为自适应参数 ns/nr为整数
ptime=0.1;     %单点暂停时间
n_par=2;       %参数个数
lowb=[-10 -10];  %下边界
upperb=[10 10];  %上边界
x_idx=1;       %显示i分量为x轴
y_idx=2;       %显示j分量为y轴
debug=0;       %是否开启调试模式,0关闭，1开启（n_par~=2时自动关闭）
method='NA';    %'MC' 蒙特卡洛算法,'NA' 邻近算法
%%---------------------------------------------------------------
%目标函数设置
%%---------------------------------------------------------------
%f=@(x)(x(1)-1.5)^2+(x(2)+3)^2+x(3)^2+1;  %目标函数
%测试函数 参照Page153 MATLAB优化算法安全分析与应用 
% f=@(x)0.5-(sin(sqrt(x(1)^2+x(2)^2))^2-0.5)/ ... 
%      (1+0.001*(x(1)^2+x(2)^2))^2;  %Schaffer函数
% f=@(x)-100*(x(1)^22-x(2))^2-(x(1)-1)^2; %Rosenbrock函数
f=@(x)1/4000*sum(x.^2)-prod(cos(x./(sqrt(1:length(x)))))+1; %Griewank函数
%  f=@(x)sum(x.^2-10*cos(2*pi*x)+10); %Rastrigin函数
% f=@(x)-20*exp(-0.2*sqrt(1/length(x)*sum(x.^2)))-exp(1/length(x)* ...
%     sum(cos(2*pi.*x)))+exp(1)+20;  %Ackley函数
%%---------------------------------------------------------------

%确保lowb、upperb与n_par一致
if length(lowb)~=n_par
    disp('lowb is not the same length as n_par');
    return;
end
if length(upperb)~=n_par
    disp('upperb is not the same length as n_par');
    return;
end

%确保上边界大于下边界
if sum(find(lowb>=upperb))>0
    disp('lowb is larger than upperb');
    return;
end

%参数个数大于2，关闭debug模式
if n_par~=2
    disp('Can''t open debug mode!');
    debug=0;
end

%生成图形，并指定坐标轴范围；
figure;
axis([lowb(x_idx) upperb(x_idx)  lowb(y_idx) upperb(y_idx)]);
grid on;
box on;
axis manual;
hold on;
xlabel([num2str(x_idx),' component']);
ylabel([num2str(y_idx),' component']);


% 参数为2且不为debug的情况下作出目标函数f(x)的等值图
if n_par==2 && debug~=1
    x=linspace(lowb(1),upperb(1),100);
    y=linspace(lowb(2),upperb(2),100);
    [X,Y]=meshgrid(x,y);
    z=X;
    for i=1:size(X,1)
        for j=1:size(X,2)
            z(i,j)=f([X(i,j) Y(i,j)]);
        end
    end
    contourf(X,Y,z);
    colorbar;
end

%Monte Carlo方法
if strcmp(method,'MC')
    data=zeros(ns*nr+nr,n_par);
    misfit=zeros(ns*nr+nr,n_par);
    for i=1:(ns*nr+nr)
    data(i,:)=lowb+(upperb-lowb).*rand([1,n_par]);
    end
    title(['Monte Carlo with point nums=',num2str(ns*nr+nr)]);
    scatter(data(:,x_idx),data(:,y_idx),'r','fill');
    [chdata, chindex]=sort(misfit);
    disp('data and misfist is:');
    disp([data(chindex(1),:),chdata(1)]);
    return
end

%%---------------------------------------------------------------
% NA算法
%%---------------------------------------------------------------
data=zeros(ns,n_par);
misfit=zeros(1,ns);
%随机生成ns个样本，并绘图
for i=1:ns
    data(i,:)=lowb+(upperb-lowb).*rand([1,n_par]);
    misfit(i)=f(data(i,:));
end
title(['neighbourhood algoritthm with ns=',num2str(ns),' nr=',num2str(nr)]);
scatter(data(:,x_idx),data(:,y_idx),'r','fill');

save parameters


if debug==1
    vdiagram=voronoi(data(:,x_idx),data(:,y_idx));
    pause;
end
%选择nr个样本点，生成ns个样本
datasum=data;
for itr=1:n_itr
    disp(['iteration number is ',num2str(itr)]);
    if  itr~=1
        data=oridata;
    end
    for i=1:size(data,1)
        misfit(i)=f(data(i,:));
    end
    %在前一次迭代生成的ns个样本中选择misfit最小的nr个
    [chdata, chindex]=sort(misfit); 
    oridata=[];
    disp(['data: ',num2str(data(chindex(1),:))]);
    disp(['misfit is ',num2str(chdata(1))]);
    for j=1:nr
        orpt=chindex(j);
        vdata=data(orpt,:);
        tmpdata=vdata;
        if debug==1
            scatter(tmpdata(x_idx),tmpdata(y_idx),140,'d');
        end 
        sca=[]; 
        for kk=1:floor(ns/nr)  
            if debug==1 
                scatmp=scatter(tmpdata(x_idx),tmpdata(y_idx),140,'p','fill');
                if j==1
                    disp(['kk=',num2str(kk),' data=',num2str(tmpdata)]);
                end
                pause;
            end
            for i=1:n_par  %进行n_par个坐标循环
%               disp(['  i=',num2str(i)]);
                range=[lowb(i) upperb(i)];
                tmp2data=tmpdata;
                tmp2data(i)=vdata(i);
                vji=vdata(i);        
                dj2=sum((vdata-tmp2data).^2); 
                for k=1:size(datasum,1)
                    if k==orpt
                        continue;
                    end
                    tmp2data(i)=datasum(k,i);
                    vki=datasum(k,i);
                    dk2=sum((datasum(k,:)-tmp2data).^2);
                    if abs(vki-vji)<1e-10
                        continue;
                    end
                    xji=1/2*(vki+vji+(dk2-dj2)/(vki-vji));%计算xj的第i个分量
                    kkk=0;
                    if xji>range(1) && xji<=tmpdata(i)
                        range(1)=xji;
                        kkk=1;
                    elseif xji<range(2) && xji>=tmpdata(i)
                        range(2)=xji;
                        kkk=1;
                    end
                    if debug==1
                        if kkk==1
                            scad=scatter(datasum(k,x_idx),datasum(k,y_idx),80,'s');
                            if i==1
                                scaa=scatter(xji,tmpdata(2),30,'r+');
                                linea=line([vdata(1) xji],[vdata(2) tmpdata(2)]);
                                lineb=line([datasum(k,1) xji],[datasum(k,2) tmpdata(2)]);
                                %disa=sum(([vdata(1) vdata(2)]-[xji tmpdata(2)]).^2);
                                %disb=sum(([data(k,1) data(k,2)]-[xji tmpdata(2)]).^2);                        
                            else
                                scaa=scatter(tmpdata(1),xji,30,'r+');
                                linea=line([vdata(1) tmpdata(1)],[vdata(2) xji]);
                                lineb=line([datasum(k,1) tmpdata(1)],[datasum(k,2) xji]);
    %                           %disa=sum(([vdata(1) vdata(2)]-[tmpdata(1) xji]).^2);
    %                           %disb=sum(([data(k,1) data(k,2)]-[tmpdata(1) xji]).^2);                        
                            end
                            sca=[sca scaa];
                            %测试||vj-xj||==||vj-kj||
                            %disp([disa disb]);
                            pause;
                            delete(linea);
                            delete(lineb);
                            delete(scad);
                        end
                    end
                end
                newx_loc=range(1)+(range(2)-range(1))*rand; %均分分布生成新的坐标
                tmpdata(i)=newx_loc;
                if debug==1
%                     disp(['range is: ',num2str(range(1)),' ',num2str(range(2))]);
%                     disp(['newx_loc: ',num2str(newx_loc)]);
                    if i==1
                        scatter(tmpdata(x_idx),tmpdata(y_idx),50,'x');
                        pause;
                    end
                end
            end
            scatter(tmpdata(x_idx),tmpdata(y_idx),40,'b','fill');
            if debug==1
                pause;
                delete(sca);
                delete(scatmp);
            else
                pause(0.01);
            end
            oridata=[oridata;tmpdata]; %该次迭代生成的ns个样本
        end
        if debug==1
             delete(findobj(gca,'marker','x'));
        end
    end
    datasum=[datasum;oridata]; %生成的总样本数
    if debug==1
        delete(findobj(gca,'marker','d'));
        delete(vdiagram);
        vdiagram=voronoi(datasum(:,x_idx),datasum(:,y_idx));
    end
end
