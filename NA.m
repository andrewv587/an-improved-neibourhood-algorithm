%%---------------------------------------------------------------
% NA���ݳ���
% ���ߣ�������
% ���ڣ�2015/12/31
% ������λ���й���ѧԺ��������������о���
% ����:huangweihang14@mails.ucas.ac.cn
% �ο����� Sambridge M. Geophysical inversion with a neighbourhood 
% algorithm��I. Searching a parameter space[J]. 
% Geophysical Journal International, 1999, 138(2): 479-494.
%%---------------------------------------------------------------

%%---------------------------------------------------------------
%��������
%%---------------------------------------------------------------
clc;clear;%close all;
n_itr=5;     %��������
ns=50;         %ns����
nr=10;         %nr���� ����Ϊ����Ӧ���� ns/nrΪ����
ptime=0.1;     %������ͣʱ��
n_par=2;       %��������
lowb=[-10 -10];  %�±߽�
upperb=[10 10];  %�ϱ߽�
x_idx=1;       %��ʾi����Ϊx��
y_idx=2;       %��ʾj����Ϊy��
debug=0;       %�Ƿ�������ģʽ,0�رգ�1������n_par~=2ʱ�Զ��رգ�
method='NA';    %'MC' ���ؿ����㷨,'NA' �ڽ��㷨
%%---------------------------------------------------------------
%Ŀ�꺯������
%%---------------------------------------------------------------
%f=@(x)(x(1)-1.5)^2+(x(2)+3)^2+x(3)^2+1;  %Ŀ�꺯��
%���Ժ��� ����Page153 MATLAB�Ż��㷨��ȫ������Ӧ�� 
% f=@(x)0.5-(sin(sqrt(x(1)^2+x(2)^2))^2-0.5)/ ... 
%      (1+0.001*(x(1)^2+x(2)^2))^2;  %Schaffer����
% f=@(x)-100*(x(1)^22-x(2))^2-(x(1)-1)^2; %Rosenbrock����
f=@(x)1/4000*sum(x.^2)-prod(cos(x./(sqrt(1:length(x)))))+1; %Griewank����
%  f=@(x)sum(x.^2-10*cos(2*pi*x)+10); %Rastrigin����
% f=@(x)-20*exp(-0.2*sqrt(1/length(x)*sum(x.^2)))-exp(1/length(x)* ...
%     sum(cos(2*pi.*x)))+exp(1)+20;  %Ackley����
%%---------------------------------------------------------------

%ȷ��lowb��upperb��n_parһ��
if length(lowb)~=n_par
    disp('lowb is not the same length as n_par');
    return;
end
if length(upperb)~=n_par
    disp('upperb is not the same length as n_par');
    return;
end

%ȷ���ϱ߽�����±߽�
if sum(find(lowb>=upperb))>0
    disp('lowb is larger than upperb');
    return;
end

%������������2���ر�debugģʽ
if n_par~=2
    disp('Can''t open debug mode!');
    debug=0;
end

%����ͼ�Σ���ָ�������᷶Χ��
figure;
axis([lowb(x_idx) upperb(x_idx)  lowb(y_idx) upperb(y_idx)]);
grid on;
box on;
axis manual;
hold on;
xlabel([num2str(x_idx),' component']);
ylabel([num2str(y_idx),' component']);


% ����Ϊ2�Ҳ�Ϊdebug�����������Ŀ�꺯��f(x)�ĵ�ֵͼ
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

%Monte Carlo����
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
% NA�㷨
%%---------------------------------------------------------------
data=zeros(ns,n_par);
misfit=zeros(1,ns);
%�������ns������������ͼ
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
%ѡ��nr�������㣬����ns������
datasum=data;
for itr=1:n_itr
    disp(['iteration number is ',num2str(itr)]);
    if  itr~=1
        data=oridata;
    end
    for i=1:size(data,1)
        misfit(i)=f(data(i,:));
    end
    %��ǰһ�ε������ɵ�ns��������ѡ��misfit��С��nr��
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
            for i=1:n_par  %����n_par������ѭ��
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
                    xji=1/2*(vki+vji+(dk2-dj2)/(vki-vji));%����xj�ĵ�i������
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
                            %����||vj-xj||==||vj-kj||
                            %disp([disa disb]);
                            pause;
                            delete(linea);
                            delete(lineb);
                            delete(scad);
                        end
                    end
                end
                newx_loc=range(1)+(range(2)-range(1))*rand; %���ֲַ������µ�����
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
            oridata=[oridata;tmpdata]; %�ôε������ɵ�ns������
        end
        if debug==1
             delete(findobj(gca,'marker','x'));
        end
    end
    datasum=[datasum;oridata]; %���ɵ���������
    if debug==1
        delete(findobj(gca,'marker','d'));
        delete(vdiagram);
        vdiagram=voronoi(datasum(:,x_idx),datasum(:,y_idx));
    end
end
