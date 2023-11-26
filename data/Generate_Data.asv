%% Generate_Data %%

clear; 
close all;
addpath 'D:\310 科研'
addpath 'D:\310 科研\000 Simulation'
load("orange.mat")
load("purple.mat")
load("matlab1109.mat")
rr = r;
clear S;
S = 24;
Time = 200;
Std_total = [];
Stability = -ones(Time,1);
Input_Abundance = zeros(Time,54);
Input_RelativeAbundance = zeros(Time,54);
Input_Presence = zeros(Time,54);
H = tiledlayout(10,10,"TileSpacing","compact");
for time = 1:Time
    Rand = randperm(54);
    species = sort(Rand(1:S));
    AA = A(species, species);
    r = rr(species);
    clear N;
    N(1,:) = rand(1,S);
    Input_Abundance(time,species) = N(1,:);
    Input_RelativeAbundance(time,species) = N(1,:)/sum(N(1,:));
    Input_Presence(time,species) = ones(1,S);
    for day = 0:(Day-1)
        for i = T*day + 2: T*day + T + 1

            for j=1:S
                k1=N(i-1,j)*(r(j)-AA(j,:)*(N(i-1,:)'))*step;
                k2=(N(i-1,j)+k1/2)*(r(j)-AA(j,:)*(N(i-1,:)'))*step;
                k3=(N(i-1,j)+k2/2)*(r(j)-AA(j,:)*(N(i-1,:)'))*step;
                k4=(N(i-1,j)+k3)*(r(j)-AA(j,:)*(N(i-1,:)'))*step;
                N(i,j)=N(i-1,j)+(1/6)*(k1+2*k2+2*k3+k4)+10^-6*step;
                if N(i,j)>1
                    N(i,j) = 1;
                end
            end

        end
        %     N(i,:) = N(i,:)/30;
    end
    Std = std(N(1800:2000,:));
    Std_total = [Std_total,max(Std)];
    if (max(Std)>0.01)
        Stability(time) = 1;
    else
        Stability(time) = 0;
    end
end
