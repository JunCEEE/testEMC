%plotFSC.m

fn = '/beegfs/desy/user/kimyoon/S14_r112_p_only_0001_cen2bst/fsc_cross_avg_S14_r112_p_only_0001_cen2bst_4series_avg20.mat';
load(fn)

figure;
p1=errorbar(fsc_avg(1:r_max_1),fsc_std(1:r_max_1));
hold on;
p2=plot(threshold_main(1:r_max_1));
xlabel('Half period resolution (nm)');
% ylabel('fsc S14 r112 poission AGIPD 0001 cen2bst 4series avg20');

% ylabel_name = ['fsc' dataset_name 'r112' noise_name bst_area '4series avg20']
% ylabel_name = ['fsc' dataset_name 'r112' noise_name bst_area '4series rand avg5']
% ylabel(ylabel_name)


p1.Color='b';p1.LineStyle='none';p1.LineWidth=1;p1.Marker='o';p1.MarkerSize=4
p2.Color='r';p2.LineStyle=':';p2.LineWidth=3;%hold off;

xticks([16 32 53 80 106 138 159])
xticklabels({'10', '5','3','2','1.5','1.15','1'})

% axis([0 138 -0.2 1])
% hold off