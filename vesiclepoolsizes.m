clearvars -except rawdata_measurement*
close all


%% User input

pathname=('V:\AG_Neurophotonik\Projekte\Fluoxetin\Messdaten\');  % Path name of recordings
commonfilenameend=('.tif'); % File type of recording
savename='results_fluox.mat'; % Save results as

yn_lzeros=1;  % Number of leading zeros in recording names
n_pics=310;     % Number of images for one recording

% Number of frames at which baseline and stimulations start and end
baseline_start=25;
baseline_stop=36;

rrpstim_start=47;
rrpstim_stop=57;

vorrecstim_start=71;
vorrecstim_stop=81;

synstim_start=126;
synstim_stop=133;

recstim_start=196;
recstim_stop=207;

vortotalstim_start=207;
vortotalstim_stop=217;

totalstim_start=227;
totalstim_stop=255;

start_measurement=1; % Number of first measurement to be analyzed
n_measurements=34;  % Number of measurements to be analyzed


type='fluox';     % Label of this experimental group for figures
show_specified_traces=0; % Do you want to have an exemplary traces plotted? 1 or 0
loaddata=0; % Should raw data be loaded from the images? 1 or 0   ___ Load raw data from the images the first time. After that, the fluorescence traces can be loaded directly from the matlab-file, which saves a lot of processing time
first_run=1; % Should fluorescence traces be loaded from the matlab-file or are they already loaded into the workspace? 1 or 0
show_control_figures=0; % Do you wish to see exemplary figures along the processing? 1 or 0
realy_show_every_synapse=0; % Do you wish to see to most detailed versions of the exemplary figures along the processing? 1 or 0  (if this is activated, the script produces a lot of images; do not activate during batch processing of several recordings)





%% Initialize variables

results=[];





%% 1. Load raw images from files
if loaddata==1;
    for measurement=start_measurement:n_measurements
        filenamestem=sprintf('vesiclepool_recording_%i_image',measurement');
        
        
        
        for i=1:n_pics
            % Attach image number to the filename stem, including the according
            % number of leading zeros
            if yn_lzeros==1
                zstr=char();
                for j=1:(length(num2str(n_pics))-length(num2str(i)))
                    zstr=[zstr num2str(0)];
                end
                zstr=[zstr num2str(i)];
            else
                zstr=num2str(i);
            end
            
            % Load images
            eval(sprintf('rawdata_measurement%i(:,:,i)=imread([pathname filenamestem zstr commonfilenameend]);',measurement));
        end
        
        fprintf('Data loaded for measurement %i \n', measurement);
        
    end
    
    
    % Save the rawdata, so that next time it doesn't take so long to
    % process anymore.
    save('rawdata_fluox.mat', 'rawdata_measurement*')
    fprintf('Data saved\n')
    
else
    if first_run==1; % If the rawdata is already saved as a matlab file deactivate loaddata and activate first_run instead to save processing time
        load('rawdata_fluox.mat')
        fprintf('Data loaded\n')
    end
end

[ x_werte ] = scale_x_axis( n_pics ); % Use this function if the recording had phases with different acquisition frequencies (e.g. low frequency to save data during equilibration phase)

if show_specified_traces==1;
    %Plot single trace
    figure
    temp=rawdata_measurement1(287,398,:);
    plot(temp(:))
    clear temp
    %Plot all single traces
    figure
    hold on
    for i=1:size(rawdata_measurement1,3)
        temp=rawdata_measurement1(:,:,i);
        temp3=temp(temp>0);
        temp_mean(i)=mean(temp3(:));
    end
    for r=1:50
        for c=1:50
            temp2=rawdata_measurement1(r,c,:);
            plot(temp2(:), 'Color', rand(1,3))
        end
    end
    plot(temp_mean, 'LineWidth', 3)
    title('All single traces')
    clear temp2 temp_mean temp temp3
end



for measurement=start_measurement:n_measurements
    
    % Get baseline
    eval(sprintf('baseline_stack=rawdata_measurement%i(:,:,baseline_start:baseline_stop);', measurement));
    baseline_stack=double(baseline_stack);
    % Average baseline
    baselineimg=mean(baseline_stack,3);
    
    % Get stimulation response
    eval(sprintf('synstim_stack=rawdata_measurement%i(:,:,synstim_start:synstim_stop);', measurement));
    synstim_stack=double(synstim_stack);
    % Average stimulation response
    synstimimg=mean(synstim_stack,3);
    
    % Subtract baseline from stimulated response to filter out unresponsive
    % artefacts
    diff_img=synstimimg-baselineimg;
    
    if show_control_figures==1
        figure
        subplot(2,2,1)
        imagesc(baselineimg)
        title('Baseline')
        subplot(2,2,2)
        imagesc(synstimimg)
        title('Stimulated response')
        subplot(2,2,3)
        imagesc(diff_img)
        title(sprintf('Diff Image of Measurement %i', measurement))
        subplot(2,2,4)
        hist(diff_img(:), 1000)
        title('Hist of diff image')
    end
    
    % Get cell detection threshold
    [ith] = thresh_img(diff_img,2,0, show_control_figures);
    
    % ROI detection with laplace filtering
    % Here a different algorithm such as feature point detection by Sbalzarini et al. can be used if
    % the results are not satisfactory.
    [region_properties, L2, bw2, num2, mF, sF, area, cent] = peakdetect_laplace_ith(diff_img,ith.threshold, show_control_figures);
    
    
    % Get fluorescence traces from the regions of interest
    numb_regions(measurement,1)=size(region_properties,1);
    c2=1;
    for i=1:size(region_properties,1)
        
        for frame=1:size(rawdata_measurement1,3)
            eval(sprintf('temp=rawdata_measurement%i(:,:,frame);', measurement));
            pixelstack(:,frame)=temp(region_properties(i,1).PixelIdxList);
            meanstack(1,frame)=mean(pixelstack(:,frame));
        end
        
        
        % For each region get the baseline level, the levels befor the
        % electrical stimulations and the fluorescence during total pool
        % activation (ammonium chloride perfusion)
        
        baseline(measurement,i)=mean(meanstack(1, baseline_start:baseline_stop));
        rrpstim(measurement,i)=mean(meanstack(1,rrpstim_start:rrpstim_stop));
        vorrecstim(measurement,i)=mean(meanstack(1,vorrecstim_start:vorrecstim_stop));
        recstim(measurement,i)=mean(meanstack(1,recstim_start:recstim_stop));
        vortotalstim(measurement,i)=mean(meanstack(1,vortotalstim_start:vortotalstim_stop));
        totalstim(measurement,i)=mean(meanstack(1,totalstim_start:totalstim_stop));
        
        
        % Normalize fluorescence traces
        norm_meanstack(i,:, measurement)=(meanstack-baseline(measurement,i))/(totalstim(measurement,i)-baseline(measurement,i));
        
        
        % Plot traces of some synapses
        if realy_show_every_synapse==1
            switch i
                case {1,3,50,64,72,90}
                    figure
                    subplot(3,1,1)
                    hold on
                    for t1=1:size(pixelstack,1)
                        plot(pixelstack(t1,:,1), 'Color', rand(1,3))
                    end
                    title(sprintf('Measurement %i - Synapse %i -  Pixelstack',measurement, i))
                    subplot(3,1,2)
                    plot(meanstack(1,:))
                    title('Mean Stack')
                    subplot(3,1,3)
                    plot(norm_meanstack(i,:,measurement), 'Color', 'red', 'LineWidth', 2)
                    title('Norm Meanstack')
            end
        end
        
        clearvars temp pixelstack meanstack
    end
    
    
    
    % Calculate relative vesicle pool sizes
    rrpoolanteil=(rrpstim-baseline)./(rrpstim-baseline+recstim-vorrecstim+totalstim-vortotalstim); %die Anteile werden für jede ROI berechnet
    recpoolanteil=(recstim-vorrecstim+rrpstim-baseline)./(rrpstim-baseline+recstim-vorrecstim+totalstim-vortotalstim);
    
    
    fprintf('Finished calculations for measurement %i\n', measurement)
    
    
end


% clean out regions of interest with relative vesicle pool sizes not
% between 0 and 1
rrpoolanteil_clean1=rrpoolanteil.*(rrpoolanteil>0);
rrpoolanteil_clean=rrpoolanteil_clean1.*(rrpoolanteil_clean1<1);
recpoolanteil_clean1=recpoolanteil.*(recpoolanteil>0);
recpoolanteil_clean=recpoolanteil_clean1.*(recpoolanteil_clean1<1);



%% Average results across multiple measurements
for measurement=start_measurement:n_measurements
    
    
    % Clean out regions of interest with relative vesicle pool sizes not
    % between 0 and 1
    rrpool_clean_row=rrpoolanteil_clean(measurement,:);
    recpool_clean_row=recpoolanteil_clean(measurement,:);
    maske_rrp=rrpool_clean_row>0;
    maske_rec=recpool_clean_row>0;
    temp=maske_rrp+maske_rec;
    maske=temp>1;
    clear temp
    
    % Average across all synapses of a recording
    rrpoolanteil_mean(measurement, 1)=mean(rrpool_clean_row(maske));
    recpoolanteil_mean(measurement,1)=mean(recpool_clean_row(maske)); % mean across all synapses
    baseline_masked(measurement,1)=mean(baseline(measurement,maske));
    rrpstim_masked(measurement,1)=mean(rrpstim(measurement,maske));
    vorrecstim_masked(measurement,1)=mean(vorrecstim(measurement,maske));
    recstim_masked(measurement,1)=mean(recstim(measurement,maske));
    vortotalstim_masked(measurement,1)=mean(vortotalstim(measurement,maske));
    totalstim_masked(measurement,1)=mean(totalstim(measurement,maske));
    
    
    temp=norm_meanstack(maske',:,measurement);
    norm_kurve(1,:,measurement)=mean(temp(:,:),1);
    norm_meanstack_clean(1:size(temp,1),1:size(temp,2),measurement)=temp;
    norm_kurve(2,:,measurement)=std(norm_meanstack(maske',:,measurement),1)./size(norm_meanstack,1);
    
    
    
    figure
    errorbar(x_werte(1:size(rawdata_measurement1,3)), norm_kurve(1,:,measurement),norm_kurve(2,:,measurement))
    hold on
    line([x_werte(1,baseline_start) x_werte(1,baseline_stop)], [(baseline_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (baseline_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'Color', 'red');
    line([x_werte(1,rrpstim_start) x_werte(1,rrpstim_stop)], [(rrpstim_masked(measurement, 1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (rrpstim_masked(measurement, 1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'Color', 'red');
    line([x_werte(1,vorrecstim_start) x_werte(1,vorrecstim_stop)], [(vorrecstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (vorrecstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'Color', 'red');
    line([x_werte(1,recstim_start) x_werte(1,recstim_stop)], [(recstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (recstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'color', 'red');
    line([x_werte(1,vortotalstim_start) x_werte(1,vortotalstim_stop)], [(vortotalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (vortotalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'Color', 'red');
    line([x_werte(1,totalstim_start) x_werte(1,totalstim_stop)], [(totalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1)) (totalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1))], 'Color', 'red');
    title(sprintf('Mean curve of measurement %i', measurement))
    
    clear  temp
    
    if realy_show_every_synapse==1
        
        
        figure
        hold on
        for t1=1:size(norm_meanstack_clean,1)
            switch t1
                case {1,3,5,6,10,15,22,31,34,40,42,43,55}
                    plot(x_werte, norm_meanstack_clean(t1,:,measurement), 'Color', rand(1,3))
            end
        end
        plot(x_werte, norm_kurve(1,:,measurement), 'Color', 'red', 'LineWidth', 2)
        title(sprintf('Measurement %i -  Norm Meanstack Clean',measurement))
        
    end
end

%%%%%%%%%% Average across multiple recordings

count=1;
for measurement=start_measurement:n_measurements
    switch measurement
        case{1,2,4,6,7,13,15,16,17,18,19,21,23,24,25,29,30,31,33,34}  % Choose which recordings should be included in the averaging
            norm_kurve_ausgewaehlt(1,:,count)=norm_kurve(1,:,measurement);
            rrpoolanteil_mean_ausgewaehlt(count,1)=rrpoolanteil_mean(measurement,1);
            recpoolanteil_mean_ausgewaehlt(count,1)=recpoolanteil_mean(measurement,1);
            baseline_ausgewaehlt(count,1)=(baseline_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            rrpstim_ausgewaehlt(count,1)=(rrpstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            vorrecstim_ausgewaehlt(count,1)=(vorrecstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            recstim_ausgewaehlt(count,1)=(recstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            vortotalstim_ausgewaehlt(count,1)=(vortotalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            totalstim_ausgewaehlt(count,1)=(totalstim_masked(measurement,1)-baseline_masked(measurement,1))/(totalstim_masked(measurement,1)-baseline_masked(measurement,1));
            count=count+1;
    end
end



mittelwert_kurve(1,:)=mean(norm_kurve_ausgewaehlt(1,:,:),3);
mittelwert_kurve(2,:)=std(norm_kurve_ausgewaehlt(1,:,:),0,3)./(size(norm_kurve_ausgewaehlt,3));
mean_rrpoolanteil=mean(rrpoolanteil_mean_ausgewaehlt);
mean_recpoolanteil=mean(recpoolanteil_mean_ausgewaehlt);
mean_baseline=mean(baseline_ausgewaehlt);
mean_rrpstim=mean(rrpstim_ausgewaehlt);
mean_vorrecstim=mean(vorrecstim_ausgewaehlt);
mean_recstim=mean(recstim_ausgewaehlt);
mean_vortotalstim=mean(vortotalstim_ausgewaehlt);
mean_totalstim=mean(totalstim_ausgewaehlt);



figure
hold on
errorbar(x_werte(1:size(rawdata_measurement1,3)), mittelwert_kurve(1,:),mittelwert_kurve(2,:), 'Color', 'red')
line([x_werte(1,baseline_start) x_werte(1,rrpstim_stop)], [mean_baseline mean_baseline], 'Color', 'blue');
line([x_werte(1,baseline_start) x_werte(1,rrpstim_stop)], [mean_rrpstim mean_rrpstim], 'Color', 'blue');
line([(x_werte(1,rrpstim_stop)-x_werte(1,baseline_start))/2+x_werte(1,baseline_start) (x_werte(1,rrpstim_stop)-x_werte(1,baseline_start))/2+x_werte(1,baseline_start)], [mean_baseline mean_rrpstim], 'Color', 'blue');
line([x_werte(1,vorrecstim_start) x_werte(1,recstim_stop)], [mean_vorrecstim mean_vorrecstim], 'Color', 'green');
line([x_werte(1,vorrecstim_start) x_werte(1,recstim_stop)], [mean_recstim mean_recstim], 'Color', 'green');
line([(x_werte(1,recstim_stop)-x_werte(1,vorrecstim_start))/2+x_werte(1,vorrecstim_start) (x_werte(1,recstim_stop)-x_werte(1,vorrecstim_start))/2+x_werte(1,vorrecstim_start)], [mean_vorrecstim mean_recstim], 'Color', 'green');
line([x_werte(1,vortotalstim_start) x_werte(1,totalstim_stop)], [mean_vortotalstim mean_vortotalstim], 'Color', [255/255;153/255;51/255]);
line([x_werte(1,vortotalstim_start) x_werte(1,totalstim_stop)], [mean_totalstim mean_totalstim], 'Color', [255/255;153/255;51/255]);
line([(x_werte(1,totalstim_stop)-x_werte(1,vortotalstim_start))/2+x_werte(1,vortotalstim_start) (x_werte(1,totalstim_stop)-x_werte(1,vortotalstim_start))/2+x_werte(1,vortotalstim_start)], [mean_vortotalstim mean_totalstim], 'Color', [255/255;153/255;51/255]);



title('mean curve over all measurements incubated with Fluoxetine');
xlabel('frames')
ylabel('fluorescence emission [a.u.]')



figure
hold on
errorbar(x_werte(1:size(rawdata_measurement1,3)), mittelwert_kurve(1,:),mittelwert_kurve(2,:), 'Color', 'red')
for measurement=start_measurement:n_measurements
    switch measurement
        case{1,2,4,6,7,13,15,16,17,18,19,21,23,24,25,29,30,31,33,34}
            fluox_frames;
            if n_pics>=380
                x_werte3=x_werte2;
            else
                x_werte3=x_werte;
            end
            norm_kurve_ausgewaehlt(2,:,measurement)=x_werte3(1:310);
            plot(x_werte3(1:size(rawdata_measurement1,3)), norm_kurve(1,:,measurement), 'Color', [0.2 measurement/34 0.4])
    end
end

for i=1:10
    count=1;
    mean_curve(1,i)=i;
    for measurement=start_measurement:n_measurements
        switch measurement
            case{1,2,4,6,7,13,15,16,17,18,19,21,23,24,25,29,30,31,33,34}
                for spalten=1:310
                    if norm_kurve_ausgewaehlt(2,spalten,measurement)>=i
                        if norm_kurve_ausgewaehlt(2,spalten,measurement)<i+1
                            tomean(count)=norm_kurve_ausgewaehlt(1,spalten,measurement);
                            count=count+1;
                        end
                    end
                end
        end
        mean_curve(2,i)=mean(tomean);
        clearvars tomean
    end
end

plot(mean_curve(1,1:size(rawdata_measurement1,3)), mean_curve(2,1:size(rawdata_measurement1,3)))


% Save results

save(savename)


