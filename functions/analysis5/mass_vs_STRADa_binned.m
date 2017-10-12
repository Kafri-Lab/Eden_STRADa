 function fun(ResultTable, Data)
  channels = Data.O.Channels;
  cell_stages = {'EG1', 'LG1', 'S', 'G2'};
  STRADa_channel = find(ismember(channels, 'STRADa'));
  cell_mass_channel = find(ismember(channels, 'SE'));
  DAPI_channel = find(ismember(channels, 'DAPI'));
  geminin_channel = find(ismember(channels, 'Geminin'));
  EG1_color = [0/255 255/255 150/255];
  LG1_color = [232/255 227/255 12/255];
  S_color = [255/255 102/255 0/255];
  G2_color = [207/255 12/255 232/255];
  DAPI = ResultTable.NInt(:,DAPI_channel);
  geminin = ResultTable.NInt(:,geminin_channel);
  cell_mass = ResultTable.CInt(:,cell_mass_channel);
  nuc_mass = ResultTable.NInt(:,cell_mass_channel);
  cyto_mass = cell_mass - nuc_mass;
  STRADa_nuc = ResultTable.NInt(:,STRADa_channel);
  STRADa_total = ResultTable.CInt(:,STRADa_channel);
  STRADa_cyto = STRADa_total - STRADa_nuc;
  STRADa_nuc = STRADa_nuc./nuc_mass;
  STRADa_total = STRADa_total./cell_mass;
  STRADa_cyto = STRADa_cyto./cyto_mass;
  STRADa_localization_ratio_nuc_cyto = STRADa_nuc./STRADa_cyto;
  STRADa_localization_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_localization_ratio_cyto_total = STRADa_cyto./STRADa_total;

  % % Eliminate outliers
  lower_prctile = 0.1;
  higher_prctile = 99.9;
  cell_mass = set_outliers_to(cell_mass, lower_prctile, higher_prctile, NaN);
  nuc_mass = set_outliers_to(nuc_mass, lower_prctile, higher_prctile, NaN);
  cyto_mass = set_outliers_to(cyto_mass, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_nuc = set_outliers_to(STRADa_nuc, lower_prctile, higher_prctile, NaN);
  STRADa_total = set_outliers_to(STRADa_total, lower_prctile, higher_prctile, NaN);
  STRADa_cyto = set_outliers_to(STRADa_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_cyto = set_outliers_to(STRADa_localization_ratio_nuc_cyto, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_nuc_total = set_outliers_to(STRADa_localization_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_localization_ratio_cyto_total = set_outliers_to(STRADa_localization_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  

  % Find cell stages
  [idx_EG1,idx_LG1,idx_G1S,idx_S,idx_G2] = FindStages_v2(DAPI,log(geminin));

  % combine LG1 and G1S into one group (LG1) because G1S wasn't detected well enough
  idx_LG1 = idx_LG1 | idx_G1S;
  idx_G1S = []; % delete so no mistake is made


  % Get stats for each cell stage
  cell_mass_EG1 = cell_mass(idx_EG1);
  cell_mass_LG1 = cell_mass(idx_LG1);
  cell_mass_S = cell_mass(idx_S);
  cell_mass_G2 = cell_mass(idx_G2);
  STRADa_total_EG1 = STRADa_total(idx_EG1);
  STRADa_total_LG1 = STRADa_total(idx_LG1);
  STRADa_total_S = STRADa_total(idx_S);
  STRADa_total_G2 = STRADa_total(idx_G2);
  STRADa_nuc_EG1 = STRADa_nuc(idx_EG1);
  STRADa_nuc_LG1 = STRADa_nuc(idx_LG1);
  STRADa_nuc_S = STRADa_nuc(idx_S);
  STRADa_nuc_G2 = STRADa_nuc(idx_G2);
  STRADa_cyto_EG1 = STRADa_cyto(idx_EG1);
  STRADa_cyto_LG1 = STRADa_cyto(idx_LG1);
  STRADa_cyto_S = STRADa_cyto(idx_S);
  STRADa_cyto_G2 = STRADa_cyto(idx_G2);
  STRADa_ratio_EG1 = STRADa_localization_ratio_nuc_cyto(idx_EG1);
  STRADa_ratio_LG1 = STRADa_localization_ratio_nuc_cyto(idx_LG1);
  STRADa_ratio_S = STRADa_localization_ratio_nuc_cyto(idx_S);
  STRADa_ratio_G2 = STRADa_localization_ratio_nuc_cyto(idx_G2);


  STRADa_nuc_EG1_mean = nanmean(STRADa_nuc_EG1);
  STRADa_nuc_LG1_mean = nanmean(STRADa_nuc_LG1);
  STRADa_nuc_S_mean = nanmean(STRADa_nuc_S);
  STRADa_nuc_G2_mean = nanmean(STRADa_nuc_G2);
  STRADa_ratio_EG1_mean = nanmean(STRADa_ratio_EG1);
  STRADa_ratio_LG1_mean = nanmean(STRADa_ratio_LG1);
  STRADa_ratio_S_mean = nanmean(STRADa_ratio_S);
  STRADa_ratio_G2_mean = nanmean(STRADa_ratio_G2);
  cell_mass_EG1_mean = nanmean(cell_mass_EG1);
  cell_mass_LG1_mean = nanmean(cell_mass_LG1);
  cell_mass_S_mean = nanmean(cell_mass_S);
  cell_mass_G2_mean = nanmean(cell_mass_G2);










  num_bins = 129;















  %% SE vs STRADa Ratio
  num_cell_stages = 4;
  EG1_means_per_bin = nan(num_bins,1);
  LG1_means_per_bin = nan(num_bins,1);
  S_means_per_bin = nan(num_bins,1);
  G2_means_per_bin = nan(num_bins,1);
  EG1_errors_per_bin = nan(num_bins,1);
  LG1_errors_per_bin = nan(num_bins,1);
  S_errors_per_bin = nan(num_bins,1);
  G2_errors_per_bin = nan(num_bins,1);
  bin_centers = nan(num_bins,1);
  % Error = std / sqrt(length(data_points))
  num_columns = num_bins*num_cell_stages;
  data_for_bars = nan(length(cell_mass),num_columns); % data for bar chart is represented as columns for each bar and rows for each datum. Since some bars will have more data than others, some columns will have a lot of NaNs (to keep the matrix shape rectangular)
  for bin_num=1:num_bins
    % Calc bin spacing
    bin_spacing = (max(cell_mass) - min(cell_mass))/num_bins;
    bin_start = min(cell_mass) + (bin_spacing * (bin_num-1));
    bin_end = min(cell_mass) + (bin_spacing * bin_num);
    bin_centers(bin_num) = bin_end - bin_spacing/2;
    % Find cells in bin by cell stage
    in_bin = bin_start < cell_mass & cell_mass < bin_end;
    EG1_in_bin = idx_EG1 & in_bin;
    LG1_in_bin = idx_LG1 & in_bin;
    S_in_bin = idx_S & in_bin;
    G2_in_bin = idx_G2 & in_bin;
    % Get cell mass
    value_EG1_in_bin = STRADa_localization_ratio_nuc_cyto(EG1_in_bin);
    value_LG1_in_bin = STRADa_localization_ratio_nuc_cyto(LG1_in_bin);
    value_S_in_bin = STRADa_localization_ratio_nuc_cyto(S_in_bin);
    value_G2_in_bin = STRADa_localization_ratio_nuc_cyto(G2_in_bin);
    % Store data points appropriately
    % Four columns will be filled with whatever amount of data is present for the cell stage and this bin
    data_for_bars(1:length(value_EG1_in_bin),((bin_num-1)*num_cell_stages)+1) = value_EG1_in_bin;
    data_for_bars(1:length(value_LG1_in_bin),((bin_num-1)*num_cell_stages)+2) = value_LG1_in_bin;
    data_for_bars(1:length(value_S_in_bin),((bin_num-1)*num_cell_stages)+3) = value_S_in_bin;
    data_for_bars(1:length(value_G2_in_bin),((bin_num-1)*num_cell_stages)+4) = value_G2_in_bin;
    % Store means
    EG1_means_per_bin(bin_num) = nanmean(value_EG1_in_bin);
    LG1_means_per_bin(bin_num) = nanmean(value_LG1_in_bin);
    S_means_per_bin(bin_num) = nanmean(value_S_in_bin);
    G2_means_per_bin(bin_num) = nanmean(value_G2_in_bin);
    % Store errors
    EG1_errors_per_bin(bin_num) = nanstd(value_EG1_in_bin)/sqrt(length(value_EG1_in_bin));
    LG1_errors_per_bin(bin_num) = nanstd(value_LG1_in_bin)/sqrt(length(value_LG1_in_bin));
    S_errors_per_bin(bin_num) = nanstd(value_S_in_bin)/sqrt(length(value_S_in_bin));
    G2_errors_per_bin(bin_num) = nanstd(value_G2_in_bin)/sqrt(length(value_G2_in_bin));
    %%% Plot to debug bins
    % x = [bin_start bin_start];
    % y = [0 1];
    % line(x,y,'Color','blue');
    % % pause
    % x = [bin_end bin_end];
    % y = [0 1];
    % line(x,y,'Color','red');
    % % pause
  end

  num_plot_columns = num_bins*(num_cell_stages+1);
  spacing = 1:num_plot_columns;
  for num=4:5:num_columns*(num_cell_stages+1)
     spacing(spacing==num+1) = [];
  end
  figure('Position',[2          42         958        1074]);
  % figure('Position',[1          41        1920        1083]);
  H=notBoxPlot(data_for_bars,spacing,'jitter',0.6);
  d=[H.data];
  hold on

  % Disable the displaying the data points
  set(d,'Visible','off')

  % Box Colors
  EG1_IND=zeros(1,num_columns);
  EG1_IND(1:4:num_columns)=1;
  LG1_IND=zeros(1,num_columns);
  LG1_IND(2:4:num_columns)=1;
  S_IND=zeros(1,num_columns);
  S_IND(3:4:num_columns)=1;
  G2_IND=zeros(1,num_columns);
  G2_IND(4:4:num_columns)=1;

  set([H(find(EG1_IND)).data],'MarkerSize',4,...
      'markerFaceColor',EG1_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(EG1_IND)).semPtch],...
      'FaceColor',EG1_color*0.25,...
      'EdgeColor','none')
  set([H(find(EG1_IND)).sdPtch],...
      'FaceColor',EG1_color*0.75,...
      'EdgeColor','none')
  set([H(find(EG1_IND)).mu],...
      'Color','b')

  set([H(find(LG1_IND)).data],'MarkerSize',4,...
      'markerFaceColor',LG1_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(LG1_IND)).semPtch],...
      'FaceColor',LG1_color*0.25,...
      'EdgeColor','none')
  set([H(find(LG1_IND)).sdPtch],...
      'FaceColor',LG1_color*0.75,...
      'EdgeColor','none')
  set([H(find(LG1_IND)).mu],...
      'Color','b')

  set([H(find(S_IND)).data],'MarkerSize',4,...
      'markerFaceColor',S_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(S_IND)).semPtch],...
      'FaceColor',S_color*0.25,...
      'EdgeColor','none')
  set([H(find(S_IND)).sdPtch],...
      'FaceColor',S_color*0.75,...
      'EdgeColor','none')
  set([H(find(S_IND)).mu],...
      'Color','b')

  set([H(find(G2_IND)).data],'MarkerSize',4,...
      'markerFaceColor',G2_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(G2_IND)).semPtch],...
      'FaceColor',G2_color*0.25,...
      'EdgeColor','none')
  set([H(find(G2_IND)).sdPtch],...
      'FaceColor',G2_color*0.75,...
      'EdgeColor','none')
  set([H(find(G2_IND)).mu],...
      'Color','b')

  % Trend lines
  bin_centers_on_plot  = (2:num_cell_stages+1:num_plot_columns)+.5;
  plot(bin_centers_on_plot,EG1_means_per_bin,'color',EG1_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,LG1_means_per_bin,'color',LG1_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,S_means_per_bin,'color',S_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,G2_means_per_bin,'color',G2_color*.9,'LineWidth', 2);

  % Mean line color
  set([H.mu],'color',[0.5 .5 .5])

  % X ticks
  xticklabels = {};
  for bin_num=1:num_bins
    % xticklabel = sprintf('%d (%0.2e)',bin_num, bin_centers(bin_num));
    xticklabel = sprintf('%d',bin_num);
    xticklabels{bin_num} = xticklabel;
  end
  set(gca,'XTick',bin_centers_on_plot) % The middle point of each bin
  set(gca,'XTickLabels',xticklabels)


  % Pretty
  box on
  title('Cell Mass (SE) vs STRADa localization ratio (STRADa Normalized by SE)')
  xlabel('Cell Mass (SE)');
  ylabel('<- More STRADa in Cytoplasm                                                Ratio                                                More STRADa in Nucleus ->');

  filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/binned_%s_mass-vs-ratio.png', num2str(num_bins));
  export_fig(filename, '-m2 -transparent -nocrop')



















  %% SE vs NUC STRADa
  num_cell_stages = 4;
  EG1_means_per_bin = nan(num_bins,1);
  LG1_means_per_bin = nan(num_bins,1);
  S_means_per_bin = nan(num_bins,1);
  G2_means_per_bin = nan(num_bins,1);
  EG1_errors_per_bin = nan(num_bins,1);
  LG1_errors_per_bin = nan(num_bins,1);
  S_errors_per_bin = nan(num_bins,1);
  G2_errors_per_bin = nan(num_bins,1);
  bin_centers = nan(num_bins,1);
  % Error = std / sqrt(length(data_points))
  num_columns = num_bins*num_cell_stages;
  data_for_bars = nan(length(cell_mass),num_columns); % data for bar chart is represented as columns for each bar and rows for each datum. Since some bars will have more data than others, some columns will have a lot of NaNs (to keep the matrix shape rectangular)
  for bin_num=1:num_bins
    % Calc bin spacing
    bin_spacing = (max(cell_mass) - min(cell_mass))/num_bins;
    bin_start = min(cell_mass) + (bin_spacing * (bin_num-1));
    bin_end = min(cell_mass) + (bin_spacing * bin_num);
    bin_centers(bin_num) = bin_end - bin_spacing/2;
    % Find cells in bin by cell stage
    in_bin = bin_start < cell_mass & cell_mass < bin_end;
    EG1_in_bin = idx_EG1 & in_bin;
    LG1_in_bin = idx_LG1 & in_bin;
    S_in_bin = idx_S & in_bin;
    G2_in_bin = idx_G2 & in_bin;
    % Get cell mass
    value_EG1_in_bin = STRADa_nuc(EG1_in_bin);
    value_LG1_in_bin = STRADa_nuc(LG1_in_bin);
    value_S_in_bin = STRADa_nuc(S_in_bin);
    value_G2_in_bin = STRADa_nuc(G2_in_bin);
    % Store data points appropriately
    % Four columns will be filled with whatever amount of data is present for the cell stage and this bin
    data_for_bars(1:length(value_EG1_in_bin),((bin_num-1)*num_cell_stages)+1) = value_EG1_in_bin;
    data_for_bars(1:length(value_LG1_in_bin),((bin_num-1)*num_cell_stages)+2) = value_LG1_in_bin;
    data_for_bars(1:length(value_S_in_bin),((bin_num-1)*num_cell_stages)+3) = value_S_in_bin;
    data_for_bars(1:length(value_G2_in_bin),((bin_num-1)*num_cell_stages)+4) = value_G2_in_bin;
    % Store means
    EG1_means_per_bin(bin_num) = nanmean(value_EG1_in_bin);
    LG1_means_per_bin(bin_num) = nanmean(value_LG1_in_bin);
    S_means_per_bin(bin_num) = nanmean(value_S_in_bin);
    G2_means_per_bin(bin_num) = nanmean(value_G2_in_bin);
    % Store errors
    EG1_errors_per_bin(bin_num) = nanstd(value_EG1_in_bin)/sqrt(length(value_EG1_in_bin));
    LG1_errors_per_bin(bin_num) = nanstd(value_LG1_in_bin)/sqrt(length(value_LG1_in_bin));
    S_errors_per_bin(bin_num) = nanstd(value_S_in_bin)/sqrt(length(value_S_in_bin));
    G2_errors_per_bin(bin_num) = nanstd(value_G2_in_bin)/sqrt(length(value_G2_in_bin));
    %%% Plot to debug bins
    % x = [bin_start bin_start];
    % y = [0 1];
    % line(x,y,'Color','blue');
    % % pause
    % x = [bin_end bin_end];
    % y = [0 1];
    % line(x,y,'Color','red');
    % % pause
  end

  spacing = 1:num_plot_columns;
  for num=4:5:num_columns*(num_cell_stages+1)
     spacing(spacing==num+1) = [];
  end
  figure('Position',[2          42         958        1074]);
  % figure('Position',[1          41        1920        1083]);

  H=notBoxPlot(data_for_bars,spacing,'jitter',0.6);
  d=[H.data];
  hold on

  % Disable the displaying the data points
  set(d,'Visible','off')

  % Box Colors
  EG1_IND=zeros(1,num_columns);
  EG1_IND(1:4:num_columns)=1;
  LG1_IND=zeros(1,num_columns);
  LG1_IND(2:4:num_columns)=1;
  S_IND=zeros(1,num_columns);
  S_IND(3:4:num_columns)=1;
  G2_IND=zeros(1,num_columns);
  G2_IND(4:4:num_columns)=1;

  set([H(find(EG1_IND)).data],'MarkerSize',4,...
      'markerFaceColor',EG1_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(EG1_IND)).semPtch],...
      'FaceColor',EG1_color*0.25,...
      'EdgeColor','none')
  set([H(find(EG1_IND)).sdPtch],...
      'FaceColor',EG1_color*0.75,...
      'EdgeColor','none')
  set([H(find(EG1_IND)).mu],...
      'Color','b')

  set([H(find(LG1_IND)).data],'MarkerSize',4,...
      'markerFaceColor',LG1_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(LG1_IND)).semPtch],...
      'FaceColor',LG1_color*0.25,...
      'EdgeColor','none')
  set([H(find(LG1_IND)).sdPtch],...
      'FaceColor',LG1_color*0.75,...
      'EdgeColor','none')
  set([H(find(LG1_IND)).mu],...
      'Color','b')

  set([H(find(S_IND)).data],'MarkerSize',4,...
      'markerFaceColor',S_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(S_IND)).semPtch],...
      'FaceColor',S_color*0.25,...
      'EdgeColor','none')
  set([H(find(S_IND)).sdPtch],...
      'FaceColor',S_color*0.75,...
      'EdgeColor','none')
  set([H(find(S_IND)).mu],...
      'Color','b')

  set([H(find(G2_IND)).data],'MarkerSize',4,...
      'markerFaceColor',G2_color*0.5,...
      'markerEdgeColor', 'none')
  set([H(find(G2_IND)).semPtch],...
      'FaceColor',G2_color*0.25,...
      'EdgeColor','none')
  set([H(find(G2_IND)).sdPtch],...
      'FaceColor',G2_color*0.75,...
      'EdgeColor','none')
  set([H(find(G2_IND)).mu],...
      'Color','b')

  % Trend lines
  bin_centers_on_plot  = (2:num_cell_stages+1:num_plot_columns)+.5;
  plot(bin_centers_on_plot,EG1_means_per_bin,'color',EG1_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,LG1_means_per_bin,'color',LG1_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,S_means_per_bin,'color',S_color*.9,'LineWidth', 2);
  plot(bin_centers_on_plot,G2_means_per_bin,'color',G2_color*.9,'LineWidth', 2);

  % Mean line color
  set([H.mu],'color',[0.5 .5 .5])

  % X ticks
  xticklabels = {};
  for bin_num=1:num_bins
    % xticklabel = sprintf('%d (%0.2e)',bin_num, bin_centers(bin_num));
    xticklabel = sprintf('%d',bin_num);
    xticklabels{bin_num} = xticklabel;
  end
  set(gca,'XTick',(2:num_cell_stages+1:num_plot_columns)+.5) % The middle point of each bin
  set(gca,'XTickLabels',xticklabels)

  % Pretty
  box on
  title('Cell Mass (SE) vs Nuclear STRADa (Normalized by SE)')
  xlabel('Cell Mass (SE)');
  ylabel('<- Less STRADa in Nucleus                                                Ratio                                                More STRADa in Nucleus ->');

  filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 5/binned_%s_mass-vs-nuc-STRADa.png', num2str(num_bins));
  export_fig(filename, '-m2 -transparent -nocrop')













  
  
end






