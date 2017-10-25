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
  STRADa_ratio = STRADa_nuc./STRADa_cyto;
  STRADa_ratio_nuc_total = STRADa_nuc./STRADa_total;
  STRADa_ratio_cyto_total = STRADa_cyto./STRADa_total;

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
  STRADa_ratio = set_outliers_to(STRADa_ratio, lower_prctile, higher_prctile, NaN);
  STRADa_ratio_nuc_total = set_outliers_to(STRADa_ratio_nuc_total, lower_prctile, higher_prctile, NaN);
  STRADa_ratio_cyto_total = set_outliers_to(STRADa_ratio_cyto_total, lower_prctile, higher_prctile, NaN);

  

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
  STRADa_ratio_EG1 = STRADa_ratio(idx_EG1);
  STRADa_ratio_LG1 = STRADa_ratio(idx_LG1);
  STRADa_ratio_S = STRADa_ratio(idx_S);
  STRADa_ratio_G2 = STRADa_ratio(idx_G2);


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


  figure('Position',[251,    91,   883,   708]);
  random_cells = round(rand(1,10000)*height(ResultTable));
  for cell_id=random_cells
    cell_stage = cell_stages(find([idx_EG1(cell_id) idx_LG1(cell_id) idx_S(cell_id) idx_G2(cell_id)]))
    if isempty(cell_stage)
        continue
    end

    clf
    subplot(3,3,5)
    hold on
    % Scatter plot by cell cycle
    X=cell_mass_EG1(1:50:end);
    Y=STRADa_ratio_EG1(1:50:end);
    h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 .6 0]);
    X=cell_mass_LG1(1:50:end);
    Y=STRADa_ratio_LG1(1:50:end);
    h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 1 .5]);
    X=cell_mass_S(1:50:end);
    Y=STRADa_ratio_S(1:50:end);
    h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.2 1 1]);
    X=cell_mass_G2(1:50:end);
    Y=STRADa_ratio_G2(1:50:end);
    h=scatter(X,Y,5,'filled', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.4 .6 1]);
    % Scatter plot by cell cycle MEAN
    h=scatter(cell_mass_EG1_mean,STRADa_ratio_EG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 .6 0]);
    h=scatter(cell_mass_LG1_mean,STRADa_ratio_LG1_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[0 1 .5]);
    h=scatter(cell_mass_S_mean,STRADa_ratio_S_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.2 1 1]);
    h=scatter(cell_mass_G2_mean,STRADa_ratio_G2_mean,150,'o','LineWidth',2,'MarkerEdgeColor','k', 'MarkerFaceAlpha',1,'MarkerFaceColor',[.4 .6 1]);
    % Plot an X for this cell
    X = cell_mass(cell_id);
    Y = STRADa_ratio(cell_id);
    text(X,Y,'x','Color','r')
    % Make plot pretty
    legend({'EG1', 'LG1', 'S', 'G2','EG1 Mean', 'LG1 Mean', 'S Mean', 'G2 Mean'});
    set(gca,'ydir','normal')
    title('STRADa Localization Ratio (average intensity) vs Cell Size (SE)')
    xlabel('Cell Size (SE)');ylabel('<- More Cyto      More Nuc ->');
    legendmarkeradjust(90);
    grid on
    axis tight

    centroid_x = ResultTable.Centroid(cell_id,1);
    centroid_y = ResultTable.Centroid(cell_id,2);

    %% Find and segment single image.
    ImageIDs=Data.O.ImageIDs;
    t1=datetime;
    t2=datetime;

    % Get image this cell is in
    I= ImageIDs.Row==ResultTable.Row(cell_id) ...
        & ImageIDs.Column==ResultTable.Column(cell_id) ...
        & ImageIDs.Field==ResultTable.Field(cell_id) ...
        & ImageIDs.Channel==1;

    addpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
    [iterTable,ImageID,O]=O_SegmentCells_v6_SingleImage(Data.O,ImageIDs(I,:),t1,t2);
    rmpath('Z:\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion')
    
    nuc_perim = bwperim(O.BW{DAPI_channel});
    cyto_perim = bwperim(O.BW{cell_mass_channel});
    
    cyto_mask = O.BW{cell_mass_channel};
    labelled_id = cyto_mask(round(centroid_y), round(centroid_x));
    cyto_mask = cyto_mask==labelled_id;

    %% Display each channel with an x to mark the cell and nuc and cyto boundries
    set(0,'DefaultFigureWindowStyle','docked');
    warning off images:initSize:adjustingMag;
    warning off images:imshow:magnificationMustBeFitForDockedFigure;

    for cid=1:length(channels)
      subplot(3,3,cid)
      % display this channel
      img=O.IM{cid};
      imshow(img, [prctile(img(:), 3), prctile(img(:), 99.9)]);
      title(channels(cid))
      hold on
      % nuc boundries
      [Y,X] = find(nuc_perim);
      h=scatter(X,Y,4,'filled','MarkerFaceAlpha',3/4);
      % cyto boundries
      [Y,X] = find(cyto_perim);
      h=scatter(X,Y,4,'filled','MarkerFaceAlpha',2/3);
      % mark the x

      text(centroid_x-1.00,centroid_y-1.00,'x','Color','r')

      xlim([centroid_x - 120 centroid_x + 120])
      ylim([centroid_y - 120 centroid_y + 120])

    end

    % Image of STRADa normalized by SE
    subplot(3,3,7)
    STRADa_img=O.IM{STRADa_channel};
    STRADa_img(find(~cyto_mask))=0;
    cell_mass_img=O.IM{cell_mass_channel};
    cell_mass_img(~cyto_mask)=0;
    img = STRADa_img ./ cell_mass_img;
    imshow(img, [prctile(img(:), 3), prctile(img(:), 97)]);
    title('STRADa Normalized by SE (STRADa./SE)')
    hold on
    % nuc boundries
    [Y,X] = find(nuc_perim);
    h=scatter(X,Y,4,'filled','MarkerFaceAlpha',3/4);
    % cyto boundries
    [Y,X] = find(cyto_perim);
    h=scatter(X,Y,4,'filled','MarkerFaceAlpha',2/3);
    xlim([centroid_x - 120 centroid_x + 120])
    ylim([centroid_y - 120 centroid_y + 120])

    %% Extra details
    subplot(3,3,8)
    hold on
    text(0,.5,['Cell Stage: ' char(cell_stage)]); 
    text(0,.4,['Cell Size Percentile (Overall): ' num2str(invprctile(cell_mass,cell_mass(cell_id)))]);
    cell_stage_cell_mass = eval(char(strcat('cell_mass_', cell_stage))); % within cell stage only
    text(0,.3,['Cell Size Percentile (within stage): ' num2str(invprctile(cell_stage_cell_mass,cell_mass(cell_id)))]);
    text(0,.2,['Average STRADa Percentile (Overall): ' num2str(invprctile(STRADa_total,STRADa_total(cell_id)))]);
    cell_stage_STRADa_total = eval(char(strcat('STRADa_total_', cell_stage))); % within cell stage only
    text(0,.1,['Average STRADa Percentile (within stage): ' num2str(invprctile(cell_stage_STRADa_total,STRADa_total(cell_id)))]);
    axis off

    %% Image of line scan from centroid to furthest point on cell boundry
    % get X any Y of each point on the perim
    [perim_y,perim_x]=find(bwperim(cyto_mask));
    % align the point so that it is relative to the nucleus's center, ie. it's origin (0,0) is the nucleus's center
    perim_x = perim_x - centroid_x;
    perim_y = perim_y - centroid_y;
    % get euclidean distance from manhatten distance
    [theta,rho] = cart2pol(perim_x,perim_y);
    [max_dist, max_dist_id] = max(rho);
    max_dist_x = perim_x(max_dist_id) + centroid_x;
    max_dist_y = perim_y(max_dist_id) + centroid_y;
    subplot(3,3,3)
    X = [centroid_x max_dist_x];
    Y = [centroid_y max_dist_y];
    line(X,Y,'Color','r')


    %% Plot line scan
    subplot(3,3,6)
    [CX,CY,C,xi,yi]=improfile(O.IM{cell_mass_channel},X,Y); 
    hold on
    plot(1:length(C),C,'-'); 
    grid on
    axis tight
    xlabel('lateral distance [px]'); 
    ylabel('intensity [a.u.]'); 
    title('STRADa linescan');

    % %% Plot line scan
    % subplot(3,3,9)
    % thickness=3;
    % [CX,CY,C_sum,C,xi,yi]=improfile_integrated(O.IM{cell_mass_channel},thickness,X,Y); 
    % if C_sum(1)<C_sum(end)
    %   C_sum=fliplr(C_sum);
    % end
    % hold on
    % plot(1:length(C),C,'.'); 
    % plot(1:length(C_sum),C_sum/(thickness+1),'-'); 
    % grid on
    % axis tight
    % xlabel('lateral distance [px]'); 
    % ylabel('intensity [a.u.]'); 
    % title('comparision of single linescan and integrated line scan');

    filename = sprintf('Z:/DanielS/Projects/Eden_STRADa/Plots/Analysis 3 by mass/%s_stage_%s_mass_id_%s.png',char(cell_stage),num2str(cell_mass(cell_id)),num2str(cell_id));
    export_fig(filename, '-m1')
  end
end