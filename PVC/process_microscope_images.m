%% Read the images and process the particle sizes
clear all;
sets = 1:2;
plasticisers = ["BPA", "DEHP"];
scales = [500/426, 500/274]; % um/pixel for each set

all_data = [];
for set = sets
    scale = scales(set);
    for p_id = 1:numel(plasticisers)
        p = plasticisers(p_id);
        
        files = dir(fullfile("Microscope Images", "Set" + set, p, "*.mask.final.tif"))';
        for f = files
            % Load the mask 
            mask = imread(fullfile(f.folder, f.name));
            assert(max(mask(:) == 1), f.name + " is not a binary image");
            
            % Compute particle sizes
            stats = regionprops(mask, "Area", "Circularity", "Eccentricity", "MajorAxisLength", "EquivDiameter", "BoundingBox", "Centroid", "PixelIdxList", "MaxFeretProperties", "MinFeretProperties");
            N = numel(stats);
            stats = struct2table(stats);
            
            Area = stats.Area * scale^2;
            MajorAxisLength = stats.MajorAxisLength * scale;
            EquivDiameter = stats.EquivDiameter * scale;
            Circularity = stats.Circularity;
            Eccentricity = stats.Eccentricity;
            MaxFeretDiameter = stats.MaxFeretDiameter * scale;
            MinFeretDiameter = stats.MinFeretDiameter * scale;
            Plasticiser = repmat(p, [N 1]);
            Filename = repmat(string(f.name), [N 1]);
            Centroid = stats.Centroid;
            
            data = table(Plasticiser, Filename, Centroid, MaxFeretDiameter, MinFeretDiameter, Area, MajorAxisLength, EquivDiameter, Circularity  , Eccentricity);
                
            if isempty(all_data)
                all_data = data;
            else
                all_data = [all_data; data];
            end
        end
    end
end

%% Display statistics
for p = plasticisers
    this_data = all_data(all_data.Plasticiser == p, :);
    fprintf("%-20s  %15s %15s %15s\n", p, "Mean (um)", "Stdev (um)", "Median (um)");
    for statistics = ["MaxFeretDiameter", "MinFeretDiameter", "MajorAxisLength", "EquivDiameter"]
        d = this_data.(statistics);
        fprintf("%20s: %15.1f %15.1f %15.1f\n", statistics, mean(d), std(d), median(d));
    end
end

%% Figure S4 
[~,~,~] = mkdir("Figures");
figh = figure();
t = tiledlayout(2,2);

label = 'a';
for p = plasticisers
    ax = nexttile();
    max_feret_diameter = all_data.MaxFeretDiameter(all_data.Plasticiser == p);
    bin_edges = 0:50:550;
    histogram(max_feret_diameter, "BinEdges",bin_edges);
    xlabel("Maximum Feret diameter (um)");
    ylabel("Count");
    ax.XTick = 0:100:500;
    ax.TickDir = "out";
    title(sprintf("(%s) %s Feret diameter", label, p));
    
    label = char(label + 1);

    ax = nexttile();
    equiv_rad = all_data.EquivDiameter(all_data.Plasticiser == p) / 2;    
    bin_edges = 0:20:160;
    histogram(equiv_rad, "BinEdges", bin_edges);
    ax.XTick = 0:40:500;
    ax.TickDir = "out";    
    xlabel("Equivalent radius (um)");
    ylabel("Count");
    title(sprintf("(%s) %s equivalent radius", label, p));
    label = char(label + 1);

end

exportgraphics(figh, fullfile("Figures", "Figure_S4.png"), "Resolution", 600);


%% Make histograms
[~,~,~] = mkdir("Figures");
for p = plasticisers
    figure();
    max_feret_diameter = all_data.MaxFeretDiameter(all_data.Plasticiser == p);
    histogram(max_feret_diameter);
    xlabel("Maximum Feret diameter (um)");
    ylabel("Count");
    sgtitle(p);
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Max Feret Diameter.pdf"));
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Max Feret Diameter.png"), "Resolution", 600);

    figure();
    equiv_rad = all_data.EquivDiameter(all_data.Plasticiser == p) / 2;
    
    h=histogram(equiv_rad, "NumBins", 15);
    bins = h.BinEdges;
    counts = h.BinCounts;

    xlabel("Equivalent radius (um)");
    ylabel("Count");
    sgtitle(p);
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Equivalent Radius.pdf"));
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Equivalent Radius.png"), "Resolution", 600);

    % Calc contribution to leaching
    radii = (bins(1:end-1) + bins(2:end))/2;
    radii_edges = bins;
    volumes = (4/3)*pi*radii.^3;
    weighting = counts .* volumes;
    weighting = weighting / sum(weighting);
    figure();
    bar(radii, weighting);
    xlabel("Equiv radius (um)");
    ylabel("Relative mass of plasticiser");
    save(fullfile("Microscope Images", "size_contribution_"+p+".mat"), "radii", "weighting", "radii_edges");
    sgtitle(p);
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Contribution.pdf"));
%     exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Contribution.png"), "Resolution", 600);
end

%% Figure S7

figh = figure();
t = tiledlayout(2,1);
label = 'a';

for p = plasticisers
    load(fullfile("Microscope Images", "size_contribution_"+p+".mat"));

    ax = nexttile();
    histogram("BinEdges",radii_edges, "BinCounts", weighting);
    xlabel("Equivalent radius (um)");
    ylabel("Relative mass of plasticiser");
    title(sprintf("(%s) %s", label, p));
    label = char(label + 1);
end
linkaxes(t.Children, 'x');

%%
exportgraphics(figh, fullfile("Figures", "Figure_S7.png"), "Resolution", 600);


