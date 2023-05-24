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

%% Make histograms
[~,~,~] = mkdir("Figures");
for p = plasticisers
    figure();
    max_feret_diameter = all_data.MaxFeretDiameter(all_data.Plasticiser == p);
    histogram(max_feret_diameter);
    xlabel("Maximum Feret diameter (um)");
    ylabel("Count");
    sgtitle(p);
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Max Feret Diameter.pdf"));
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Max Feret Diameter.png"), "Resolution", 600);

    figure();
    equiv_rad = all_data.EquivDiameter(all_data.Plasticiser == p) / 2;
    
    h=histogram(equiv_rad, "NumBins", 15);
    bins = h.BinEdges;
    counts = h.BinCounts;

    xlabel("Equivalent radius (um)");
    ylabel("Count");
    sgtitle(p);
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Equivalent Radius.pdf"));
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Equivalent Radius.png"), "Resolution", 600);

    % Calc contribution to leaching
    radii = (bins(1:end-1) + bins(2:end))/2;
    volumes = (4/3)*pi*radii.^3;
    weighting = counts .* volumes;
    weighting = weighting / sum(weighting);
    figure();
    bar(radii, weighting);
    xlabel("Equiv radius (um)");
    ylabel("Relative mass of plasticiser");
    save(fullfile("Microscope Images", "size_contribution_"+p+".mat"), "radii", "weighting");
    sgtitle(p);
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Contribution.pdf"));
    exportgraphics(gcf(), fullfile("Figures", "Size distribution " + p + " Contribution.png"), "Resolution", 600);
end
