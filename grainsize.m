function [grain_ID, grain_comp, grain_D, isize, w, bound] = grainsize(model, listDB, currentFolder, phase_name, phase_ID)
% Loads a phase map of either monophase or composite aggregates composed of
% N phases. Unknown phases are returned as not-indexed pixels.
%
% RETURNS:
%   grain_ID   = [1-by-Ngrains] cell array. Each cell contains the GLOBAL
%                pixel indices belonging to that grain.
%   grain_comp = [1-by-Ngrains] vector. Each element is the mineral index
%                (from listDB) of that grain's phase.
%   grain_D    = [1-by-Ngrains] vector of equivalent-circle diameters [mm].
%   isize      = [rows cols] size of the indexed image.
%   w          = scalar width of image [mm].
%   bound      = binary mask: 1 = grain interior, 0 = boundary.
%--------------------------------------------------------------------------

%% --- Number of phases ---
if model == 0
    N = 1;   % monophase
else
    N = numel(phase_name);   % composite
end

%% --- User instructions ---
msg1 = msgbox({[char(963), ' and ', char(941), ' colormaps require uploading']; ...
    'a color-coded phase map composed of at least 1 mineral species.'; ...
    'The selected phase(s) must be included in the mineral database.'; ...
    'Make phase boundaries BLACK to optimize grain-size analysis!'}, ...
    'Phase Map Setup');
uiwait(msg1);

%% --- Load image ---
[I, path] = uigetfile(fullfile(currentFolder, 'SiiEii_map', 'microstructural_map', '*.png'), ...
    'Load Phase Map');
if isequal(I, 0)
    error('grainsize:noFile', 'No file selected. Aborting.');
end

%% --- Physical width ---
options.Resize      = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
w = str2double(inputdlg({'Set width of phasemap [mm]:'}, ...
    ['Set phasemap ', I, ': '], 1, {'5'}, options));
if isempty(w) || isnan(w)
    error('grainsize:noWidth', 'No valid width entered. Aborting.');
end

%% --- Read image and build boundary mask from TRUE black pixels ---
% FIX 1: Do NOT use rgb2ind for boundary detection — it re-quantizes colors
%         and can reassign near-black pixels, corrupting the boundary mask.
%         Instead, detect boundaries directly from the original RGB image
%         by looking for pixels that are truly black (all channels < threshold).
im    = imread(fullfile(path, I));   % original RGB image
isize = size(im);
isize = isize(1:2);                  % [rows, cols]

BLACK_THRESH = 30;   % pixels with R,G,B all below this are treated as boundaries
isBlack = (im(:,:,1) < BLACK_THRESH) & ...
          (im(:,:,2) < BLACK_THRESH) & ...
          (im(:,:,3) < BLACK_THRESH);

bound = ~isBlack;   % 1 = interior, 0 = boundary

%% --- Connected-component labelling (find grains) ---
cc      = bwconncomp(bound, 4);
ngrains = cc.NumObjects;

%% --- Phase identification ---
if model == 0
    %% Monophase: no color selection needed
    phaseID = phase_ID;
    phind   = NaN;   % not used in monophase path

    f1 = figure('Name', 'Phase Map', 'NumberTitle', 'off');
    imshow(im);
    msg = msgbox(['MYflow detected ', num2str(ngrains), ' grains in phasemap ', I]);
    uiwait(msg);

else
    %% Composite: user clicks one grain per phase to sample its color
    % FIX 2: Keep rgb2ind ONLY for color-index sampling of user-clicked
    %         points (N+1 colors), NOT for boundary detection.
    imIndexed = double(rgb2ind(im, N + 1));   % used only for phind lookup

    phaseID = zeros(N, 1);
    phind   = nan(N, 1);

    msg2 = msgbox({['Set the composition of phase map ', I]; ...
        '1. Left-click at the center of a grain on the phase map.'; ...
        '2. Select the corresponding mineral from the database.'; ...
        '3. Repeat until all phases are classified.'; ...
        ['This map requires ', num2str(N), ' selection(s).']}, ...
        'Phase Map Setup');
    uiwait(msg2);

    f1 = figure('Name','Phase Map',...
            'NumberTitle','off',...
            'MenuBar','none',...
            'ToolBar','none');
imshow(im);
drawnow;
figure(f1);

    for j = 1:N
        while true   % loop until a valid (non-duplicate) phase is chosen
            title(sprintf('Click inside a grain of phase %d / %d', j, N));
            waitforbuttonpress;
            pt = get(gca,'CurrentPoint');
            x = pt(1,1);
            y = pt(1,2);

            % --- Get user click ---
            [x, y] = ginput(1);

            % FIX 3: Replace the original bisection-based coordinate lookup
            % with direct rounding and clamping of pixel coordinates.
            %         ginput returns floating-point figure coords; round them
            %         directly to get the nearest pixel.
            xi = min(max(round(x), 1), isize(2));
            yi = min(max(round(y), 1), isize(1));

            % Warn if the user clicked on a boundary pixel
            if isBlack(yi, xi)
                msg = msgbox('WARNING: You clicked on a boundary pixel. Please click inside a grain.');
                uiwait(msg);
                continue;
            end

            % --- Select mineral from database ---
            [local_ind, ok] = listdlg('ListString', phase_name, ...
                'Name', 'Set mineral', 'SelectionMode', 'single');
            if ~ok
                msg = msgbox('No mineral selected. Please try again.');
                uiwait(msg);
                continue;
            end

            ind        = 1:numel(listDB);
            global_ind = ind(strcmp(phase_name{local_ind}, listDB));

            % FIX 4: Simplify duplicate-phase detection using a direct
            % comparison against previously assigned phases.
            %         Original condition "sum(phaseID~=global_ind)==N" is always
            %         true on first iteration and does not correctly detect
            %         duplicates. Use "any(phaseID==global_ind)" instead.
            if any(phaseID == global_ind)
                msg = msgbox({['WARNING: "', listDB{global_ind}, '" is already assigned.']; ...
                    'Please choose a different phase.'});
                uiwait(msg);
                continue;
            end

            % Record phase ID and the color-index at the clicked pixel
            phaseID(j) = global_ind;
            phind(j)   = imIndexed(yi, xi);
            break;
        end
    end
end

%% --- Validate all phases were assigned ---
if any(phaseID == 0)
    msgbox('WARNING: Not all phases were assigned. Aborting.');
    grain_ID   = [];
    grain_comp = [];
    grain_D    = [];
    return;
end

close(f1);

%% --- Map color indices back to mineral database IDs ---
% FIX 5: Work on a copy of imIndexed (composite) or a simple label array
%         (monophase) so the original RGB image is untouched.
if model == 0
    % Monophase: every interior pixel belongs to the single phase
    imLabelled = ones(isize) .* phaseID;
else
    imLabelled = -double(rgb2ind(im, N + 1));   % negative to avoid clash with positive DB IDs
    for i = 1:N
        imLabelled(imLabelled == -phind(i)) = phaseID(i);
    end
    % Any pixel not matched to a known phase → 0 (unclassified)
    imLabelled(imLabelled < 0) = 0;
end

%% --- Assign composition to each grain ---
grain_comp = zeros(1, ngrains);
for i = 1:ngrains
    grain_comp(i) = mode(imLabelled(cc.PixelIdxList{i}));
end

% FIX 6: Exclude components whose dominant label is 0
% from grain statistics. (unclassified / boundary
%         fragments) — these are the spurious extra grains in composite mode.
validGrain = grain_comp > 0;
grain_ID   = cc.PixelIdxList(validGrain);
grain_comp = grain_comp(validGrain);
ngrains    = sum(validGrain);

%% --- Grain-size computation ---
% Pixel area in mm^2 (assumes square pixels)
pixres = (w / isize(2))^2;

graindata = regionprops(cc, 'Area');
allAreas  = [graindata.Area];
grainArea = allAreas(validGrain) .* pixres;   % keep only valid grains
grain_D   = 2 * sqrt(grainArea ./ pi);        % equivalent-circle diameter [mm]

%% --- Grain-size statistics plots ---
f2 = figure('Name', 'Grain size statistics', 'NumberTitle', 'off');
colorpalette = {'k', 'b', 'm', 'g', 'r', 'c', 'y'};

phase_area  = zeros(1, N);
valid_phase = 0;
for i = 1:N
    if ~strcmp(listDB(phaseID(i)), 'not defined')
        valid_phase   = valid_phase + 1;
        phase_area(i) = sum(grainArea(grain_comp == phaseID(i)));
    end
end
phase_total = sum(phase_area);

subplot_idx = 0;
for i = 1:N
    if ~strcmp(listDB(phaseID(i)), 'not defined')
        subplot_idx    = subplot_idx + 1;
        phase_diameter = grain_D(grain_comp == phaseID(i));
        area_fraction  = phase_area(i) / phase_total;

        subplot(ceil(valid_phase / 3), 3, subplot_idx);

        % FIX 7: Ensure nbins >= 1 to avoid histogram crash for small datasets.
        nbins = max(1, round(numel(phase_diameter) / 5));

        % FIX 8: Plot raw diameters (not divided by area_fraction) —
        %         dividing distorted the x-axis scale.
        h1 = histogram(phase_diameter, nbins, 'Normalization', 'Probability');
        h1.FaceColor = colorpalette{i};

        title([listDB(phaseID(i)), ', N = ', num2str(numel(phase_diameter)), ...
            sprintf(' (\\phi = %.2f)', area_fraction)], 'FontSize', 7);
        xlabel('Diameter [mm]');
        ylabel('Probability');
        hold on;

        xav = median(phase_diameter);
        yav = max(h1.Values);
        text(xav, yav * 0.95, ...
            ['d_{50} = ', num2str(median(phase_diameter) * 1000, '%.1f'), ' \mum'], ...
            'BackgroundColor', 'w', 'Color', colorpalette{i}, 'FontSize', 7, ...
            'VerticalAlignment', 'top');
    end
end

msg = msgbox(['Phase Map successfully created! ', num2str(ngrains), ' valid grains detected.']);
uiwait(msg);
if isvalid(f2)
    close(f2);
end

end