function loadAllFile(ssdot, root, group)
    % Rebase path and load each node in the hierarchy
    parentPath = filesepStandard(root, 'dir');
    rebasePaths(group, parentPath);

    % Minimal popup logger
    fig = uifigure('Name','Status','Position',[100 100 520 360]);
    txa = uitextarea(fig, 'Position',[12 12 496 336], 'Editable','off', ...
                          'Value', {'📂 Loading data...'});
    drawnow;

    for c = 1:numel(group)
        msg = sprintf('Loading Group %d of %d...', c, numel(group));
        txa.Value(end+1) = {msg};

        err = safeLoad(group(c));
        checkLoad(txa, err, 0);  % <-- pass txa

        for i = 1:numel(group(c).subjs)
            msg = sprintf('\tSubject %d of %d...', i, numel(group(c).subjs));
            txa.Value(end+1) = {msg};

            err = safeLoad(group(c).subjs(i));
            checkLoad(txa, err, 1);

            for j = 1:numel(group(c).subjs(i).sess)
                msg = sprintf('\t\tSession %d of %d...', j, numel(group(c).subjs(i).sess));
                txa.Value(end+1) = {msg};

                err = safeLoad(group(c).subjs(i).sess(j));
                checkLoad(txa, err, 2);

                for k = 1:numel(group(c).subjs(i).sess(j).runs)
                    msg = sprintf('\t\t\tRun %d of %d...', k, numel(group(c).subjs(i).sess(j).runs));
                    txa.Value(end+1) = {msg};

                    err = safeLoad(group(c).subjs(i).sess(j).runs(k));
                    checkLoad(txa, err, 3);

                    % save all runs to ssdot
                    ssdot.setVar('Run', group(c).subjs(i).sess(j).runs(k));
                end
            end

            drawnow limitrate
        end
    end

    txa.Value(end+1) = {'✅ Done loading all data!'};
    drawnow;
end

% Set .path for group/subj/sess/run to the dataset root
function rebasePaths(group, root)
    group.path = filesepStandard(root, 'dir');
    for i = 1:numel(group.subjs)
        group.subjs(i).path = group.path;
        for j = 1:numel(group.subjs(i).sess)
            group.subjs(i).sess(j).path = group.path;
            for k = 1:numel(group.subjs(i).sess(j).runs)
                group.subjs(i).sess(j).runs(k).path = group.path;
            end
        end
    end
end

% -------- local helpers --------
function ok = safeLoad(node)
    try
        err = node.Load();
        fprintf('[DEBUG] %s.Load() returned: ', class(node));
        disp(err);
        if islogical(err)
            ok = err;
        else
            ok = (err == 0);
        end
    catch
        ok = false;
    end
end

function checkLoad(txa, result, level)
    if islogical(result)
        ok = result;
    else
        ok = (result == 0);
    end
    indent = repmat(' ', 1, 3*(level+1));  % 3 spaces per level
    if ok
        txa.Value(end+1) = {sprintf('%s✅ Loaded', indent)};
    else
        txa.Value(end+1) = {sprintf('%s⚠️ Failed / skipped', indent)};
    end
    drawnow limitrate
end
