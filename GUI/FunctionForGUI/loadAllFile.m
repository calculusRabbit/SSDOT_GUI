function loadAllFile(app, root, group)
    % Rebase path and load each node in the hierarchy

    parentPath = filesepStandard(root, 'dir');
    rebasePaths(group, parentPath);

    app.StatusTextArea.Value = {'📂 Loading data...'};
    drawnow;

    for c = 1:numel(group)
        msg = sprintf('Loading Group %d of %d...', c, numel(group));
        app.StatusTextArea.Value(end+1) = {msg};

        err = safeLoad(group(c));
        checkLoad(app, err, 0);

        for i = 1:numel(group(c).subjs)
            msg = sprintf('\tSubject %d of %d...', i, numel(group(c).subjs));
            app.StatusTextArea.Value(end+1) = {msg};

            err = safeLoad(group(c).subjs(i));
            checkLoad(app, err, 1); % level 1 = subject

            for j = 1:numel(group(c).subjs(i).sess)
                msg = sprintf('\t\tSession %d of %d...', j, numel(group(c).subjs(i).sess));
                app.StatusTextArea.Value(end+1) = {msg};

                err = safeLoad(group(c).subjs(i).sess(j));
                checkLoad(app, err, 2); % level 2 = session

                for k = 1:numel(group(c).subjs(i).sess(j).runs)
                    msg = sprintf('\t\t\tRun %d of %d...', k, numel(group(c).subjs(i).sess(j).runs));
                    app.StatusTextArea.Value(end+1) = {msg};

                    err = safeLoad(group(c).subjs(i).sess(j).runs(k));
                    checkLoad(app, err, 3);
                end
            end

            % Optionally flush less often to speed UI:
            drawnow limitrate
        end
    end

    app.StatusTextArea.Value(end+1) = {'✅ Done loading all data!'};
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


% -------- local helpers (same file / methods block) --------
function ok = safeLoad(node)
    % Wrap Load() to return true/false regardless of errcode/exception.
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

function checkLoad(app, result, level)
    % result can be boolean ok or errcode (0=ok)
    if islogical(result)
        ok = result;
    else
        ok = (result == 0);
    end
    indent = repmat(' ', 1, 3*(level+1));  % 3 spaces per level
    if ok
        app.StatusTextArea.Value(end+1) = {sprintf('%s✅ Loaded', indent)};
    else
        app.StatusTextArea.Value(end+1) = {sprintf('%s⚠️ Failed / skipped', indent)};
    end
    % Use limitrate to keep UI responsive without over-updating
    drawnow limitrate
end
