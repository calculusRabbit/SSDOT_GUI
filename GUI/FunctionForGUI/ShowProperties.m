function showObjectInspector(rootValue, titleStr)
% showObjectInspector - interactive viewer for nested structs/objects/cells
% Usage:
%   showObjectInspector(app.ssdot)
%   showObjectInspector(app.ssdot, 'SS-DOT Inspector')

    if nargin < 2
        try
            titleStr = sprintf('%s Inspector', class(rootValue));
        catch
            titleStr = 'Inspector';
        end
    end

    % ---- Figure & layout ----
    fig = uifigure('Name', titleStr, 'Position', [100 100 900 520]);
    gl = uigridlayout(fig,[1 3]);
    gl.ColumnWidth = {280, '1x', 120};   % Tree | Table | Buttons
    gl.RowHeight   = {'1x'};

    % Tree
    tree = uitree(gl,'multiselect','off');
    tree.FontName = 'Consolas';
    tree.NodeExpandedFcn   = @(src,evt)onExpand(evt.Node);
    tree.SelectionChangedFcn = @(src,evt)showDetails(evt, tbl);

    % Details table
    tbl = uitable(gl,'ColumnName',{'Name','Size','Class','Preview'}, ...
                     'ColumnEditable',[false false false false]);
    tbl.FontName = 'Consolas';

    % Controls
    p = uipanel(gl,'BorderType','none');
    vb = uigridlayout(p,[6 1]); vb.RowHeight = {30,30,30,30,'1x',30};
    uilabel(vb,'Text','Controls','FontWeight','bold');

    uibutton(vb,'Text','Refresh', ...
        'ButtonPushedFcn', @(~,~)rebuildTree());

    uibutton(vb,'Text','Expand All', ...
        'ButtonPushedFcn', @(~,~)expandAll(tree));

    uibutton(vb,'Text','Collapse All', ...
        'ButtonPushedFcn', @(~,~)collapseAll(tree));

    pathField = uieditfield(vb,'text','Placeholder','Selected path...');
    pathField.Editable = false;

    uibutton(vb,'Text','Copy Path', ...
        'ButtonPushedFcn', @(~,~)copy(pathField.Value));

    % ---- Build root ----
    rebuildTree();

    % ===== nested helpers =====
    function rebuildTree()
        delete(tree.Children);
        makeNode(tree, rootValue, "<root>", []);
        tree.SelectedNodes = tree.Children(1);
        showDetails(struct('SelectedNodes',tree.SelectedNodes), tbl);
    end

    function node = makeNode(parent, val, name, path)
        % path is a cell array of subsrefs to reach val from the root
        node = uitreenode(parent);
        node.Text  = sprintf('%s  :  %s', name, shortPreview(val));
        node.NodeData = struct('path',{path}, 'name',name, 'class',classSafe(val));
        node.Icon = matlab.internal.webwindow.icon.getIconURI('folder'); %#ok<NASGU>
        if isExpandable(val)
            % add a placeholder so the chevron shows
            uitreenode(node,'Text','(double-click to expand)','NodeData',[]);
        end
        node.ContextMenu = uicontextmenu(fig);
        uimenu(node.ContextMenu,'Text','Expand here','MenuSelectedFcn',@(~,~)onExpand(node));
        uimenu(node.ContextMenu,'Text','Copy value preview','MenuSelectedFcn',@(~,~)copy(shortPreview(val)));
    end

    function onExpand(node)
        % Expand lazily. If already populated, skip.
        if ~isscalar(node.Children) || ~strcmp(node.Children(1).Text,'(double-click to expand)')
            return
        end
        % Clear placeholder
        delete(node.Children);

        val = getByPath(rootValue, node.NodeData.path);
        kids = enumerateChildren(val);
        for k = 1:numel(kids)
            child = kids(k);
            makeNode(node, child.value, child.name, [node.NodeData.path, child.subs]);
        end
    end

    function showDetails(evt, tableHandle)
        nd = evt.SelectedNodes;
        if isempty(nd); return; end
        nd = nd(1);
        pathField.Value = pathToString(nd.NodeData.path);

        val = getByPath(rootValue, nd.NodeData.path);
        rows = detailsRows(val);
        tableHandle.Data = rows;
        tableHandle.ColumnName = {'Name','Size','Class','Preview'};
        tableHandle.ColumnWidth = {220,100,160,'1x'};
    end

end

% ---------- utility functions (outside main) ----------

function cls = classSafe(v)
    try, cls = class(v); catch, cls = 'unknown'; end
end

function tf = isExpandable(v)
    tf = isstruct(v) || isobject(v) || iscell(v) || isa(v,'containers.Map') ...
        || (isnumeric(v) && ~isscalar(v)) || (islogical(v) && ~isscalar(v)) ...
        || istable(v);
end

function kids = enumerateChildren(v)
% Return struct array with fields: name, value, subs
    kids = struct('name',{},'value',{},'subs',{});
    if isstruct(v)
        fn = fieldnames(v);
        for i = 1:numel(fn)
            name = fn{i};
            value = v.(name);
            kids(end+1) = struct('name',name,'value',value,'subs',{{'.',name}}); %#ok<AGROW>
        end
    elseif isobject(v)
        try
            pn = properties(v);
        catch
            pn = {};
        end
        for i = 1:numel(pn)
            name = pn{i};
            try
                value = v.(name);
                kids(end+1) = struct('name',name,'value',value,'subs',{{'.',name}}); %#ok<AGROW>
            catch
                % unreadable property
            end
        end
    elseif istable(v)
        for i = 1:width(v)
            name  = v.Properties.VariableNames{i};
            value = v.(name);
            kids(end+1) = struct('name',name,'value',value,'subs',{{'.',name}}); %#ok<AGROW>
        end
    elseif iscell(v)
        for i = 1:numel(v)
            name = sprintf('{%d}', i);
            value = v{i};
            kids(end+1) = struct('name',name,'value',value,'subs',{{'{}',i}}); %#ok<AGROW>
        end
    elseif isa(v,'containers.Map')
        k = keys(v);
        for i = 1:numel(k)
            name = sprintf('("%s")',string(k{i}));
            value = v(k{i});
            kids(end+1) = struct('name',name,'value',value,'subs',{{'()',k{i}}}); %#ok<AGROW>
        end
    elseif (isnumeric(v) || islogical(v)) && ~isscalar(v)
        % show first level elements as linear indices
        n = min(numel(v), 200); % cap to prevent huge trees
        for i = 1:n
            name = sprintf('(%d)', i);
            value = v(i);
            kids(end+1) = struct('name',name,'value',value,'subs',{{'()',i}}); %#ok<AGROW>
        end
        if numel(v) > n
            kids(end+1) = struct('name',sprintf('... (%d more)', numel(v)-n), ...
                                 'value',[], 'subs',{{}});
        end
    end
end

function s = sizeStr(v)
    try
        sz = size(v);
        s = sprintf('%dx%d', sz(1), max(1, numel(sz) > 1 && sz(2) || 1));
        if numel(sz) > 2
            s = sprintf('%s x%s', s, sprintf('x%d', sz(3:end)));
        end
    catch
        s = '?';
    end
end

function pv = shortPreview(v)
    cls = classSafe(v);
    if isstruct(v)
        pv = sprintf('%s %s', sizeStr(v), 'struct');
    elseif istable(v)
        pv = sprintf('%dx%d table', size(v,1), size(v,2));
    elseif iscell(v)
        pv = sprintf('%s cell', sizeStr(v));
    elseif isa(v,'containers.Map')
        pv = sprintf('%d-key map', v.Count);
    elseif isstring(v)
        if isscalar(v), pv = '"' + v + '"'; else, pv = sprintf('%s string', sizeStr(v)); end
    elseif ischar(v)
        pv = ['''' v ''''];
    elseif isnumeric(v) || islogical(v)
        if isempty(v), pv = '[]';
        elseif isscalar(v), pv = num2str(v);
        else, pv = sprintf('%s %s', sizeStr(v), cls);
        end
    elseif isobject(v)
        pv = sprintf('%s object', class(v));
    else
        pv = cls;
    end
end

function rows = detailsRows(v)
% Convert a value into table rows of {Name, Size, Class, Preview}
    if isstruct(v) || isobject(v) || istable(v) || iscell(v) || isa(v,'containers.Map')
        kids = enumerateChildren(v);
        rows = cell(numel(kids),4);
        for i = 1:numel(kids)
            rows{i,1} = kids(i).name;
            rows{i,2} = sizeStr(kids(i).value);
            rows{i,3} = classSafe(kids(i).value);
            rows{i,4} = shortPreview(kids(i).value);
        end
    else
        rows = {
            '(value)', sizeStr(v), classSafe(v), shortPreview(v)
        };
    end
end

function out = getByPath(root, path)
% path is like {'.','A','.', 'Amatrix_brain','.', 'w1_HbO'} etc.
    if isempty(path)
        out = root; return
    end
    s = struct('type',{},'subs',{});
    for i = 1:2:numel(path)
        s(end+1).type = path{i}; %#ok<AGROW>
        s(end).subs   = path{i+1};
    end
    out = subsref(root, s);
end

function s = pathToString(path)
    if isempty(path), s = '(root)'; return; end
    buf = "";
    i = 1;
    while i <= numel(path)
        t = path{i}; a = path{i+1};
        if t=='.'
            if buf==""; buf = string(a); else; buf = buf + "." + string(a); end
        elseif t=='{}'
            buf = buf + sprintf('{%d}', a);
        elseif t=='()'
            if ischar(a) || isstring(a)
                buf = buf + sprintf('("%s")', string(a));
            else
                buf = buf + sprintf('(%d)', a);
            end
        end
        i = i+2;
    end
    s = char(buf);
end

function expandAll(tree)
    stack = tree.Children;
    while ~isempty(stack)
        n = stack(1); stack(1) = [];
        try, n.expand; catch, end %#ok<TRYNC>
        % trigger lazy expansion
        try, feval(tree.NodeExpandedFcn, tree, struct('Node',n)); catch, end %#ok<TRYNC>
        stack = [stack; n.Children]; %#ok<AGROW>
    end
end

function collapseAll(tree)
    stack = tree.Children;
    while ~isempty(stack)
        n = stack(1); stack(1) = [];
        try, n.collapse; catch, end %#ok<TRYNC>
        stack = [stack; n.Children]; %#ok<AGROW>
    end
end
