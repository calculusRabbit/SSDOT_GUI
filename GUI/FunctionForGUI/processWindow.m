function processWindow(app)
    loadH = @() loadAV(app, []);
    dispH = @() plotSD(app, []);   

    win = LoadDisplayWindow();
    win.addStep('AV', loadH, dispH);
end
